program prova

  use functions, only : dimSpin1_Sz0, init_random_seed
  use matrices, only: buildSz0_HSwap, buildSz0_HMBL, &
    & print_hamiltonian_Sz0, print_unitary_Sz0
  use states, only: buildstate => buildNeelState_Sz0, printstate_Sz0
  use observables, only: sigmaz_Sz0
  use printing
  use exponentiate, only: diagSYM, expSYM, diagUN
  implicit none

  complex (c_double_complex), parameter :: C_UNIT = dcmplx(0._c_double, 1._c_double)
  real (c_double), parameter :: pi = 4.d0 * datan(1.d0)
  integer (c_int), parameter :: dimSpin1 = 3

  integer (c_int) :: nspin, dim, dim_Sz0, Sz
  integer (c_int), allocatable :: config(:), idxSz0(:)
  integer (c_int) :: i, j, k, n 

  complex (c_double_complex) :: alpha(dimSpin1)

  real (c_double) :: norm
  complex (c_double_complex), allocatable :: init_state(:), state(:), psi_swap(:), psi_mbl(:)
  integer (c_int) :: n_disorder, steps
  real (c_double) :: hz_coupling, J_coupling, V_coupling, T0, T1, kick
  real (c_double), allocatable :: H(:,:), hz(:), Jxy(:), Vz(:)
  real (c_double), allocatable :: E(:), W_r(:,:)
  complex (c_double_complex), allocatable :: USwap(:,:), PH(:), W(:,:), U(:,:), UMBL(:,:)

  real (c_double), allocatable :: imb_avg(:), imb_sq(:)
  real (c_double), allocatable :: sigmaz_avg(:,:,:), sigmaz_sq(:,:,:)

  integer(c_int) :: count_beginning, count_end, count_rate

  character(len=200) :: filestring_Neel, filestring_imb, filestring_sigmaz
  character(len=:), allocatable :: columns
  integer (c_int) :: unit_sigmaz

  logical :: SELECT
  EXTERNAL SELECT


  !------------ Input ----------------!

  write (*,*) "Number of Spins"
  read (*,*) nspin
  print *, ""

  write (*,*) "Number of Disorder Realizations"
  read (*,*) n_disorder
  print*,""

  write (*,*) "Number of Steps"
  read (*,*) steps
  print*,""

  write (*,*) "Period T0"
  read (*,*) T0
  print*,""
 
  write (*,*) "Transverse Interaction Constant -J * (XX + YY)"
  read (*,*) J_coupling
  print*,""

  write (*,*) "Longitudinal Interaction Constant -V * ZZ"
  read (*,*) V_coupling
  print*,""

  write (*,*) "Longitudinal Field h_z * Z"
  read (*,*) hz_coupling
  print*,""

  write (*,*) "Perturbation on Kick T1 = pi/4 + kick"
  read (*,*) kick
  print*,""

  call system_clock(count_beginning, count_rate)

  dim = dimSpin1**nspin
  dim_Sz0 = dimSpin1_Sz0(nspin)
  T1 = pi/2 + kick

  !------------------- File names ------------------------!

  write(filestring_Neel,93) "data/dynamics/sigmaz_Swap_vs_MBL_Neel_nspin", &
    & nspin, "_steps", steps, "_period", T0, "_n_disorder", n_disorder, &
    & "_J", J_coupling, "_V", int(V_coupling), V_coupling-int(V_coupling), &
    & "_hz", int(hz_coupling), hz_coupling-int(hz_coupling), &
    & "_kick", kick, ".txt"

  93  format(A,I0, A,I0, A,F4.2, A,I0, A,F4.2, A,I0,F0.2, A,I0,F0.2, A,F4.2, A)

  open(newunit=unit_sigmaz,file=filestring_Neel)

  call write_info(unit_sigmaz, "Neel")


  !-------- Swap Operator --------------------!
  allocate(H(dim_Sz0,dim_Sz0))

  call buildSz0_HSwap(nspin, dim_Sz0, H)
  !call print_hamiltonian_Sz0(nspin, dim_Sz0, H)

  allocate(E(dim_Sz0), W_r(dim_Sz0,dim_Sz0), USwap(dim_Sz0,dim_Sz0))
  call diagSYM( 'V', dim_Sz0, H, E, W_r)
  deallocate(H)

  call expSYM( dim_Sz0, -C_UNIT*T1, E, W_r, USwap )
  deallocate(E, W_r)
  !call print_unitary_Sz0(nspin, dim_Sz0, USwap)

  !------------- Initial State ------------------!
  allocate(init_state(dim_Sz0))
  call buildstate(nspin, dim_Sz0, init_state)
  !call printstate_Sz0(nspin, dim_Sz0, init_state, "Neel state:")


  !Allocate state, evolution operators, observables
  allocate(psi_mbl(dim_Sz0), psi_swap(dim_Sz0))

  allocate(hz(nspin), Jxy(nspin-1), Vz(nspin-1))
  allocate(U(dim_Sz0,dim_Sz0), UMBL(dim_Sz0,dim_Sz0), H(dim_Sz0,dim_Sz0), E(dim_Sz0), W_r(dim_Sz0,dim_Sz0))

  !allocate(imb_avg(steps), imb_sq(steps))
  allocate(sigmaz_avg(steps,nspin,2), sigmaz_sq(steps,nspin,2))
  
  !imb_avg = 0
  !imb_sq = 0
  sigmaz_avg = 0
  sigmaz_sq = 0


  !Time evolution over various realizations of disorder
  !$OMP PARALLEL
  call init_random_seed() 
  !print *, "Size of Thread team: ", omp_get_num_threads()
  !print *, "Verify if current code segment is in parallel: ", omp_in_parallel()
  !$OMP do reduction(+: sigmaz_avg, sigmaz_sq) private(i, hz, Vz, norm, j, &
  !$OMP & psi_mbl, psi_swap, H, E, W_r, U, UMBL )
  do i = 1, n_disorder

    if (n_disorder < 10) then
     print *, "Current disorder realization = ", i
    else
      if (mod(i, n_disorder/10)==0) then 
        print *, "Current disorder realization = ", i
      endif
    endif

    !--------------- H_MBL Parameters ------------!
    !call random_number(Jxy)
    call random_number(Vz)
    call random_number(hz)

    hz = hz_coupling * 2*(hz-0.5) !hz in [-h, h]
    Jxy = -J_coupling
    Vz = -V_coupling + V_coupling * (Vz - 0.5) !Vz in [-3/2V, -V/2] <=> -V \pm V/2

    !------------ Floquet Operator(s) --------------!
    call buildSz0_HMBL(nspin, dim_Sz0, Jxy, Vz, hz, H)
    call diagSYM( 'V', dim_Sz0, H, E, W_r)
    call expSYM( dim_Sz0, -C_UNIT*T0, E, W_r, UMBL )
    U = matmul(USwap,UMBL)


    !allocate(PH(dim_Sz0),W(dim_Sz0,dim_Sz0))
    !call diagUN( SELECT, dim_Sz0, U, PH, W)

    psi_swap = init_state
    psi_mbl = init_state
    j = 1
    sigmaz_avg(j,:,1) = sigmaz_avg(j,:,1) + sigmaz_Sz0(nspin, dim_Sz0, psi_mbl)
    sigmaz_sq(j,:,1) = sigmaz_sq(j,:,1) + sigmaz_Sz0(nspin, dim_Sz0, psi_mbl) ** 2
    sigmaz_avg(j,:,2) = sigmaz_avg(j,:,2) + sigmaz_Sz0(nspin, dim_Sz0, psi_swap)
    sigmaz_sq(j,:,2) = sigmaz_sq(j,:,2) + sigmaz_Sz0(nspin, dim_Sz0, psi_swap) ** 2
    !print *, i, j, sigmaz_avg(j,:,:)
    do j = 2, steps

      if(mod(j,10)==0) then
        norm = real(dot_product(psi_mbl,psi_mbl))
        psi_mbl = psi_mbl/sqrt(norm)
        norm = real(dot_product(psi_swap,psi_swap))
        psi_swap = psi_swap/sqrt(norm)
      endif
      psi_mbl = matmul(UMBL, psi_mbl)
      psi_swap = matmul(U, psi_swap)
      sigmaz_avg(j,:,1) = sigmaz_avg(j,:,1) + sigmaz_Sz0(nspin, dim_Sz0, psi_mbl)
      sigmaz_sq(j,:,1) = sigmaz_sq(j,:,1) + sigmaz_Sz0(nspin, dim_Sz0, psi_mbl) ** 2
      sigmaz_avg(j,:,2) = sigmaz_avg(j,:,2) + sigmaz_Sz0(nspin, dim_Sz0, psi_swap)
      sigmaz_sq(j,:,2) = sigmaz_sq(j,:,2) + sigmaz_Sz0(nspin, dim_Sz0, psi_swap) ** 2
      !print *, i, j, sigmaz_avg(j,:,:)

    enddo

  enddo
  !$OMP END DO
  !$OMP END PARALLEL 

  sigmaz_avg = sigmaz_avg / n_disorder
  sigmaz_sq = sqrt( (sigmaz_sq/n_disorder - sigmaz_avg**2) / n_disorder )
  !write (unit_sigmaz,"(*(A))") 
  !write(string,"(A26)") "sigma_k^z"
  !write(unit_ent, "(A12,A,*(A24,2X))") "Iteration" , repeat(trim(string),nspin), "LI", "IPR"!, "MBE", "I^2", "CORR_Z", "E"
  !columns = column_titles(nspin)
  !print *, columns
  do j = 1, steps
    !write (*,*) j*T0, sigmaz_avg(j,:,1), sigmaz_sq(j,:,1)
    write (unit_sigmaz,*) j*T0, sigmaz_avg(j,:,:), sigmaz_sq(j,:,:)
  enddo


  close(unit_sigmaz)

  call take_time(count_rate, count_beginning, count_end, 'T', "Program")

end program

logical function SELECT(z)

  implicit none

  complex(kind=8), intent(in) :: z

  SELECT = .TRUE.
  RETURN

end

function column_titles(nspin) result(columns)

  use iso_c_binding
  implicit none

  integer (c_int), intent(in) :: nspin
  character (26*(4*nspin+3)):: columns
  character (nspin) :: s_index
  integer (c_int) :: i, j1, j2, j3, j4


  columns = ""
  do i = 1, nspin
    write(s_index, "(I0)") i
    j1 = 26*(i-1)+2
    j2 = 26*i+2
    j3 = 26*(i+nspin-1)+3
    j4 = 26*(i+nspin)+3
    print *, j1, j2, j3, j4
    columns(j1:j2) = trim("sigma_z^")//trim(s_index)
    columns(j3:j4) = trim("Var(sigma_z^")//trim(s_index)//trim(")")
    !print *, string
  enddo
  write(columns,"(3X,A4,20X,A)") "j*T0", trim(columns)

end function


subroutine write_info(unit_file, state_name)

  integer, intent(in) :: unit_file
  character(len=*) :: state_name

  write (unit_file,*) "Some info: "
  write (unit_file,*) "Dynamics of average of magnetization at integer multiples of the period."
  write (unit_file,*) "Spin-1 chain with hamiltonian H = sum hz * Z + V * ZZ + J * (XX + YY)."
  write (unit_file,*) "Periodic perturbed swap, U_swap = exp(-i(pi/2 + kick) * sum (sigma*sigma)^2 - (sigma*sigma)."
  write (unit_file,*) "V is taken in [-3V/2, V/2];  hz is taken in [-hz, hz];  J is uniform."
  write (unit_file,*) "Exact diagonalization of the dense matrix has been used to compute U(t)."
  write (unit_file,*) "Initial state is "//trim(state_name)//trim(".")

end subroutine
