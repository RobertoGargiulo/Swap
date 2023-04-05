program spectrum_Swap

  use functions, only: dimSpin1_Sz0, init_random_seed
  use matrices, only: buildSz0_HSwap, buildSz0_HMBL, &
    & print_hamiltonian_Sz0, print_unitary_Sz0
  use observables, only: gap_ratio
  use printing
  use exponentiate, only: diagSYM, expSYM, diagUN
  use sorts, only: dpquicksort
  implicit none

  complex (c_double_complex), parameter :: C_UNIT = dcmplx(0._c_double, 1._c_double)
  real (c_double), parameter :: pi = 4.d0 * datan(1.d0)
  integer (c_int), parameter :: dimSpin1 = 3

  integer (c_int) :: nspin, dim, dim_Sz0, Sz
  integer (c_int), allocatable :: config(:), idxSz0(:)
  integer (c_int) :: i, j, k, n, l

  complex (c_double_complex) :: alpha(dimSpin1)

  integer (c_int) :: n_disorder
  real (c_double) :: hz_coupling, J_coupling, V_coupling, T0, T1, kick
  real (c_double), allocatable :: H(:,:), hz(:), Jxy(:), Vz(:)
  real (c_double), allocatable :: E(:), W_r(:,:), E_MBL(:,:), QE(:,:)
  complex (c_double_complex), allocatable :: USwap(:,:), PH(:), W(:,:), U(:,:)

  real (c_double), dimension(:), allocatable :: r_avg, r_sq, r_avg2, r_sq2
  real(c_double) :: r_dis_avg, r_dis_sigma, r_dis_avg2, r_dis_sigma2

  integer(c_int) :: count_beginning, count_end, count_rate

  character(len=200) :: filestring_qe
  character(len=:), allocatable :: columns
  integer (c_int) :: unit_qe

  logical :: SELECT
  EXTERNAL SELECT

  !--------------- (Explicit) Interfaces ------------!
  interface
    function column_titles(nspin) result(columns)
      integer, intent(in) :: nspin
      character(26*(4*nspin+3)) :: columns
    end function
  end interface

  !------------ Input ----------------!

  write (*,*) "Number of Spins"
  read (*,*) nspin
  print *, ""

  write (*,*) "Number of Disorder Realizations"
  read (*,*) n_disorder
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

  write(filestring_qe,93) "data/dynamics/quasienergies_Swap_nspin", &
    & nspin, "_period", T0, "_n_disorder", n_disorder, &
    & "_J", J_coupling, "_V", int(V_coupling), V_coupling-int(V_coupling), &
    & "_hz", int(hz_coupling), hz_coupling-int(hz_coupling), &
    & "_kick", kick, ".txt"

  93  format(A,I0, A,F4.2, A,I0, A,F4.2, A,I0,F0.2, A,I0,F0.2, A,F4.2, A)

  !open(newunit=unit_qe,file=filestring_qe)
  !call write_info(unit_qe)


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

  allocate(hz(nspin), Jxy(nspin-1), Vz(nspin-1))
  allocate(U(dim_Sz0,dim_Sz0), H(dim_Sz0,dim_Sz0), E(dim_Sz0), W_r(dim_Sz0,dim_Sz0))

  allocate(PH(dim_Sz0),W(dim_Sz0,dim_Sz0))
  allocate(E_MBL(n_disorder,dim_Sz0), QE(n_disorder,dim_Sz0))

  allocate( r_avg(n_disorder), r_sq(n_disorder))
  allocate( r_avg2(n_disorder), r_sq2(n_disorder))

  r_avg = 0
  r_sq = 0
  r_avg2 = 0
  r_sq2 = 0

  !Time evolution over various realizations of disorder
  !$OMP PARALLEL
  call init_random_seed() 
  !print *, "Size of Thread team: ", omp_get_num_threads()
  !print *, "Verify if current code segment is in parallel: ", omp_in_parallel()
  !$OMP do private(i, hz, Vz, &
  !$OMP & H, E, W_r, U, W, PH )
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

    E_MBL(i,:) = E
    call gap_ratio(dim_Sz0, E, r_avg2(i), r_sq2(i))

    call expSYM( dim_Sz0, -C_UNIT*T0, E, W_r, U )
    U = matmul(USwap,U)

    call diagUN( SELECT, dim_Sz0, U, PH, W)
    E = real(C_UNIT*log(PH))
    call dpquicksort(E)
    call gap_ratio(dim_Sz0, E, r_avg(i), r_sq(i))
    QE(i,:) = E


  enddo
  !$OMP END DO
  !$OMP END PARALLEL 

  !write (unit_qe,"(*(A))") 
  !write(unit_ent, "(A12,A,*(A26))") ""Disorder Realization" , repeat(trim(string),nspin), "LI", "IPR"!, "MBE", "I^2", "CORR_Z", "E"

  !print "(A12,*(A26))", "Sample_Dis", "Quasienergy", "Energy MBL"
  !do i = 1, n_disorder
  !  do l = 1, dim_Sz0
  !    write (*,*) i, QE(i,l), E_MBL(i,l)
  !    !write (unit_qe,*) i, QE(i,l), E_MBL(i,l)
  !  enddo
  !enddo

  r_dis_avg = sum(r_avg) / n_disorder
  r_dis_sigma = sqrt( ( sum(r_sq)/n_disorder - r_dis_avg**2 ) / n_disorder )
  r_dis_avg2 = sum(r_avg2) / n_disorder
  r_dis_sigma2 = sqrt( ( sum(r_sq2)/n_disorder - r_dis_avg2**2 ) / n_disorder )

  print *, "Average and Variance of Gap Ratio (over the spectrum and then disorder)"
  print "(*(A26))", "<r>_Swap", "sigma(r)_Swap", "<r>_MBL", "sigma(r)_MBL"
  print *, r_dis_avg, r_dis_sigma, r_dis_avg2, r_dis_sigma2

  !close(unit_qe)

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


subroutine write_info(unit_file)

  integer, intent(in) :: unit_file

  write (unit_file,*) "Some info: "
  write (unit_file,*) "Quasi-Energies of Floquet Operator U_F = e^(-i H) U_swap."
  write (unit_file,*) "Spin-1 chain with hamiltonian H = sum hz * Z + V * ZZ + J * (XX + YY)."
  write (unit_file,*) "Periodic perturbed swap, U_swap = exp(-i(pi/2 + kick) * sum (sigma*sigma)^2 - (sigma*sigma)."
  write (unit_file,*) "V is taken in [-3V/2, V/2];  hz is taken in [-hz, hz];  J is uniform."
  write (unit_file,*) "Exact diagonalization of the dense matrix has been used to compute U_F and diagonalize it."

end subroutine
