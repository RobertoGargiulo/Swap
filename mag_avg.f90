program average_magnetization

  use iso_c_binding

  implicit none

  integer (c_int) :: i, j, k, u_int, u_fluct, iteration, steps, nspin, n_iterations
  real (c_double) :: J_coupling, V_coupling, h_coupling, hz_coupling, kick
  real (c_double) :: mag, T0, time
  real (c_double), allocatable, dimension(:) :: avg, fluct
  character(len=200) :: filestring, string, commandstring


  !---------------------------
  write (*,*) "Number of Spins"
  read (*,*) nspin
  print*,""

  write (*,*) "Number of Steps"
  read (*,*) steps
  print*,""

  write (*,*) "Number of Iterations"
  read (*,*) n_iterations
  print*,""

  write (*,*) "Time Step/Period"
  read (*,*) T0
  print*,""

  write (*,*) "Longitudinal Interaction Constant J * ZZ"
  read (*,*) J_coupling
  print*,""

  write (*,*) "Transverse Interaction Constant V * (XX + YY)"
  read (*,*) V_coupling
  print*,""

  write (*,*) "Transverse Field h_x * X"
  read (*,*) h_coupling
  print*,""

  write (*,*) "Longitudinal Field h_z * Z"
  read (*,*) hz_coupling
  print*,""

  write (*,*) "Perturbation on Kick, epsilon = T1 - pi/4"
  read (*,*) kick
  print*,""


  !----------------------------
  !write(filestring,91) "data/magnetizations/Swap_Sz_nspin", nspin, "_steps", steps    , &
  !&  "_iterations", n_iterations, ".txt"
  write(filestring,92) "data/magnetizations/Clean_MBL_Imbalance_nspin", nspin, "_steps", steps, &
    &  "_iterations", n_iterations, "_J", J_coupling, "_V", V_coupling,"_h", h_coupling, "_hz", hz_coupling, &
   & "_no_kick", kick, ".txt"

  
  !  commandstring = "awk 'BEGIN{m=0}{m=$2}END{print m}' " // trim(filestring) // " > steps"

!  call execute_command_line(commandstring)
!
!  open (10, file="steps")
!  read (10, *) steps
!  close (10)
!  
!  call execute_command_line("rm steps")

!  commandstring = "grep iter " // trim(filestring) // " | tail -1 | awk '{print $3}' > iteration"
!  call execute_command_line(commandstring)
!  open (10, file="iteration")
!  read (10,*) iteration
!  close (10)
!  call execute_command_line("rm iteration")

  allocate(avg(steps), fluct(steps))
  !------------------------------

  !Compute vector of averages 'avg(j=1, ..., steps)
  open(newunit=u_int, file=filestring)

  avg = 0
  fluct = 0
  do i = 1, n_iterations
    read(u_int, *) string
    do j = 1, steps
      read(u_int, *) mag, time
      avg(j) = avg(j) + mag
      fluct(j) = fluct(j) + mag*mag
    enddo
  enddo
  avg = avg/real(n_iterations)
  fluct = fluct/real(n_iterations)
  fluct = fluct - avg**2
  close(u_int)
  !-------------------------

  !Write data to file
  !write(filestring,91) "data/magnetizations/Swap_Sz_AVG_nspin", nspin, "_steps", steps    , &
  !&  "_iterations", n_iterations, ".txt"

  write(filestring,92) "data/magnetizations/Clean_MBL_Imbalance_AVG_nspin", nspin, "_steps", steps, &
    &  "_iterations", n_iterations, "_J", J_coupling, "_V", V_coupling,"_h", h_coupling, "_hz", hz_coupling, &
   & "_no_kick", kick, ".txt"

!  write(filestring,92) "data/magnetizations/Swap_Sz_AVG_nspin", nspin, "_steps", steps, &
!    &  "_iterations", n_iterations, "_J", J_coupling, "_V", V_coupling,"_h", h_coupling, "_hz", hz_coupling, "_kick", kick, ".txt"
  open(newunit=u_int, file=filestring)

  write(filestring,92) "data/magnetizations/Clean_MBL_Imbalance_FLUCT_nspin", nspin, "_steps", steps, &
    &  "_iterations", n_iterations, "_J", J_coupling, "_V", V_coupling,"_h", h_coupling, "_hz", hz_coupling, &
   & "_no_kick", kick, ".txt"

!  write(filestring,92) "data/magnetizations/Swap_Sz_FLUCT_nspin", nspin, "_steps", steps, &
!    &  "_iterations", n_iterations, "_J", J_coupling, "_V", V_coupling,"_h", h_coupling, "_hz", hz_coupling, "_kick", kick, ".txt"
  open(newunit=u_fluct, file=filestring)

  do j = 1, steps
    !print *, avg(j)
    write (u_int,*) avg(j), j*T0
    write (u_fluct,*) fluct(j), j*T0
  enddo
  close(u_int)
  close(u_fluct)
  !--------------------


  91  format(A,I0,A,I0,A,I0,A)
  92  format(A,I0, A,I0, A,I0, A,F4.2, A,F4.2, A,F4.2, A,F4.2, A,F4.2, A)
end program
