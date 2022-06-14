program average_magnetization

  use iso_c_binding

  implicit none

  integer (c_int) :: i, j, k, u_int, u_fluct, iteration, steps, nspin, n_iterations
  real (c_double) :: mag, h_coupling, kick
  real (c_double), allocatable, dimension(:) :: avg, fluct
  character(len=100) :: filestring, string, commandstring


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

  write (*,*) "Coupling Constant, h = h_x / J"
  read (*,*) h_coupling
  print*,""
  
  write (*,*) "Perturbation on Kick, epsilon = T1 - pi/2"
  read (*,*) kick
  print*,""


  !----------------------------
  !write(filestring,91) "data/magnetizations/Swap_Sz_nspin", nspin, "_steps", steps    , &
  !&  "_iterations", n_iterations, ".txt"
  write(filestring,92) "data/magnetizations/Swap_Sz_nspin", nspin, "_steps", steps, &
    &  "_iterations", n_iterations, "_h", h_coupling, "_kick", kick, ".txt"
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
      read(u_int, *) mag, k
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
  write(filestring,92) "data/magnetizations/Swap_Sz_AVG_nspin", nspin, "_steps", steps, &
    &  "_iterations", n_iterations, "_h", h_coupling, "_kick", kick, ".txt"
  open(newunit=u_int, file=filestring)


  write(filestring,92) "data/magnetizations/Swap_Sz_FLUCT_nspin", nspin, "_steps", steps, &
    &  "_iterations", n_iterations, "_h", h_coupling, "_kick", kick, ".txt"
  open(newunit=u_fluct, file=filestring)

  do j = 1, steps
    !print *, avg(j)
    write (u_int,*) avg(j), j
    write (u_fluct,*) fluct(j), j
  enddo
  close(u_int)
  close(u_fluct)
  !--------------------


  91  format(A,I0,A,I0,A,I0,A)
  92  format(A,I0, A,I0, A,I0, A,F4.2, A,F4.2, A)

end program