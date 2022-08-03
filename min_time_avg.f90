program min_time_avg

  use genmat
  use iso_c_binding

  implicit none

  integer (c_int)     ::  nspin, iteration, steps, n_iterations, start
  integer (c_int)     ::  i, j, k, p
  integer (c_int)     ::  unit_mag, unit_ph, unit_w, unit_avg

  real(c_double) :: T0, T1, J_coupling, V_coupling, hz_coupling, kick 
  real (c_double) :: norm, time_j, t_avg, t_sigma
  real (c_double), dimension(:), allocatable :: avg, sigma
  character(len=200) :: filestring



  read (*,*) nspin

  read (*,*) n_iterations

  read (*,*) steps

  read (*,*) krylov_dim

  read (*,*) T0
  
  read (*,*) J_coupling

  read (*,*) V_coupling

  read (*,*) hz_coupling



  write(filestring,91) "data/magnetizations/Sz0_DENSE_MBL_hz_Disorder_AVG_FLUCT_Imbalance_nspin", &
    & nspin, "_steps", steps, "_time_step", T0, &
    &  "_iterations", n_iterations, "_J", J_coupling, "_V", V_coupling, "_hz", hz_coupling, ".txt"
  open(newunit=unit_avg,file=filestring)
  91  format(A,I0, A,I0, A,I0, A,F4.2, A,F4.2, A,F4.2, A,F4.2, A)

  !print *, filestring
  
  !Allocate observables and averages
  !allocate( avg(steps), sigma(steps) )

  do j = 1, steps
    read(unit_avg,*) 
  enddo
  
  start = int(100/T0)
  call time_avg(steps, start, avg, sigma, t_avg, t_sigma)
  print "(3(F4.2,2X), 2(F6.3,2X), I0)", J_coupling, V_coupling, hz_coupling, t_avg, t_sigma, start

  close(unit_avg)



end program min_time_avg
