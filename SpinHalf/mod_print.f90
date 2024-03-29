module printing_subrtns

  use iso_c_binding

  implicit none

contains

  subroutine printmat_C(dim,M,t)

    integer, intent(in) :: dim
    complex(c_double_complex), intent(in), dimension(dim,dim) :: M
    character, intent(in) :: t*1

    integer :: i

    if (t == 'R') then
      do i = 1,dim
        write (*,"('|',*(f5.1))",advance='no') real(M(i,:))
        print *, "|"
      enddo
    else if (t == 'C') then
      do i = 1,dim
        write (*,97,advance='no') M(i,:)
        print *, "|"
        97 format('|',(*(sf8.2,spf8.2,1x'i':1x)))
      enddo
    else if (t == 'I') then
      do i = 1, dim
        write (*,96,advance='no') int(M(i,:))
        print *,"|"
        96 format('|',(*(1X,I0)))
        !print*,""
      enddo
    else if (t == 'A') then
      do i = 1,dim
        write (*,"('|',*(f5.1))",advance='no') abs(M(i,:))
        print *, "|"
      enddo
    end if
    print *,""

  end subroutine printmat_C

  subroutine printmat_R(dim,M,t)

    integer, intent(in) :: dim
    real(c_double), intent(in), dimension(dim,dim) :: M
    character, intent(in) :: t*1

    integer :: i

    if (t == 'R') then
      do i = 1,dim
        write (*,"('|',*(f5.1))",advance='no') M(i,:)
        print *, "|"
      enddo
    else if (t == 'I') then
      do i = 1, dim
        write (*,96,advance='no') int(M(i,:))
        print *,"|"
        96 format('|',(*(1X,I0)))
        !print*,""
      enddo
    end if
    print *,""

  end subroutine printmat_R

  subroutine writemat_C(u, dim, M, t)

    integer, intent(in) :: dim, u
    complex(c_double_complex), intent(in), dimension(dim,dim) :: M
    character, intent(in) :: t*1

    integer :: i

    if (t == 'C') then
      do i = 1,dim
        write (u,*) M(i,:)
      enddo
    else if (t == 'R') then
      do i = 1,dim
        write (u,*) real(M(i,:))
      enddo
    else if (t == 'I') then
      do i = 1, dim
        write (u,*) int(M(i,:))
      enddo
    end if
    !print *,""

  end subroutine writemat_C

  subroutine writemat_R(u, dim, M, t)
    integer, intent(in) :: dim, u
    real(c_double), intent(in), dimension(dim,dim) :: M
    character, intent(in) :: t*1

    integer :: i

    if (t == 'R') then
      do i = 1,dim
        write (u,*) M(i,:)
      enddo
    else if (t == 'I') then
      do i = 1, dim
        write (u,*) int(M(i,:))
      enddo
    end if
    !print *,""

  end subroutine writemat_R

  subroutine writevec_C(u, dim, V, t)

    integer, intent(in) :: dim, u
    complex(c_double_complex), intent(in), dimension(dim) :: V
    character, intent(in) :: t*1

    integer :: i

    if (t == 'I') then
      do i = 1,dim
        write (u,*) int(V(i))
      enddo
    else if (t == 'R') then
      do i = 1, dim
        write (u,*) real(V(i))
      enddo
    else if (t == 'C') then
      do i = 1,dim
        write (u,*) V(i)
      enddo
    else if (t == 'G') then
      write (u,'(f15.10,f15.10)') V(:)
    end if
    !print *,""

  end subroutine writevec_C


  subroutine writevec_R(u, dim, V, t)

    integer, intent(in) :: dim, u
    real (c_double), intent(in), dimension(dim) :: V
    character, intent(in) :: t*1

    integer :: i

    if (t == 'I') then
      do i = 1,dim
        write (u,*) int(V(i))
      enddo
    else if (t == 'R') then
      do i = 1, dim
        write (u,*) V(i)  
      enddo
    else if (t == 'G') then
      do i = 1, dim
        write (u,"(f15.10)") V(i)
      enddo
    end if
    !print *,""

  end subroutine writevec_R

  subroutine printvec_C(dim, V, t)

    integer, intent(in) :: dim
    complex(c_double_complex), intent(in), dimension(dim) :: V
    character, intent(in) :: t*1

    integer :: i

    if (t == 'R') then
      do i = 1,dim
        write (*,"('|',*(f6.2))",advance='no') real(V(i))
        print *, "|"
      enddo
    else if (t == 'C') then
      do i = 1,dim
        write (*,97,advance='no') V(i)
        print *, "|"
        97 format('|',(*(sf8.2spf8.2,1x'i':1x)))
      enddo
    else if (t == 'I') then
      do i = 1, dim
        write (*,96,advance='no') int(V(i))
        print *,"|"
        96 format('|',(*(1X,I0)))
        !print*,""
      enddo
    else if (t == 'A') then
      do i = 1,dim
        write (*,"('|',*(f6.2))",advance='no') abs(V(i))
        print *, "|"
      enddo
    end if
    print *,""

  end subroutine printvec_C

  subroutine printvec_R(dim, V, t)

    integer, intent(in) :: dim
    real (c_double), intent(in), dimension(dim) :: V
    character, intent(in) :: t*1

    integer :: i

    if (t == 'R') then
      do i = 1,dim
        write (*,"('|',*(f6.2))",advance='no') V(i)
        print *, "|"
      enddo
    else if (t == 'I') then
      do i = 1, dim
        write (*,96,advance='no') int(V(i))
        print *,"|"
        96 format('|',(*(1X,I0)))
        !print*,""
      enddo
    end if
    print *,""

  end subroutine printvec_R

  subroutine printvec_I(dim, V, t)

    integer, intent(in) :: dim
    integer (c_int), intent(in), dimension(dim) :: V
    character, intent(in) :: t*1

    integer :: i

    if (t == 'I') then
      do i = 1, dim
        write (*,96,advance='no') V(i)
        print *,"|"
        96 format('|',(*(1X,I0)))
        !print*,""
      enddo
    end if
    print *,""

  end subroutine printvec_I

  subroutine take_time(count_rate, count_start, count_end, opt, filestring)
    implicit none
    integer (c_long), intent(in) :: count_rate, count_start
    integer (c_long), intent(out) :: count_end
    character :: opt*1
    character (*) :: filestring
  
    real (c_double) :: time_s
    integer(c_long) :: time_min, time_hour, time_day
  
    call system_clock(count_end)
  
    time_s = real(count_end - count_start) / real(count_rate)
    !print *, "time_s = ", time_s
    time_day = int(time_s/(60*60*24),kind(time_day))
    time_s = time_s - 60*60*24*time_day
    time_hour = int(time_s/(60*60),kind(time_day))
    time_s = time_s - 60*60*time_hour
    time_min = int(time_s/(60),kind(time_day))
    time_s = time_s - 60*time_min
    !print *, "time_s = ", time_s

    if(opt == 'T') then
      print "(A,A,A,3(1X,I0,A,1X),1X,F14.10,A)", "Elapsed Time for ", filestring, ": ", time_day, " days", time_hour, " hours", time_min, " min", time_s, " sec"
    else if(opt == 'F') then
      print *, ""
    endif
    !print *, ""
  end subroutine take_time

end module printing_subrtns

module printing
  use printing_subrtns

  interface printmat
    module procedure printmat_R, printmat_C
  end interface

  interface writemat
    module procedure writemat_R, writemat_C
  end interface

  interface writevec
    module procedure writevec_R, writevec_C
  end interface

  interface printvec
    module procedure printvec_R, printvec_C, printvec_I
  end interface
  



end module printing
