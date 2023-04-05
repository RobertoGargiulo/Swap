program integerkinds
    use iso_fortran_env
    use iso_c_binding

    implicit none

    integer :: i
    integer(kind=int8)  :: i8
    integer(kind=int16) :: i16
    integer(kind=int32) :: i32
    integer(kind=int64) :: i64
    integer (c_int) :: i_c
    integer (c_long) :: i_c_long
    real(c_double) :: t_double
    real (c_long_double) :: t_long_double

    integer(kind=selected_int_kind(6)) :: j6
    integer(kind=selected_int_kind(15)):: j15

    t_double = 1
    t_long_double = 1

    print *,'Default:'
    print *, huge(i)
    print *,'Int8:'
    print *, huge(i8)
    print *,'Int16:'
    print *, huge(i16)
    print *,'Int32:'
    print *, huge(i32)
    print *,'Int64:'
    print *, huge(i64)
    print *, "c_int:"
    print *, huge(i_c)
    print *, "c_long_int:"
    print *, huge(i_c_long)

    print *, "real double"
    print *, huge(i_c_long)*t_double
    print *, "real long double"
    print *, huge(i_c_long)*t_long_double

    print *, "max real double"
    print *, huge(t_double)
    print *, "max real long double"
    print *, huge(t_long_double)

    print *,''

    print *,'Selected Integer Kind 6:'
    print *, huge(j6)

    print *,'Selected Integer Kind 15:'
    print *, huge(j15)

end program integerkinds

