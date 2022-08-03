module MBL_subrtns

  use iso_c_binding
  implicit none

contains

  subroutine gap_ratio_R(dim, energies, r_avg, r_sigma)

    integer, intent(in) :: dim
    real (c_double), intent(in) :: energies(dim)

    integer :: n
    real (c_double) :: r_avg, r_sigma, gap1, gap2

    r_avg = 0
    r_sigma = 0
    do n = 1, dim - 2

      gap1 = energies(n+1) - energies(n)
      gap2 = energies(n+2) - energies(n+1)

      r_avg = r_avg + min(gap1,gap2)/max(gap1,gap2)
      r_sigma = r_sigma + (min(gap1,gap2)/max(gap1,gap2))**2

    enddo

    r_avg = r_avg/(dim-2)
    r_sigma = sqrt( (dim-2)/(dim-3) * ( r_sigma/(dim-2) - r_avg**2  ) )

  end subroutine gap_ratio_R

  !real function avg_gap_ratio_C(nspin, quasi_energies)

    !

  !end function avg_gap_ratio_C


  
end module MBL_subrtns


module MBL

  use MBL_subrtns

  interface gap_ratio
    module procedure gap_ratio_R!, level_statistic_C
  end interface

end module MBL



