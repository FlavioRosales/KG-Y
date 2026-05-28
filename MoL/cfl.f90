module cfl
  use iso_fortran_env, only: real64
  implicit none
  private
  public :: dt_cfl

contains
  real(real64) function dt_cfl(dr, alpha, beta, guu, cfl) result(dt)
    real(real64), intent(in) :: dr, alpha(:), beta(:), guu(:), cfl
    real(real64) :: smax
    integer :: i, n
    n = size(alpha); smax = 0.0_real64
    do i=1,n
      smax = max( smax, abs(-beta(i) - alpha(i)*sqrt(max(guu(i),0.0_real64))), &
                        abs(-beta(i) + alpha(i)*sqrt(max(guu(i),0.0_real64))) )
    end do
    if (smax <= 0.0_real64) then
      dt = 1.0e99_real64
    else
      dt = cfl * dr / smax
    end if
  end function
end module
