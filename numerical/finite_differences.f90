module finite_differences
  use iso_fortran_env, only: real64
  implicit none
  private

  !===============================
  !  API pública
  !===============================
  public :: first_derivative_x_2
  public :: advec_x
  public :: first_x_2_u_sub, first_x_2_n_sub
  public :: advec_x_u_sub, advec_x_n_sub

  !--- interfaces de usuario ---
  interface first_derivative_x_2
    module procedure first_x_2_u
    module procedure first_x_2_n
  end interface


  interface advec_x
    module procedure advec_x_u
    module procedure advec_x_n
  end interface


contains
!===============================================================================
! Núcleos 2º ORDEN — derivada primera
!===============================================================================
subroutine first_x_2_u_sub(f, h, df)
  use iso_fortran_env, only: real64
  implicit none

  complex(real64), intent(in),  contiguous :: f(:)
  real(real64),    intent(in)              :: h
  complex(real64), intent(out), contiguous :: df(:)

  integer       :: n, i
  real(real64)  :: res

  n   = size(f)
  res = 0.5_real64 / h


  df(1) = (-3.0_real64*f(1) + 4.0_real64*f(2) - f(3)) * res
  df(n) = ( 3.0_real64*f(n) - 4.0_real64*f(n-1) + f(n-2)) * res

  !do i = 2, n-1
  !  df(i) = (f(i+1) - f(i-1)) * res
  !end do

  df(2:n-1) = (f(3:n) - f(1:n-2)) * res

end subroutine first_x_2_u_sub


subroutine first_x_2_n_sub(f, x, df)
    complex(real64), intent(in),  contiguous :: f(:)
    real(real64),    intent(in),  contiguous :: x(:)
    complex(real64), intent(out), contiguous :: df(:)
    integer :: i, n
    real(real64) :: h0, h1, denom

    n = size(f)

    ! i=1 (forward 2º orden)
    h0 = x(2) - x(1)
    h1 = x(3) - x(2)
    denom = h0 * (h0 + h1)
    df(1) = (-(2.0_real64*h0 + h1)*f(1) + (h0 + h1)*f(2) - h0*f(3)) / denom

    ! interior
    do i = 2, n-1
      h0 = x(i)   - x(i-1)
      h1 = x(i+1) - x(i)
      df(i) = (-h1*f(i-1) + (h1 - h0)*f(i) + h0*f(i+1)) / &
              (0.5_real64*(h0 + h1)*h0*h1)
    end do

    ! i=n (backward 2º orden)
    h0 = x(n)   - x(n-1)
    h1 = x(n-1) - x(n-2)
    denom = h0 * (h0 + h1)
    df(n) = ((2.0_real64*h0 + h1)*f(n) - (h0 + h1)*f(n-1) + h0*f(n-2)) / denom

  end subroutine first_x_2_n_sub

!===============================================================================
! Advección df/dx — 2º ORDEN
!===============================================================================

subroutine advec_x_u_sub(f, beta, dx, df)
  use iso_fortran_env, only: real64
  implicit none

  complex(real64), intent(in),  contiguous :: f(:)
  real(real64), intent(in),  contiguous :: beta(:)
  real(real64),    intent(in)              :: dx
  complex(real64), intent(out), contiguous :: df(:)

  real(real64)  :: res
  integer       :: n

  n   = size(f)
  res = 0.5_real64 / dx

  df(n-1) = (( f(n) - f(n-1))/dx )*beta(n-1)
  df(1:n-2) = ((-3.0_real64*f(1:n-2) + 4.0_real64*f(2:n-1) - f(3:n)) * res ) * beta(1:n-2)
  df(n) = ((3.0_real64*f(n) - 4.0_real64*f(n-1) + f(n-2)) * res) * beta(n)
end subroutine advec_x_u_sub

subroutine advec_x_n_sub(f, beta, x, df)
    complex(real64), intent(in),  contiguous :: f(:)
    real(real64),    intent(in),  contiguous :: x(:)
    real(real64),    intent(in),  contiguous :: beta(:)
    complex(real64), intent(out), contiguous :: df(:)
    integer :: n, i
    real(real64) :: h0, h1, invden

    n = size(f)

    do i = 1, n-2
        h0 = x(i+1)-x(i)
        h1 = x(i+2)-x(i+1)
        invden = 1/(h0*(h0+h1))
        df(i) = ((-(2*h0+h1)*f(i) + (h0+h1)*f(i+1) - h0*f(i+2)) * invden ) * beta(i)
    end do

    ! i = n-1
    h1 = x(n)-x(n-1)
    df(n-1) = ((f(n)-f(n-1))/(h1) )*beta(n-1)

end subroutine advec_x_n_sub


!===============================================================================
! WRAPPERS PÚBLICOS
!===============================================================================

function first_x_2_u(f, h) result(df)
    complex(real64), intent(in) :: f(:)
    real(real64),    intent(in) :: h
    complex(real64) :: df(size(f))
    call first_x_2_u_sub(f, h, df)
  end function first_x_2_u

function first_x_2_n(f, x) result(df)
    complex(real64), intent(in) :: f(:)
    real(real64),    intent(in) :: x(:)
    complex(real64) :: df(size(f))
    call first_x_2_n_sub(f, x, df)
  end function first_x_2_n

function advec_x_u(f, beta, dx) result(df)
    complex(real64), intent(in) :: f(:)
    real(real64), intent(in) :: beta(:)
    real(real64),    intent(in) :: dx
    complex(real64) :: df(size(f))
    call advec_x_u_sub(f, beta, dx, df)
  end function advec_x_u

function advec_x_n(f, beta, x) result(df)
    complex(real64), intent(in) :: f(:)
    real(real64), intent(in) :: beta(:)
    real(real64),    intent(in) :: x(:)
    complex(real64) :: df(size(f))
    call advec_x_n_sub(f, beta, x, df)
end function advec_x_n

end module finite_differences
