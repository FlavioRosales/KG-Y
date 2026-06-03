module finite_differences
  use geometry,        only: geometry_t
  use mesh,            only: mesh_t
  implicit none
  private

  !===============================
  !  API pública
  !===============================
  public :: first_derivative_x_2
  public :: advec_x
  public :: first_x_2_u_sub, first_x_2_n_sub, first_x_2_hybrid_sub
  public :: advec_x_u_sub, advec_x_n_sub, advec_x_hybrid_sub
  public :: first_x_2_hybrid, advec_x_hybrid

  !--- interfaces de usuario ---
  interface first_derivative_x_2
    module procedure first_x_2_u
    module procedure first_x_2_n
    module procedure first_x_2_hybrid_fn
  end interface


  interface advec_x
    module procedure advec_x_u
    module procedure advec_x_n
    module procedure advec_x_hybrid_fn
  end interface


contains
!===============================================================================
! Núcleos 2º ORDEN — derivada primera
!===============================================================================
subroutine first_x_2_u_sub(f, h, df)
  implicit none

  complex(kind=8), intent(in),  contiguous :: f(:)
  real(kind=8),    intent(in)              :: h
  complex(kind=8), intent(out), contiguous :: df(:)

  integer       :: n, i
  real(kind=8)  :: res

  n   = size(f)
  res = 0.5d0 / h


  df(1) = (-3.0d0*f(1) + 4.0d0*f(2) - f(3)) * res
  df(n) = ( 3.0d0*f(n) - 4.0d0*f(n-1) + f(n-2)) * res

  !do i = 2, n-1
  !  df(i) = (f(i+1) - f(i-1)) * res
  !end do

  df(2:n-1) = (f(3:n) - f(1:n-2)) * res

end subroutine first_x_2_u_sub

subroutine first_x_2_n_sub(f, x, df)
  implicit none

  complex(kind=8), intent(in),  contiguous :: f(:)
  real(kind=8),    intent(in),  contiguous :: x(:)
  complex(kind=8), intent(out), contiguous :: df(:)

  integer :: i, n
  real(kind=8) :: h0, h1
  real(kind=8) :: a, b, c

  n = size(f)

  if (size(x) /= n) error stop "first_x_2_n_sub: size(x) /= size(f)"
  if (size(df) /= n) error stop "first_x_2_n_sub: size(df) /= size(f)"
  if (n < 3) error stop "first_x_2_n_sub: n must be at least 3"

  ! Forward en i = 1
  h0 = x(2) - x(1)
  h1 = x(3) - x(2)

  a = -(2.0d0*h0 + h1)/(h0*(h0 + h1))
  b =  (h0 + h1)/(h0*h1)
  c = -h0/(h1*(h0 + h1))

  df(1) = a*f(1) + b*f(2) + c*f(3)

  ! Centrada en puntos interiores
  do i = 2, n-1
    h0 = x(i)   - x(i-1)
    h1 = x(i+1) - x(i)

    a = -h1/(h0*(h0 + h1))
    b =  (h1 - h0)/(h0*h1)
    c =  h0/(h1*(h0 + h1))

    df(i) = a*f(i-1) + b*f(i) + c*f(i+1)
  end do

  ! Backward en i = n
  h0 = x(n)   - x(n-1)
  h1 = x(n-1) - x(n-2)

  a =  (2.0d0*h0 + h1)/(h0*(h0 + h1))
  b = -(h0 + h1)/(h0*h1)
  c =  h0/(h1*(h0 + h1))

  df(n) = a*f(n) + b*f(n-1) + c*f(n-2)

end subroutine first_x_2_n_sub

subroutine first_x_2_hybrid_sub(f, M, df)
  implicit none

  complex(kind=8), intent(in),  contiguous :: f(:)
  type(mesh_t),    intent(in)              :: M
  complex(kind=8), intent(out), contiguous :: df(:)

  integer :: i, n
  real(kind=8) :: h0, h1, res
  real(kind=8) :: a, b, c

  n = size(f)


  res = 0.5d0 / M%dr_max

  ! Región no uniforme: usamos los nodos reales de M%r hasta antes de i_flat.
  h0 = M%r(2) - M%r(1)
  h1 = M%r(3) - M%r(2)

  a = -(2.0d0*h0 + h1)/(h0*(h0 + h1))
  b =  (h0 + h1)/(h0*h1)
  c = -h0/(h1*(h0 + h1))

  df(1) = a*f(1) + b*f(2) + c*f(3)

  do i = 2, M%i_flat - 1
    h0 = M%r(i)   - M%r(i-1)
    h1 = M%r(i+1) - M%r(i)

    a = -h1/(h0*(h0 + h1))
    b =  (h1 - h0)/(h0*h1)
    c =  h0/(h1*(h0 + h1))

    df(i) = a*f(i-1) + b*f(i) + c*f(i+1)
  end do

  df(M%i_flat:n-1) = (f(M%i_flat+1:n) - f(M%i_flat-1:n-2)) * res

  df(n) = (3.0d0*f(n) - 4.0d0*f(n-1) + f(n-2)) * res

end subroutine first_x_2_hybrid_sub

subroutine first_x_2_hybrid(f, M, df)
  implicit none

  complex(kind=8), intent(in),  contiguous :: f(:)
  type(mesh_t),    intent(in)              :: M
  complex(kind=8), intent(out), contiguous :: df(:)

  call first_x_2_hybrid_sub(f, M, df)
end subroutine first_x_2_hybrid


!===============================================================================
! Advección df/dx — 2º ORDEN
!===============================================================================

subroutine advec_x_u_sub(f, beta, dx, df)
  implicit none

  complex(kind=8), intent(in),  contiguous :: f(:)
  real(kind=8), intent(in),  contiguous :: beta(:)
  real(kind=8),    intent(in)              :: dx
  complex(kind=8), intent(out), contiguous :: df(:)

  real(kind=8)  :: res
  integer       :: n

  n   = size(f)
  res = 0.5d0 / dx

  df(n-1) = (( f(n) - f(n-1))/dx )*beta(n-1)
  df(1:n-2) = ((-3.0d0*f(1:n-2) + 4.0d0*f(2:n-1) - f(3:n)) * res ) * beta(1:n-2)
  df(n) = ((3.0d0*f(n) - 4.0d0*f(n-1) + f(n-2)) * res) * beta(n)
end subroutine advec_x_u_sub

subroutine advec_x_n_sub(f, beta, x, df)
  implicit none

  complex(kind=8), intent(in),  contiguous :: f(:)
  real(kind=8),    intent(in),  contiguous :: x(:)
  real(kind=8),    intent(in),  contiguous :: beta(:)
  complex(kind=8), intent(out), contiguous :: df(:)

  integer :: n, i
  real(kind=8) :: h0, h1, denom

  n = size(f)

  ! Forward de segundo orden para i = 1, ..., n-2
  do i = 1, n-2

    h0 = x(i+1) - x(i)
    h1 = x(i+2) - x(i+1)

    denom = h0*h1*(h0 + h1)

    df(i) = beta(i) * ( &
         -(2.0d0*h0 + h1)*h1*f(i) &
         + (h0 + h1)**2*f(i+1)          &
         - h0**2*f(i+2)                 &
         ) / denom

  end do

  ! i = n-1
  !
  ! Aquí ya no existe f(i+2). Con los puntos disponibles,
  ! una opción de segundo orden es usar stencil centrado:
  !
  ! h0 = x(n-1) - x(n-2)
  ! h1 = x(n)   - x(n-1)

  h0 = x(n-1) - x(n-2)
  h1 = x(n)   - x(n-1)

  denom = h0*h1*(h0 + h1)

  df(n-1) = beta(n-1) * ( &
       - h1**2*f(n-2)              &
       + (h1**2 - h0**2)*f(n-1)    &
       + h0**2*f(n)                &
       ) / denom

  ! i = n
  !
  ! Para evitar dejar df(n) sin inicializar, usamos backward
  ! de segundo orden.

  h0 = x(n)   - x(n-1)
  h1 = x(n-1) - x(n-2)

  denom = h0*h1*(h0 + h1)

  df(n) = beta(n) * ( &
       (2.0d0*h0 + h1)*h1*f(n) &
       - (h0 + h1)**2*f(n-1)        &
       + h0**2*f(n-2)               &
       ) / denom

end subroutine advec_x_n_sub

subroutine advec_x_hybrid_sub(f, G, M, df)
  implicit none

  complex(kind=8), intent(in),  contiguous :: f(:)
  type(geometry_t), intent(in)             :: G
  type(mesh_t),     intent(in)             :: M
  complex(kind=8), intent(out), contiguous :: df(:)

  integer :: n, i
  real(kind=8) :: h0, h1, denom, res

  n = size(f)
  
  res = 0.5d0 / M%dr_max

  do i = 1, M%i_flat - 1
    h0 = M%r(i+1) - M%r(i)
    h1 = M%r(i+2) - M%r(i+1)

    denom = h0*h1*(h0 + h1)

    df(i) = G%beta(i) * ( &
         -(2.0d0*h0 + h1)*h1*f(i) &
         + (h0 + h1)**2*f(i+1)          &
         - h0**2*f(i+2)                 &
         ) / denom
  end do

  df(M%i_flat:n-2) = G%beta(M%i_flat:n-2) * ( &
        -3.0d0*f(M%i_flat:n-2) &
        + 4.0d0*f(M%i_flat+1:n-1) &
        - f(M%i_flat+2:n) ) * res

  df(n-1) = G%beta(n-1) * (f(n) - f(n-2)) * res
  df(n) = G%beta(n) * ( &
       3.0d0*f(n) - 4.0d0*f(n-1) + f(n-2) &
       ) * res

end subroutine advec_x_hybrid_sub

subroutine advec_x_hybrid(f, G, M, df)
  implicit none

  complex(kind=8), intent(in),  contiguous :: f(:)
  type(geometry_t), intent(in)             :: G
  type(mesh_t),     intent(in)             :: M
  complex(kind=8), intent(out), contiguous :: df(:)

  call advec_x_hybrid_sub(f, G, M, df)
end subroutine advec_x_hybrid


!===============================================================================
! WRAPPERS PÚBLICOS
!===============================================================================

function first_x_2_u(f, h) result(df)
    complex(kind=8), intent(in) :: f(:)
    real(kind=8),    intent(in) :: h
    complex(kind=8) :: df(size(f))
    call first_x_2_u_sub(f, h, df)
  end function first_x_2_u

function first_x_2_n(f, x) result(df)
    complex(kind=8), intent(in) :: f(:)
    real(kind=8),    intent(in) :: x(:)
    complex(kind=8) :: df(size(f))
    call first_x_2_n_sub(f, x, df)
  end function first_x_2_n

function first_x_2_hybrid_fn(f, M) result(df)
    complex(kind=8), intent(in) :: f(:)
    type(mesh_t),    intent(in) :: M
    complex(kind=8) :: df(size(f))
    call first_x_2_hybrid_sub(f, M, df)
  end function first_x_2_hybrid_fn

function advec_x_u(f, beta, dx) result(df)
    complex(kind=8), intent(in) :: f(:)
    real(kind=8), intent(in) :: beta(:)
    real(kind=8),    intent(in) :: dx
    complex(kind=8) :: df(size(f))
    call advec_x_u_sub(f, beta, dx, df)
  end function advec_x_u

function advec_x_n(f, beta, x) result(df)
    complex(kind=8), intent(in) :: f(:)
    real(kind=8), intent(in) :: beta(:)
    real(kind=8),    intent(in) :: x(:)
    complex(kind=8) :: df(size(f))
    call advec_x_n_sub(f, beta, x, df)
end function advec_x_n

function advec_x_hybrid_fn(f, G, M) result(df)
    complex(kind=8), intent(in) :: f(:)
    type(geometry_t), intent(in) :: G
    type(mesh_t),     intent(in) :: M
    complex(kind=8) :: df(size(f))
    call advec_x_hybrid_sub(f, G, M, df)
end function advec_x_hybrid_fn


end module finite_differences
