module finite_differences
  use geometry, only: geometry_t
  use mesh,     only: mesh_t
  implicit none
  private
  integer, parameter :: Re = 1, Im = 2
  public :: first_derivative_x_2, advec_x, advec_full, first_derivative_escaled

  interface first_derivative_x_2
     module procedure first_x_2_u 
     module procedure first_x_2_n
     module procedure first_x_2_hybrid
  end interface

  interface advec_x
     module procedure advec_x_u
     module procedure advec_x_n
     module procedure advec_x_hybrid
  end interface

  interface advec_full
     module procedure advec_full_hybrid
  end interface

  interface first_derivative_escaled
      module procedure first_derivative_escaled_hybrid
      module procedure derivative_escaled_hybrid_product
  end interface

contains

!===============================================================================
! Derivada primera uniforme — complejo 1D
!===============================================================================
function first_x_2_u(f, h) result(df)
  implicit none

  real(kind=8), intent(in),  contiguous :: f(:)
  real(kind=8),    intent(in)              :: h
  real(kind=8) :: df(size(f))

  integer      :: n
  real(kind=8) :: res

  n   = size(f)
  res = 0.5d0 / h

  df(1) = (-3.0d0*f(1) + 4.0d0*f(2) - f(3)) * res

  df(2:n-1) = (f(3:n) - f(1:n-2)) * res

  df(n) = (3.0d0*f(n) - 4.0d0*f(n-1) + f(n-2)) * res

end function first_x_2_u


!===============================================================================
! Derivada primera no uniforme — complejo 1D
!===============================================================================
function first_x_2_n(f, x) result(df)
  implicit none

  real(kind=8), intent(in),  contiguous :: f(:)
  real(kind=8),    intent(in),  contiguous :: x(:)
  real(kind=8) :: df(size(f))

  integer      :: i, n
  real(kind=8) :: h0, h1
  real(kind=8) :: a, b, c

  n = size(f)

  h0 = x(2) - x(1)
  h1 = x(3) - x(2)

  a = -(2.0d0*h0 + h1)/(h0*(h0 + h1))
  b =  (h0 + h1)/(h0*h1)
  c = -h0/(h1*(h0 + h1))

  df(1) = a*f(1) + b*f(2) + c*f(3)

  do i = 2, n-1
     h0 = x(i)   - x(i-1)
     h1 = x(i+1) - x(i)

     a = -h1/(h0*(h0 + h1))
     b =  (h1 - h0)/(h0*h1)
     c =  h0/(h1*(h0 + h1))

     df(i) = a*f(i-1) + b*f(i) + c*f(i+1)
  end do

  h0 = x(n)   - x(n-1)
  h1 = x(n-1) - x(n-2)

  a =  (2.0d0*h0 + h1)/(h0*(h0 + h1))
  b = -(h0 + h1)/(h0*h1)
  c =  h0/(h1*(h0 + h1))

  df(n) = a*f(n) + b*f(n-1) + c*f(n-2)

end function first_x_2_n


!===============================================================================
! Derivada primera híbrida
!
! Requiere en mesh_t:
!   M%d1_a(:), M%d1_b(:), M%d1_c(:)
!   M%inv_2dr_max
!===============================================================================
function first_x_2_hybrid(f, M) result(df)
  implicit none

  real(kind=8), intent(in),  contiguous :: f(:,:)
  type(mesh_t), intent(in)              :: M
  real(kind=8) :: df(size(f,1),size(f,2))

  integer      :: i, n, i_flat
  real(kind=8) :: res

  n      = size(f,1)
  i_flat = M%i_flat
  res    = M%inv_2dr_max

  df(1,Re) = M%d1_a(1)*f(1,Re) + M%d1_b(1)*f(2,Re) + M%d1_c(1)*f(3,Re)
  df(1,Im) = M%d1_a(1)*f(1,Im) + M%d1_b(1)*f(2,Im) + M%d1_c(1)*f(3,Im)

  do concurrent (i = 2:i_flat-1)
     df(i,Re) = M%d1_a(i)*f(i-1,Re) + M%d1_b(i)*f(i,Re) + M%d1_c(i)*f(i+1,Re)
     df(i,Im) = M%d1_a(i)*f(i-1,Im) + M%d1_b(i)*f(i,Im) + M%d1_c(i)*f(i+1,Im)
  end do

  do concurrent (i = i_flat:n-1)
     df(i,Re) = (f(i+1,Re) - f(i-1,Re)) * res
     df(i,Im) = (f(i+1,Im) - f(i-1,Im)) * res
  end do

  df(n,Re) = (3.0d0*f(n,Re) - 4.0d0*f(n-1,Re) + f(n-2,Re)) * res
  df(n,Im) = (3.0d0*f(n,Im) - 4.0d0*f(n-1,Im) + f(n-2,Im)) * res

end function first_x_2_hybrid


!===============================================================================
! Advección uniforme — complejo 1D
!===============================================================================
function advec_x_u(f, beta, dx) result(df)
  implicit none

  real(kind=8), intent(in),  contiguous :: f(:)
  real(kind=8), intent(in),  contiguous :: beta(:)
  real(kind=8), intent(in)              :: dx
  real(kind=8) :: df(size(f))

  integer      :: n
  real(kind=8) :: res

  n   = size(f)
  res = 0.5d0 / dx

  df(1:n-2) = beta(1:n-2) * ( &
       -3.0d0*f(1:n-2)       &
     +  4.0d0*f(2:n-1)       &
     -        f(3:n) ) * res

  df(n-1) = beta(n-1) * (f(n) - f(n-2)) * res

  df(n) = beta(n) * ( &
       3.0d0*f(n) - 4.0d0*f(n-1) + f(n-2) ) * res

end function advec_x_u


!===============================================================================
! Advección no uniforme — complejo 1D
!===============================================================================
function advec_x_n(f, beta, x) result(df)
  implicit none

  real(kind=8), intent(in),  contiguous :: f(:)
  real(kind=8),    intent(in),  contiguous :: beta(:)
  real(kind=8),    intent(in),  contiguous :: x(:)
  real(kind=8) :: df(size(f))

  integer      :: i, n
  real(kind=8) :: h0, h1, denom

  n = size(f)

  do i = 1, n-2
     h0 = x(i+1) - x(i)
     h1 = x(i+2) - x(i+1)

     denom = h0*h1*(h0 + h1)

     df(i) = beta(i) * ( &
          -(2.0d0*h0 + h1)*h1*f(i) &
        + (h0 + h1)*(h0 + h1)*f(i+1) &
        - h0*h0*f(i+2) ) / denom
  end do

  h0 = x(n-1) - x(n-2)
  h1 = x(n)   - x(n-1)

  denom = h0*h1*(h0 + h1)

  df(n-1) = beta(n-1) * ( &
       - h1*h1*f(n-2)              &
       + (h1*h1 - h0*h0)*f(n-1)    &
       + h0*h0*f(n) ) / denom

  h0 = x(n)   - x(n-1)
  h1 = x(n-1) - x(n-2)

  denom = h0*h1*(h0 + h1)

  df(n) = beta(n) * ( &
       (2.0d0*h0 + h1)*h1*f(n)       &
     - (h0 + h1)*(h0 + h1)*f(n-1)    &
     + h0*h0*f(n-2) ) / denom

end function advec_x_n


!===============================================================================
! Advección híbrida 
!
! Requiere en geometry_t:
!   G%adv_a(:),  G%adv_b(:),  G%adv_c(:)
!   G%adv_u0(:), G%adv_u1(:), G%adv_u2(:)
!
!===============================================================================
function advec_x_hybrid(f, G, M) result(df)
  implicit none

  real(kind=8), intent(in),  contiguous :: f(:,:)
  type(geometry_t), intent(in)          :: G
  type(mesh_t),     intent(in)          :: M
  real(kind=8)  :: df(size(f,1),size(f,2))

  integer :: i, n, i_flat

  n      = size(f,1)
  i_flat = M%i_flat

  do concurrent (i = 1:i_flat-1)
     df(i,Re) = G%adv_a(i)*f(i,Re)   &
             + G%adv_b(i)*f(i+1,Re) &
             + G%adv_c(i)*f(i+2,Re)
     df(i,Im) = G%adv_a(i)*f(i,Im)   &
             + G%adv_b(i)*f(i+1,Im) &
             + G%adv_c(i)*f(i+2,Im)

  end do

  do concurrent (i = i_flat:n-2)
     df(i,Re) = G%adv_u0(i)*f(i,Re)   &
             + G%adv_u1(i)*f(i+1,Re) &
             + G%adv_u2(i)*f(i+2,Re)
     df(i,Im) = G%adv_u0(i)*f(i,Im)   &
             + G%adv_u1(i)*f(i+1,Im) &
             + G%adv_u2(i)*f(i+2,Im)
  end do

  df(n-1,Re) = G%adv_u0(n-1)*f(n-2,Re) + G%adv_u2(n-1)*f(n,Re)
  df(n-1,Im) = G%adv_u0(n-1)*f(n-2,Im) + G%adv_u2(n-1)*f(n,Im)

  df(n,Re) = G%adv_u0(n)*f(n-2,Re) &
          + G%adv_u1(n)*f(n-1,Re) &
          + G%adv_u2(n)*f(n,Re)
  df(n,Im) = G%adv_u0(n)*f(n-2,Im) &
          + G%adv_u1(n)*f(n-1,Im) &
          + G%adv_u2(n)*f(n,Im)


end function advec_x_hybrid

!===============================================================================!
! Advección híbrida completa: incluye el término de advección con la derivada de beta
! Requiere en geometry_t:
!   G%adv_a(:),  G%adv_b(:),  G%adv_c(:)
!   G%adv_u0(:), G%adv_u1(:), G%adv_u2(:)
!   G%dbeta(:)
!===============================================================================!

function advec_full_hybrid(f, G, M) result(df)
  implicit none

  real(kind=8), intent(in),  contiguous :: f(:,:)
  type(geometry_t), intent(in)          :: G
  type(mesh_t),     intent(in)          :: M
  real(kind=8)  :: df(size(f,1),size(f,2))

  integer :: i, n, i_flat

  n      = size(f,1)
  i_flat = M%i_flat

  do concurrent (i = 1:i_flat-1)
     df(i,Re) = G%adv_a(i)*f(i,Re)   &
             + G%adv_b(i)*f(i+1,Re) &
             + G%adv_c(i)*f(i+2,Re) & 
             + G%dbeta(i) * f(i,Re)
     df(i,Im) = G%adv_a(i)*f(i,Im)   &
             + G%adv_b(i)*f(i+1,Im) &
             + G%adv_c(i)*f(i+2,Im) & 
             + G%dbeta(i) * f(i,Im)

  end do

  do concurrent (i = i_flat:n-2)
     df(i,Re) = G%adv_u0(i)*f(i,Re)   &
             + G%adv_u1(i)*f(i+1,Re) &
             + G%adv_u2(i)*f(i+2,Re) & 
             + G%dbeta(i) * f(i,Re)
     df(i,Im) = G%adv_u0(i)*f(i,Im)   &
             + G%adv_u1(i)*f(i+1,Im) &
             + G%adv_u2(i)*f(i+2,Im) & 
             + G%dbeta(i) * f(i,Im)
  end do

  df(n-1,Re) = G%adv_u0(n-1)*f(n-2,Re) + G%adv_u2(n-1)*f(n,Re) + G%dbeta(n-1) * f(n-1,Re)
  df(n-1,Im) = G%adv_u0(n-1)*f(n-2,Im) + G%adv_u2(n-1)*f(n,Im) + G%dbeta(n-1) * f(n-1,Im)

  df(n,Re) = G%adv_u0(n)*f(n-2,Re) &
          + G%adv_u1(n)*f(n-1,Re) &
          + G%adv_u2(n)*f(n,Re) & 
          + G%dbeta(n) * f(n,Re)
  df(n,Im) = G%adv_u0(n)*f(n-2,Im) &
          + G%adv_u1(n)*f(n-1,Im) &
          + G%adv_u2(n)*f(n,Im) & 
          + G%dbeta(n) * f(n,Im)

end function advec_full_hybrid

function first_derivative_escaled_hybrid(a,f,M) result(df)
  implicit none

  real(kind=8), intent(in),  contiguous :: f(:,:)
  real(kind=8), intent(in),  contiguous :: a(:)
  type(mesh_t), intent(in)              :: M
  real(kind=8) :: df(size(f,1),size(f,2))

  integer      :: i, n, i_flat
  real(kind=8) :: res

  n      = size(f,1)
  i_flat = M%i_flat
  res    = M%inv_2dr_max


  df(1,Re) = M%d1_a(1)*f(1,Re)*a(1) + M%d1_b(1)*f(2,Re)*a(2) + M%d1_c(1)*f(3,Re)*a(3)
  df(1,Im) = M%d1_a(1)*f(1,Im)*a(1) + M%d1_b(1)*f(2,Im)*a(2) + M%d1_c(1)*f(3,Im)*a(3)

  do concurrent (i = 2:i_flat-1)
     df(i,Re) = M%d1_a(i)*f(i-1,Re)*a(i-1) + M%d1_b(i)*f(i,Re)*a(i) + M%d1_c(i)*f(i+1,Re)*a(i+1)
     df(i,Im) = M%d1_a(i)*f(i-1,Im)*a(i-1) + M%d1_b(i)*f(i,Im)*a(i) + M%d1_c(i)*f(i+1,Im)*a(i+1)
  end do

  do concurrent (i = i_flat:n-1)
     df(i,Re) = (f(i+1,Re)*a(i+1) - f(i-1,Re)*a(i-1)) * res 
     df(i,Im) = (f(i+1,Im)*a(i+1) - f(i-1,Im)*a(i-1)) * res 
  end do

  df(n,Re) = (3.0d0*f(n,Re)*a(n) - 4.0d0*f(n-1,Re)*a(n-1) + f(n-2,Re)*a(n-2)) * res
  df(n,Im) = (3.0d0*f(n,Im)*a(n) - 4.0d0*f(n-1,Im)*a(n-1) + f(n-2,Im)*a(n-2)) * res

end function first_derivative_escaled_hybrid

function derivative_escaled_hybrid_product(b,a,f,M) result(df)
  implicit none

  real(kind=8), intent(in),  contiguous :: f(:,:)
  real(kind=8), intent(in),  contiguous :: a(:), b(:)
  type(mesh_t), intent(in)              :: M
  real(kind=8) :: df(size(f,1),size(f,2))

  integer      :: i, n, i_flat
  real(kind=8) :: res

  n      = size(f,1)
  i_flat = M%i_flat
  res    = M%inv_2dr_max


  df(1,Re) = (M%d1_a(1)*f(1,Re)*a(1) + M%d1_b(1)*f(2,Re)*a(2) + M%d1_c(1)*f(3,Re)*a(3)) * b(1)
  df(1,Im) = (M%d1_a(1)*f(1,Im)*a(1) + M%d1_b(1)*f(2,Im)*a(2) + M%d1_c(1)*f(3,Im)*a(3)) * b(1)

  do concurrent (i = 2:i_flat-1)
     df(i,Re) = (M%d1_a(i)*f(i-1,Re)*a(i-1) + M%d1_b(i)*f(i,Re)*a(i) + M%d1_c(i)*f(i+1,Re)*a(i+1)) * b(i)
     df(i,Im) = (M%d1_a(i)*f(i-1,Im)*a(i-1) + M%d1_b(i)*f(i,Im)*a(i) + M%d1_c(i)*f(i+1,Im)*a(i+1)) * b(i)
  end do

  do concurrent (i = i_flat:n-1)
     df(i,Re) = (f(i+1,Re)*a(i+1) - f(i-1,Re)*a(i-1)) * res * b(i)
     df(i,Im) = (f(i+1,Im)*a(i+1) - f(i-1,Im)*a(i-1)) * res * b(i)
  end do

  df(n,Re) = (3.0d0*f(n,Re)*a(n) - 4.0d0*f(n-1,Re)*a(n-1) + f(n-2,Re)*a(n-2)) * res * b(n)
  df(n,Im) = (3.0d0*f(n,Im)*a(n) - 4.0d0*f(n-1,Im)*a(n-1) + f(n-2,Im)*a(n-2)) * res * b(n)

end function derivative_escaled_hybrid_product

end module finite_differences