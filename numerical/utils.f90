module utils
  use mesh, only: mesh_t
  implicit none
  private

  public :: trapezium
  public :: simpson

contains

function trapezium(f, Mesh) result(integral)
  implicit none

  real(kind=8), intent(in), contiguous :: f(:)
  type(mesh_t), intent(in)             :: Mesh
  real(kind=8)                         :: integral
  integer                              :: n, i_flat

  n      = size(f)
  i_flat = Mesh%i_flat

  ! Parte no uniforme: intervalos [1,2], ..., [i_flat-1,i_flat]
  integral = sum( 0.5d0 * (f(1:i_flat-1) + f(2:i_flat)) * &
                  (Mesh%r(2:i_flat) - Mesh%r(1:i_flat-1)) )

  ! Parte uniforme: intervalos [i_flat,i_flat+1], ..., [n-1,n]
  integral = integral + Mesh%dr_max * ( &
       0.5d0 * (f(i_flat) + f(n)) + sum(f(i_flat+1:n-1)) )

end function trapezium

  function simpson(f, dx) result(integral)
  !-----------------------------------------------------------------------
  ! Integración 1D (mayor precisión) para datos reales (vectorizada)
  !
  ! - Si n es impar (m = n-1 subintervalos par): Simpson 1/3 compuesto (O(dx^4))
  ! - Si n es par: Simpson 1/3 en [1..n-1] + trapecio en el último intervalo
  !
  ! integral ≈ (dx/3)[ f1 + fn + 4*sum(f2,f4,...) + 2*sum(f3,f5,...) ]
  !-----------------------------------------------------------------------
  implicit none
  real(kind=8), intent(in) :: dx
  real(kind=8), intent(in) :: f(:)
  real(kind=8)             :: integral
  integer                  :: n, nm1

  n = size(f)
  if (n < 2) then
     integral = 0.0d0
     return
  end if

  ! Caso ideal: n impar -> Simpson 1/3 compuesto en todo el dominio
  if (mod(n,2) == 1) then
     integral = (dx/3.0d0) * ( f(1) + f(n) &
                 + 4.0d0*sum(f(2:n-1:2)) &
                 + 2.0d0*sum(f(3:n-2:2)) )
     return
  end if

  ! Caso n par: aplica Simpson 1/3 hasta n-1 (que es impar), y el último tramo con trapecio
  nm1 = n - 1
  integral = (dx/3.0d0) * ( f(1) + f(nm1) &
              + 4.0d0*sum(f(2:nm1-1:2)) &
              + 2.0d0*sum(f(3:nm1-2:2)) )

  integral = integral + dx * 0.5d0 * ( f(nm1) + f(n) )

  end function simpson


end module
