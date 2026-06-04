module initial_data
  !!
  !! Condiciones iniciales usando la base radial de Sturm–Liouville.
  !!
  !! Para cada modo (ell, m):
  !!   phi(r,0) = sum_n C_{ell m n} R_{ell n}(r)
  !!   psi      = d phi / dr
  !!   pi       = (dtphi0 - G%beta psi)/G%alpha, con dtphi0 = psi
  !!
  use state,            only: state_t
  use geometry,            only: geometry_t
  use mesh,            only: mesh_t
  use finite_differences, only: dx => first_derivative_x_2
  use sl_spectrum,      only: Clmn_amp
  implicit none
  private
  public :: set_IC_mode_from_SL
  public :: mode_from_SL_all
  public :: set_IC_incoming_gaussian_positive_charge

contains

!-------------------------------------------------------------------
  !  IC: paquete gaussiano entrante (aprox) con carga de Noether positiva
!-------------------------------------------------------------------
  subroutine set_IC_incoming_gaussian_positive_charge(k0, S, G, Mesh)
    real(kind=8), intent(in)    :: k0
    type(state_t),    intent(inout) :: S
    type(geometry_t), intent(in)    :: G
    type(mesh_t),     intent(in)    :: Mesh
    integer :: Nr, i
    real(kind=8) :: mu0, sigma, A, Mass_field
    real(kind=8) :: x, env, phase
    complex(kind=8) :: ci

    Nr  = size(Mesh%r)
    ci  = dcmplx(0.0d0, 1.0d0)
    mu0 = 80.0d0
    sigma = 10.0d0



    do i = 1, Nr
      x     = (Mesh%r(i) - mu0) / sigma
      env   =  dexp(-0.5d0 * x**2)
      phase = -k0 * Mesh%r(i) 
      
      S%phi(i) = dcmplx(env * dcos(phase), env * dsin(phase))
    end do

    S%psi = dx(S%phi, Mesh)

    S%pi  = dsqrt(G%guu) * S%psi

    Mass_field = S%NC(G,Mesh)

    A = G%M_bh / Mass_field

    S%phi = dsqrt(A) * S%phi
    S%psi = dsqrt(A) * S%psi
    S%pi  = dsqrt(A) * S%pi

    

  end subroutine set_IC_incoming_gaussian_positive_charge

subroutine set_IC_mode_from_SL(lambda, eigenR, S, G, Mesh)
  !!
  !! IC complejas para un modo (ell,m):
  !!   phi(r_i) = sum_n C_n * R(i,n)
  !!   psi = d phi / dr
  !!   pi  = -i * sum_n Omega_n * C_n * R(i,n)
  !
  use sl_spectrum,     only: Clmn_amp
  implicit none

  real(kind=8), intent(in)    :: lambda(:), eigenR(:,:)
  type(state_t),        intent(inout) :: S
  type(geometry_t),     intent(in) :: G
  type(mesh_t),     intent(in) :: Mesh
  integer :: n



  do n = 1, size(lambda)
        S%phi = S%phi + Clmn_amp(lambda(n)) * eigenR(:,n)
  end do


  S%phi = S%phi * exp(- (Mesh%r - 80.0d0)**2 / (2.0d0 * 20.0d0**2) ) 
  S%psi = dx(S%phi, Mesh)
  S%pi = dcmplx( dsqrt(G%guu), 0.0d0 ) * S%psi



end subroutine set_IC_mode_from_SL

subroutine mode_from_SL_all(mu, lambda, eigenR, S, Mesh)
  !!
  !! IC complejas para un modo (ell,m):
  !!   phi(r_i) = sum_n C_n * R(i,n)
  !!   psi = d phi / dr
  !!   pi  = sqrt(guu) * psi   (como en tu elección actual)
  !
  use sl_spectrum,     only: Clmn_amp
  implicit none

  real(kind=8), intent(in)    :: lambda(:), eigenR(:,:)
  real(kind=8), intent(in)    :: mu
  type(state_t),        intent(inout) :: S
  type(mesh_t),     intent(in) :: Mesh

  integer :: n
  complex(kind=8) :: Omega_ln
  complex(kind=8), parameter :: imag_unit = dcmplx(0.0d0, 1.0d0)

  do n = 1, size(lambda)
        Omega_ln = dcmplx(dsqrt(mu**2 + lambda(n)), 0.0d0)
        S%phi = S%phi + Clmn_amp(lambda(n)) * eigenR(:,n)
        S%pi = S%pi - imag_unit * Omega_ln * Clmn_amp(lambda(n)) * eigenR(:,n)
  end do

  S%psi = dx(S%phi, Mesh)

end subroutine mode_from_SL_all


end module initial_data
