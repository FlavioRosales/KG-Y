module initial_data
  !!
  !! Condiciones iniciales usando la base radial de Sturm–Liouville.
  !!
  !! Para cada modo (ell, m):
  !!   phi(r,0) = sum_n C_{ell m n} R_{ell n}(r)
  !!   psi      = d phi / dr
  !!   pi       = (dtphi0 - G%beta psi)/G%alpha, con dtphi0 = psi
  !!
  use iso_fortran_env,  only: real64
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
  subroutine set_IC_incoming_gaussian_positive_charge(ell, m, k0, mu_field, S, G, Mesh)
    use iso_fortran_env, only: real64
    integer,      intent(in)    :: ell, m
    real(real64), intent(in)    :: k0
    real(real64), intent(in)    :: mu_field
    type(state_t),    intent(inout) :: S
    type(geometry_t), intent(in)    :: G
    type(mesh_t),     intent(in)    :: Mesh
    integer :: Nr, i
    real(real64) :: mu0, sigma, A, Mass_field
    real(real64) :: x, env, phase
    complex(real64) :: ci

    Nr  = size(Mesh%r)
    ci  = cmplx(0.0_real64, 1.0_real64, kind=real64)
    mu0 = 80.0_real64
    sigma = 10.0_real64



    do i = 1, Nr
      x     = (Mesh%r(i) - mu0) / sigma
      env   =  exp(-0.5_real64 * x**2)
      phase = -k0 * Mesh%r(i) 
      
      S%phi(i) = cmplx(env * cos(phase), env * sin(phase), kind=real64)
    end do

    S%psi = dx(S%phi, Mesh)

    S%pi  = dsqrt(G%guu) * S%psi

    !Mass_field = mu_field*S%NC(G,Mesh)
    Mass_field = S%NC(G,Mesh)

    A = G%M_bh / Mass_field

    S%phi = sqrt(A) * S%phi
    S%psi = sqrt(A) * S%psi
    S%pi  = sqrt(A) * S%pi

    

  end subroutine set_IC_incoming_gaussian_positive_charge

subroutine set_IC_mode_from_SL(ell, m, lambda, eigenR, S, G, Mesh)
  !!
  !! IC complejas para un modo (ell,m):
  !!   phi(r_i) = sum_n C_n * R(i,n)
  !!   psi = d phi / dr
  !!   pi  = -i * sum_n Omega_n * C_n * R(i,n)
  !!
  use iso_fortran_env, only: real64
  use sl_spectrum,     only: Clmn_amp
  implicit none

  integer,      intent(in)    :: ell, m
  real(real64), intent(in)    :: lambda(:), eigenR(:,:)
  real(real64)  :: Mass_field, A
  type(state_t),        intent(inout) :: S
  type(geometry_t),     intent(in) :: G
  type(mesh_t),     intent(in) :: Mesh

  integer :: n, i
  complex(real64), allocatable :: C(:)
  real(real64), parameter :: alpha_eps = 1.0e-14_real64

  allocate(C(size(lambda)))

  ! --- coeficientes C_n (complejos gaussianos con varianza P(k))
  do n = 1, size(lambda)
    C(n) = Clmn_amp(ell, m, n, lambda(n))
  end do

  do n = 1, size(lambda)
        S%phi = S%phi + C(n) * eigenR(:,n)
  end do

  deallocate(C)

  S%phi = S%phi * exp(- (Mesh%r - 80.0_real64)**2 / (2.0_real64 * 20.0_real64**2) ) 
  S%psi = dx(S%phi, Mesh)
  S%pi = cmplx( sqrt(G%guu), 0.0_real64, kind=real64 ) * S%psi


  ! Mass_field = S%NC(G,Mesh)
  ! A = G%M_bh / Mass_field

  ! S%phi = sqrt(A) * S%phi
  ! S%psi = sqrt(A) * S%psi
  ! S%pi  = sqrt(A) * S%pi


end subroutine set_IC_mode_from_SL

subroutine mode_from_SL_all(ell, m, mu, lambda, eigenR, S, G, Mesh)
  !!
  !! IC complejas para un modo (ell,m):
  !!   phi(r_i) = sum_n C_n * R(i,n)
  !!   psi = d phi / dr
  !!   pi  = sqrt(guu) * psi   (como en tu elección actual)
  !!
  use iso_fortran_env, only: real64
  use sl_spectrum,     only: Clmn_amp
  implicit none

  integer,      intent(in)    :: ell, m
  real(real64), intent(in)    :: lambda(:), eigenR(:,:)
  real(real64), intent(in)    :: mu
  type(state_t),        intent(inout) :: S
  type(geometry_t),     intent(in) :: G
  type(mesh_t),     intent(in) :: Mesh

  integer :: n
  complex(real64), allocatable :: C(:)
  complex(real64) :: Omega_ln
  complex(real64), parameter :: imag_unit = cmplx(0.0_real64, 1.0_real64, kind=real64)

  allocate(C(size(lambda)))

  ! --- coeficientes C_n (complejos gaussianos con varianza P(k))
  do n = 1, size(lambda)
    C(n) = Clmn_amp(ell, m, n, lambda(n))
  end do

  do n = 1, size(lambda)
        Omega_ln = cmplx(dsqrt(mu**2 + lambda(n)), 0.0_real64, kind=real64)
        S%phi = S%phi + C(n) * eigenR(:,n)
        S%pi = S%pi - imag_unit * Omega_ln * C(n) * eigenR(:,n)
  end do

  deallocate(C)

  S%psi = dx(S%phi, Mesh)
  ! S%pi = - G%beta / G%alpha * S%psi


  ! Mass_field = S%NC(G,Mesh)
  ! A = G%M_bh / Mass_field

  ! S%phi = sqrt(A) * S%phi
  ! S%psi = sqrt(A) * S%psi
  ! S%pi  = sqrt(A) * S%pi


end subroutine mode_from_SL_all




! module initial_data
!   use iso_fortran_env, only: real64
!   use state,           only: state_t
!   use finite_differences, only: dx => first_derivative_x_2
!   implicit none
!   private
!   public :: set_IC_from_projection

! contains

!   !------------------------------------------------------------
!   ! Amplitud determinista "pseudo-aleatoria" para cada (ell,m)
!   ! Distinta para cada par, fija entre corridas, sin RNG.
!   !------------------------------------------------------------
!   pure complex(real64) function Clm_amp(ell, m) result(C)
!     integer, intent(in) :: ell, m
!     real(real64) :: a, b

!     a = sin( 1.234567_real64 * real(ell+1,real64) + &
!              2.345678_real64 * real(m+37,real64) )

!     b = cos( 3.456789_real64 * real(ell+1,real64) + &
!              4.567891_real64 * real(m+73,real64) )

!     C = cmplx(a, b, kind=real64)
! !    C = 1.0d0
!   end function Clm_amp

!   !------------------------------------------------------------
!   ! IC para un modo (ell,m):
!   !   phi(r,0) = C_{ell m} * exp( - (r-r0)^2 / (2 sigma^2) )
!   !   psi      = d phi / dr
!   !   pi       = (dtphi0 - G%beta psi)/G%alpha, con dtphi0 = psi
!   !------------------------------------------------------------
!   subroutine set_IC_from_projection(r, dr, G%alpha, G%beta, ell, m, Nth, Nph, S)
!     use iso_fortran_env, only: real64
!     real(real64), intent(in)     :: r(:), dr
!     real(real64), intent(in)     :: G%alpha(:), G%beta(:)
!     integer,      intent(in)     :: ell, m, Nth, Nph
!     type(state_t), intent(inout) :: S

!     integer        :: n, i
!     real(real64)   :: r0, sigma, alpha_eps
!     complex(real64):: C

!     n = size(r)
!     if (S%Nr /= n) stop "set_IC_from_projection: tamaños inconsistentes (S%Nr vs r)"
!     if (size(G%alpha) /= n .or. size(G%beta) /= n) stop "set_IC_from_projection: coef geom"

!     ! Amplitud fija por modo (ell,m)
!     C     = Clm_amp(ell, m)
!     r0    = 80.0_real64
!     sigma = 1.0_real64

!     ! phi(r)
!     do i = 1, n
!       S%phi(i) = C * exp( - (r(i) - r0)**2 / (2.0_real64 * sigma**2) )
!     end do

!     ! psi = d(phi)/dr
!     S%psi = dx(S%phi, dr)

!     ! pi = (dtphi0 - G%beta * psi)/G%alpha, con dtphi0 = psi
!     alpha_eps = 1.0e-14_real64
!     do i = 1, n
!       S%pi(i) = ( S%psi(i) - cmplx(G%beta(i),0.0_real64,kind=real64)*S%psi(i) ) / &
!                 cmplx( max(G%alpha(i),alpha_eps), 0.0_real64, kind=real64 )
!     end do

!   end subroutine set_IC_from_projection

end module initial_data
