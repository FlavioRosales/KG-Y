module time_integrators
  use geometry,        only: geometry_t
  use mesh,            only: mesh_t
  use state,           only: state_t
  use rhs,             only: compute_rhs_state
  implicit none
  private

  public :: ssprk3_step_state
  public :: bc_init_from_geometry, bc_set_rgrid

  ! ---- Almacenamiento opcional de coeficientes de BC / malla ----
  real(kind=8), allocatable, save :: rgrid(:), alpha(:), bta(:), guu(:)

contains

  !============================
  ! Inicialización (opcional) de BC desde la geometría
  !============================
  subroutine bc_init_from_geometry(G)
    type(geometry_t), intent(in) :: G
    call assign_or_copy(alpha, G%alpha)
    call assign_or_copy(bta,   G%beta)
    call assign_or_copy(guu,   G%guu)
  end subroutine bc_init_from_geometry

  subroutine bc_set_rgrid(r)
    real(kind=8), intent(in) :: r(:)
    call assign_or_copy(rgrid, r)
  end subroutine bc_set_rgrid

  !============================
  ! SSP-RK3 actuando directamente sobre state_t
  !============================
  subroutine ssprk3_step_state(dt, S, G, M, ell, mu, &
                              phi_p, psi_p, pi_p)
    use rhs,             only: compute_rhs_state
    implicit none
    real(kind=8),  intent(in)    :: dt
    type(state_t), intent(inout) :: S
    type(geometry_t), intent(in) :: G
    type(mesh_t),     intent(in) :: M
    integer,          intent(in) :: ell
    real(kind=8),     intent(in) :: mu

    complex(kind=8), intent(inout) :: phi_p(:), psi_p(:), pi_p(:)
    complex(kind=8) :: dphi(S%Nr), dpsi(S%Nr), dpi(S%Nr)

    integer :: N, rk

    N = S%Nr

    ! Guardamos estado de referencia en t^n
    phi_p = S%phi
    psi_p = S%psi
    pi_p  = S%pi

    do rk = 1, 3

      call compute_rhs_state(S, G, M, ell, mu, dphi, dpsi, dpi)

      select case (rk)
      case (1)
        S%phi = phi_p + dt*dphi
        S%psi = psi_p + dt*dpsi
        S%pi  = pi_p  + dt*dpi

      case (2)
        S%phi = (3.0d0/4.0d0)*phi_p + &
                (1.0d0/4.0d0)*(S%phi + dt*dphi)
        S%psi = (3.0d0/4.0d0)*psi_p + &
                (1.0d0/4.0d0)*(S%psi + dt*dpsi)
        S%pi  = (3.0d0/4.0d0)*pi_p  + &
                (1.0d0/4.0d0)*(S%pi  + dt*dpi)

      case (3)
        S%phi = (1.0d0/3.0d0)*phi_p + &
                (2.0d0/3.0d0)*(S%phi + dt*dphi)
        S%psi = (1.0d0/3.0d0)*psi_p + &
                (2.0d0/3.0d0)*(S%psi + dt*dpsi)
        S%pi  = (1.0d0/3.0d0)*pi_p  + &
                (2.0d0/3.0d0)*(S%pi  + dt*dpi)
      end select

      !call bc_apply_state(S) 
      call bc_apply_state(S, G, M)
    end do

  end subroutine ssprk3_step_state


subroutine bc_apply_state(S, G, M)
  implicit none
  class(state_t),    intent(inout) :: S
  class(geometry_t), intent(in)    :: G
  class(mesh_t),     intent(in)    :: M

  integer :: N
  real(kind=8) :: sN
  complex(kind=8) :: wplus_bnd, wminus_bnd

  N = size(S%phi)
  if (N < 4) stop "bc_apply_state: N < 4 no soportado"

  !===========================================================
  ! 1) Frontera interna (i = 1): extrapolación 2º orden
  !===========================================================
  !S%phi(1) = S%phi(2) !3.0d0*S%phi(2) - 3.0d0*S%phi(3) + S%phi(4)
  !S%psi(1) = S%psi(2) !3.0d0*S%psi(2) - 3.0d0*S%psi(3) + S%psi(4)
  !S%pi(1)  = S%pi(2)  !3.0d0*S%pi(2)  - 3.0d0*S%pi(3)  + S%pi(4)
S%phi(1) = 2.0d0*S%phi(2) - S%phi(3)
S%psi(1) = 2.0d0*S%psi(2) - S%psi(3)
S%pi(1)  = 2.0d0*S%pi(2)  - S%pi(3)
  S%phi(N) = S%phi(N-1)
  S%psi(N) = S%psi(N-1)
  S%pi(N)  = S%pi(N-1)

end subroutine bc_apply_state





  !============================
  ! Utilidad para copiar/asignar arrays reales allocatables
  !============================
  subroutine assign_or_copy(a, b)
    real(kind=8), allocatable, intent(inout) :: a(:)
    real(kind=8),           intent(in)       :: b(:)
    if (allocated(a)) then
      if (size(a) /= size(b)) then
        deallocate(a)
        allocate(a(size(b)))
      end if
    else
      allocate(a(size(b)))
    end if
    a = b
  end subroutine assign_or_copy

end module time_integrators
