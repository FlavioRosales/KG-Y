module time_integrators
  use geometry, only: geometry_t
  use mesh,     only: mesh_t
  use state,    only: state_t
  use rhs,      only: compute_rhs_state
  implicit none
  private

  integer, parameter :: Re = 1, Im = 2

  public :: ssprk3_step_state

contains

  subroutine ssprk3_step_state(dt, S, G, M, ell, mu)
    implicit none

    real(kind=8),     intent(in)    :: dt
    type(state_t),    intent(inout) :: S
    type(geometry_t), intent(in)    :: G
    type(mesh_t),     intent(in)    :: M
    integer,          intent(in)    :: ell
    real(kind=8),     intent(in)    :: mu

    real(kind=8) :: phi_p(S%Nr,2), psi_p(S%Nr,2), pi_p(S%Nr,2)
    real(kind=8) :: dphi(S%Nr,2),  dpsi(S%Nr,2),  dpi(S%Nr,2)

    integer :: i, rk, n

    n = S%Nr

    ! Garantiza fronteras consistentes antes del primer RHS.
    call bc_apply_state(S, G, M)

    ! Estado de referencia en t^n.
    phi_p = S%phi
    psi_p = S%psi
    pi_p  = S%pi

    do rk = 1, 3

      call compute_rhs_state(S, G, M, ell, mu, dphi, dpsi, dpi)

      select case (rk)

      case (1)

        do concurrent (i = 1:n)
          S%phi(i,Re) = phi_p(i,Re) + dt*dphi(i,Re)
          S%phi(i,Im) = phi_p(i,Im) + dt*dphi(i,Im)

          S%psi(i,Re) = psi_p(i,Re) + dt*dpsi(i,Re)
          S%psi(i,Im) = psi_p(i,Im) + dt*dpsi(i,Im)

          S%pi(i,Re)  = pi_p(i,Re)  + dt*dpi(i,Re)
          S%pi(i,Im)  = pi_p(i,Im)  + dt*dpi(i,Im)
        end do

      case (2)

        do concurrent (i = 1:n)
          S%phi(i,Re) = 0.75d0*phi_p(i,Re) + 0.25d0*(S%phi(i,Re) + dt*dphi(i,Re))
          S%phi(i,Im) = 0.75d0*phi_p(i,Im) + 0.25d0*(S%phi(i,Im) + dt*dphi(i,Im))

          S%psi(i,Re) = 0.75d0*psi_p(i,Re) + 0.25d0*(S%psi(i,Re) + dt*dpsi(i,Re))
          S%psi(i,Im) = 0.75d0*psi_p(i,Im) + 0.25d0*(S%psi(i,Im) + dt*dpsi(i,Im))

          S%pi(i,Re)  = 0.75d0*pi_p(i,Re)  + 0.25d0*(S%pi(i,Re)  + dt*dpi(i,Re))
          S%pi(i,Im)  = 0.75d0*pi_p(i,Im)  + 0.25d0*(S%pi(i,Im)  + dt*dpi(i,Im))
        end do

      case (3)

        do concurrent (i = 1:n)
          S%phi(i,Re) = (1.0d0/3.0d0)*phi_p(i,Re) &
                      + (2.0d0/3.0d0)*(S%phi(i,Re) + dt*dphi(i,Re))
          S%phi(i,Im) = (1.0d0/3.0d0)*phi_p(i,Im) &
                      + (2.0d0/3.0d0)*(S%phi(i,Im) + dt*dphi(i,Im))

          S%psi(i,Re) = (1.0d0/3.0d0)*psi_p(i,Re) &
                      + (2.0d0/3.0d0)*(S%psi(i,Re) + dt*dpsi(i,Re))
          S%psi(i,Im) = (1.0d0/3.0d0)*psi_p(i,Im) &
                      + (2.0d0/3.0d0)*(S%psi(i,Im) + dt*dpsi(i,Im))

          S%pi(i,Re)  = (1.0d0/3.0d0)*pi_p(i,Re) &
                      + (2.0d0/3.0d0)*(S%pi(i,Re) + dt*dpi(i,Re))
          S%pi(i,Im)  = (1.0d0/3.0d0)*pi_p(i,Im) &
                      + (2.0d0/3.0d0)*(S%pi(i,Im) + dt*dpi(i,Im))
        end do

      end select

      call bc_apply_state(S, G, M)

    end do

  end subroutine ssprk3_step_state


  subroutine bc_apply_state(S, G, M)
    implicit none

    class(state_t),    intent(inout) :: S
    class(geometry_t), intent(in)    :: G
    class(mesh_t),     intent(in)    :: M

    integer :: n

    n = size(S%phi, dim=1)

    ! Frontera interna: extrapolación lineal.
    S%phi(1,Re) = 2.0d0*S%phi(2,Re) - S%phi(3,Re)
    S%phi(1,Im) = 2.0d0*S%phi(2,Im) - S%phi(3,Im)

    S%psi(1,Re) = 2.0d0*S%psi(2,Re) - S%psi(3,Re)
    S%psi(1,Im) = 2.0d0*S%psi(2,Im) - S%psi(3,Im)

    S%pi(1,Re)  = 2.0d0*S%pi(2,Re)  - S%pi(3,Re)
    S%pi(1,Im)  = 2.0d0*S%pi(2,Im)  - S%pi(3,Im)

    ! Frontera externa: gradiente cero.
    S%phi(n,Re) = S%phi(n-1,Re)
    S%phi(n,Im) = S%phi(n-1,Im)

    S%psi(n,Re) = S%psi(n-1,Re)
    S%psi(n,Im) = S%psi(n-1,Im)

    S%pi(n,Re)  = S%pi(n-1,Re)
    S%pi(n,Im)  = S%pi(n-1,Im)

  end subroutine bc_apply_state

end module time_integrators