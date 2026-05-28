module rhs
  use iso_fortran_env, only: real64
  use geometry,        only: geometry_t
  use mesh,            only: mesh_t
  use state,           only: state_t
  implicit none
  private
  public :: compute_rhs_state

contains

  subroutine compute_rhs_state(S, G, M, ell, mu, dphi_dt, dpsi_dt, dpi_dt)
    use finite_differences, &
      only: dx    => first_derivative_x_2, &
            advec => advec_x
    implicit none
    ! ---- Entradas ----
    class(state_t),    intent(in)  :: S
    class(geometry_t), intent(in)  :: G
    class(mesh_t),     intent(in)  :: M
    integer,           intent(in)  :: ell
    real(real64),      intent(in)  :: mu
    ! ---- Salidas ----
    complex(real64),   intent(out) :: dphi_dt(:)
    complex(real64),   intent(out) :: dpsi_dt(:)
    complex(real64),   intent(out) :: dpi_dt(:)

    ! --- Ecuación de φ:  ∂_t φ = α π + β ψ ---
    dphi_dt = G%alpha * S%pi + G%beta * S%psi

    ! --- Ecuación de ψ ---
    ! dψ/dt = ∂_r(α π) + (dβ/dr) ψ + advección(ψ; β)
    dpsi_dt = dx(G%alpha * S%pi, M)    &
            + G%dbeta  * S%psi        &
            + advec(S%psi, G, M)

    ! --- Ecuación de π ---
    dpi_dt =  advec(S%pi, G, M)                                             &
            +  G%alpha / M%r**2 * dx( M%r**2 * G%guu * S%psi, M )            &
            +  G%alpha * ( G%K_c * S%pi                                      &
            -  ( mu**2 + dble(ell*(ell + 1)) / M%r**2 ) * S%phi )

  end subroutine compute_rhs_state

end module rhs
