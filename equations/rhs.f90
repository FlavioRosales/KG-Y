module rhs
  use geometry,        only: geometry_t
  use mesh,            only: mesh_t
  use state,           only: state_t
  implicit none
  private
  public :: compute_rhs_state

  integer, parameter :: Re = 1, Im = 2

contains

  subroutine compute_rhs_state(S, G, M, ell, mu, dphi_dt, dpsi_dt, dpi_dt)
    use finite_differences, &
      only: dx    => first_derivative_x_2, &
            advec => advec_x, & 
            advec_full, & 
            dx_esc => first_derivative_escaled
    implicit none
    ! ---- Entradas ----
    class(state_t),    intent(in)  :: S
    class(geometry_t), intent(in)  :: G
    class(mesh_t),     intent(in)  :: M
    integer,           intent(in)  :: ell
    real(kind=8),      intent(in)  :: mu
    ! ---- Salidas ----
    real(kind=8),   intent(out) :: dphi_dt(:,:)
    real(kind=8),   intent(out) :: dpsi_dt(:,:)
    real(kind=8),   intent(out) :: dpi_dt(:,:)
    ! ----- Variables locales -----
    integer :: i, n
    real(kind = 8) :: mu2, ell_factor

    n = size(dphi_dt, dim=1)



    ! --- Ecuación de φ:  ∂_t φ = α π + β ψ ---
    do concurrent (i = 1:n)
      dphi_dt(i,Re) = G%alpha(i)*S%pi(i,Re) + G%beta(i)*S%psi(i,Re)
      dphi_dt(i,Im) = G%alpha(i)*S%pi(i,Im) + G%beta(i)*S%psi(i,Im)
    end do
    
    ! --- Ecuación de ψ ---
    ! dψ/dt = ∂_r(α π) + (dβ/dr) ψ + advección(ψ; β)

    dpsi_dt =  dx_esc(G%alpha, S%pi, M) + advec_full(S%psi, G, M)

    ! --- Ecuación de π ---
    ! dπ/dt = advección(π; β) + α/M/r^2 ∂_r(r^2 guu ψ) + α (K_c π - (μ^2 + l(l+1)/r^2) φ)
    mu2 = mu**2
    ell_factor = dble(ell*(ell + 1)) 

    
    dpi_dt = advec(S%pi, G, M) + dx_esc(G%alpha * G%r2_inv, M%r**2 * G%guu, S%psi, M) 
    do concurrent (i = 1:n)
      dpi_dt(i,Re) = dpi_dt(i,Re)  & 
            + G%alpha(i)*(G%K_c(i)*S%pi(i,Re) - (mu2 + ell_factor*G%r2_inv(i))*S%phi(i,Re) )
      dpi_dt(i,Im) = dpi_dt(i,Im)  & 
            + G%alpha(i)*(G%K_c(i)*S%pi(i,Im) - (mu2 + ell_factor*G%r2_inv(i))*S%phi(i,Im) )
    end do

  end subroutine compute_rhs_state



end module rhs
