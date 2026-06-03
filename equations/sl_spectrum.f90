module sl_spectrum
  !! Coeficientes complejos C_{ell m n} ~ CN(0, P(k))
  !! Replica Python: z = (N(0,1) + i N(0,1))/sqrt(2)
  !! C = sqrt(P(k)) * z, con P(k)=exp(-k^2/(2*sigma_k^2))
  !!
  implicit none
  private

  real(kind=8), parameter :: PI = acos(-1.0d0)
  real(kind=8) :: sigma_k = 0.1d0

  public :: set_sigma_k, P_of_k, Clmn_amp

contains

  subroutine set_sigma_k(s)
    real(kind=8), intent(in) :: s
    sigma_k = max(s, 1.0d-12)
  end subroutine set_sigma_k

  pure real(kind=8) function P_of_k(k) result(P)
    real(kind=8), intent(in) :: k
    real(kind=8) :: kk
    kk = max(k, 0.0d0)
    P  = dexp( - (kk**2) / (2.0d0 * sigma_k**2) )
  end function P_of_k

  subroutine randn2(z1, z2)
    !! Devuelve dos N(0,1) independientes (Box–Muller)
    real(kind=8), intent(out) :: z1, z2
    real(kind=8) :: u1, u2, s
    call random_number(u1); if (u1 <= 0.0d0) u1 = 1.0d-16
    call random_number(u2)
    s  = sqrt(-2.0d0*log(u1))
    z1 = s*cos(2.0d0*PI*u2)
    z2 = s*sin(2.0d0*PI*u2)
  end subroutine randn2

  function Clmn_amp(ell, m, n, lambda_n) result(C)
    integer,      intent(in) :: ell, m, n
    real(kind=8), intent(in) :: lambda_n
    complex(kind=8)          :: C
    real(kind=8) :: k, std, z1, z2

    k   = dsqrt(max(lambda_n, 0.0d0))
    std = dsqrt(max(P_of_k(k), 0.0d0))

    call randn2(z1, z2)  ! z1,z2 ~ N(0,1)
    C = std * cmplx(z1, z2) / sqrt(2.0d0)
  end function Clmn_amp

end module sl_spectrum

