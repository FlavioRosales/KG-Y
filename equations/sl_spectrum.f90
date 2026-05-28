module sl_spectrum
  !! Coeficientes complejos C_{ell m n} ~ CN(0, P(k))
  !! Replica Python: z = (N(0,1) + i N(0,1))/sqrt(2)
  !! C = sqrt(P(k)) * z, con P(k)=exp(-k^2/(2*sigma_k^2))
  !!
  use iso_fortran_env, only: real64
  implicit none
  private

  real(real64), parameter :: PI = acos(-1.0_real64)
  real(real64) :: sigma_k = 0.1_real64

  public :: set_sigma_k, P_of_k, Clmn_amp

contains

  subroutine set_sigma_k(s)
    real(real64), intent(in) :: s
    sigma_k = max(s, 1.0e-12_real64)
  end subroutine set_sigma_k

  pure real(real64) function P_of_k(k) result(P)
    real(real64), intent(in) :: k
    real(real64) :: kk
    kk = max(k, 0.0_real64)
    P  = exp( - (kk**2) / (2.0_real64 * sigma_k**2) )
  end function P_of_k

  subroutine randn2(z1, z2)
    !! Devuelve dos N(0,1) independientes (Box–Muller)
    real(real64), intent(out) :: z1, z2
    real(real64) :: u1, u2, s
    call random_number(u1); if (u1 <= 0.0_real64) u1 = 1.0e-16_real64
    call random_number(u2)
    s  = sqrt(-2.0_real64*log(u1))
    z1 = s*cos(2.0_real64*PI*u2)
    z2 = s*sin(2.0_real64*PI*u2)
  end subroutine randn2

  function Clmn_amp(ell, m, n, lambda_n) result(C)
    integer,      intent(in) :: ell, m, n
    real(real64), intent(in) :: lambda_n
    complex(real64)          :: C
    real(real64) :: k, std, z1, z2

    k   = sqrt(max(lambda_n, 0.0_real64))
    std = sqrt(max(P_of_k(k), 0.0_real64))

    call randn2(z1, z2)  ! z1,z2 ~ N(0,1)
    C = std * cmplx(z1, z2, kind=real64) / sqrt(2.0_real64)
  end function Clmn_amp

end module sl_spectrum

