module utils
  use iso_fortran_env, only: real64, int32
  ! SHTOOLS: Legendre asociados ortonormalizados y util de indexado
  !use SHTOOLS,       only: PlmON, PlmIndex
  implicit none
  private
 ! public :: Ylm_complex, project_phi_lm_uniform
 ! public :: project_phi_lm_GL
 ! public :: gl_nodes_weights
 ! public :: phi_sampler
  public :: trapezium
  public :: simpson

  ! --- Interfaz de muestreador de φ(r,θ,φ) complejo ---
  ! abstract interface
  !   complex(real64) function phi_sampler(r, th, ph)
  !     use iso_fortran_env, only: real64
  !     real(real64), intent(in) :: r, th, ph
  !   end function phi_sampler
  ! end interface


contains
  ! !========================================================
  ! !  Y_lm complejo, ortonormal en S^2: ∫|Y|^2 dΩ = 1
  ! !  Implementado con SHTOOLS: PlmON (cnorm=1, csphase=-1)
  ! !  Y_lm(θ,φ) = \bar P_l^m(cosθ) e^{i m φ}
  ! !  Para m<0: Y_{l,-m} = (-1)^m conj(Y_{l,m})
  ! !========================================================
  ! complex(real64) function Ylm_complex(ell, m, th, ph) result(Y)
  !   integer, intent(in) :: ell, m
  !   real(real64), intent(in) :: th, ph
  !   integer(int32) :: lmax32, csphase, cnorm, idx
  !   real(real64), allocatable :: P(:)
  !   integer :: mm
  !   complex(real64) :: eimp

  !   if (abs(m) > ell) then
  !     Y = (0.0_real64,0.0_real64); return
  !   end if

  !   lmax32 = int(ell, int32)
  !   csphase = -1_int32     ! incluye fase de Condon–Shortley (-1)^m
  !   cnorm   =  1_int32     ! normalización "complex": ∫_{-1}^1 P^2 dx = 1/(2π)

  !   allocate(P( (ell+1)*(ell+2)/2 ))
  !   call PlmON(P, lmax32, cos(th), csphase, cnorm)  ! SHTOOLS, estable hasta l~2800
  !   mm  = abs(m)
  !   idx = ell*(ell+1)/2 + mm + 1

  !   eimp = cmplx(cos(real(mm,real64)*ph), sin(real(mm,real64)*ph), kind=real64)
  !   if (m >= 0) then
  !     Y = P(idx) * eimp
  !   else
  !     ! Y_{l,-m} = (-1)^m conj(Y_{l,m})
  !     Y = ((-1.0_real64)**mm) * conjg( P(idx) * eimp )
  !   end if
  !   deallocate(P)
  ! end function

  ! !========================================================
  ! !  Proyección uniforme (midpoint) con peso sinθ
  ! !  φ_lm(r_i) = ∫ φ(r_i,θ,φ) conj(Y_lm) sinθ dθ dφ
  ! !  Aceleraciones:
  ! !    - e^{-i m φ_k} precomputado
  ! !    - \bar P_l^{|m|}(cosθ_j) desde SHTOOLS (P es real)
  ! !========================================================
  ! subroutine project_phi_lm_uniform(r, ell, m, Nth, Nph, sampler, phi_lm)
  !   real(real64), intent(in)     :: r(:)
  !   integer,      intent(in)     :: ell, m, Nth, Nph
  !   procedure(phi_sampler)       :: sampler
  !   complex(real64), intent(out) :: phi_lm(:)

  !   integer :: i, j, k, Nr, mm
  !   real(real64) :: th, ph, dth, dph, w, ri
  !   complex(real64) :: acc
  !   complex(real64), allocatable :: epos(:), eneg(:)
  !   ! Legendre ortonormalizados en el θ actual:
  !   real(real64), allocatable :: P(:)
  !   integer(int32) :: csphase, cnorm, lmax32, idx

  !   Nr = size(r)
  !   if (size(phi_lm) /= Nr) stop "project_phi_lm_uniform: tamaño incompatible"
  !   if (Nth < 1 .or. Nph < 1) stop "project_phi_lm_uniform: Nth/Nph >= 1"

  !   dth = acos(-1.0_real64) / real(Nth,real64)
  !   dph = 2.0_real64*acos(-1.0_real64) / real(Nph,real64)
  !   mm  = abs(m)

  !   allocate(epos(Nph), eneg(Nph))
  !   do k = 1, Nph
  !     ph   = (real(k,real64)-0.5_real64) * dph
  !     epos(k) = cmplx(cos(real(mm,real64)*ph),  sin(real(mm,real64)*ph),  kind=real64)  ! e^{+i m φ}
  !     eneg(k) = cmplx(cos(real(mm,real64)*ph), -sin(real(mm,real64)*ph),  kind=real64)  ! e^{-i m φ}
  !   end do

  !   lmax32 = int(ell, int32); csphase = -1_int32; cnorm = 1_int32
  !   allocate(P( (ell+1)*(ell+2)/2 ))

  !   do i = 1, Nr
  !     acc = (0.0_real64,0.0_real64)
  !     ri  = r(i)
  !     do j = 0, Nth-1
  !       th = (j + 0.5_real64) * dth
  !       w  = sin(th) * dth
  !       call PlmON(P, lmax32, cos(th), csphase, cnorm)
  !       idx = PlmIndex(int(ell,int32), int(mm,int32))
  !       if (m >= 0) then
  !         do k = 1, Nph
  !           acc = acc + sampler(ri, th, (real(k,real64)-0.5_real64)*dph) * (P(idx) * eneg(k)) * w * dph
  !         end do
  !       else
  !         ! conj(Y_{l,-m}) = (-1)^m Y_{l,m} = (-1)^mm P * e^{+i mm φ}
  !         do k = 1, Nph
  !           acc = acc + sampler(ri, th, (real(k,real64)-0.5_real64)*dph) * (((-1.0_real64)**mm) * P(idx) * epos(k)) * w * dph
  !         end do
  !       end if
  !     end do
  !     phi_lm(i) = acc
  !   end do

  !   deallocate(P, epos, eneg)
  ! end subroutine

  ! !========================================================
  ! !  Gauss–Legendre nodes/weights on [-1,1] (orden N>=1)
  ! !========================================================
  ! subroutine gl_nodes_weights(N, x, w)
  !   integer,      intent(in)  :: N
  !   real(real64), intent(out) :: x(N), w(N)
  !   integer :: i, j, maxit
  !   real(real64) :: xi, dx, Pn, dPn, tol
  !   if (N < 1) stop "gl_nodes_weights: N must be >= 1"
  !   tol = 1.0e-14_real64; maxit = 50

  !   do i = 1, N
  !     xi = cos( (real(2*i-1,real64) * acos(-1.0_real64)) / real(2*N,real64) )
  !     do j = 1, maxit
  !       call legendre_P_and_dP(N, xi, Pn, dPn)
  !       dx = -Pn / dPn;  xi = xi + dx
  !       if (abs(dx) <= tol) exit
  !     end do
  !     x(i) = xi
  !     w(i) = 2.0_real64 / ( (1.0_real64 - xi*xi) * dPn*dPn )
  !   end do
  !   call sort_increasing(x, w)
  ! end subroutine

  ! subroutine legendre_P_and_dP(N, x, Pn, dPn)
  !   integer,      intent(in)  :: N
  !   real(real64), intent(in)  :: x
  !   real(real64), intent(out) :: Pn, dPn
  !   real(real64) :: Pnm1, Pnm2, dPnm1, dPnm2
  !   integer :: nn
  !   if (N == 0) then
  !     Pn  = 1.0_real64; dPn = 0.0_real64; return
  !   else if (N == 1) then
  !     Pn  = x; dPn = 1.0_real64; return
  !   end if
  !   Pnm2=1.0_real64; Pnm1=x; dPnm2=0.0_real64; dPnm1=1.0_real64
  !   do nn = 2, N
  !     Pn  = ( (2.0_real64*nn-1.0_real64)*x*Pnm1 - (nn-1.0_real64)*Pnm2 ) / real(nn,real64)
  !     dPn = ( (2.0_real64*nn-1.0_real64)*(Pnm1 + x*dPnm1) - (nn-1.0_real64)*dPnm2 ) / real(nn,real64)
  !     Pnm2=Pnm1; Pnm1=Pn; dPnm2=dPnm1; dPnm1=dPn
  !   end do
  ! end subroutine

  ! subroutine sort_increasing(x, w)
  !   real(real64), intent(inout) :: x(:), w(:)
  !   integer :: i, j, n; real(real64) :: tx, tw
  !   n = size(x)
  !   do i = 1, n-1
  !     do j = i+1, n
  !       if (x(j) < x(i)) then
  !         tx=x(i); x(i)=x(j); x(j)=tx
  !         tw=w(i); w(i)=w(j); w(j)=tw
  !       end if
  !     end do
  !   end do
  ! end subroutine

  ! !========================================================
  ! !  Proyección GL en θ y trapecio en φ (midpoints)
  ! !  Acelerada con SHTOOLS y exponenciales precomputados
  ! !========================================================
  ! subroutine project_phi_lm_GL(r, ell, m, Nth, Nph, sampler, phi_lm)
  !   real(real64), intent(in)     :: r(:)
  !   integer,      intent(in)     :: ell, m, Nth, Nph
  !   procedure(phi_sampler)       :: sampler
  !   complex(real64), intent(out) :: phi_lm(:)

  !   integer :: i, j, k, Nr, mm
  !   real(real64) :: dph, ph, ri, th
  !   real(real64), allocatable :: x(:), w(:)
  !   complex(real64) :: acc
  !   complex(real64), allocatable :: epos(:), eneg(:)
  !   real(real64), allocatable :: P(:)
  !   integer(int32) :: lmax32, csphase, cnorm, idx

  !   Nr = size(r)
  !   if (size(phi_lm) /= Nr) stop "project_phi_lm_GL: tamaño incompatible"
  !   if (Nth < 1 .or. Nph < 1) stop "project_phi_lm_GL: Nth/Nph >= 1"

  !   allocate(x(Nth), w(Nth))
  !   call gl_nodes_weights(Nth, x, w)           ! en [-1,1]
  !   dph = 2.0_real64*acos(-1.0_real64) / real(Nph,real64)

  !   mm = abs(m)
  !   allocate(epos(Nph), eneg(Nph))
  !   do k = 1, Nph
  !     ph = (real(k,real64)-0.5_real64) * dph
  !     epos(k) = cmplx(cos(real(mm,real64)*ph),  sin(real(mm,real64)*ph),  kind=real64) ! e^{+i m φ}
  !     eneg(k) = cmplx(cos(real(mm,real64)*ph), -sin(real(mm,real64)*ph),  kind=real64) ! e^{-i m φ}
  !   end do

  !   lmax32 = int(ell, int32); csphase = -1_int32; cnorm = 1_int32
  !   allocate(P( (ell+1)*(ell+2)/2 ))

  !   do i = 1, Nr
  !     acc = (0.0_real64,0.0_real64)
  !     ri  = r(i)
  !     do j = 1, Nth
  !       th = acos( x(j) )                       ! θ ∈ [0,π]
  !       call PlmON(P, lmax32, x(j), csphase, cnorm)
  !       idx = PlmIndex(int(ell,int32), int(mm,int32))
  !       if (m >= 0) then
  !         do k = 1, Nph
  !           acc = acc + sampler(ri, th, (real(k,real64)-0.5_real64)*dph) * (P(idx) * eneg(k)) * w(j) * dph
  !         end do
  !       else
  !         do k = 1, Nph
  !           acc = acc + sampler(ri, th, (real(k,real64)-0.5_real64)*dph) * (((-1.0_real64)**mm) * P(idx) * epos(k)) * w(j) * dph
  !         end do
  !       end if
  !     end do
  !     phi_lm(i) = acc
  !   end do

  !   deallocate(P, epos, eneg, x, w)
  ! end subroutine

  function trapezium(f, dx) result(integral)
  !-----------------------------------------------------------------------
  ! Integración trapezoidal 1D para datos complejos (vectorizada)
  ! integral ≈ dx * [ ½(f₁ + f_N) + sum_{i=2}^{N-1} f_i ]
  !-----------------------------------------------------------------------
  implicit none
  real(kind=8),    intent(in) :: dx
  real(kind=8), intent(in) :: f(:)
  real(kind=8)             :: integral
  integer                     :: n

  n = size(f)

  integral = dx * ( 0.5d0*(f(1) + f(n)) + sum(f(2:n-1)) )

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
