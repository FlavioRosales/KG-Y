module plane_projection
  use mesh, only: mesh_t
  implicit none
  private

  integer, parameter, public :: PLANE_XY = 1
  integer, parameter, public :: PLANE_XZ = 2
  integer, parameter, public :: PART_RE  = 1
  integer, parameter, public :: PART_IM  = 2

  real(kind=8), parameter :: PI = acos(-1.0d0)

  type, public :: plane_projector_t
    integer :: nr     = 0
    integer :: nlocal = 0
    integer :: nphi   = 0
    integer :: ntheta = 0

    real(kind=8), allocatable :: r(:)
    real(kind=8), allocatable :: varphi(:)
    real(kind=8), allocatable :: theta(:)

    ! Y_xy(angle, local_mode), theta = pi/2.
    real(kind=8), allocatable :: yxy_re(:,:), yxy_im(:,:)

    ! Y_xz(angle, local_mode), phi = 0.
    real(kind=8), allocatable :: yxz_re(:,:), yxz_im(:,:)
  end type plane_projector_t

  public :: plane_projector_init
  public :: plane_projector_free
  public :: plane_projector_add_xy
  public :: plane_projector_add_xz

  interface
    subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
      character(len=1), intent(in) :: transa, transb
      integer,          intent(in) :: m, n, k, lda, ldb, ldc
      real(kind=8),     intent(in) :: alpha, beta
      real(kind=8),     intent(in) :: a(lda,*), b(ldb,*)
      real(kind=8),     intent(inout) :: c(ldc,*)
    end subroutine dgemm
  end interface

contains

  subroutine plane_projector_init(P, M, ell_local, mm_local, r_plane_max, nphi, ntheta)
    type(plane_projector_t), intent(inout) :: P
    type(mesh_t),            intent(in)    :: M
    integer,                 intent(in)    :: ell_local(:), mm_local(:)
    real(kind=8),            intent(in)    :: r_plane_max
    integer,                 intent(in)    :: nphi, ntheta

    integer :: q, ia

    call plane_projector_free(P)

    P%nlocal = size(ell_local)
    P%nr     = size(M%r)
    if (r_plane_max > 0.0d0) P%nr = count(M%r <= r_plane_max)

    P%nphi   = nphi
    P%ntheta = ntheta

    allocate(P%r(P%nr))
    allocate(P%varphi(P%nphi))
    allocate(P%theta(P%ntheta))

    allocate(P%yxy_re(P%nphi,   P%nlocal))
    allocate(P%yxy_im(P%nphi,   P%nlocal))
    allocate(P%yxz_re(P%ntheta, P%nlocal))
    allocate(P%yxz_im(P%ntheta, P%nlocal))

    P%r = M%r(1:P%nr)

    do ia = 1, P%nphi
      P%varphi(ia) = 2.0d0*PI*dble(ia - 1)/dble(P%nphi)
    end do

    if (P%ntheta == 1) then
      P%theta(1) = 0.0d0
    else
      do ia = 1, P%ntheta
        P%theta(ia) = PI*dble(ia - 1)/dble(P%ntheta - 1)
      end do
    end if

    do q = 1, P%nlocal

      do ia = 1, P%nphi
        call spherical_harmonic(ell_local(q), mm_local(q), 0.5d0*PI, P%varphi(ia), &
                                P%yxy_re(ia,q), P%yxy_im(ia,q))
      end do

      do ia = 1, P%ntheta
        call spherical_harmonic(ell_local(q), mm_local(q), P%theta(ia), 0.0d0, &
                                P%yxz_re(ia,q), P%yxz_im(ia,q))
      end do

    end do
  end subroutine plane_projector_init


  subroutine plane_projector_free(P)
    type(plane_projector_t), intent(inout) :: P

    if (allocated(P%r))       deallocate(P%r)
    if (allocated(P%varphi))  deallocate(P%varphi)
    if (allocated(P%theta))   deallocate(P%theta)
    if (allocated(P%yxy_re))  deallocate(P%yxy_re)
    if (allocated(P%yxy_im))  deallocate(P%yxy_im)
    if (allocated(P%yxz_re))  deallocate(P%yxz_re)
    if (allocated(P%yxz_im))  deallocate(P%yxz_im)

    P%nr     = 0
    P%nlocal = 0
    P%nphi   = 0
    P%ntheta = 0
  end subroutine plane_projector_free


  ! Accumulates the selected real or imaginary component of
  ! sum_q u_q(r) Y_q(pi/2,varphi) into plane(nr,nphi).
  subroutine plane_projector_add_xy(P, first_mode, u_re, u_im, part, plane)
    type(plane_projector_t), intent(in)    :: P
    integer,                 intent(in)    :: first_mode
    real(kind=8),            intent(in)    :: u_re(:,:), u_im(:,:)
    integer,                 intent(in)    :: part
    real(kind=8),            intent(inout) :: plane(:,:)

    integer :: nb

    nb = size(u_re, 2)

    if (part == PART_RE) then
      call dgemm('N', 'T', P%nr, P%nphi, nb,  1.0d0, u_re, P%nr, &
                 P%yxy_re(:,first_mode:first_mode+nb-1), P%nphi, 1.0d0, plane, P%nr)

      call dgemm('N', 'T', P%nr, P%nphi, nb, -1.0d0, u_im, P%nr, &
                 P%yxy_im(:,first_mode:first_mode+nb-1), P%nphi, 1.0d0, plane, P%nr)
    else
      call dgemm('N', 'T', P%nr, P%nphi, nb,  1.0d0, u_re, P%nr, &
                 P%yxy_im(:,first_mode:first_mode+nb-1), P%nphi, 1.0d0, plane, P%nr)

      call dgemm('N', 'T', P%nr, P%nphi, nb,  1.0d0, u_im, P%nr, &
                 P%yxy_re(:,first_mode:first_mode+nb-1), P%nphi, 1.0d0, plane, P%nr)
    end if
  end subroutine plane_projector_add_xy


  ! Accumulates the selected real or imaginary component of
  ! sum_q u_q(r) Y_q(theta,0) into plane(nr,ntheta).
  ! For the x<0 half-plane, main passes coefficients multiplied by (-1)^m.
  subroutine plane_projector_add_xz(P, first_mode, u_re, u_im, part, plane)
    type(plane_projector_t), intent(in)    :: P
    integer,                 intent(in)    :: first_mode
    real(kind=8),            intent(in)    :: u_re(:,:), u_im(:,:)
    integer,                 intent(in)    :: part
    real(kind=8),            intent(inout) :: plane(:,:)

    integer :: nb

    nb = size(u_re, 2)

    if (part == PART_RE) then
      call dgemm('N', 'T', P%nr, P%ntheta, nb,  1.0d0, u_re, P%nr, &
                 P%yxz_re(:,first_mode:first_mode+nb-1), P%ntheta, 1.0d0, plane, P%nr)

      call dgemm('N', 'T', P%nr, P%ntheta, nb, -1.0d0, u_im, P%nr, &
                 P%yxz_im(:,first_mode:first_mode+nb-1), P%ntheta, 1.0d0, plane, P%nr)
    else
      call dgemm('N', 'T', P%nr, P%ntheta, nb,  1.0d0, u_re, P%nr, &
                 P%yxz_im(:,first_mode:first_mode+nb-1), P%ntheta, 1.0d0, plane, P%nr)

      call dgemm('N', 'T', P%nr, P%ntheta, nb,  1.0d0, u_im, P%nr, &
                 P%yxz_re(:,first_mode:first_mode+nb-1), P%ntheta, 1.0d0, plane, P%nr)
    end if
  end subroutine plane_projector_add_xz


  subroutine spherical_harmonic(ell, mm, theta, varphi, yre, yim)
    integer,      intent(in)  :: ell, mm
    real(kind=8), intent(in)  :: theta, varphi
    real(kind=8), intent(out) :: yre, yim

    integer      :: mabs, sign_m
    real(kind=8) :: qlm, phase

    mabs = abs(mm)

    call normalized_legendre(ell, mabs, dcos(theta), qlm)

    phase = dble(mabs)*varphi
    yre   = qlm*dcos(phase)
    yim   = qlm*dsin(phase)

    if (mm < 0) then
      if (mod(mabs,2) == 0) then
        sign_m = 1
      else
        sign_m = -1
      end if
      yre = dble(sign_m)*yre
      yim = -dble(sign_m)*yim
    end if
  end subroutine spherical_harmonic


  ! qlm = sqrt((2l+1)/(4pi) * (l-m)!/(l+m)!) P_l^m(x)
  ! with the Condon-Shortley phase included in P_l^m.
  subroutine normalized_legendre(ell, mm, x_in, qlm)
    integer,      intent(in)  :: ell, mm
    real(kind=8), intent(in)  :: x_in
    real(kind=8), intent(out) :: qlm

    integer      :: j
    real(kind=8) :: x, s, qmm, qlm1, qlm2, qnew, a, b

    x = max(-1.0d0, min(1.0d0, x_in))
    s = dsqrt(max(0.0d0, 1.0d0 - x*x))

    qmm = 1.0d0/dsqrt(4.0d0*PI)

    do j = 1, mm
      qmm = -dsqrt(dble(2*j + 1)/dble(2*j))*s*qmm
    end do

    if (ell == mm) then
      qlm = qmm
      return
    end if

    qlm2 = qmm
    qlm1 = dsqrt(dble(2*mm + 3))*x*qmm

    if (ell == mm + 1) then
      qlm = qlm1
      return
    end if

    do j = mm + 2, ell
      a = dsqrt(dble(4*j*j - 1)/dble(j*j - mm*mm))
      b = dsqrt(dble((2*j + 1)*(j + mm - 1)*(j - mm - 1)) / &
                dble((2*j - 3)*(j + mm)*(j - mm)))

      qnew = a*x*qlm1 - b*qlm2
      qlm2 = qlm1
      qlm1 = qnew
    end do

    qlm = qlm1
  end subroutine normalized_legendre

end module plane_projection
