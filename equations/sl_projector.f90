module sl_projection
  !!
  !! Transformada espectral SL sobre el dominio activo de evolución.
  !!
  !!   a_lmn = <R_ln, phi_lm>_active
  !!
  !!   P_ln^local = sum_{m local} |a_lmn|^2
  !!
  !! No se invierte la Gram de la base restringida. La excisión rompe
  !! la ortogonalidad exacta de la base SL completa y esa inversión es
  !! numéricamente inestable.
  !!
  use mesh,     only: mesh_t
  use geometry, only: geometry_t
  use utils,    only: trapezium
  use hdf5_lib, only: read_sl_modes_h5

  implicit none
  private

  integer, parameter :: dp = kind(1.0d0)

  public :: sl_projector_t
  public :: sl_projector_load
  public :: sl_projector_project
  public :: sl_projector_power_phi
  public :: sl_projector_free

  type :: sl_projector_t
    integer :: ell   = -1
    integer :: nr    = 0
    integer :: nspec = 0

    real(dp), allocatable :: lambda(:)
    real(dp), allocatable :: k(:)

    ! R(i,n) = R_{ell n}(r_i), restringida al dominio activo.
    real(dp), allocatable :: R(:,:)
  end type sl_projector_t

contains

    subroutine sl_projector_load(P, filename, ell, M)

    type(sl_projector_t), intent(inout) :: P
    character(len=*),     intent(in)    :: filename
    integer,              intent(in)    :: ell
    type(mesh_t),         intent(in)    :: M

    real(kind=8), allocatable :: lambda_all(:)
    real(kind=8), allocatable :: eigenR_all(:,:)

    call read_sl_modes_h5(trim(filename), ell, lambda_all, eigenR_all, M)

    call sl_projector_free(P)

    P%ell   = ell
    P%nr    = M%Nr
    P%nspec = size(lambda_all)

    call move_alloc(lambda_all, P%lambda)
    call move_alloc(eigenR_all, P%R)

    allocate(P%k(P%nspec))
    P%k = sqrt(max(P%lambda, 0.0d0))

  end subroutine sl_projector_load


  subroutine sl_projector_project(P, phi_re, phi_im, M, G, a_re, a_im)
    !!
    !! phi_re(i,m), phi_im(i,m) -> a_re(n,m), a_im(n,m)
    !!
    type(sl_projector_t), intent(in)  :: P
    real(dp),             intent(in)  :: phi_re(:,:), phi_im(:,:)
    type(mesh_t),         intent(in)  :: M
    type(geometry_t),     intent(in)  :: G
    real(dp),             intent(out) :: a_re(:,:), a_im(:,:)

    integer :: n, jm, nm_local

    nm_local = size(phi_re,2)

    do jm = 1, nm_local
      do n = 1, P%nspec
        a_re(n,jm) = sl_inner(P%R(:,n), phi_re(:,jm), M, G)
        a_im(n,jm) = sl_inner(P%R(:,n), phi_im(:,jm), M, G)
      end do
    end do

  end subroutine sl_projector_project


  subroutine sl_projector_power_phi(P, phi_re, phi_im, M, G, power_ln_local)
    !!
    !! P_ln^local = sum_{m local} |a_lmn|^2
    !!
    type(sl_projector_t), intent(in)  :: P
    real(dp),             intent(in)  :: phi_re(:,:), phi_im(:,:)
    type(mesh_t),         intent(in)  :: M
    type(geometry_t),     intent(in)  :: G
    real(dp),             intent(out) :: power_ln_local(:)

    real(dp), allocatable :: a_re(:,:), a_im(:,:)
    integer :: nm_local

    nm_local = size(phi_re,2)

    allocate(a_re(P%nspec,nm_local))
    allocate(a_im(P%nspec,nm_local))

    call sl_projector_project(P, phi_re, phi_im, M, G, a_re, a_im)

    power_ln_local = sum(a_re*a_re + a_im*a_im, dim=2)

    deallocate(a_re)
    deallocate(a_im)

  end subroutine sl_projector_power_phi


  subroutine sl_projector_free(P)
    type(sl_projector_t), intent(inout) :: P

    if (allocated(P%lambda)) deallocate(P%lambda)
    if (allocated(P%k))      deallocate(P%k)
    if (allocated(P%R))      deallocate(P%R)

    P%ell   = -1
    P%nr    = 0
    P%nspec = 0

  end subroutine sl_projector_free


  function sl_inner(f, h, M, G) result(value)
    !!
    !! <f,h>_active =
    !!
    !! ∫_{rmin}^{rmax} f(r) h(r) r^2/sqrt(gamma_rr) dr
    !!
    real(dp),         intent(in) :: f(:), h(:)
    type(mesh_t),     intent(in) :: M
    type(geometry_t), intent(in) :: G
    real(dp)                     :: value

    real(dp) :: integrand(size(f))

    integrand = f * h * M%r**2 / sqrt(G%grr)

    value = trapezium(integrand, M)

  end function sl_inner

end module sl_projection