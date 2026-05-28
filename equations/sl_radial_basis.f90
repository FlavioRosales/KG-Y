module sl_radial_basis
  !! Base radial de Sturm–Liouville para cada ell.
  !!
  !! Operador continuo:
  !!   P(r) = r^2 / sqrt(gamma_rr)
  !!   q(r) = ell(ell+1) * sqrt(gamma_rr)
  !!   w(r) = r^2 * sqrt(gamma_rr)
  !!
  !! Discretización 2ª orden + condiciones Robin en rmin y rmax.
  !!
  use iso_fortran_env, only: real64
  implicit none
  private

  type :: radial_modes_t
    integer :: Nm = 0
    real(real64), allocatable :: lambda(:) ! eigenvalores seleccionados
    real(real64), allocatable :: R(:,:)    ! modos radiales (Nr, Nm)
  end type radial_modes_t

  public :: radial_modes_t
  public :: build_radial_modes_for_ell
  public :: check_radial_modes   ! <-- nuevo: checos

contains

  pure real(real64) function sigma_in(rmin) result(s)
    real(real64), intent(in) :: rmin
    s = -1.0_real64 / max(rmin, 1.0e-14_real64)
  end function sigma_in

  pure real(real64) function sigma_out(rmax) result(s)
    real(real64), intent(in) :: rmax
    s = -1.0_real64 / max(rmax, 1.0e-14_real64)
  end function sigma_out

  !------------------------------------------------------------
  ! Criterio alias-free:
  !
  ! lambda <= min_r [ ell(ell+1)/r^2 + (pi/dr)^2 safety^2 / gamma_rr(r) ] * margin
  !------------------------------------------------------------
  subroutine select_modes_l_gamma(lambdas, ell, r, grr, dr, safety, margin, &
                                  idx_keep, Nm_keep, lam_max)
    real(real64), intent(in)  :: lambdas(:)
    integer,      intent(in)  :: ell
    real(real64), intent(in)  :: r(:), grr(:)
    real(real64), intent(in)  :: dr, safety, margin
    integer,      allocatable, intent(out) :: idx_keep(:)
    integer,      intent(out) :: Nm_keep
    real(real64), intent(out) :: lam_max

    integer      :: Nr, j, count_keep
    real(real64) :: kappa_max_sq, rhs_min, rhs, ellfac

    Nr = size(r)
    if (size(grr) /= Nr) then
      Nm_keep = 0
      lam_max = 0.0_real64
      allocate(idx_keep(0))
      return
    end if

    ellfac       = real(ell,real64)*real(ell+1,real64)
    kappa_max_sq = (acos(-1.0_real64)/dr)**2 * safety * safety
    rhs_min      = huge(1.0_real64)

    do j = 1, Nr
      rhs = ellfac / max(r(j)*r(j), 1.0e-30_real64) + &
            kappa_max_sq / max(grr(j), 1.0e-30_real64)
      if (rhs < rhs_min) rhs_min = rhs
    end do

    lam_max = rhs_min * margin

    count_keep = 0
    do j = 1, size(lambdas)
      if (lambdas(j) <= lam_max) count_keep = count_keep + 1
    end do

    Nm_keep = count_keep
    if (Nm_keep == 0) then
      allocate(idx_keep(0))
      return
    end if

    allocate(idx_keep(Nm_keep))
    count_keep = 0
    do j = 1, size(lambdas)
      if (lambdas(j) <= lam_max) then
        count_keep = count_keep + 1
        idx_keep(count_keep) = j
      end if
    end do

  end subroutine select_modes_l_gamma

  !------------------------------------------------------------
  ! Construye matrices L, W y el peso w_sl(r) para un ell dado.
  ! L y W deben venir ya con tamaño (Nr,Nr).
  !------------------------------------------------------------
  subroutine build_LW_for_ell(ell, r, grr, L, W, w_sl, dr)
    integer,      intent(in)  :: ell
    real(real64), intent(in)  :: r(:), grr(:)
    real(real64), intent(out) :: L(:,:), W(:,:)
    real(real64), allocatable, intent(out) :: w_sl(:)
    real(real64), intent(out) :: dr

    integer :: Nr, i
    real(real64), allocatable :: sA(:), p(:), q(:), pface(:)
    real(real64) :: inv_dr2, sigL, sigR, Fleft, ellfac

    Nr = size(r)
    if (Nr < 3) stop 'build_LW_for_ell: Nr < 3'
    if (size(grr) /= Nr) stop 'build_LW_for_ell: size(grr) /= size(r)'

    dr       = r(2) - r(1)
    inv_dr2  = 1.0_real64 / (dr*dr)
    ellfac   = real(ell,real64)*real(ell+1,real64)

    allocate(sA(Nr))
    allocate(p(Nr))
    allocate(q(Nr))
    allocate(pface(0:Nr))
    allocate(w_sl(Nr))

    do i = 1, Nr
      if (grr(i) <= 0.0_real64) stop 'build_LW_for_ell: grr <= 0'
      sA(i)   = sqrt(grr(i))
      p(i)    = r(i)*r(i) / sA(i)
      q(i)    = ellfac * sA(i)
      w_sl(i) = r(i)*r(i) * sA(i)
    end do

    ! p en caras: promedio armónico
    do i = 1, Nr-1
      pface(i) = 2.0_real64 * p(i)*p(i+1) / (p(i) + p(i+1) + 1.0e-300_real64)
    end do
    pface(0)  = p(1)
    pface(Nr) = p(Nr)

    ! Inicializar matrices
    L = 0.0_real64
    W = 0.0_real64

    ! Interior: 2..Nr-1
    do i = 2, Nr-1
      L(i,i-1) = -pface(i-1)*inv_dr2
      L(i,i)   = (pface(i-1)+pface(i))*inv_dr2 + q(i)
      L(i,i+1) = -pface(i)*inv_dr2
    end do

    ! BC Robin en rmin
    sigL  = sigma_in(r(1))
    Fleft = 2.0_real64 * pface(0) * sigL
    L(1,1) = pface(1)*inv_dr2 + q(1) + Fleft/dr
    L(1,2) = -pface(1)*inv_dr2

    ! BC Robin en rmax
    sigR       = sigma_out(r(Nr))
    L(Nr,Nr-1) = -pface(Nr-1)*inv_dr2
    L(Nr,Nr)   = pface(Nr-1)*inv_dr2 + q(Nr) + (2.0_real64*pface(Nr)*sigR)/dr

    ! W = diag(w_sl * dr)
    do i = 1, Nr
      W(i,i) = w_sl(i) * dr
    end do

    deallocate(sA)
    deallocate(p)
    deallocate(q)
    deallocate(pface)
  end subroutine build_LW_for_ell

  !------------------------------------------------------------
  ! Resuelve L v = lambda W v y construye modos_l para un ell.
  !------------------------------------------------------------
  subroutine build_radial_modes_for_ell(ell, r, grr, safety, margin, &
                                        nmodes_cap, modes)
    integer,      intent(in)  :: ell
    real(real64), intent(in)  :: r(:), grr(:)
    real(real64), intent(in)  :: safety, margin
    integer,      intent(in)  :: nmodes_cap   ! <=0 => sin tope
    type(radial_modes_t), intent(out) :: modes

    integer :: Nr, info, Nm_keep
    real(real64), allocatable :: L(:,:), W(:,:), eig(:), work(:), w_sl(:)
    integer,      allocatable :: idx_keep(:)
    real(real64) :: dr, lam_max, norm2
    integer :: i, j, lwork
    real(real64), parameter :: tiny_n = 1.0e-30_real64

    external :: dsygv

    Nr = size(r)
    if (Nr < 3) then
      modes%Nm = 0
      return
    end if

    allocate(L(Nr,Nr))
    allocate(W(Nr,Nr))
    allocate(eig(Nr))

    call build_LW_for_ell(ell, r, grr, L, W, w_sl, dr)

    ! Workspace query
    lwork = -1
    allocate(work(1))
    call dsygv(1, 'V', 'U', Nr, L, Nr, W, Nr, eig, work, lwork, info)
    if (info /= 0) then
      modes%Nm = 0
      deallocate(L, W, eig, w_sl, work)
      return
    end if

    lwork = int(work(1))
    deallocate(work)
    allocate(work(lwork))

    ! Solve
    call dsygv(1, 'V', 'U', Nr, L, Nr, W, Nr, eig, work, lwork, info)
    deallocate(work)
    if (info /= 0) then
      modes%Nm = 0
      deallocate(L, W, eig, w_sl)
      return
    end if

    ! Selección alias-free
    call select_modes_l_gamma(eig, ell, r, grr, dr, safety, margin, &
                              idx_keep, Nm_keep, lam_max)

    if (Nm_keep <= 0) then
      modes%Nm = 0
      deallocate(L, W, eig, w_sl, idx_keep)
      return
    end if

    if (nmodes_cap > 0 .and. Nm_keep > nmodes_cap) Nm_keep = nmodes_cap

    allocate(modes%lambda(Nm_keep))
    allocate(modes%R(Nr, Nm_keep))

    do j = 1, Nm_keep
      modes%lambda(j) = eig(idx_keep(j))
      modes%R(:,j)    = L(:, idx_keep(j))   ! columnas de L son eigenvectores
    end do
    modes%Nm = Nm_keep

    ! Normalización con peso w_sl
    do j = 1, Nm_keep
      norm2 = 0.0_real64
      do i = 1, Nr
        norm2 = norm2 + modes%R(i,j)**2 * w_sl(i) * dr
      end do
      if (norm2 > tiny_n) then
        modes%R(:,j) = modes%R(:,j) / sqrt(norm2)
      else
        modes%R(:,j) = 0.0_real64
      end if
    end do

    deallocate(L)
    deallocate(W)
    deallocate(eig)
    deallocate(w_sl)
    deallocate(idx_keep)

  end subroutine build_radial_modes_for_ell

  !------------------------------------------------------------
  ! CHEQUEO: ortogonalidad + residuo del eigenproblema
  !------------------------------------------------------------
  subroutine check_radial_modes(ell, r, grr, modes, tol_orth, tol_res, verbose)
    integer,      intent(in) :: ell
    real(real64), intent(in) :: r(:), grr(:)
    type(radial_modes_t), intent(in) :: modes
    real(real64), intent(in) :: tol_orth   ! tolerancia ortogonalidad (ej. 1.0e-6)
    real(real64), intent(in) :: tol_res    ! tolerancia residuo (ej. 1.0e-6)
    logical,      intent(in) :: verbose
    real(real64) :: LRn, WRn


    integer :: Nr, Nm
    real(real64), allocatable :: L(:,:), W(:,:), w_sl(:)
    real(real64) :: dr
    integer :: i, j, n
    real(real64) :: inner, max_offdiag, min_diag, max_diag
    real(real64) :: num, den, res_n, max_res

    if (modes%Nm <= 0) then
      if (verbose) then
        write(*,'(A,I3,A)') 'check_radial_modes: ell=', ell, ' sin modos (Nm=0)'
      end if
      return
    end if

    Nm = modes%Nm
    Nr = size(r)

    if (.not. allocated(modes%R)) then
      if (verbose) write(*,*) 'check_radial_modes: modes%R no asignado'
      return
    end if
    if (.not. allocated(modes%lambda)) then
      if (verbose) write(*,*) 'check_radial_modes: modes%lambda no asignado'
      return
    end if

    if (size(modes%R,1) /= Nr .or. size(modes%R,2) /= Nm) then
      if (verbose) write(*,*) 'check_radial_modes: tamaños inconsistenes en R'
      return
    end if

    ! Reconstruir L, W y w_sl para este ell
    allocate(L(Nr,Nr))
    allocate(W(Nr,Nr))
    call build_LW_for_ell(ell, r, grr, L, W, w_sl, dr)

    ! ===== Chequeo de ortogonalidad (R^T W R ~ I) =====
    max_offdiag = 0.0_real64
    min_diag    = huge(1.0_real64)
    max_diag    = -huge(1.0_real64)

    do i = 1, Nm
      do j = 1, Nm
        inner = 0.0_real64
        do n = 1, Nr
          inner = inner + modes%R(n,i) * modes%R(n,j) * w_sl(n) * dr
        end do
        if (i == j) then
          if (inner < min_diag) min_diag = inner
          if (inner > max_diag) max_diag = inner
        else
          if (abs(inner) > max_offdiag) max_offdiag = abs(inner)
        end if
      end do
    end do

    ! ===== Chequeo de residuo ||L R_n - λ_n W R_n|| / ||R_n|| =====
    max_res = 0.0_real64
    do j = 1, Nm
      num = 0.0_real64
      den = 0.0_real64
      do n = 1, Nr
        ! (L R)_n
        LRn = 0.0_real64
        WRn = 0.0_real64
        ! producto fila n de L con columna j
        LRn = sum( L(n,1:Nr) * modes%R(1:Nr,j) )
        ! producto fila n de W con columna j (W es diagonal)
        WRn = W(n,n) * modes%R(n,j)

        num = num + ( LRn - modes%lambda(j)*WRn )**2
        den = den + modes%R(n,j)**2
      end do
      if (den > 0.0_real64) then
        res_n = sqrt(num/den)
        if (res_n > max_res) max_res = res_n
      end if
    end do

    if (verbose) then
      write(*,'(A,I3)')    'check_radial_modes: ell = ', ell
      write(*,'(A,I6)')    '  Nm            = ', Nm
      write(*,'(A,1P,E10.3)') '  diag min(I)   ~ ', min_diag
      write(*,'(A,1P,E10.3)') '  diag max(I)   ~ ', max_diag
      write(*,'(A,1P,E10.3)') '  max offdiag   = ', max_offdiag
      write(*,'(A,1P,E10.3)') '  max residuo   = ', max_res
      write(*,'(A,1P,E10.3)') '  tol_orth      = ', tol_orth
      write(*,'(A,1P,E10.3)') '  tol_res       = ', tol_res
    end if

    if (max_offdiag > tol_orth) then
      write(*,'(A,I3,A,1P,E10.3)') &
        'WARNING: ell=', ell, ' ortogonalidad deficiente, max_offdiag=', max_offdiag
    end if

    if (abs(min_diag-1.0_real64) > tol_orth .or. abs(max_diag-1.0_real64) > tol_orth) then
      write(*,'(A,I3,A,1P,E10.3,1X,E10.3)') &
        'WARNING: ell=', ell, ' normas diag alejadas de 1: [min,max]=', min_diag, max_diag
    end if

    if (max_res > tol_res) then
      write(*,'(A,I3,A,1P,E10.3)') &
        'WARNING: ell=', ell, ' residuo grande en eigenpares, max_res=', max_res
    end if

    deallocate(L)
    deallocate(W)
    deallocate(w_sl)

  end subroutine check_radial_modes

end module sl_radial_basis
