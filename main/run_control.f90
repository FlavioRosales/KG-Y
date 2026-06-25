module run_control
  implicit none
  private
  public :: sim_setup, sim_finalize
  public :: rmin_p, rmax_p, nr_p, t_end_p, cfl_p, nmax_p, ell_min_p, ell_max_p
  public :: dt_out_r_p
  public :: out_t_p
  public :: every_0D_p, every_1D_p
  public :: dt_out_0D_p, dt_out_1D_p
  public :: save_planes_p, every_planes_p, dt_out_planes_p
  public :: r_plane_max_p, nphi_plane_p, ntheta_plane_p
  public :: p0_p, k0_p
  public :: M_bh_p, mu_p
  public :: only_m0_p
  public :: initial_conditions_p
  public :: r_flat_p, dr_min_p, dr_max_p

  ! -------- Defaults --------
  real(kind=8) :: rmin_p = 1.5d0
  real(kind=8) :: rmax_p = 120.0d0
  real(kind=8) :: dr_min_p = 0.1d0
  real(kind=8) :: dr_max_p = 1.0d0
  real(kind=8) :: r_flat_p = 1000.0d0
  integer      :: nr_p   = 256
  integer      :: every_0D_p = 1
  integer      :: every_1D_p = 1
  real(kind=8) :: dt_out_0D_p = -1.0d0
  real(kind=8) :: dt_out_1D_p = -1.0d0

  logical      :: save_planes_p   = .false.
  integer      :: every_planes_p  = 100
  real(kind=8) :: dt_out_planes_p = -1.0d0
  real(kind=8) :: r_plane_max_p   = 150.0d0
  integer      :: nphi_plane_p    = 0
  integer      :: ntheta_plane_p  = 0
  integer      :: ell_min_p   = 0
  integer      :: ell_max_p   = 1

  real(kind=8) :: t_end_p = 1.0d0
  real(kind=8) :: cfl_p   = 0.5d0
  real(kind=8) :: p0_p   = 1.0d0
  real(kind=8) :: k0_p   = 0.0d0
  integer      :: nmax_p  = 1000

  real(kind=8) :: M_bh_p  = 1.0d0
  real(kind=8) :: mu_p  = 0.0d0
  real(kind=8) :: dt_out_r_p = 1.0d0
  integer :: out_t_p = 2

  logical :: only_m0_p = .false.
  character(len=32) :: initial_conditions_p = "gaussian_positive_charge"


contains

  subroutine sim_setup()
    integer :: u, ios
    logical :: fexist
    namelist /grid/    rmin_p, rmax_p, nr_p, ell_min_p, ell_max_p, dr_min_p, dr_max_p, r_flat_p
    namelist /physics/ M_bh_p, mu_p, p0_p, k0_p, initial_conditions_p
    namelist /run/     t_end_p, cfl_p, nmax_p, dt_out_r_p, out_t_p, every_0D_p, every_1D_p, &
                        dt_out_0D_p, dt_out_1D_p, only_m0_p, save_planes_p, every_planes_p, &
                        dt_out_planes_p, r_plane_max_p, nphi_plane_p, ntheta_plane_p

    inquire(file="params.nml", exist=fexist)
    if (.not. fexist) then
      call sanity_check()
      return
    end if

    ! --- Intento 1: NAMELIST estándar ---
    open(newunit=u, file="params.nml", status="old", action="read", iostat=ios)
    if (ios == 0) then
      ios = 0; rewind(u); read(u, nml=grid,    iostat=ios)
      if (ios == 0) then
        ios = 0; rewind(u); read(u, nml=physics, iostat=ios)
        ios = 0; rewind(u); read(u, nml=run,     iostat=ios)
        close(u)
        call post_process()
        call sanity_check()
        return
      end if
      close(u)
    end if

    ! --- Intento 2: Parser manual tolerante ---
    call parse_params_manual("params.nml")
    call post_process()
    call sanity_check()
  end subroutine

  subroutine sim_finalize()
  end subroutine



  !==============================!
  !===  Parser manual simple  ===!
  !==============================!
  subroutine parse_params_manual(path)
    character(len=*), intent(in) :: path
    integer :: u, ios
    character(len=:), allocatable :: line, sec, key, val
    character(len=512) :: buf
    logical :: fexist

    inquire(file=path, exist=fexist)
    if (.not. fexist) then
      write(*,*) "[run_control] No existe ", trim(path), " (usando defaults)"
      return
    end if

    open(newunit=u, file=path, status="old", action="read", iostat=ios, encoding="UTF-8")
    if (ios /= 0) then
      write(*,*) "[run_control] No se pudo abrir params.nml (usando defaults)."
      return
    end if

    sec = ""; ios = 0
    do
      read(u,'(A)', iostat=ios) buf
      if (ios /= 0) exit
      call strip_bom_inplace(buf)
      line = to_lower(trim(buf))
      call strip_inline_comment(line)
      if (len_trim(line) == 0) cycle

      if (line(1:1) == "&") then
        sec = trim(line(2:))
        cycle
      else if (line == "/") then
        sec = ""
        cycle
      end if

      if (len_trim(sec) == 0) cycle
      call split_kv(line, key, val)
      if (len_trim(key) == 0) cycle

      select case (trim(sec))
      case ("grid")
        select case (trim(key))
        case ("rmin", "rmin_p");       call parse_real(val, rmin_p)
        case ("rmax", "rmax_p");       call parse_real(val, rmax_p)
        case ("nr", "nr_p");           call parse_int (val, nr_p)
        case ("ell_min", "ell_min_p"); call parse_int (val, ell_min_p)
        case ("ell_max", "ell_max_p"); call parse_int (val, ell_max_p)
        case ("dr_min", "dr_min_p");   call parse_real(val, dr_min_p)
        case ("dr_max", "dr_max_p");   call parse_real(val, dr_max_p)
        case ("r_flat", "r_flat_p");   call parse_real(val, r_flat_p)
        end select
      case ("physics")
        select case (trim(key))
        case ("m_bh", "m_bh_p"); call parse_real(val, M_bh_p)
        case ("p0", "p0_p");     call parse_real(val, p0_p)
        case ("k0", "k0_p");     call parse_real(val, k0_p)
        case ("mu", "mu_p");     call parse_real(val, mu_p)
        case ("initial_conditions", "initial_conditions_p"); call parse_str(val, initial_conditions_p)
        end select
      case ("run")
        select case (trim(key))
        case ("t_end", "t_end_p");       call parse_real(val, t_end_p)
        case ("cfl", "cfl_p");           call parse_real(val, cfl_p)
        case ("nmax", "nmax_p");         call parse_int (val, nmax_p)
        case ("dt_out_r", "dt_out_r_p"); call parse_real(val, dt_out_r_p)
        case ("out_t", "out_t_p");       call parse_int (val, out_t_p)
        case ("every_0d", "every_0d_p"); call parse_int (val, every_0D_p)
        case ("every_1d", "every_1d_p");     call parse_int (val, every_1D_p)
        case ("dt_out_0d", "dt_out_0d_p");     call parse_real(val, dt_out_0D_p)
        case ("dt_out_1d", "dt_out_1d_p");     call parse_real(val, dt_out_1D_p)
        case ("save_planes", "save_planes_p"); call parse_logical(val, save_planes_p)
        case ("every_planes", "every_planes_p"); call parse_int(val, every_planes_p)
        case ("dt_out_planes", "dt_out_planes_p"); call parse_real(val, dt_out_planes_p)
        case ("r_plane_max", "r_plane_max_p"); call parse_real(val, r_plane_max_p)
        case ("nphi_plane", "nphi_plane_p");   call parse_int(val, nphi_plane_p)
        case ("ntheta_plane", "ntheta_plane_p"); call parse_int(val, ntheta_plane_p)
        case ("only_m0", "only_m0_p");         call parse_logical(val, only_m0_p)
        end select
      end select
    end do
    close(u)
  end subroutine

  !--- helpers de parsing ---!
  pure function to_lower(s) result(out)
    character(len=*), intent(in) :: s
    character(len=len(s))        :: out
    integer :: i, ia
    do i=1,len(s)
      ia = iachar(s(i:i))
      if (ia >= iachar('A') .and. ia <= iachar('Z')) then
        out(i:i) = achar(ia + 32)
      else
        out(i:i) = s(i:i)
      end if
    end do
  end function

  subroutine strip_inline_comment(s)
    character(len=:), allocatable, intent(inout) :: s
    integer :: p1
    p1 = index(s, "!")
    if (p1==0) p1 = index(s, "#")
    if (p1>0) s = s(:p1-1)
    s = trim(s)
  end subroutine

  subroutine split_kv(line, key, val)
    character(len=*), intent(in)  :: line
    character(len=:), allocatable, intent(out) :: key, val
    integer :: peq
    peq = index(line, "=")
    if (peq == 0) then
      key = ""; val = ""; return
    end if
    key = trim(adjustl(line(:peq-1)))
    val = trim(adjustl(line(peq+1:)))
  end subroutine

  subroutine parse_real(s, x)
    character(len=*), intent(in) :: s
    real(kind=8),     intent(out):: x
    character(len=:), allocatable :: t
    integer :: ios
    t = dequote(s)
    read(t,*,iostat=ios) x
    if (ios /= 0) then
      ! deja default si falla
    end if
  end subroutine

  subroutine parse_logical(s, x)
    character(len=*), intent(in) :: s
    logical,         intent(out):: x
    character(len=:), allocatable :: t
    t = to_lower(dequote(s))
    select case (t)
    case ("true", ".true.", "t", "1")
      x = .true.
    case ("false", ".false.", "f", "0")
      x = .false.
    end select
  end subroutine

  subroutine parse_int(s, x)
    character(len=*), intent(in) :: s
    integer,          intent(out):: x
    character(len=:), allocatable :: t
    integer :: ios
    t = dequote(s)
    read(t,*,iostat=ios) x
    if (ios /= 0) then
      ! deja default si falla
    end if
  end subroutine

  subroutine parse_str(s, x)
    character(len=*), intent(in) :: s
    character(len=*), intent(out):: x
    character(len=:), allocatable :: t
    t = trim(adjustl(dequote(s)))
    x = t(1:min(len(x),len(t)))
  end subroutine

  pure function dequote(s) result(t)
    character(len=*), intent(in) :: s
    character(len=len_trim(s))   :: t
    integer :: L
    t = trim(s)
    L = len_trim(t)
    if (L >= 2) then
      if ( (t(1:1) == '"' .and. t(L:L) == '"') .or. (t(1:1) == "'" .and. t(L:L) == "'") ) then
        t = t(2:L-1)
      end if
    end if
  end function

  subroutine strip_bom_inplace(buf)
    character(len=*), intent(inout) :: buf
    if (len_trim(buf) >= 3) then
      if (iachar(buf(1:1))==239 .and. iachar(buf(2:2))==187 .and. iachar(buf(3:3))==191) then
        buf = buf(4:)
      end if
    end if
  end subroutine

  subroutine post_process()
    initial_conditions_p = trim(adjustl(initial_conditions_p))
  end subroutine post_process

  subroutine sanity_check()
    if (rmin_p < 0.0d0) error stop "run_control: rmin must be non-negative"
    if (rmax_p <= rmin_p) error stop "run_control: rmax must be greater than rmin"
    if (dr_min_p <= 0.0d0) error stop "run_control: dr_min must be positive"
    if (dr_max_p <= 0.0d0) error stop "run_control: dr_max must be positive"
    if (dr_max_p < dr_min_p) error stop "run_control: dr_max must be >= dr_min"
    if (r_flat_p <= 0.0d0) error stop "run_control: r_flat must be positive"
    if (ell_min_p < 0) error stop "run_control: ell_min must be non-negative"
    if (ell_max_p < ell_min_p) error stop "run_control: ell_max must be >= ell_min"
    if (ell_min_p /= 0 .and. .not. only_m0_p) &
      error stop "run_control: ell_min /= 0 currently requires only_m0 = .true."
    if (t_end_p <= 0.0d0) error stop "run_control: t_end must be positive"
    if (cfl_p <= 0.0d0) error stop "run_control: cfl must be positive"
  end subroutine sanity_check

end module