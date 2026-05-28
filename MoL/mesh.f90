module mesh
  use iso_fortran_env, only: real64
  use run_control
  implicit none
  private
  public :: mesh_t, build_mesh, print_mesh_info

  type :: mesh_t
     integer :: Nr = 0
     real(real64) :: rmin = 0.0_real64, rmax = 0.0_real64, dr = 0.0_real64
     real(real64), allocatable :: r(:)
     real(real64) :: dr_min = 0.0_real64, dr_max = 0.0_real64
     integer :: i_flat = 0
     real(real64) :: r_flat = 0.0_real64
     integer :: i_sl_first = 0, i_sl_last = 0
  end type mesh_t

contains

  subroutine build_mesh(M)
    type(mesh_t), intent(inout) :: M
    integer :: i, n_intervals, nr_sl, i_first, i_last
    real(real64), allocatable :: r_sl(:)
    real(real64) :: r_outer, rmin_req, rmax_req
    real(real64), parameter :: eps = 1.0e-12_real64

    rmin_req = rmin_p
    rmax_req = rmax_p
    M%dr_min = dr_min_p
    M%dr_max = dr_max_p
    M%r_flat = r_flat_p

    call count_sl_grid(rmax_req, n_intervals, r_outer)
    nr_sl = n_intervals - 1
    if (nr_sl < 1) error stop "build_mesh: SL grid has no interior points"

    allocate(r_sl(0:n_intervals))
    call build_sl_grid(n_intervals, r_sl)

    i_first = 0
    do i = 1, nr_sl
      if (r_sl(i) >= rmin_req - eps) then
        i_first = i
        exit
      end if
    end do
    if (i_first == 0) error stop "build_mesh: rmin is outside the SL interior grid"

    i_last = nr_sl
    if (i_last < i_first) error stop "build_mesh: empty evolution grid after excision"

    if (allocated(M%r)) deallocate(M%r)
    M%Nr = i_last - i_first + 1
    allocate(M%r(M%Nr))
    M%r = r_sl(i_first:i_last)

    M%i_sl_first = i_first
    M%i_sl_last  = i_last
    M%rmin = M%r(1)
    M%rmax = M%r(M%Nr)

    M%i_flat = 0
    do i = 1, M%Nr
        if (M%r(i) >= M%r_flat) then
            M%i_flat = i + 1
            exit
        end if
    end do

    if (M%Nr > 1) then
      M%dr_min = minval(M%r(2:M%Nr) - M%r(1:M%Nr-1))
      M%dr_max = maxval(M%r(2:M%Nr) - M%r(1:M%Nr-1))
      M%dr = M%dr_min
    else
      M%dr = 0.0_real64
    end if

    deallocate(r_sl)

  end subroutine build_mesh

  subroutine print_mesh_info(M)
    type(mesh_t), intent(in) :: M
    write(*,'(A)')          "---- Mesh info ----"
    write(*,'(A,I0)')       "Nr       = ", M%Nr
    write(*,'(A,ES14.6)')   "rmin     = ", M%rmin
    write(*,'(A,ES14.6)')   "rmax     = ", M%rmax
    write(*,'(A,ES14.6)')   "dr       = ", M%dr
    write(*,'(A,ES14.6)')   "dr_min   = ", M%dr_min
    write(*,'(A,ES14.6)')   "dr_max   = ", M%dr_max
    write(*,'(A,I0)')       "i_flat   = ", M%i_flat
    write(*,'(A,I0)')       "i_sl_ini = ", M%i_sl_first
    write(*,'(A,I0)')       "i_sl_end = ", M%i_sl_last
    write(*,'(A,ES14.6)')   "r(1)     = ", M%r(1)
    write(*,'(A,ES14.6)')   "r(N)     = ", M%r(M%Nr)
    write(*,'(A,I0)')       "ell_max  = ", ell_max_p
    write(*,'(A)')          "-------------------"
  end subroutine print_mesh_info

  subroutine count_sl_grid(r_target, n_intervals, r_outer)
    real(real64), intent(in)  :: r_target
    integer,      intent(out) :: n_intervals
    real(real64), intent(out) :: r_outer

    real(real64) :: x, dx
    real(real64), parameter :: eps = 1.0e-12_real64

    x = 0.0_real64
    n_intervals = 0
    do while (x < r_target - eps)
      dx = next_dr(x, r_flat_p, dr_min_p, dr_max_p)
      if (dx <= 0.0_real64) error stop "count_sl_grid: dx <= 0"
      x = x + dx
      n_intervals = n_intervals + 1
    end do
    r_outer = x
  end subroutine count_sl_grid

  subroutine build_sl_grid(n_intervals, r)
    integer,      intent(in)  :: n_intervals
    real(real64), intent(out) :: r(0:n_intervals)

    integer :: i
    real(real64) :: dx

    r(0) = 0.0_real64
    do i = 0, n_intervals - 1
      dx = next_dr(r(i), r_flat_p, dr_min_p, dr_max_p)
      if (dx <= 0.0_real64) error stop "build_sl_grid: dx <= 0"
      r(i+1) = r(i) + dx
    end do
  end subroutine build_sl_grid

  real(kind=8) function local_dr(x, r_flat, dr_min, dr_max)
     real(kind=8), intent(in) :: x
     real(kind=8), intent(in) :: r_flat
     real(kind=8), intent(in) :: dr_min
     real(kind=8), intent(in) :: dr_max
     real(kind=8) :: xi

    if (x <= 1.0d0) then
        local_dr = dr_min
    else if (x >= r_flat) then
        local_dr = dr_max
    else
        xi = dlog(0.5d0*x*(x + 1.0d0)) / dlog(0.5d0*r_flat*(r_flat + 1.0d0))
        local_dr = dr_min + (dr_max - dr_min)*xi
    end if
end function local_dr

real(kind=8) function next_dr(x,r_flat,dr_min, dr_max)
    real(kind=8), intent(in) :: x
    real(kind=8), intent(in) :: r_flat
    real(kind=8), intent(in) :: dr_min
    real(kind=8), intent(in) :: dr_max
    real(kind=8), parameter  :: eps = 1.0d-12

    next_dr = local_dr(x, r_flat, dr_min, dr_max)

    ! Keep r=1 as an exact node. Do not force r_target/r_max.
    if (x < 1.0d0 - eps .and. x + next_dr > 1.0d0) then
        next_dr = 1.0d0 - x
    end if
end function next_dr

end module mesh
