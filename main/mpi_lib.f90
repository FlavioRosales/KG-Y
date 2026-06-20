module kg_mpi
  use mpi
  implicit none
    !
  integer :: ierr
  !
  ! rank of the calling process in the group of comm
  !
  integer :: rank
  !
  ! number of processes in the group of comm
  !
  integer :: nproc
  !
  integer :: master = 0

contains

  subroutine MPI_CREATE
    implicit none
    logical :: inited
    call MPI_Initialized(inited, ierr)
    if (.not. inited) call MPI_Init(ierr)    ! <-- INICIALIZA MPI AQUÍ
    call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank,  ierr)
  end subroutine MPI_CREATE


  pure function modes_total(lmax) result(n_total)
    integer,        intent(in) :: lmax
    integer             :: n_total
    integer :: S
    S       = lmax + 1
    n_total = S**2
  end function

  pure subroutine lm_of_n(n, l, m)
    integer,        intent(in)  :: n
    integer,        intent(out) :: l, m
    integer                     :: ll, r
    ll = int(floor( dsqrt( dble(n) )))
    l  = int(ll)
    r  = n - ll**2
    m  = int(r) - l
  end subroutine

  pure function n_of_lm(l, m) result(n)
    integer,        intent(in) :: l, m
    integer                    :: n
    n = l**2 + m + l
  end function

  pure subroutine n_range_for_rank(rank, nprocs, n_total, n_first, n_last)
    implicit none

    integer,        intent(in)  :: rank, nprocs
    integer,        intent(in)  :: n_total
    integer,        intent(out) :: n_first, n_last

    integer                     :: base, extra
    integer                     :: r, p

    r = int(rank)
    p = int(nprocs)

    ! Aquí n_total debe ser Nell = ell_max + 1 (ell=0..ell_max)
    if (n_total <= 0 .or. p <= 0) then
      n_first = 0
      n_last  = -1
      return
    end if

    base  = n_total / p
    extra = mod(n_total, p)

    if (r < extra) then
      n_first = r * (base + 1)
      n_last  = n_first + (base + 1) - 1
    else
      n_first = extra * (base + 1) + (r - extra) * base
      n_last  = n_first + base - 1
    end if
  end subroutine n_range_for_rank



  subroutine n_iter_cyclic(rank, nprocs, n_total, n, has_more)
    integer,        intent(in)    :: rank, nprocs
    integer,        intent(in)    :: n_total
    integer,        intent(inout) :: n
    logical,        intent(out)   :: has_more
    if (n < 0) then
      n = int(rank)
    else
      n = n + int(nprocs)
    end if
    has_more = (n < n_total)
  end subroutine

  subroutine n_iter_block_cyclic(rank, nprocs, n_total, B, n0, n1, has_more)
    integer,        intent(in)    :: rank, nprocs, B
    integer,        intent(in)    :: n_total
    integer,        intent(inout) :: n0, n1
    logical,        intent(out)   :: has_more
    integer                       :: stride, start
    stride = nprocs*B
    if (n0 < 0) then
      start = rank*B; n0 = start
    else
      n0 = n0 + stride
    end if
    if (n0 >= n_total) then
      has_more = .false.; n1 = -1; return
    end if
    n1 = min(n0 + B-1, n_total-1); has_more = .true.
  end subroutine

  pure function format_mode_dir(base, l, m, rank) result(path)
    character(len=*), intent(in) :: base
    integer,          intent(in) :: l, m, rank
    character(len=256)           :: path
    write(path,'(A,"/l",I2.2,"_m",SP,I3.3,"/r",I0)') trim(base), l, m, rank
  end function

end module kg_mpi
