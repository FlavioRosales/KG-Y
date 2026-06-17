module sl_spectrum
  !! Coeficientes complejos C_n ~ CN(0, P(k))
  !!
  !! Mantiene la misma llamada:
  !!
  !!   C = Clmn_amp(lambda_n)
  !!
  !! Cambios:
  !!   - Inicializa random_seed una sola vez.
  !!   - Usa el rank MPI desde variables de entorno:
  !!       OMPI_COMM_WORLD_RANK, PMI_RANK, SLURM_PROCID
  !!   - Evita que distintos procesos MPI arranquen con la misma secuencia.
  !!   - Permite diagnóstico con:
  !!       export SL_RNG_DEBUG=1
  !!
  implicit none
  private

  real(kind=8), parameter :: PI = acos(-1.0d0)

  real(kind=8) :: sigma_k = 0.1d0

  logical :: rng_initialized = .false.
  logical :: debug_rng       = .false.

  integer :: rng_rank        = 0
  integer :: rng_call_count  = 0

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


  subroutine ensure_rng_initialized()
    integer :: nseed, i
    integer, allocatable :: seed(:)
    integer(kind=8) :: clock_count, clock_rate
    integer(kind=8) :: base, mixed
    character(len=64) :: env
    integer :: status

    if (rng_initialized) return

    rng_rank = 0

    call get_environment_variable("OMPI_COMM_WORLD_RANK", env, status=status)
    if (status == 0) read(env, *, iostat=status) rng_rank

    if (status /= 0) then
       call get_environment_variable("PMI_RANK", env, status=status)
       if (status == 0) read(env, *, iostat=status) rng_rank
    end if

    if (status /= 0) then
       call get_environment_variable("SLURM_PROCID", env, status=status)
       if (status == 0) read(env, *, iostat=status) rng_rank
    end if

    call get_environment_variable("SL_RNG_DEBUG", env, status=status)
    if (status == 0) then
       if (trim(env) == "1" .or. trim(env) == "true" .or. trim(env) == "TRUE") then
          debug_rng = .true.
       end if
    end if

    call random_seed(size=nseed)
    allocate(seed(nseed))

    call system_clock(count=clock_count, count_rate=clock_rate)

    base = clock_count                                      &
         + 104729_8   * int(rng_rank + 1, kind=8)           &
         + 15485863_8 * int(nseed + 1, kind=8)

    do i = 1, nseed

       mixed = base                                      &
             + 32452843_8 * int(i, kind=8)               &
             + 49979687_8 * int(rng_rank + 1, kind=8)

       mixed = modulo(mixed, int(huge(0), kind=8) - 1_8)

       if (mixed <= 0_8) mixed = mixed + 13579_8

       seed(i) = int(mixed)

    end do

    call random_seed(put=seed)

    deallocate(seed)

    rng_initialized = .true.

    if (debug_rng) then
       write(*,'(A,I8,A,I8,A,I8)') &
            "[sl_spectrum] RNG initialized. rank=", rng_rank, &
            " nseed=", nseed, " clock=", int(modulo(clock_count, 100000000_8))
    end if

  end subroutine ensure_rng_initialized


  subroutine randn2(z1, z2)
    !! Devuelve dos N(0,1) independientes con Box-Muller.
    real(kind=8), intent(out) :: z1, z2
    real(kind=8) :: u1, u2, s

    call ensure_rng_initialized()

    call random_number(u1)
    call random_number(u2)

    u1 = max(u1, 1.0d-16)

    s  = dsqrt(-2.0d0*dlog(u1))

    z1 = s*dcos(2.0d0*PI*u2)
    z2 = s*dsin(2.0d0*PI*u2)

  end subroutine randn2


  function Clmn_amp(lambda_n) result(C)
    real(kind=8), intent(in) :: lambda_n
    real(kind=8)             :: C(2)

    real(kind=8) :: k, std, z1, z2

    rng_call_count = rng_call_count + 1

    k   = dsqrt(max(lambda_n, 0.0d0))
    std = dsqrt(max(P_of_k(k), 0.0d0))

    call randn2(z1, z2)

    C(1) = std * z1 / dsqrt(2.0d0)
    C(2) = std * z2 / dsqrt(2.0d0)

    if (debug_rng .and. rng_call_count <= 12) then
       write(*,'(A,I6,A,I8,A,ES14.6,A,2ES16.8)') &
            "[sl_spectrum] rank=", rng_rank, &
            " call=", rng_call_count, &
            " lambda=", lambda_n, &
            " C=", C(1), C(2)
    end if

  end function Clmn_amp

end module sl_spectrum