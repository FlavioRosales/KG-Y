module geometry
  !! Geometría OOP: arreglos métricos (centros e interfaces) y sus derivadas.
  !! Unidades geométricas: G = c = 1  → Rs = 2*M_bh
  implicit none
  private
  public :: geometry_t

  type :: geometry_t
     ! Parámetros
     integer :: ell = 0          ! multipolo ℓ
     integer :: Nr  = 0          ! número de celdas / puntos
     real(kind=8) :: M_bh = 0.5d0  ! masa del BH (unidades geométricas)

     ! Centros (tamaño Nr)
     real(kind=8), allocatable :: alpha(:), beta(:), grr(:), guu(:), Vell_c(:)
     real(kind=8), allocatable :: dalpha(:), dbeta(:), dgrr(:), dguu(:)
     real(kind=8), allocatable :: K_c(:)

   contains
     procedure :: init => geom_init
     procedure :: build => geom_build
     procedure :: info => geom_info
     procedure :: free => geom_free
  end type geometry_t

contains

  subroutine geom_init(this, r, M_bh)
    !! Inicializa y construye geometría en centros e interfaces.
    class(geometry_t), intent(inout) :: this
    real(kind=8),      intent(in)    :: r(:)
    real(kind=8),      intent(in),   optional :: M_bh

    this%Nr = size(r)
    if (present(M_bh)) this%M_bh = M_bh

    call this%build(r)
  end subroutine geom_init

  subroutine geom_build(this, r)
    class(geometry_t), intent(inout) :: this
    real(kind=8),      intent(in)    :: r(:)

    integer :: i, n
    real(kind=8) :: rr, Rs, dr_loc

    n        = size(r)
    this%Nr  = n
    Rs       = 2.0d0 * this%M_bh           ! radio de Schwarzschild

    ! Paso local (asumimos malla uniforme)
    if (n > 1) then
       dr_loc = r(2) - r(1)
    else
       dr_loc = 0.0d0
    end if

    ! (Re)asigna arreglos de centros
    call safe_resize(this%alpha,  n)
    call safe_resize(this%beta,   n)
    call safe_resize(this%grr,    n)
    call safe_resize(this%guu,    n)
    call safe_resize(this%dalpha, n)
    call safe_resize(this%dbeta,  n)
    call safe_resize(this%dgrr,   n)
    call safe_resize(this%dguu,   n)
    call safe_resize(this%K_c,   n)

    ! -------- Centros --------
    do i = 1, n
      rr = r(i)

      ! Métrica en centros
      this%grr(i)    = 1.0d0 + Rs/rr                 ! gamma_rr
      this%alpha(i)  = (this%grr(i))**(-0.5d0)       ! alpha(r) = 1 / sqrt(1 + Rs/r)
      this%beta(i)   = Rs / (rr * this%grr(i))            ! beta^r
      this%guu(i)    = 1.0d0 / this%grr(i)           ! gamma^rr

      ! Derivadas analíticas
      this%dalpha(i) = 0.5d0 * this%alpha(i)**3 * Rs / (rr**2)
      this%dgrr(i)   = -Rs / rr**2
      this%dguu(i)   = -1.0d0 / (this%grr(i)**2) * this%dgrr(i)
      this%dbeta(i)  = -this%beta(i) * ( 1.0d0/rr + this%dgrr(i)/this%grr(i) )
      ! Curvatura extrínseca
      this%K_c(i) = ( this%dbeta(i)                                   &
              + 0.5d0 * this%beta(i) * this%dgrr(i) / this%grr(i) &
              + 2.0d0/rr * this%beta(i) ) / this%alpha(i)

    end do


  end subroutine geom_build

  subroutine geom_info(this)
    class(geometry_t), intent(in) :: this
    write(*,'(A)')       "---- Geometry info ----"
    write(*,'(A,I0)')    "ell = ", this%ell
    write(*,'(A,I0)')    "Nr  = ", this%Nr
    write(*,'(A,ES14.6)')"M_bh = ", this%M_bh
    write(*,'(A)')       "------------------------"
  end subroutine geom_info

  subroutine geom_free(this)
    class(geometry_t), intent(inout) :: this

    ! Centros
    call safe_dealloc(this%alpha)
    call safe_dealloc(this%beta)
    call safe_dealloc(this%grr)
    call safe_dealloc(this%guu)
    call safe_dealloc(this%dalpha)
    call safe_dealloc(this%dbeta)
    call safe_dealloc(this%dgrr)
    call safe_dealloc(this%dguu)

    this%Nr  = 0
    this%ell = 0
    this%M_bh = 0.5d0
  end subroutine geom_free

  ! ===== Helpers =====
  subroutine safe_resize(a, n)
    real(kind=8), allocatable, intent(inout) :: a(:)
    integer,                intent(in)       :: n
    if (allocated(a)) then
       if (size(a) /= n) then
          deallocate(a); allocate(a(n))
       end if
    else
       allocate(a(n))
    end if
  end subroutine safe_resize

  subroutine safe_dealloc(a)
    real(kind=8), allocatable, intent(inout) :: a(:)
    if (allocated(a)) deallocate(a)
  end subroutine safe_dealloc

end module geometry
