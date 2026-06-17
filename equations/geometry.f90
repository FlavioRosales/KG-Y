module geometry
  !! Geometría OOP: arreglos métricos (centros e interfaces) y sus derivadas.
  !! Unidades geométricas: G = c = 1  → Rs = 2*M_bh
  use mesh, only: mesh_t
  implicit none
  public :: geometry_t

  type :: geometry_t
     ! Parámetros
     integer :: ell = 0          ! multipolo ℓ
     integer :: Nr  = 0          ! número de celdas / puntos
     real(kind=8) :: M_bh = 0.5d0  ! masa del BH (unidades geométricas)

     ! Centros (tamaño Nr)
     real(kind=8), allocatable :: alpha(:), beta(:), grr(:), guu(:)
     !real(kind=8), allocatable :: dalpha(:), dgrr(:), dguu(:)
     real(kind=8), allocatable :: K_c(:), dbeta(:)
     real(kind=8), allocatable :: r2_inv(:)
     real(kind=8), allocatable :: adv_a(:), adv_b(:), adv_c(:)
     real(kind=8), allocatable :: adv_u0(:), adv_u1(:), adv_u2(:)

   contains
     procedure :: init => geom_init
     procedure :: build => geom_build
     procedure :: build_advec => geom_build_advec
     procedure :: info => geom_info
     procedure :: free => geom_free
  end type geometry_t

contains

  subroutine geom_init(this,M_bh,Mesh)
    !! Inicializa y construye geometría en centros e interfaces.
    class(geometry_t), intent(inout) :: this
    real(kind=8),      intent(in),   optional :: M_bh
    type(mesh_t),      intent(in)    :: Mesh

    this%Nr = Mesh%Nr
    if (present(M_bh)) this%M_bh = M_bh

    call this%build(Mesh%r)
    call this%build_advec(Mesh)
  end subroutine geom_init

  subroutine geom_build(this, r)
    class(geometry_t), intent(inout) :: this
    real(kind=8),      intent(in)    :: r(:)

    integer :: i, n
    real(kind=8) :: rr, Rs
    real(kind=8) :: dgrr

    n        = size(r)
    this%Nr  = n
    Rs       = 2.0d0 * this%M_bh           ! radio de Schwarzschild



    allocate(this%alpha(n), this%beta(n), this%grr(n), this%guu(n), &
             this%dbeta(n), this%K_c(n) )

!    allocate(this%dalpha(n), this%dgrr(n), this%dguu(n), this%dbeta(n))

    allocate(this%adv_a(n), this%adv_b(n), this%adv_c(n), &
             this%adv_u0(n), this%adv_u1(n), this%adv_u2(n) )

    allocate(this%r2_inv(n))

     do concurrent (i = 1:n)
       this%r2_inv(i) = 1.0d0 / r(i)**2
     end do

    do i = 1, n
      rr = r(i)

      this%grr(i)    = 1.0d0 + Rs/rr                 ! gamma_rr
      this%alpha(i)  = (this%grr(i))**(-0.5d0)       ! alpha(r) = 1 / sqrt(1 + Rs/r)
      this%beta(i)   = Rs / (rr * this%grr(i))            ! beta^r
      this%guu(i)    = 1.0d0 / this%grr(i)           ! gamma^rr

      !this%dalpha(i) = 0.5d0 * this%alpha(i)**3 * Rs / (rr**2)
      dgrr   = -Rs / rr**2
      !this%dguu(i)   = -1.0d0 / (this%grr(i)**2) * this%dgrr(i)
      this%dbeta(i)  = -this%beta(i) * ( 1.0d0/rr + dgrr/this%grr(i) )


      this%K_c(i) = ( this%dbeta(i)                                   &
              + 0.5d0 * this%beta(i) * dgrr / this%grr(i) &
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

    call safe_dealloc(this%alpha)
    call safe_dealloc(this%beta)
    call safe_dealloc(this%grr)
    call safe_dealloc(this%guu)
    !call safe_dealloc(this%dalpha)
    call safe_dealloc(this%dbeta)
    !call safe_dealloc(this%dgrr)
    !call safe_dealloc(this%dguu)
    call safe_dealloc(this%K_c)
    call safe_dealloc(this%adv_a)
    call safe_dealloc(this%adv_b)
    call safe_dealloc(this%adv_c)
    call safe_dealloc(this%adv_u0)
    call safe_dealloc(this%adv_u1)
    call safe_dealloc(this%adv_u2)
    call safe_dealloc(this%r2_inv)

    this%Nr  = 0
    this%ell = 0
    this%M_bh = 0.5d0
  end subroutine geom_free

  subroutine safe_dealloc(a)
    real(kind=8), allocatable, intent(inout) :: a(:)
    if (allocated(a)) deallocate(a)
  end subroutine safe_dealloc

  subroutine geom_build_advec(this, M)
  implicit none

  class(geometry_t), intent(inout) :: this
  type(mesh_t),      intent(in)    :: M

  integer      :: i, n, i_flat
  real(kind=8) :: h0, h1, denom, res
  real(kind=8) :: a, b, c

  n      = M%Nr
  i_flat = M%i_flat
  res    = M%inv_2dr_max

  !-----------------------------------------------------------------------
  ! Región no uniforme: beta(i) * stencil forward no uniforme
  !-----------------------------------------------------------------------
  do i = 1, i_flat - 1
     h0 = M%r(i+1) - M%r(i)
     h1 = M%r(i+2) - M%r(i+1)

     denom = h0*h1*(h0 + h1)

     a = -(2.0d0*h0 + h1)*h1 / denom
     b =  (h0 + h1)*(h0 + h1) / denom
     c = -h0*h0 / denom

     this%adv_a(i) = this%beta(i) * a
     this%adv_b(i) = this%beta(i) * b
     this%adv_c(i) = this%beta(i) * c
  end do

  !-----------------------------------------------------------------------
  ! Región uniforme: beta(i) * forward de segundo orden
  !
  ! d f_i = (-3 f_i + 4 f_{i+1} - f_{i+2})/(2 dr)
  !-----------------------------------------------------------------------
  do i = i_flat, n - 2
     this%adv_u0(i) = -3.0d0 * this%beta(i) * res
     this%adv_u1(i) =  4.0d0 * this%beta(i) * res
     this%adv_u2(i) = -1.0d0 * this%beta(i) * res
  end do

  !-----------------------------------------------------------------------
  ! Penúltimo punto: beta(i) * centrado
  !
  ! d f_{n-1} = (f_n - f_{n-2})/(2 dr)
  !-----------------------------------------------------------------------
  this%adv_u0(n-1) = -this%beta(n-1) * res
  this%adv_u1(n-1) =  0.0d0
  this%adv_u2(n-1) =  this%beta(n-1) * res

  !-----------------------------------------------------------------------
  ! Último punto: beta(n) * backward de segundo orden
  !
  ! d f_n = (3 f_n - 4 f_{n-1} + f_{n-2})/(2 dr)
  !-----------------------------------------------------------------------
  this%adv_u0(n) =  this%beta(n) * res
  this%adv_u1(n) = -4.0d0 * this%beta(n) * res
  this%adv_u2(n) =  3.0d0 * this%beta(n) * res

end subroutine geom_build_advec

end module geometry
