module state
    use geometry,   only: geometry_t        ! <-- módulo que define TYPE :: geometry
    use mesh,       only: mesh_t          ! <-- módulo que define TYPE :: mesh_t

  implicit none
  private
  public :: state_t

  real(kind=8), parameter :: pii = acos(-1.0d0)
  real(kind=8), parameter :: EPS = 1.0d-8
  integer,      parameter :: Re = 1, Im = 2

  type :: state_t
     integer :: Nr = 0
     integer :: Nt = 0
     real(kind=8) :: mu_field
     real(kind=8), allocatable :: phi(:,:), psi(:,:), pi(:,:)
     real(kind=8), allocatable :: Noether(:), Accretion_Noether(:)
     !real(kind=8), allocatable :: energy_density(:), klm(:), Omegalm(:), omglm(:), vr(:) 
   contains
     procedure :: init  => state_init
     procedure :: zero  => state_zero
     procedure :: info  => state_info
     procedure :: free  => state_free
     procedure :: diagnostic => diagnostic_state
     procedure :: NC => Noether_charge 
     procedure :: Noether_rate 
     procedure :: radial_wavenumber
     procedure :: local_frequency


  end type

contains

  subroutine state_init(this, Nr, Nt, mu_field)
    class(state_t), intent(inout) :: this
    integer,        intent(in)    :: Nr
    integer,        intent(in)    :: Nt
    real(kind=8),   intent(in)    :: mu_field

    call this%free()

    this%Nr = Nr
    this%Nt = Nt
    this%mu_field = mu_field

    allocate(this%phi(Nr,2), this%psi(Nr,2), this%pi(Nr,2))

    call this%zero()
  end subroutine

subroutine diagnostic_state(this, G, M, idx,hit)
  implicit none

  class(state_t),    intent(inout) :: this
  class(geometry_t), intent(in)    :: G
  class(mesh_t),     intent(in)    :: M
  integer,           intent(in)    :: idx
  logical,           intent(in)    :: hit

  !this%Accretion_Noether(idx) = this%Noether_rate(G, M)
  !this%Noether(idx) = this%NC(G, M)


end subroutine diagnostic_state

real(kind=8) function Noether_rate(this, G, M) result(dQdt)
  implicit none

  class(state_t),    intent(in) :: this
  class(geometry_t), intent(in) :: G
  class(mesh_t),     intent(in) :: M

  integer :: i_hor, Nr
  real(kind=8) :: flux
  real(kind=8) :: r_hor

  Nr = size(M%r)
  r_hor = 1.0d0

  ! Índices
  i_hor = nint( (r_hor - M%r(1)) / M%dr ) + 1
  i_hor = max(i_hor, 1)
  

  !flux = conjg(this%phi(i_hor)) * ( G%guu(i_hor) * this%psi(i_hor) + &
  !                                 (G%beta(i_hor) / G%alpha(i_hor)) * this%pi(i_hor) )
  flux = this%phi(i_hor, Re) * (G%guu(i_hor)*this%psi(i_hor, Im) + G%beta(i_hor)/G%alpha(i_hor)*this%pi(i_hor, Im)) &
        - this%phi(i_hor, Im) * (G%guu(i_hor)*this%psi(i_hor, Re) + G%beta(i_hor)/G%alpha(i_hor)*this%pi(i_hor, Re))
  dQdt = - G%alpha(i_hor) * M%r(i_hor)**2 * sqrt(G%grr(i_hor)) * flux


end function Noether_rate


real(kind=8) function Noether_charge(this, G, M) result(Q)
  use utils,          only: trapezium
  implicit none

  class(state_t),    intent(in) :: this
  class(geometry_t), intent(in) :: G
  class(mesh_t),     intent(in) :: M

  integer :: i_hor, Nr, i_inf

  Nr = size(M%r)

  i_hor = nint( (1.0d0 - M%r(1)) / M%dr ) + 1
  i_inf = nint( (300.0d0 - M%r(1)) / M%dr ) + 1
  i_hor = max(i_hor, 1)

  Q = - trapezium( M%r(i_hor:)**2 * sqrt(G%grr(i_hor:)) * &
                  this%phi(i_hor:,Re)*this%pi(i_hor:,Im) - this%phi(i_hor:,Im)*this%pi(i_hor:,Re), M)

end function Noether_charge



function radial_wavenumber(this,G,M) result(k_r)
  implicit none

  class(state_t),    intent(in) :: this
  class(geometry_t), intent(in) :: G
  class(mesh_t),     intent(in) :: M
  real(kind=8) :: k_r(M%Nr)


  k_r = &
    1.0d0 / max(dsqrt(this%phi(:,Re)**2 + this%phi(:,Im)**2),EPS) * & 
     (this%phi(:,Re)*this%psi(:,Im) - this%phi(:,Im)*this%psi(:,Re))
end function

function local_frequency(this,G,M) result(Omega)
  implicit none

  class(state_t),    intent(in) :: this
  class(geometry_t), intent(in) :: G
  class(mesh_t),     intent(in) :: M
  real(kind=8) :: Omega(M%Nr)

  Omega = &
    -1.0d0 / max(dsqrt(this%phi(:,Re)**2 + this%phi(:,Im)**2),EPS) * & 
     (this%phi(:,Re)*this%pi(:,Im) - this%phi(:,Im)*this%pi(:,Re))
end function


  subroutine state_zero(this)
    class(state_t), intent(inout) :: this
    this%phi =  0.0d0
    this%psi =  0.0d0
    this%pi  =  0.0d0
!    this%energy_density = 0.0d0
  end subroutine

  subroutine state_info(this)
    class(state_t), intent(in) :: this
    write(*,'(A)') "---- State info ----"
    write(*,'(A,I0)') "Nr = ", this%Nr
    write(*,'(A,ES14.6)') "max|phi| = ", maxval(dsqrt(this%phi(:,Re)**2 + this%phi(:,Im)**2))
    write(*,'(A,ES14.6)') "max|psi| = ", maxval(dsqrt(this%psi(:,Re)**2 + this%psi(:,Im)**2))
    write(*,'(A,ES14.6)') "max|pi|  = ", maxval(dsqrt(this%pi(:,Re)**2 + this%pi(:,Im)**2))
    write(*,'(A)') "--------------------"
  end subroutine

  subroutine state_free(this)
    class(state_t), intent(inout) :: this
    if (allocated(this%phi)) deallocate(this%phi)
    if (allocated(this%psi)) deallocate(this%psi)
    if (allocated(this%pi))  deallocate(this%pi)
    this%Nr = 0
    this%Nt = 0
  end subroutine

end module
