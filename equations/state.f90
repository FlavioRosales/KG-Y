module state
  use iso_fortran_env, only: real64
    use geometry,   only: geometry_t        ! <-- módulo que define TYPE :: geometry
    use mesh,       only: mesh_t          ! <-- módulo que define TYPE :: mesh_t

  implicit none
  private
  public :: state_t

  real(real64), parameter :: pii = acos(-1.0_real64)
  real(real64), parameter :: EPS = 1.0e-8_real64

  type :: state_t
     integer :: Nr = 0
     integer :: Nt = 0
     real(real64) :: mu_field
     integer :: ell
     complex(real64), allocatable :: phi(:), psi(:), pi(:)
     real(real64), allocatable :: Noether(:), Accretion_Noether(:)
     real(real64), allocatable :: energy_density(:), klm(:), Omegalm(:), omglm(:), vr(:) 
   contains
     procedure :: init  => state_init
     procedure :: zero  => state_zero
     procedure :: info  => state_info
     procedure :: free  => state_free
     procedure :: diagnostic => diagnostic_state
     procedure :: NC => Noether_charge 
     procedure :: Noether_rate 
     procedure :: energy
     procedure :: radial_wavenumber
     procedure :: local_frequency


  end type

contains

  subroutine state_init(this, Nr, Nt, ell, mu_field)
    class(state_t), intent(inout) :: this
    integer,        intent(in)    :: Nr
    integer,        intent(in)    :: Nt
    integer,        intent(in)    :: ell
    real(real64),   intent(in)    :: mu_field
    
    this%Nr = Nr
    this%Nt = Nt
    this%ell = ell
    this%mu_field = mu_field

    allocate(this%phi(Nr), this%psi(Nr), this%pi(Nr), &
    this%Noether(0:Nt), this%Accretion_Noether(0:Nt), &
    this%energy_density(Nr), this%klm(Nr), this%Omegalm(Nr), this%omglm(Nr), this%vr(Nr) )

    call this%zero()
  end subroutine

subroutine diagnostic_state(this, G, M, idx,hit)
  use iso_fortran_env, only: real64
  implicit none

  class(state_t),    intent(inout) :: this
  class(geometry_t), intent(in)    :: G
  class(mesh_t),     intent(in)    :: M
  integer,           intent(in)    :: idx
  logical,           intent(in)    :: hit
  real(real64) :: den(M%Nr)

  if(hit) then
    this%energy_density = this%energy(G, M)
    this%klm = this%radial_wavenumber(G, M)
    this%Omegalm = this%local_frequency(G, M)
    this%omglm = G%alpha * this%Omegalm - G%beta * this%klm
    den = sign(max(abs(this%Omegalm), EPS), this%Omegalm)
    this%vr = -G%beta + G%alpha * G%guu * this%klm / den
  end if
    ! Carga Noether
  this%Accretion_Noether(idx) = this%Noether_rate(G, M)
  this%Noether(idx) = this%NC(G, M)


end subroutine diagnostic_state

real(real64) function Noether_rate(this, G, M) result(dQdt)
  use iso_fortran_env, only: real64
  use utils,          only: trapezium
  implicit none

  class(state_t),    intent(in) :: this
  class(geometry_t), intent(in) :: G
  class(mesh_t),     intent(in) :: M

  integer :: i_hor, Nr
  complex(real64) :: flux
  real(real64) :: r_hor, f_hor

  Nr = size(M%r)
  r_hor = 1.0_real64

  ! Índices
  i_hor = nint( (r_hor - M%r(1)) / M%dr ) + 1
  i_hor = max(i_hor, 1)
  

  flux = conjg(this%phi(i_hor)) * ( G%guu(i_hor) * this%psi(i_hor) + &
                                   (G%beta(i_hor) / G%alpha(i_hor)) * this%pi(i_hor) )

  f_hor = G%alpha(i_hor) * M%r(i_hor)**2 * sqrt(G%grr(i_hor)) * aimag(flux)

  dQdt = - f_hor 

end function Noether_rate


real(real64) function Noether_charge(this, G, M) result(Q)
  use iso_fortran_env, only: real64
  use utils,          only: trapezium
  implicit none

  class(state_t),    intent(in) :: this
  class(geometry_t), intent(in) :: G
  class(mesh_t),     intent(in) :: M

  integer :: i_hor, Nr, i_inf

  Nr = size(M%r)

  i_hor = nint( (1.0_real64 - M%r(1)) / M%dr ) + 1
  i_inf = nint( (300.0d0 - M%r(1)) / M%dr ) + 1
  i_hor = max(i_hor, 1)

  Q = - trapezium( M%r(i_hor:)**2 * sqrt(G%grr(i_hor:)) * &
                   aimag( conjg(this%phi(i_hor:)) * this%pi(i_hor:) ), M%dr )

end function Noether_charge

function energy(this,G,M) result(rho_bar)
  use iso_fortran_env, only: real64
  use utils,          only: trapezium
  implicit none

  class(state_t),    intent(in) :: this
  class(geometry_t), intent(in) :: G
  class(mesh_t),     intent(in) :: M
  real(real64) :: rho_bar(M%Nr)

  rho_bar = &
    abs(this%pi)**2 + G%guu * abs(this%psi)**2 + &
    this%mu_field**2  +  (this%ell * (this%ell + 1) / M%r*2 ) * abs(this%phi)**2

end function

function radial_wavenumber(this,G,M) result(k_r)
  use iso_fortran_env, only: real64
  use utils,          only: trapezium
  implicit none

  class(state_t),    intent(in) :: this
  class(geometry_t), intent(in) :: G
  class(mesh_t),     intent(in) :: M
  real(real64) :: k_r(M%Nr)


  k_r = &
    1.0_real64 / max(abs(this%phi)**2,EPS) * aimag(conjg(this%phi) * this%psi)
end function

function local_frequency(this,G,M) result(Omega)
  use iso_fortran_env, only: real64
  use utils,          only: trapezium
  implicit none

  class(state_t),    intent(in) :: this
  class(geometry_t), intent(in) :: G
  class(mesh_t),     intent(in) :: M
  real(real64) :: Omega(M%Nr)

  Omega = &
    -1.0_real64 / max(abs(this%phi)**2,EPS) * aimag(conjg(this%phi) * this%pi)
end function


  subroutine state_zero(this)
    class(state_t), intent(inout) :: this
    this%phi = (0.0_real64, 0.0_real64)
    this%psi = (0.0_real64, 0.0_real64)
    this%pi  = (0.0_real64, 0.0_real64)
    this%Noether = 0.0_real64
    this%Accretion_Noether = 0.0_real64
    this%energy_density = 0.0_real64
  end subroutine

  subroutine state_info(this)
    class(state_t), intent(in) :: this
    write(*,'(A)') "---- State info ----"
    write(*,'(A,I0)') "Nr = ", this%Nr
    write(*,'(A,ES14.6)') "max|phi| = ", maxval(abs(this%phi))
    write(*,'(A,ES14.6)') "max|psi| = ", maxval(abs(this%psi))
    write(*,'(A,ES14.6)') "max|pi|  = ", maxval(abs(this%pi))
    write(*,'(A)') "--------------------"
  end subroutine

  subroutine state_free(this)
    class(state_t), intent(inout) :: this
    if (allocated(this%phi)) deallocate(this%phi)
    if (allocated(this%psi)) deallocate(this%psi)
    if (allocated(this%pi))  deallocate(this%pi)
    if (allocated(this%Noether)) deallocate(this%Noether)
    if (allocated(this%Accretion_Noether)) deallocate(this%Accretion_Noether)
    this%Nr = 0
  end subroutine

end module
