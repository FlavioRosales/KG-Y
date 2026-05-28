module mesh
  use iso_fortran_env, only: real64
  use run_control,    only: rmin_p, rmax_p, nr_p, ell_max_p
  implicit none
  private
  public :: mesh_t, build_mesh, print_mesh_info

  type :: mesh_t
     integer :: Nr = 0
     real(real64) :: rmin = 0.0_real64, rmax = 0.0_real64, dr = 0.0_real64
     real(real64), allocatable :: r(:)      ! puntos "de celda" (como ya tenías)
  end type mesh_t

contains

  subroutine build_mesh(M)
    type(mesh_t), intent(inout) :: M
    integer :: i

    M%Nr   = nr_p
    M%rmin = rmin_p
    M%rmax = rmax_p

    M%dr = (M%rmax - M%rmin) / real(M%Nr-1,real64)

    ! M%rmin = log(M%rmin)
    ! M%rmax = log(M%rmax)

    ! M%dr = (M%rmax - M%rmin) / real(M%Nr-1,real64)


    allocate(M%r(M%Nr))
    do i = 1, M%Nr
       M%r(i) = M%rmin + real(i-1,real64)*M%dr
    end do
    !allocate(M%r(M%Nr))
    !do i = 1, M%Nr
    !   M%r(i) = exp(M%rmin + real(i-1,real64)*M%dr)
    !end do



  end subroutine build_mesh

  subroutine print_mesh_info(M)
    type(mesh_t), intent(in) :: M
    write(*,'(A)')          "---- Mesh info ----"
    write(*,'(A,I0)')       "Nr       = ", M%Nr
    write(*,'(A,ES14.6)')   "rmin     = ", M%rmin
    write(*,'(A,ES14.6)')   "rmax     = ", M%rmax
    write(*,'(A,ES14.6)')   "dr       = ", M%dr
    write(*,'(A,ES14.6)')   "r(1)     = ", M%r(1)
    write(*,'(A,ES14.6)')   "r(N)     = ", M%r(M%Nr)
    write(*,'(A,I0)')       "ell_max  = ", ell_max_p
    write(*,'(A)')          "-------------------"
  end subroutine print_mesh_info

end module mesh
