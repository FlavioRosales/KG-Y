module hdf5_lib
  use hdf5
  implicit none
  private

  public :: h5_init, h5_finalize
  public :: h5_open_new, h5_close_file
  public :: h5_create_group, h5_close_group
  public :: h5_close_dataset
  public :: h5_write_real_1d, h5_write_int_1d
  public :: h5_create_real_time_series, h5_write_real_time_sample
  public :: h5_create_real_0d_series, h5_write_real_0d_sample
  public :: h5_create_real_1d_series, h5_write_real_1d_block
  public :: h5_create_real_2d_series, h5_write_real_2d_sample
  public :: read_sl_modes_h5

contains

  subroutine h5_init(ierr)
    integer, intent(out) :: ierr
    call h5open_f(ierr)
  end subroutine h5_init


  subroutine h5_finalize(ierr)
    integer, intent(out) :: ierr
    call h5close_f(ierr)
  end subroutine h5_finalize


  subroutine h5_open_new(filename, fid, ierr)
    character(len=*), intent(in)  :: filename
    integer(HID_T),   intent(out) :: fid
    integer,          intent(out) :: ierr

    call h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, fid, ierr)
  end subroutine h5_open_new


  subroutine h5_close_file(fid, ierr)
    integer(HID_T), intent(inout) :: fid
    integer,        intent(out)   :: ierr

    call h5fclose_f(fid, ierr)
    fid = -1
  end subroutine h5_close_file


  subroutine h5_create_group(loc_id, name, gid, ierr)
    integer(HID_T),   intent(in)  :: loc_id
    character(len=*), intent(in)  :: name
    integer(HID_T),   intent(out) :: gid
    integer,          intent(out) :: ierr

    call h5gcreate_f(loc_id, trim(name), gid, ierr)
  end subroutine h5_create_group


  subroutine h5_close_group(gid, ierr)
    integer(HID_T), intent(inout) :: gid
    integer,        intent(out)   :: ierr

    call h5gclose_f(gid, ierr)
    gid = -1
  end subroutine h5_close_group


  subroutine h5_close_dataset(did, ierr)
    integer(HID_T), intent(inout) :: did
    integer,        intent(out)   :: ierr

    call h5dclose_f(did, ierr)
    did = -1
  end subroutine h5_close_dataset


  subroutine h5_write_real_1d(loc_id, name, x, ierr)
    integer(HID_T),   intent(in)  :: loc_id
    character(len=*), intent(in)  :: name
    real(kind=8),     intent(in)  :: x(:)
    integer,          intent(out) :: ierr

    integer(HID_T)   :: sid, did
    integer(HSIZE_T) :: dims(1)
    integer          :: ierr2

    dims(1) = size(x, kind=HSIZE_T)

    call h5screate_simple_f(1, dims, sid, ierr)
    call h5dcreate_f(loc_id, trim(name), H5T_NATIVE_DOUBLE, sid, did, ierr)
    call h5dwrite_f(did, H5T_NATIVE_DOUBLE, x, dims, ierr)

    call h5dclose_f(did, ierr2)
    call h5sclose_f(sid, ierr2)
  end subroutine h5_write_real_1d


  subroutine h5_write_int_1d(loc_id, name, x, ierr)
    integer(HID_T),   intent(in)  :: loc_id
    character(len=*), intent(in)  :: name
    integer,          intent(in)  :: x(:)
    integer,          intent(out) :: ierr

    integer(HID_T)   :: sid, did
    integer(HSIZE_T) :: dims(1)
    integer          :: ierr2

    dims(1) = size(x, kind=HSIZE_T)

    call h5screate_simple_f(1, dims, sid, ierr)
    call h5dcreate_f(loc_id, trim(name), H5T_NATIVE_INTEGER, sid, did, ierr)
    call h5dwrite_f(did, H5T_NATIVE_INTEGER, x, dims, ierr)

    call h5dclose_f(did, ierr2)
    call h5sclose_f(sid, ierr2)
  end subroutine h5_write_int_1d


  subroutine h5_create_real_time_series(loc_id, name, nt, did, ierr)
    integer(HID_T),   intent(in)  :: loc_id
    character(len=*), intent(in)  :: name
    integer,          intent(in)  :: nt
    integer(HID_T),   intent(out) :: did
    integer,          intent(out) :: ierr

    integer(HID_T)   :: sid
    integer(HSIZE_T) :: dims(1)
    integer          :: ierr2

    dims(1) = int(nt, HSIZE_T)

    call h5screate_simple_f(1, dims, sid, ierr)
    call h5dcreate_f(loc_id, trim(name), H5T_NATIVE_DOUBLE, sid, did, ierr)
    call h5sclose_f(sid, ierr2)
  end subroutine h5_create_real_time_series


  subroutine h5_write_real_time_sample(did, it, value, ierr)
    integer(HID_T), intent(in)  :: did
    integer,        intent(in)  :: it
    real(kind=8),   intent(in)  :: value
    integer,        intent(out) :: ierr

    integer(HID_T)   :: fspace, mspace
    integer(HSIZE_T) :: start(1), count(1)
    integer          :: ierr2
    real(kind=8)     :: x(1)

    x(1) = value

    start(1) = int(it - 1, HSIZE_T)
    count(1) = 1_HSIZE_T

    call h5dget_space_f(did, fspace, ierr)
    call h5sselect_hyperslab_f(fspace, H5S_SELECT_SET_F, start, count, ierr)
    call h5screate_simple_f(1, count, mspace, ierr)

    call h5dwrite_f(did, H5T_NATIVE_DOUBLE, x, count, ierr, &
                    mem_space_id=mspace, file_space_id=fspace)

    call h5sclose_f(mspace, ierr2)
    call h5sclose_f(fspace, ierr2)
  end subroutine h5_write_real_time_sample


  ! dataset(nmodes, nt)
  subroutine h5_create_real_0d_series(loc_id, name, nmodes, nt, did, ierr)
    integer(HID_T),   intent(in)  :: loc_id
    character(len=*), intent(in)  :: name
    integer,          intent(in)  :: nmodes, nt
    integer(HID_T),   intent(out) :: did
    integer,          intent(out) :: ierr

    integer(HID_T)   :: sid
    integer(HSIZE_T) :: dims(2)
    integer          :: ierr2

    dims = [int(nmodes, HSIZE_T), int(nt, HSIZE_T)]

    call h5screate_simple_f(2, dims, sid, ierr)
    call h5dcreate_f(loc_id, trim(name), H5T_NATIVE_DOUBLE, sid, did, ierr)
    call h5sclose_f(sid, ierr2)
  end subroutine h5_create_real_0d_series


  ! dataset(:,it)
  subroutine h5_write_real_0d_sample(did, it, x, ierr)
    integer(HID_T), intent(in)  :: did
    integer,        intent(in)  :: it
    real(kind=8),   intent(in)  :: x(:)
    integer,        intent(out) :: ierr

    integer(HID_T)   :: fspace, mspace
    integer(HSIZE_T) :: start(2), count(2), memdims(1)
    integer          :: ierr2

    start   = [0_HSIZE_T, int(it - 1, HSIZE_T)]
    count   = [size(x, kind=HSIZE_T), 1_HSIZE_T]
    memdims = [size(x, kind=HSIZE_T)]

    call h5dget_space_f(did, fspace, ierr)
    call h5sselect_hyperslab_f(fspace, H5S_SELECT_SET_F, start, count, ierr)
    call h5screate_simple_f(1, memdims, mspace, ierr)

    call h5dwrite_f(did, H5T_NATIVE_DOUBLE, x, memdims, ierr, &
                    mem_space_id=mspace, file_space_id=fspace)

    call h5sclose_f(mspace, ierr2)
    call h5sclose_f(fspace, ierr2)
  end subroutine h5_write_real_0d_sample


  ! dataset(nr, nmodes, nt)
  subroutine h5_create_real_1d_series(loc_id, name, nr, nmodes, nt, did, ierr)
    integer(HID_T),   intent(in)  :: loc_id
    character(len=*), intent(in)  :: name
    integer,          intent(in)  :: nr, nmodes, nt
    integer(HID_T),   intent(out) :: did
    integer,          intent(out) :: ierr

    integer(HID_T)   :: sid
    integer(HSIZE_T) :: dims(3)
    integer          :: ierr2

    dims = [int(nr, HSIZE_T), int(nmodes, HSIZE_T), int(nt, HSIZE_T)]

    call h5screate_simple_f(3, dims, sid, ierr)
    call h5dcreate_f(loc_id, trim(name), H5T_NATIVE_DOUBLE, sid, did, ierr)
    call h5sclose_f(sid, ierr2)
  end subroutine h5_create_real_1d_series


  ! dataset(:, first_mode:first_mode+nb-1, it)
  subroutine h5_write_real_1d_block(did, it, first_mode, x, ierr)
    integer(HID_T), intent(in)  :: did
    integer,        intent(in)  :: it, first_mode
    real(kind=8),   intent(in)  :: x(:,:)
    integer,        intent(out) :: ierr

    integer(HID_T)   :: fspace, mspace
    integer(HSIZE_T) :: start(3), count(3), memdims(2)
    integer          :: ierr2

    start = [0_HSIZE_T, int(first_mode - 1, HSIZE_T), int(it - 1, HSIZE_T)]
    count = [size(x,1,kind=HSIZE_T), size(x,2,kind=HSIZE_T), 1_HSIZE_T]
    memdims = [size(x,1,kind=HSIZE_T), size(x,2,kind=HSIZE_T)]

    call h5dget_space_f(did, fspace, ierr)
    call h5sselect_hyperslab_f(fspace, H5S_SELECT_SET_F, start, count, ierr)
    call h5screate_simple_f(2, memdims, mspace, ierr)

    call h5dwrite_f(did, H5T_NATIVE_DOUBLE, x, memdims, ierr, &
                    mem_space_id=mspace, file_space_id=fspace)

    call h5sclose_f(mspace, ierr2)
    call h5sclose_f(fspace, ierr2)
  end subroutine h5_write_real_1d_block


  ! dataset(nr, nangle, nt)
  subroutine h5_create_real_2d_series(loc_id, name, nr, nangle, nt, did, ierr)
    integer(HID_T),   intent(in)  :: loc_id
    character(len=*), intent(in)  :: name
    integer,          intent(in)  :: nr, nangle, nt
    integer(HID_T),   intent(out) :: did
    integer,          intent(out) :: ierr

    integer(HID_T)   :: sid
    integer(HSIZE_T) :: dims(3)
    integer          :: ierr2

    dims = [int(nr, HSIZE_T), int(nangle, HSIZE_T), int(nt, HSIZE_T)]

    call h5screate_simple_f(3, dims, sid, ierr)
    call h5dcreate_f(loc_id, trim(name), H5T_NATIVE_DOUBLE, sid, did, ierr)
    call h5sclose_f(sid, ierr2)
  end subroutine h5_create_real_2d_series


  ! dataset(:,:,it)
  subroutine h5_write_real_2d_sample(did, it, x, ierr)
    integer(HID_T), intent(in)  :: did
    integer,        intent(in)  :: it
    real(kind=8),   intent(in)  :: x(:,:)
    integer,        intent(out) :: ierr

    integer(HID_T)   :: fspace, mspace
    integer(HSIZE_T) :: start(3), count(3), memdims(2)
    integer          :: ierr2

    start = [0_HSIZE_T, 0_HSIZE_T, int(it - 1, HSIZE_T)]
    count = [size(x,1,kind=HSIZE_T), size(x,2,kind=HSIZE_T), 1_HSIZE_T]
    memdims = [size(x,1,kind=HSIZE_T), size(x,2,kind=HSIZE_T)]

    call h5dget_space_f(did, fspace, ierr)
    call h5sselect_hyperslab_f(fspace, H5S_SELECT_SET_F, start, count, ierr)
    call h5screate_simple_f(2, memdims, mspace, ierr)

    call h5dwrite_f(did, H5T_NATIVE_DOUBLE, x, memdims, ierr, &
                    mem_space_id=mspace, file_space_id=fspace)

    call h5sclose_f(mspace, ierr2)
    call h5sclose_f(fspace, ierr2)
  end subroutine h5_write_real_2d_sample


  subroutine read_sl_modes_h5(filename, ell, lambda, eigenR, M)
    use mesh, only: mesh_t
    implicit none

    character(len=*), intent(in) :: filename
    integer,          intent(in) :: ell
    type(mesh_t),     intent(in), optional :: M

    real(kind=8), allocatable, intent(out) :: lambda(:)
    real(kind=8), allocatable, intent(out) :: eigenR(:,:)

    integer(HID_T) :: fid, gid, did, sid, mspace
    integer        :: ierr, ierr2, nm
    integer(HSIZE_T) :: dims1(1), maxdims1(1)
    integer(HSIZE_T) :: dims2(2), maxdims2(2)
    integer(HSIZE_T) :: start(2), count(2)
    character(len=32) :: gname
    real(kind=8), allocatable :: k(:)

    write(gname, '("ell_", I4.4)') ell

    call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, fid, ierr)
    call h5gopen_f(fid, trim(gname), gid, ierr)

    call h5dopen_f(gid, "k", did, ierr)
    call h5dget_space_f(did, sid, ierr)
    call h5sget_simple_extent_dims_f(sid, dims1, maxdims1, ierr)

    nm = int(dims1(1))

    allocate(k(nm), lambda(nm))
    call h5dread_f(did, H5T_NATIVE_DOUBLE, k, dims1, ierr)
    lambda = k**2
    deallocate(k)

    call h5sclose_f(sid, ierr2)
    call h5dclose_f(did, ierr2)

    call h5dopen_f(gid, "modes", did, ierr)
    call h5dget_space_f(did, sid, ierr)

    if (present(M)) then
      allocate(eigenR(M%Nr, nm))

      start = [int(M%i_sl_first - 1, HSIZE_T), 0_HSIZE_T]
      count = [int(M%Nr, HSIZE_T), int(nm, HSIZE_T)]

      call h5sselect_hyperslab_f(sid, H5S_SELECT_SET_F, start, count, ierr)
      call h5screate_simple_f(2, count, mspace, ierr)

      call h5dread_f(did, H5T_NATIVE_DOUBLE, eigenR, count, ierr, &
                     mem_space_id=mspace, file_space_id=sid)

      call h5sclose_f(mspace, ierr2)
    else
      call h5sget_simple_extent_dims_f(sid, dims2, maxdims2, ierr)
      allocate(eigenR(int(dims2(1)), int(dims2(2))))
      call h5dread_f(did, H5T_NATIVE_DOUBLE, eigenR, dims2, ierr)
    end if

    call h5sclose_f(sid, ierr2)
    call h5dclose_f(did, ierr2)
    call h5gclose_f(gid, ierr2)
    call h5fclose_f(fid, ierr2)
  end subroutine read_sl_modes_h5

end module hdf5_lib