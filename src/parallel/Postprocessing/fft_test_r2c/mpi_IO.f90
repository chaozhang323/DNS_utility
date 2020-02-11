!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This example calculates the divergency of a random field using
!   (1) global transposition
!   (2) halo-cell exchange
! The two method should give identical results
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program ddy_fft

  use decomp_2d
  use decomp_2d_io
  use decomp_2d_fft
  use HDF5

  implicit none

  include 'mpif.h'

!  integer, parameter :: imaxo=400, jmaxo=300, kmaxo=130
  integer, parameter :: imax=130, jmax=400, kmax=300
  integer, parameter :: p_row=3, p_col=1

  character(50) :: fname = "flowdata_00000000.h5"
  character(10) :: fname2 = "rgrid.h5"

  character(5) :: dname1 = "u"
  character(5) :: dname2 = "v"
  character(5) :: dname3 = "w"

  integer(HID_T) :: file_id
  integer(HID_T) :: dset_id
  integer(HID_T) :: filespace
  integer(HID_T) :: memspace
  integer(HID_T) :: plist_id

  integer(HSIZE_T), dimension(3) :: dimsf
  integer(HSIZE_T), dimension(3) :: dimsf_sub

  integer(HSIZE_T), dimension(3) :: count
  integer(HSIZE_T), dimension(3) :: offset
  integer(HSIZE_T), dimension(3) :: offset2
  integer :: rank = 3

  integer :: i,j,k, ierror,hdferr
  real(8), dimension(:,:,:), allocatable :: u, v, w
  real(8), dimension(:,:,:), allocatable :: u2
  complex(8), dimension(:,:,:), allocatable :: out
  integer, dimension(3) :: fft_start, fft_end, fft_size
!  real(8), dimension(:,:,:), allocatable :: xx, yy, zz

  call MPI_INIT(ierror)
  call decomp_2d_init(imax,jmax,kmax,p_row,p_col)
!  call decomp_2d_fft_init
     call decomp_2d_fft_init(PHYSICAL_IN_Z)
     call decomp_2d_fft_get_size(fft_start,fft_end,fft_size)

!  allocate(X(xsize(1),xsize(2),xsize(3)),Y(xsize(1),xsize(2),xsize(3)),Z(xsize(1),xsize(2),xsize(3)))
  allocate(u(zsize(1),zsize(2),zsize(3)))
!print *, 'nrank = ', nrank
!print *, 'xsize(1) = ', xsize(1)
!print *, 'xsize(2) = ', xsize(2)
!print *, 'xsize(3) = ', xsize(3)
!print *, 'ysize(1) = ', ysize(1)
!print *, 'ysize(2) = ', ysize(2)
!print *, 'ysize(3) = ', ysize(3)
!print *, 'zsize(1) = ', zsize(1)
!print *, 'zsize(2) = ', zsize(2)
!print *, 'zsize(3) = ', zsize(3)
!print *, 'xstart(1) = ', xstart(1), 'xend(1) = ', xend(1)
!print *, 'xstart(2) = ', xstart(2), 'xend(2) = ', xend(2)
!print *, 'xstart(3) = ', xstart(3), 'xend(3) = ', xend(3)
!print *, ystart
!print *, zstart

     dimsf(1) = zsize(1) !imax
     dimsf(2) = zsize(2) !jmax
     dimsf(3) = zsize(3) !kmax

     count(1) = zsize(1)
     count(2) = zsize(2)
     count(3) = zsize(3)
     offset(1) = zstart(1) - 1
     offset(2) = zstart(2) - 1
     offset(3) = zstart(3) - 1

     dimsf_sub(1) = count(1)
     dimsf_sub(2) = count(2)
     dimsf_sub(3) = count(3)

     offset2(1) = 0
     offset2(2) = 0
     offset2(3) = 0

     CALL h5open_f(hdferr)
 !    print *, 'hdferr1 = ', hdferr
     CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdferr)
 !    print *, 'hdferr2 = ', hdferr
     CALL h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, hdferr)
 !    print *, 'hdferr3 = ', hdferr
     CALL h5fopen_f(fname,H5F_ACC_RDONLY_F, file_id, hdferr, access_prp = plist_id )
 !    print *, 'hdferr4 = ', hdferr
     CALL h5pclose_f(plist_id, hdferr)

     CALL h5dopen_f(file_id, dname1,dset_id, hdferr)
     CALL h5dget_space_f(dset_id, filespace, hdferr)
     CALL h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, hdferr)
     CALL h5screate_simple_f(rank, dimsf_sub, memspace, hdferr)
     CALL h5sselect_hyperslab_f (memspace, H5S_SELECT_SET_F, offset2, count, hdferr)
     CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdferr)
     CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdferr)
     CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE,u, dimsf, hdferr, &
                     file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
     CALL h5sclose_f(filespace, hdferr)
     CALL h5sclose_f(memspace, hdferr)
     CALL h5dclose_f(dset_id, hdferr)
     CALL h5pclose_f(plist_id, hdferr)
     CALL h5fclose_f(file_id, hdferr)

     ! using Z-pencil input
!     call decomp_2d_fft_init(PHYSICAL_IN_Z)
!     call decomp_2d_fft_get_size(fft_start,fft_end,fft_size)
     allocate(out(fft_start(1):fft_end(1), &
                  fft_start(2):fft_end(2), &
                  fft_start(3):fft_end(3) ))

     call decomp_2d_fft_3d(u,out)

     print *, 'nrank = ', nrank
     print *, 'zsize(1) = ', zsize(1)
     print *, 'zsize(2) = ', zsize(2)
     print *, 'zsize(3) = ', zsize(3)
     print *, 'fft_start(1) = ', fft_start(1)
     print *, 'fft_end(1) = ', fft_end(1)
     print *, 'fft_start(2) = ', fft_start(2)
     print *, 'fft_end(2) = ', fft_end(2)
     print *, 'fft_start(3) = ', fft_start(3)
     print *, 'fft_end(3) = ', fft_end(3)
     print *, 'fft_size(1) = ', fft_size(1)
     print *, 'fft_size(2) = ', fft_size(2)
     print *, 'fft_size(3) = ', fft_size(3)

 !    do i=fft_start(1), fft_end(1)
 !      out(i,:,:) = out(i,:,:)*(0,i-1)
 !    enddo







  
  call decomp_2d_finalize 
  call MPI_FINALIZE(ierror)

end program ddy_fft
