  !> paraller version for read/write HDF5 files
  !>    InitGridHDF5
  !> 1. ReadHDF5_scalar         : read scalar data
  !> 2. WriteHDF5_scalar        : write scalar data
  !> 3. ReadtsGridMetrics_P     : read timeseries GridMetrics file "timeseries_GridMetrics.h5"
  !> 4. Readtsflow_P            : read timeseries flow data
  !> 5. Writetsplanes_P         : write timeseries plane data (not complete)
  !> 6. WriteHDF5grid_P         : write grid.h5
  !> 7. WriteHDF5flow_P         : write flowdata_xxxxxxxx.h5
  !> 8. ReadHDF5grid_P          : read grid.h5
  !> 9. ReadHDF5sol_P           : read flowdata_xxxxxxxx.h5
  !>10. Write3dHDF5_svariable_P : write 3d single variables
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    module MPRWHDF5
      use HDF5
      implicit none
      include 'mpif.h'

      integer(HID_T), private :: file_id, group_id, dset_id, dspace_id
      integer(HID_T), private :: memspace, plist_id, plist_file, plist_cdrw
      integer, private :: hdferr, errcode, ierr
      logical :: IsHDF5Initialized = .false.

      type tp_rdwt_hdf5
        integer :: rank, dnum, comm, info
        character(300) :: fname
        character(100) :: gname
        character(50)  :: sname  ! scalar name
        character(50), dimension(:), allocatable :: dname
        integer(HSIZE_T), dimension(:), allocatable :: dimsf
        integer(HSIZE_T), dimension(:), allocatable :: dimsm
        integer(HSIZE_T), dimension(:), allocatable :: count, offset, block, stride
        logical :: IsHSInitialized = .false.
        logical :: IsMultiGroup = .false. !this was disabled
      end type tp_rdwt_hdf5

     type tp_DNSIndex
        integer :: ibe, iend, iskip
        integer :: jbe, jend, jskip
        integer :: kbe, kend, kskip
        integer :: ibuffer, jbuffer
        integer :: nzloc, njloc, niloc
        integer :: nkplane, njplane, niplane
        integer, dimension(:), allocatable :: kplane, jplane, iplane
        real(8), dimension(:), allocatable :: zloc, dzdk, jloc, iloc
        real(8) :: dx, dy
     end type tp_DNSIndex

     type fprop
       character(400) :: filepath
       integer :: ntpoint, nskip
       integer :: file_be, file_end, file_skip
     end type fprop

  contains

     subroutine InitHDF5()
      ! Initialize HDF5 library and Fortran interfaces.
       CALL h5open_f(hdferr)
       IsHDF5Initialized = .True.
     end subroutine InitHDF5

     subroutine FinalizeHDF5()
      ! Close FORTRAN interfaces and HDF5 library.
       CALL h5close_f(hdferr)
     end subroutine FinalizeHDF5

     subroutine InitGridHDF5(hslab)
       type(tp_rdwt_hdf5), intent(out) :: hslab

       hslab%comm = MPI_COMM_WORLD
       hslab%info = MPI_INFO_NULL
       hslab%gname = '/'
       hslab%dnum = 3
       hslab%rank = 3

       allocate(hslab%dname(hslab%dnum), hslab%dimsf(hslab%rank), hslab%dimsm(hslab%rank))
       allocate(hslab%count(hslab%rank), hslab%offset(hslab%rank), hslab%block(hslab%rank), hslab%stride(hslab%rank))

       hslab%dname(1) = 'x'
       hslab%dname(2) = 'y'
       hslab%dname(3) = 'z'
       hslab%IsHSInitialized = .true.
       ! default value
       hslab%count = 1
       hslab%stride = 1
       hslab%offset = 0

     end subroutine InitGridHDF5


     subroutine InitFlowHDF5(hslab)
       type(tp_rdwt_hdf5), intent(out) :: hslab

       hslab%comm = MPI_COMM_WORLD
       hslab%info = MPI_INFO_NULL
       hslab%gname = '/'
       hslab%dnum = 5
       hslab%rank = 3

       allocate(hslab%dname(hslab%dnum), hslab%dimsf(hslab%rank), hslab%dimsm(hslab%rank))
       allocate(hslab%count(hslab%rank), hslab%offset(hslab%rank), hslab%block(hslab%rank), hslab%stride(hslab%rank))

       hslab%dname(1) = 'u'
       hslab%dname(2) = 'v'
       hslab%dname(3) = 'w'
       hslab%dname(4) = 'p'
       hslab%dname(5) = 'T'
       hslab%sname = 'time'
       hslab%IsHSInitialized = .true.
       ! default value
       hslab%count = 1
       hslab%stride = 1
       hslab%offset = 0
     end subroutine InitFlowHDF5

    subroutine ReadDNS_index_kplane(nt, dx, dy, dt, DNSIndex, filepath)
       integer, intent(out) :: nt
       real(8), intent(out) :: dx, dy, dt
       type(tp_DNSIndex), intent(out) :: DNSIndex
       character(*), intent(in) :: filepath
       type(tp_rdwt_hdf5) :: Index_kplane
       real(8), dimension(:,:), allocatable :: buffer_real
       integer, dimension(:,:), allocatable :: buffer_integer
       integer, dimension(:,:), allocatable :: buffer_klocs
       real(8) :: time1, time2

       open(11,file=trim(filepath)//'series_time_ascii.dat',status='old')
         read(11,*) time1
         read(11,*) time2
       close(11)
       dt = time2 - time1

       Index_kplane%gname = '/kplane'
       Index_kplane%rank  = 1

       Index_kplane%dnum  = 2
       allocate(Index_kplane%dname(Index_kplane%dnum),Index_kplane%dimsf(Index_kplane%rank))
       allocate(Index_kplane%dimsm(Index_kplane%rank),Index_kplane%count(Index_kplane%rank))
       allocate(Index_kplane%offset(Index_kplane%rank),Index_kplane%block(Index_kplane%rank))
       allocate(Index_kplane%stride(Index_kplane%rank))
       Index_kplane%dname(1)  = 'dx'
       Index_kplane%dname(2)  = 'dy'
       Index_kplane%dimsf = (/1/)
       Index_kplane%dimsm = Index_kplane%dimsf
       Index_kplane%block = Index_kplane%dimsm
       Index_kplane%offset = 0
       Index_kplane%count  = 1
       Index_kplane%stride = 1
       Index_kplane%IsHSInitialized = .true.
       allocate(buffer_real(1,2))

!       Index_kplane%fname = 'DNS_index.h5'
       Index_kplane%fname = trim(filepath)//'DNS_index.h5'

       call ReadTSHDF5_1D(Index_kplane,buffer_real)

       dx = buffer_real(1,1)
       dy = buffer_real(1,2)

       deallocate(Index_kplane%dname)
       Index_kplane%dnum  = 8
       allocate(Index_kplane%dname(Index_kplane%dnum))

       Index_kplane%dname(1)  = 'istart'
       Index_kplane%dname(2)  = 'iend'
       Index_kplane%dname(3)  = 'iskip'
       Index_kplane%dname(4)  = 'jstart'
       Index_kplane%dname(5)  = 'jend'
       Index_kplane%dname(6)  = 'jskip'
       Index_kplane%dname(7)  = 'nkplane'
       Index_kplane%dname(8)  = 'ntpoint'

       allocate(buffer_integer(1,8))

       call ReadTSHDF5_1D_integer(Index_kplane,buffer_integer)

       DNSIndex%ibe     = buffer_integer(1,1)
       DNSIndex%iend    = buffer_integer(1,2)
       DNSIndex%iskip   = buffer_integer(1,3)
       DNSIndex%jbe     = buffer_integer(1,4)
       DNSIndex%jend    = buffer_integer(1,5)
       DNSIndex%jskip   = buffer_integer(1,6)
       DNSIndex%nkplane = buffer_integer(1,7)
       nt = buffer_integer(1,8)

       deallocate(Index_kplane%dname)
       Index_kplane%dnum  = 1
       allocate(Index_kplane%dname(Index_kplane%dnum))
       Index_kplane%dname(1) = 'klocs'
       Index_kplane%dimsf = (/DNSIndex%nkplane/)
       Index_kplane%dimsm = Index_kplane%dimsf
       Index_kplane%block = Index_kplane%dimsm
       allocate(buffer_klocs(DNSIndex%nkplane,1))

       call ReadTSHDF5_1D_integer(Index_kplane,buffer_klocs)

       allocate(DNSIndex%kplane(DNSIndex%nkplane))
       DNSIndex%kplane = buffer_klocs(:,1)

     end subroutine ReadDNS_index_kplane

   subroutine ReadTSHDF5_1D_integer(TSplane, buffer)
      type(tp_rdwt_hdf5), intent(in) :: TSplane
      integer, intent(out) :: buffer(TSplane%dimsf(1),TSplane%dnum)

      integer :: n
      integer(HID_T) :: file_id, group_id, memspace
      integer(HID_T) :: dspace_id, dset_id

      if(.not.(IsHDF5Initialized.and.TSplane%IsHSInitialized)) then
          print *, 'HDF5 or ty_Hyperslab NOT intialized !!! STOP'
          stop
      endif

      CALL h5fopen_f(trim(TSplane%fname), H5F_ACC_RDONLY_F, file_id, hdferr)
      CALL h5gopen_f(file_id, trim(TSplane%gname), group_id, hdferr)
      CALL h5screate_simple_f(TSplane%rank, TSplane%dimsm, memspace, hdferr)
      do n = 1, TSplane%dnum
         CALL h5dopen_f(group_id, trim(TSplane%dname(n)), dset_id, hdferr)
         CALL h5dget_space_f(dset_id, dspace_id, hdferr)
         CALL h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, TSplane%offset, &
                                    TSplane%count, hdferr, TSplane%stride, TSplane%block)
         CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, buffer(:,n),  &
                             TSplane%dimsf, hdferr, file_space_id=dspace_id, mem_space_id=memspace)
         CALL h5sclose_f(dspace_id, hdferr)
         CALL h5dclose_f(dset_id, hdferr)
     enddo
     CALL h5sclose_f(memspace, hdferr)
     CALL h5gclose_f(group_id, hdferr)
     CALL h5fclose_f(file_id, hdferr)

   end subroutine ReadTSHDF5_1D_integer

   subroutine ReadTSHDF5_1D(TSplane, buffer)
      type(tp_rdwt_hdf5), intent(in) :: TSplane
      real(8), intent(out) :: buffer(TSplane%dimsf(1),TSplane%dnum)

      integer :: n
      integer(HID_T) :: file_id, group_id, memspace
      integer(HID_T) :: dspace_id, dset_id

      if(.not.(IsHDF5Initialized.and.TSplane%IsHSInitialized)) then
          print *, 'HDF5 or ty_Hyperslab NOT intialized !!! STOP'
          stop
      endif

      CALL h5fopen_f(trim(TSplane%fname), H5F_ACC_RDONLY_F, file_id, hdferr)
      CALL h5gopen_f(file_id, trim(TSplane%gname), group_id, hdferr)
      CALL h5screate_simple_f(TSplane%rank, TSplane%dimsm, memspace, hdferr)
      do n = 1, TSplane%dnum
         CALL h5dopen_f(group_id, trim(TSplane%dname(n)), dset_id, hdferr)
         CALL h5dget_space_f(dset_id, dspace_id, hdferr)
         CALL h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, TSplane%offset, &
                                    TSplane%count, hdferr, TSplane%stride, TSplane%block)
         CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, buffer(:,n),  &
                             TSplane%dimsf, hdferr, file_space_id=dspace_id, mem_space_id=memspace)
         CALL h5sclose_f(dspace_id, hdferr)
         CALL h5dclose_f(dset_id, hdferr)
     enddo
     CALL h5sclose_f(memspace, hdferr)
     CALL h5gclose_f(group_id, hdferr)
     CALL h5fclose_f(file_id, hdferr)

   end subroutine ReadTSHDF5_1D

    ! Read scalar in the file
    subroutine ReadHDF5_scalar(hslab, scalar)
      type(tp_rdwt_hdf5), intent(in) :: hslab
      real(8), intent(out) :: scalar
      integer(HSIZE_T), dimension(3) :: dimsf
      dimsf = 1

      if(.not.(IsHDF5Initialized)) then
          print *, 'HDF5 or type hslab NOT intialized !!! STOP'
          stop
      endif

      CALL h5fopen_f(trim(hslab%fname), H5F_ACC_RDONLY_F, file_id, hdferr)
      CALL h5dopen_f(file_id, trim(hslab%sname), dset_id, hdferr)
      CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, scalar, dimsf, hdferr)
      CALL h5dclose_f(dset_id, hdferr)
      CALL h5fclose_f(file_id, hdferr)

    end subroutine ReadHDF5_scalar


    ! Write scalar to the file
    subroutine WriteHDF5_scalar(hslab, scalar)
      type(tp_rdwt_hdf5), intent(in) :: hslab
      real(8), intent(in) :: scalar
      integer(HSIZE_T), dimension(3) :: dimsf
      dimsf = 1

      if(.not.(IsHDF5Initialized)) then
          print *, 'HDF5 or type hslab NOT intialized !!! STOP'
          stop
      endif

      CALL h5fopen_f(trim(hslab%fname),H5F_ACC_RDWR_F, file_id, hdferr)
      CALL h5screate_f(H5S_SCALAR_F, dspace_id, hdferr)
      CALL h5dcreate_f(file_id, trim(hslab%sname),H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, scalar, dimsf, hdferr)
      CALL h5dclose_f(dset_id, hdferr)
      CALL h5sclose_f(dspace_id, hdferr)
      CALL h5fclose_f(file_id, hdferr)

    end subroutine WriteHDF5_scalar

   subroutine ReadHDF5_4D(hslab, buffer)
      type(tp_rdwt_hdf5), intent(in) :: hslab
      real(8), intent(out) :: buffer(hslab%dimsm(1),hslab%dimsm(2),hslab%dimsm(3),hslab%dimsm(4),hslab%dnum)
      integer :: n

      if(.not.(IsHDF5Initialized.and.hslab%IsHSInitialized)) then
          print *, 'HDF5 or type hslab NOT intialized !!! STOP'
          errcode = 314
          call  MPI_Abort(MPI_COMM_WORLD,errcode,ierr)
      endif

      ! Setup file access property list with parallel I/O access.
      CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_file, hdferr)
      CALL h5pset_fapl_mpio_f(plist_file, hslab%comm, hslab%info, hdferr)

     ! Create property list for collective dataset read/write
      CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_cdrw, hdferr)
      CALL h5pset_dxpl_mpio_f(plist_cdrw, H5FD_MPIO_COLLECTIVE_F, hdferr)

      CALL h5fopen_f(trim(hslab%fname), H5F_ACC_RDONLY_F, file_id, hdferr, access_prp = plist_file)
      CALL h5gopen_f(file_id, trim(hslab%gname), group_id, hdferr)
      CALL h5screate_simple_f(hslab%rank, hslab%dimsm, memspace, hdferr)
      do n = 1, hslab%dnum
         CALL h5dopen_f(group_id, trim(hslab%dname(n)), dset_id, hdferr)
         CALL h5dget_space_f(dset_id, dspace_id, hdferr)
         CALL h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, hslab%offset, &
                                    hslab%count, hdferr, hslab%stride, hslab%block)
         CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, buffer(:,:,:,:,n),  &
              hslab%dimsm, hdferr, file_space_id=dspace_id, mem_space_id=memspace, xfer_prp = plist_cdrw)
         CALL h5sclose_f(dspace_id, hdferr)
         CALL h5dclose_f(dset_id, hdferr)
     enddo
     CALL h5sclose_f(memspace, hdferr)
     CALL h5gclose_f(group_id, hdferr)
     CALL h5fclose_f(file_id, hdferr)
     CALL h5pclose_f(plist_file, hdferr)
     CALL h5pclose_f(plist_cdrw, hdferr)

   end subroutine ReadHDF5_4D

   subroutine ReadHDF5_3D(hslab, buffer)
      type(tp_rdwt_hdf5), intent(in) :: hslab
      real(8), intent(out) :: buffer(hslab%dimsm(1),hslab%dimsm(2),hslab%dimsm(3),hslab%dnum)
      integer :: n
!      integer(HID_T) :: file_id, group_id, memspace
!      integer(HID_T) :: dspace_id, dset_id
!      integer(HID_T) :: plist_file, plist_cdrw

      if(.not.(IsHDF5Initialized.and.hslab%IsHSInitialized)) then
          print *, 'HDF5 or type hslab NOT intialized !!! STOP'
          errcode = 410
          call  MPI_Abort(MPI_COMM_WORLD,errcode,ierr)
      endif

      ! Setup file access property list with parallel I/O access.
      CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_file, hdferr)
      CALL h5pset_fapl_mpio_f(plist_file, hslab%comm, hslab%info, hdferr)

     ! Create property list for collective dataset read/write
      CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_cdrw, hdferr)
      CALL h5pset_dxpl_mpio_f(plist_cdrw, H5FD_MPIO_COLLECTIVE_F, hdferr)

      CALL h5fopen_f(trim(hslab%fname), H5F_ACC_RDONLY_F, file_id, hdferr, access_prp = plist_file)
      CALL h5gopen_f(file_id, trim(hslab%gname), group_id, hdferr)
      CALL h5screate_simple_f(hslab%rank, hslab%dimsm, memspace, hdferr)
      do n = 1, hslab%dnum
         CALL h5dopen_f(group_id, trim(hslab%dname(n)), dset_id, hdferr)
         CALL h5dget_space_f(dset_id, dspace_id, hdferr)
         CALL h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, hslab%offset, &
                                    hslab%count, hdferr, hslab%stride, hslab%block)
         CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, buffer(:,:,:,n),  &
              hslab%dimsm, hdferr, file_space_id=dspace_id, mem_space_id=memspace, xfer_prp = plist_cdrw)
         CALL h5sclose_f(dspace_id, hdferr)
         CALL h5dclose_f(dset_id, hdferr)
     enddo
     CALL h5sclose_f(memspace, hdferr)
     CALL h5gclose_f(group_id, hdferr)
     CALL h5fclose_f(file_id, hdferr)
     CALL h5pclose_f(plist_file, hdferr)
     CALL h5pclose_f(plist_cdrw, hdferr)

   end subroutine ReadHDF5_3D

  subroutine WriteHDF5_3D(hslab, buffer)
      type(tp_rdwt_hdf5), intent(in) :: hslab
      real(8), intent(in) :: buffer(hslab%dimsm(1),hslab%dimsm(2),hslab%dimsm(3),hslab%dnum)
      integer :: n
!      integer(HID_T) :: file_id, group_id, memspace
!      integer(HID_T) :: dspace_id, dset_id
!      integer(HID_T) :: plist_file, plist_cdrw

      if(.not.(IsHDF5Initialized.and.hslab%IsHSInitialized)) then
          print *, 'HDF5 or type hslab NOT intialized !!! STOP'
          errcode = 493
          call  MPI_Abort(MPI_COMM_WORLD,errcode,ierr)
      endif

      ! Setup file access property list with parallel I/O access.
      CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_file, hdferr)
      CALL h5pset_fapl_mpio_f(plist_file, hslab%comm, hslab%info, hdferr)

     ! Create property list for collective dataset read/write
      CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_cdrw, hdferr)
      CALL h5pset_dxpl_mpio_f(plist_cdrw, H5FD_MPIO_COLLECTIVE_F, hdferr)

      CALL h5fcreate_f(trim(hslab%fname), H5F_ACC_TRUNC_F, file_id, hdferr,access_prp = plist_file)
!      CALL h5gcreate_f(file_id,trim(hslab%gname), group_id, hdferr)
      CALL h5screate_simple_f(hslab%rank, hslab%dimsm, memspace, hdferr)
      CALL h5screate_simple_f(hslab%rank, hslab%dimsf, dspace_id, hdferr)
!      print *, 'hanging at 4'
      do n = 1, hslab%dnum
!         CALL h5dcreate_f(group_id, trim(hslab%dname(n)), &
!                          H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)
         CALL h5dcreate_f(file_id, trim(hslab%dname(n)), &
                          H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)

         CALL h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, hslab%offset, &
                                    hslab%count, hdferr, hslab%stride, hslab%block)
         CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, buffer(:,:,:,n),  &
                         hslab%dimsf, hdferr, file_space_id=dspace_id, mem_space_id=memspace, xfer_prp = plist_cdrw)
         CALL h5dclose_f(dset_id, hdferr)
     enddo
!     print *, 'hanging at 5'
     CALL h5sclose_f(dspace_id, hdferr)
!     print *, 'hanging at 6'
     CALL h5sclose_f(memspace, hdferr)
!     CALL h5gclose_f(group_id, hdferr)
     CALL h5fclose_f(file_id, hdferr)
!     print *, 'hanging at 7'
     CALL h5pclose_f(plist_file, hdferr)
!     print *, 'hanging at 8'
     CALL h5pclose_f(plist_cdrw, hdferr)

   end subroutine WriteHDF5_3D

    ! Read 3D data
    ! hslab%dimsf and hslab%offset (optional) should be provided
    ! by default offset = 0
    subroutine ReadHDF5_3D_serial(hslab, buffer)
      type(tp_rdwt_hdf5), intent(in) :: hslab
      real(8), intent(out) :: buffer(hslab%dimsf(1),hslab%dimsf(2),hslab%dimsf(3),hslab%dnum)
      integer :: n

      if(.not.(IsHDF5Initialized.and.hslab%IsHSInitialized)) then
          print *, 'HDF5 or type hslab NOT intialized !!! STOP'
          stop
      endif
      ! hslab%dimsm = hslab%dimsf
      ! hslab%block = hslab%dimsm

      CALL h5fopen_f(trim(hslab%fname), H5F_ACC_RDONLY_F, file_id, hdferr)
      CALL h5gopen_f(file_id, trim(hslab%gname), group_id, hdferr)
      !CALL h5screate_simple_f(hslab%rank, dimsm, memspace_id, hdferr)
      CALL h5screate_simple_f(hslab%rank, hslab%dimsf, memspace, hdferr)
      do n=1, hslab%dnum
        CALL h5dopen_f(group_id, trim(hslab%dname(n)), dset_id, hdferr)
        CALL h5dget_space_f(dset_id, dspace_id, hdferr)
        CALL h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, hslab%offset, &
                                   hslab%count, hdferr, hslab%stride, hslab%dimsf)
        CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, buffer(:,:,:,n), &
                       hslab%dimsf, hdferr, file_space_id=dspace_id, mem_space_id=memspace)

        ! CALL h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, hslab%offset, &
        !                             hslab%count, hdferr, hslab%stride, hslab%block)
        ! CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, buffer(:,:,:,n), &
        !                hslab%dimsm, hdferr, file_space_id=dspace_id, mem_space_id=memspace_id)
        CALL h5sclose_f(dspace_id, hdferr)
        CALL h5dclose_f(dset_id, hdferr)
      enddo
      CALL h5sclose_f(memspace, hdferr)
      CALL h5gclose_f(group_id, hdferr)
      CALL h5fclose_f(file_id, hdferr)

    end subroutine ReadHDF5_3D_serial

    subroutine WriteHDF5_3D_serial(hslab, buffer)
      type(tp_rdwt_hdf5), intent(in) :: hslab
      real(8), intent(in) :: buffer(hslab%dimsf(1),hslab%dimsf(2),hslab%dimsf(3),hslab%dnum)
      integer :: n

      if(.not.(IsHDF5Initialized.and.hslab%IsHSInitialized)) then
          print *, 'HDF5 or type hslab NOT intialized !!! STOP'
          stop
      endif

      if(.not.hslab%IsMultiGroup) then
        CALL h5fcreate_f(trim(hslab%fname), H5F_ACC_TRUNC_F, file_id, hdferr)
      else
        CALL h5fopen_f(trim(hslab%fname),H5F_ACC_RDWR_F, file_id, hdferr)
      endif

      if(trim(hslab%gname).ne.'/') then
        CALL h5gcreate_f(file_id, trim(hslab%gname),group_id, hdferr)
      endif
      CALL h5screate_simple_f(hslab%rank, hslab%dimsf, dspace_id, hdferr)
      do n=1, hslab%dnum
        if(trim(hslab%gname).ne.'/') then
          CALL h5dcreate_f(group_id, trim(hslab%dname(n)), H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr )
        else
          CALL h5dcreate_f(file_id, trim(hslab%dname(n)), H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr )
        endif
        CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, buffer(:,:,:,n),hslab%dimsf, hdferr)
        CALL h5dclose_f(dset_id, hdferr)
      enddo
      CALL h5sclose_f(dspace_id, hdferr)
      if(trim(hslab%gname).ne.'/') then
        CALL h5gclose_f(group_id, hdferr)
      endif
      CALL h5fclose_f(file_id, hdferr)

    end subroutine WriteHDF5_3D_serial

   !---------------------------------------------------
    subroutine ReadHDF5_3D_1V(hslab, buffer, nv)
      type(tp_rdwt_hdf5), intent(in) :: hslab
      integer,intent(in) :: nv
      real(8), intent(out) :: buffer(hslab%dimsm(1),hslab%dimsm(2),hslab%dimsm(3))
      integer :: n
!      integer(HID_T) :: file_id, group_id, memspace
!      integer(HID_T) :: dspace_id, dset_id
!      integer(HID_T) :: plist_file, plist_cdrw

      if(.not.(IsHDF5Initialized.and.hslab%IsHSInitialized)) then
          print *, 'HDF5 or type hslab NOT intialized !!! STOP'
          errcode = 410
          call  MPI_Abort(MPI_COMM_WORLD,errcode,ierr)
      endif

      ! Setup file access property list with parallel I/O access.
      CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_file, hdferr)
      CALL h5pset_fapl_mpio_f(plist_file, hslab%comm, hslab%info, hdferr)

     ! Create property list for collective dataset read/write
      CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_cdrw, hdferr)
      CALL h5pset_dxpl_mpio_f(plist_cdrw, H5FD_MPIO_COLLECTIVE_F, hdferr)

      CALL h5fopen_f(trim(hslab%fname), H5F_ACC_RDONLY_F, file_id, hdferr, access_prp = plist_file)
      CALL h5gopen_f(file_id, trim(hslab%gname), group_id, hdferr)
      CALL h5screate_simple_f(hslab%rank, hslab%dimsm, memspace, hdferr)
         CALL h5dopen_f(group_id, trim(hslab%dname(nv)), dset_id, hdferr)
         CALL h5dget_space_f(dset_id, dspace_id, hdferr)
         CALL h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, hslab%offset, &
                                    hslab%count, hdferr, hslab%stride, hslab%block)
         CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, buffer,  &
              hslab%dimsm, hdferr, file_space_id=dspace_id, mem_space_id=memspace, xfer_prp = plist_cdrw)
         CALL h5sclose_f(dspace_id, hdferr)
         CALL h5dclose_f(dset_id, hdferr)
      CALL h5sclose_f(memspace, hdferr)
      CALL h5gclose_f(group_id, hdferr)
      CALL h5fclose_f(file_id, hdferr)
      CALL h5pclose_f(plist_file, hdferr)
      CALL h5pclose_f(plist_cdrw, hdferr)

    end subroutine ReadHDF5_3D_1V

   subroutine WriteHDF5_3D_1V(hslab, buffer, nv)
      type(tp_rdwt_hdf5), intent(in) :: hslab
      integer,intent(in) :: nv
      real(8), intent(in) :: buffer(hslab%dimsm(1),hslab%dimsm(2),hslab%dimsm(3))
!      integer(HID_T) :: file_id, group_id, memspace
!      integer(HID_T) :: dspace_id, dset_id
!      integer(HID_T) :: plist_file, plist_cdrw

      if(.not.(IsHDF5Initialized.and.hslab%IsHSInitialized)) then
          print *, 'HDF5 or type hslab NOT intialized !!! STOP'
          errcode = 493
          call  MPI_Abort(MPI_COMM_WORLD,errcode,ierr)
      endif

      ! Setup file access property list with parallel I/O access.
      CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_file, hdferr)
      CALL h5pset_fapl_mpio_f(plist_file, hslab%comm, hslab%info, hdferr)

      ! Create property list for collective dataset read/write
      CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_cdrw, hdferr)
      CALL h5pset_dxpl_mpio_f(plist_cdrw, H5FD_MPIO_COLLECTIVE_F, hdferr)

      if(nv.eq.1) then
          CALL h5fcreate_f(trim(hslab%fname), H5F_ACC_TRUNC_F, file_id, hdferr,access_prp = plist_file)
      else
          call h5fopen_f(trim(hslab%fname),H5F_ACC_RDWR_F, file_id, hdferr, access_prp = plist_file)
      endif     ! endif nv.eq.1
!      CALL h5gcreate_f(file_id,trim(hslab%gname), group_id, hdferr)
      CALL h5screate_simple_f(hslab%rank, hslab%dimsm, memspace, hdferr)
      CALL h5screate_simple_f(hslab%rank, hslab%dimsf, dspace_id, hdferr)
!      print *, 'hanging at 4'
!         CALL h5dcreate_f(group_id, trim(hslab%dname(n)), &
!                          H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)
         CALL h5dcreate_f(file_id, trim(hslab%dname(nv)), &
                          H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)

         CALL h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, hslab%offset, &
                                    hslab%count, hdferr, hslab%stride, hslab%block)
         CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, buffer,  &
                         hslab%dimsf, hdferr, file_space_id=dspace_id, mem_space_id=memspace, xfer_prp = plist_cdrw)
         CALL h5dclose_f(dset_id, hdferr)
!     print *, 'hanging at 5'
     CALL h5sclose_f(dspace_id, hdferr)
!     print *, 'hanging at 6'
     CALL h5sclose_f(memspace, hdferr)
!     CALL h5gclose_f(group_id, hdferr)
     CALL h5fclose_f(file_id, hdferr)
!     print *, 'hanging at 7'
     CALL h5pclose_f(plist_file, hdferr)
!     print *, 'hanging at 8'
     CALL h5pclose_f(plist_cdrw, hdferr)

   end subroutine WriteHDF5_3D_1V
   !---------------------------------------------------


     ! parallel version
     ! read timeseries GridMetrics file "timeseries_GridMetrics.h5"
     ! file name: 'fname'. group name: 'gname'.
     ! dim1, dim2, dim3 corresponding to the dimf in each partition
     ! dim1st, dim2st, dim3st corresponding to the begin index number in each direction

     subroutine ReadtsGridMetrics_P(fname,gname,dim1,dim2,dim3,dim1st,dim2st,dim3st,vars)
       character(*), intent(in) :: fname
       character(*), intent(in) :: gname
       integer, intent(in) :: dim1, dim2, dim3
       integer, intent(in) :: dim1st, dim2st, dim3st
       real(8), intent(out) :: vars(:,:,:,:)
       character(10), dimension(12) :: dname
       integer(HSIZE_T), dimension(3) :: dimsf
       integer(HSIZE_T), dimension(3) :: count, offset, stride, block
       integer :: i, j, k
       integer :: rank, num_dset

       if(.not.IsHDF5Initialized) then
         print *, 'HDF5 should be Initialized'
         stop
       endif
       if(dim1.ne.size(vars,dim=1).or.dim2.ne.size(vars,dim=2).or.dim3.ne.size(vars,dim=3)) then
         print *, 'The input dimension (dim1,dim2,dim3) should be the same with the dimension of vars'
         stop
       endif

       rank = 3
       num_dset = 12

       dimsf(1) = dim1
       dimsf(2) = dim2
       dimsf(3) = dim3

       dname(1) = "x"
       dname(2) = "y"
       dname(3) = "z"
       dname(4) = "didx"
       dname(5) = "djdx"
       dname(6) = "dkdx"
       dname(7) = "didy"
       dname(8) = "djdy"
       dname(9) = "dkdy"
       dname(10) = "didz"
       dname(11) = "djdz"
       dname(12) = "dkdz"

       count(1) = dimsf(1)
       count(2) = dimsf(2)
       count(3) = dimsf(3)
       offset(1) = dim1st - 1
       offset(2) = dim2st - 1
       offset(3) = dim3st - 1
!       stride(1) = 1
!       stride(2) = 1
!       stride(3) = 1
!       block(1) = dim1
!       block(2) = dim2
!       block(3) = dim3

       CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdferr)
       CALL h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, hdferr)
       CALL h5fopen_f(fname,H5F_ACC_RDONLY_F, file_id, hdferr, access_prp = plist_id )
       CALL h5pclose_f(plist_id, hdferr)

       CALL h5gopen_f(file_id,trim(gname),group_id,hdferr)
       CALL h5screate_simple_f(rank, dimsf, memspace, hdferr)
       do i=1, num_dset
         CALL h5dopen_f(group_id, dname(i), dset_id, hdferr)
         CALL h5dget_space_f(dset_id, dspace_id, hdferr)
         CALL h5sselect_hyperslab_f (dspace_id, H5S_SELECT_SET_F, offset, count, hdferr) !, stride, block)

         CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdferr)
         CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdferr)
         CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE,vars(:,:,:,i), dimsf, hdferr, &
                        file_space_id = dspace_id, mem_space_id = memspace, xfer_prp = plist_id)
         CALL h5sclose_f(dspace_id, hdferr)
         CALL h5dclose_f(dset_id, hdferr)
         CALL h5pclose_f(plist_id, hdferr)
       enddo
       CALL h5sclose_f(memspace, hdferr)
       CALL h5gclose_f(group_id, hdferr)
       CALL h5fclose_f(file_id, hdferr)


     end subroutine ReadtsGridMetrics_P

     ! parallel version
     ! read timeseries flow data.
     ! file name: fname. group name: gname
     ! nt for # of time steps
     ! dim1, dim2, dim3 for
     ! n: 1 for kplane, 2 for iplane, 3 for jplane
     subroutine Readtsflow_P(fname,gname,nt,dim1,dim2,dim3,dim1st,dim2st,dim3st,vars,n)
       character(*), intent(in) :: fname
       character(*), intent(in) :: gname
       integer, intent(in) :: nt, dim1, dim2, dim3, n
       integer, intent(in) :: dim1st, dim2st, dim3st
       real(8), intent(out) :: vars(:,:,:,:,:)
       character(10), dimension(10) :: dname
       integer(HSIZE_T), dimension(4) :: dimsf
       integer(HSIZE_T), dimension(4) :: count, offset, stride, block
       integer :: i, j, k
       integer :: rank, num_dset

       if(nt.ne.size(vars,dim=1).or.dim1.ne.size(vars,dim=2).or.dim2.ne.size(vars,dim=3).or.dim3.ne.size(vars,dim=4)) then
         print *, 'The input dimension (nt,dim1,dim2,dim3) should be the same with the dimension of vars'
         stop
       endif
       rank = 4
       num_dset = 10

       dimsf(1) = nt
       dimsf(2) = dim1
       dimsf(3) = dim2
       dimsf(4) = dim3
       count(1) = dimsf(1)
       count(2) = dimsf(2)
       count(3) = dimsf(3)
       count(4) = dimsf(4)
       offset(1) = 0
       offset(2) = dim1st - 1
       offset(3) = dim2st - 1
       offset(4) = dim3st - 1
       dname(1) = "u"
       dname(2) = "v"
       dname(3) = "w"
       dname(4) = "p"
       dname(5) = "T"
       if(n.eq.1) then
         dname(6) = "uk"
         dname(7) = "vk"
         dname(8) = "wk"
         dname(9) = "pk"
         dname(10) = "Tk"
       elseif(n.eq.2) then
         dname(6) = "ui"
         dname(7) = "vi"
         dname(8) = "wi"
         dname(9) = "pi"
         dname(10) = "Ti"
       elseif(n.eq.3) then
         dname(6) = "uj"
         dname(7) = "vj"
         dname(8) = "wj"
         dname(9) = "pj"
         dname(10) = "Tj"
       else
         print *, ' n should be 1: kplane, 2: iplane, 3: iplane'
         stop
       endif


       CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdferr)
       CALL h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, hdferr)
       CALL h5fopen_f(fname,H5F_ACC_RDONLY_F, file_id, hdferr, access_prp = plist_id )
       CALL h5pclose_f(plist_id, hdferr)

       CALL h5gopen_f(file_id,trim(gname),group_id,hdferr)
       CALL h5screate_simple_f(rank, dimsf, memspace, hdferr)
       do i=1, num_dset
         CALL h5dopen_f(group_id, dname(i), dset_id, hdferr)
         CALL h5dget_space_f(dset_id, dspace_id, hdferr)
         CALL h5sselect_hyperslab_f (dspace_id, H5S_SELECT_SET_F, offset, count, hdferr) !, stride, block)

         CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdferr)
         CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdferr)
         CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE,vars(:,:,:,:,i), dimsf, hdferr, &
                        file_space_id = dspace_id, mem_space_id = memspace, xfer_prp = plist_id)
         CALL h5sclose_f(dspace_id, hdferr)
         CALL h5dclose_f(dset_id, hdferr)
         CALL h5pclose_f(plist_id, hdferr)
       enddo
       CALL h5sclose_f(memspace, hdferr)
       CALL h5gclose_f(group_id, hdferr)
       CALL h5fclose_f(file_id, hdferr)

     end subroutine Readtsflow_P

     subroutine Writetsplanes_P(fname,dim1,dim2,dim3,dim1st,dim2st,dim3st,vars)
       character(*), intent(in) :: fname
       integer, intent(in) :: dim1, dim2, dim3
       integer, intent(in) :: dim1st, dim2st, dim3st
       real(8), intent(in) :: vars(:,:,:)
       character(10), dimension(10) :: dname
       integer(HSIZE_T), dimension(3) :: dimsf, dimsf_sub
       integer(HSIZE_T), dimension(3) :: count, offset, stride, block
       integer :: i, j, k
       integer :: rank, num_dset

       rank = 3
       num_dset = 1
       dname(1) = 'grad_rho'

       dimsf(1) = dim1
       dimsf(2) = dim2
       dimsf(3) = dim3
       dimsf_sub(1) = size(vars,1)
       dimsf_sub(2) = size(vars,2)
       dimsf_sub(3) = size(vars,3)
       count(1) = dimsf_sub(1)
       count(2) = dimsf_sub(2)
       count(3) = dimsf_sub(3)
       offset(1) = dim1st - 1
       offset(2) = dim2st - 1
       offset(3) = dim3st - 1

!       print *, 'dim1 = ', dimsf_sub(1)
!       print *, 'dim2 = ', dimsf_sub(2)
!       print *, 'dim3 = ', dimsf_sub(3)
!       print *, 'count(1) = ', count(1)
!       print *, 'count(2) = ', count(2)
!       print *, 'count(3) = ', count(3)
!       print *, 'offset(1) = ', offset(1)
!       print *, 'offset(2) = ', offset(2)
!       print *, 'offset(3) = ', offset(3)

       CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdferr)
       CALL h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, hdferr)
       CALL h5fcreate_f(fname, H5F_ACC_TRUNC_F, file_id, hdferr, &
                      access_prp = plist_id)
       CALL h5pclose_f(plist_id, hdferr)
       CALL h5screate_simple_f(rank, dimsf_sub, memspace, hdferr)

       do i=1, num_dset
         CALL h5screate_simple_f(rank, dimsf, dspace_id, hdferr)
         CALL h5dcreate_f(file_id, dname(i), H5T_NATIVE_DOUBLE, dspace_id, &
                      dset_id, hdferr)
         CALL h5sclose_f(dspace_id, hdferr)
         CALL h5dget_space_f(dset_id, dspace_id, hdferr)
         CALL h5sselect_hyperslab_f (dspace_id, H5S_SELECT_SET_F, offset, count, hdferr)
         CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdferr)
         CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdferr)
         CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, vars, dimsf_sub, hdferr, &
                         file_space_id = dspace_id, mem_space_id = memspace, xfer_prp = plist_id)
         CALL h5sclose_f(dspace_id, hdferr)
         CALL h5pclose_f(plist_id, hdferr)
         CALL h5dclose_f(dset_id, hdferr)
       enddo

       CALL h5sclose_f(memspace, hdferr)
       CALL h5fclose_f(file_id, hdferr)

     end subroutine Writetsplanes_P


     ! parallel version
     ! write grid.h5
     ! idim, jdim, kdim for the dimension of x, y, z
     ! ist, jst, kst for the offset
     subroutine WriteHDF5grid_P(fname,kdim,idim,jdim,kst,ist,jst,x,y,z)
       character(*), intent(in) :: fname
       integer, intent(in) :: kdim, idim, jdim
       integer, intent(in) :: kst, ist, jst
       !real, intent(out) :: x(1:kdim,1:idim,1:jdim), y(1:kdim,1:idim,1:jdim), z(1:kdim,1:idim,1:jdim)
       real, intent(out) :: x(:,:,:), y(:,:,:), z(:,:,:)
       integer(HSIZE_T), dimension(3) :: dimsf, dimsf_sub
       integer(HSIZE_T), dimension(3) :: count, offset, stride, block
       character(10), dimension(:), allocatable :: dname
       integer :: rank, num_dset
       integer :: i, j, k

       if(.not.IsHDF5Initialized) then
         print *, 'HDF5 should be Initialized'
         stop
       endif

       rank = 3
       num_dset = 3
       allocate(dname(num_dset))

       dname(1) = "x"
       dname(2) = "y"
       dname(3) = "z"
       dimsf(1) = kdim
       dimsf(2) = idim
       dimsf(3) = jdim

       dimsf_sub(1) = size(x,1)
       dimsf_sub(2) = size(x,2)
       dimsf_sub(3) = size(x,3)

       !block(1) = dimsf(1)
       !block(2) = dimsf(2)
       !block(3) = dimsf(3)
       !stride(1) = dimsf_sub(1)
       !stride(2) = dimsf_sub(1)
       !stride(3) = dimsf_sub(1)
       count(1)  = dimsf_sub(1)
       count(2)  = dimsf_sub(2)
       count(3)  = dimsf_sub(3)
       offset(1) = kst - 1
       offset(2) = ist - 1
       offset(3) = jst - 1


       CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdferr)
       CALL h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, hdferr)
       CALL h5fcreate_f(trim(fname),H5F_ACC_TRUNC_F, file_id, hdferr, access_prp=plist_id)
       !CALL h5fopen_f(trim(fname), H5F_ACC_RDONLY_F, file_id, hdferr, access_prp=plist_id)
       CALL h5pclose_f(plist_id, hdferr)
       CALL h5screate_simple_f(rank, dimsf_sub, memspace, hdferr)

         ! write x
       !do i=1, num_dset
         CALL h5screate_simple_f(rank, dimsf, dspace_id, hdferr)
         CALL h5dcreate_f(file_id, dname(1), H5T_NATIVE_DOUBLE, dspace_id, &
                      dset_id, hdferr)
         CALL h5sclose_f(dspace_id, hdferr)
         CALL h5dget_space_f(dset_id, dspace_id, hdferr)
         CALL h5sselect_hyperslab_f (dspace_id, H5S_SELECT_SET_F, offset, count, hdferr)
         CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdferr)
         CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdferr)
         CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, x, dimsf_sub, hdferr, &
                         file_space_id = dspace_id, mem_space_id = memspace, xfer_prp = plist_id)
         CALL h5sclose_f(dspace_id, hdferr)
         CALL h5pclose_f(plist_id, hdferr)
         CALL h5dclose_f(dset_id, hdferr)
       !enddo
         CALL h5screate_simple_f(rank, dimsf, dspace_id, hdferr)
         CALL h5dcreate_f(file_id, dname(2), H5T_NATIVE_DOUBLE, dspace_id, &
                      dset_id, hdferr)
         CALL h5sclose_f(dspace_id, hdferr)
         CALL h5dget_space_f(dset_id, dspace_id, hdferr)
         CALL h5sselect_hyperslab_f (dspace_id, H5S_SELECT_SET_F, offset, count, hdferr)
         CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdferr)
         CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdferr)
         CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, y, dimsf_sub, hdferr, &
                         file_space_id = dspace_id, mem_space_id = memspace, xfer_prp = plist_id)
         CALL h5sclose_f(dspace_id, hdferr)
         CALL h5pclose_f(plist_id, hdferr)
         CALL h5dclose_f(dset_id, hdferr)

         CALL h5screate_simple_f(rank, dimsf, dspace_id, hdferr)
         CALL h5dcreate_f(file_id, dname(3), H5T_NATIVE_DOUBLE, dspace_id, &
                      dset_id, hdferr)
         CALL h5sclose_f(dspace_id, hdferr)
         CALL h5dget_space_f(dset_id, dspace_id, hdferr)
         CALL h5sselect_hyperslab_f (dspace_id, H5S_SELECT_SET_F, offset, count, hdferr)
         CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdferr)
         CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdferr)
         CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, z, dimsf_sub, hdferr, &
                         file_space_id = dspace_id, mem_space_id = memspace, xfer_prp = plist_id)
         CALL h5sclose_f(dspace_id, hdferr)
         CALL h5pclose_f(plist_id, hdferr)
         CALL h5dclose_f(dset_id, hdferr)

       CALL h5sclose_f(memspace, hdferr)
       CALL h5fclose_f(file_id, hdferr)

     end subroutine WriteHDF5grid_P


    ! parallel version
     ! write flowdata_xxxxxxxx.h5
     ! idim, jdim, kdim for the dimension of u, v, w, p, T
     ! ist, jst, kst for the offset
     subroutine WriteHDF5flow_P(fname,kdim,idim,jdim,kst,ist,jst,buffer)
       character(*), intent(in) :: fname
       integer, intent(in) :: kdim, idim, jdim
       integer, intent(in) :: kst, ist, jst
       real, intent(in) :: buffer(:,:,:,:)
       integer(HSIZE_T), dimension(3) :: dimsf, dimsf_sub
       integer(HSIZE_T), dimension(3) :: count, offset, stride, block
       character(10), dimension(:), allocatable :: dname
       integer :: rank, num_dset
       integer :: i, j, k

       !real(8) :: time = 0.d0

       if(.not.IsHDF5Initialized) then
         print *, 'HDF5 should be Initialized'
         stop
       endif

       rank = 3
       num_dset = 5
       allocate(dname(num_dset))

       dname(1) = "u"
       dname(2) = "v"
       dname(3) = "w"
       dname(4) = "p"
       dname(5) = "T"
       dimsf(1) = kdim
       dimsf(2) = idim
       dimsf(3) = jdim

       dimsf_sub(1) = size(buffer,1)
       dimsf_sub(2) = size(buffer,2)
       dimsf_sub(3) = size(buffer,3)

       !block(1) = dimsf(1)
       !block(2) = dimsf(2)
       !block(3) = dimsf(3)
       !stride(1) = dimsf_sub(1)
       !stride(2) = dimsf_sub(1)
       !stride(3) = dimsf_sub(1)
       count(1)  = dimsf_sub(1)
       count(2)  = dimsf_sub(2)
       count(3)  = dimsf_sub(3)
       offset(1) = kst - 1
       offset(2) = ist - 1
       offset(3) = jst - 1


       CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdferr)
       CALL h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, hdferr)
       CALL h5fcreate_f(trim(fname),H5F_ACC_TRUNC_F, file_id, hdferr, access_prp=plist_id)
       CALL h5pclose_f(plist_id, hdferr)
       CALL h5screate_simple_f(rank, dimsf_sub, memspace, hdferr)

       ! write variables
       do i=1, num_dset
         CALL h5screate_simple_f(rank, dimsf, dspace_id, hdferr)
         CALL h5dcreate_f(file_id, dname(i), H5T_NATIVE_DOUBLE, dspace_id, &
                      dset_id, hdferr)
         CALL h5sclose_f(dspace_id, hdferr)
         CALL h5dget_space_f(dset_id, dspace_id, hdferr)
         CALL h5sselect_hyperslab_f (dspace_id, H5S_SELECT_SET_F, offset, count, hdferr)
         CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdferr)
         CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdferr)
         CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, buffer(:,:,:,i), dimsf_sub, hdferr, &
                         file_space_id = dspace_id, mem_space_id = memspace, xfer_prp = plist_id)
         CALL h5sclose_f(dspace_id, hdferr)
         CALL h5pclose_f(plist_id, hdferr)
         CALL h5dclose_f(dset_id, hdferr)
       enddo

       !CALL h5screate_f(H5S_SCALAR_F, dspace_id, hdferr)
       !CALL h5dcreate_f(file_id, trim("time"),H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)
       !CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, time, dimsf, hdferr)
       !CALL h5dclose_f(dset_id, hdferr)
       !CALL h5sclose_f(dspace_id, hdferr)

       CALL h5sclose_f(memspace, hdferr)
       CALL h5fclose_f(file_id, hdferr)

     end subroutine WriteHDF5flow_P

     ! parallel version
     ! read grid.h5
     ! idim, jdim, kdim for the dimension of x, y, z
     ! ist, jst, kst for the offset
     subroutine ReadHDF5grid_P(fname,kdim,idim,jdim,kst,ist,jst,x,y,z)
       character(*), intent(in) :: fname
       integer, intent(in) :: kdim, idim, jdim
       integer, intent(in) :: kst, ist, jst
       real, intent(out) :: x(1:kdim,1:idim,1:jdim), y(1:kdim,1:idim,1:jdim), z(1:kdim,1:idim,1:jdim)
       integer(HSIZE_T), dimension(3) :: dimsf
       integer(HSIZE_T), dimension(3) :: count, offset, stride, block
       character(10), dimension(:), allocatable :: dname
       integer :: rank, num_dset
       integer :: i, j, k

       if(.not.IsHDF5Initialized) then
         print *, 'HDF5 should be Initialized'
         stop
       endif

       rank = 3
       num_dset = 3
       allocate(dname(num_dset))

       dname(1) = "x"
       dname(2) = "y"
       dname(3) = "z"
       dimsf(1) = kdim
       dimsf(2) = idim
       dimsf(3) = jdim

       block(1) = dimsf(1)
       block(2) = dimsf(2)
       block(3) = dimsf(3)
       stride(1) = 1
       stride(2) = 1
       stride(3) = 1
       count(1)  = 1
       count(2)  = 1
       count(3)  = 1
       offset(1) = kst - 1
       offset(2) = ist - 1
       offset(3) = jst - 1


       CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdferr)
       CALL h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, hdferr)
       CALL h5fopen_f(trim(fname), H5F_ACC_RDONLY_F, file_id, hdferr, access_prp=plist_id)
       CALL h5pclose_f(plist_id, hdferr)
       CALL h5screate_simple_f(rank, dimsf, memspace, hdferr)

         ! read x
         CALL h5dopen_f(file_id, dname(1), dset_id, hdferr)
         CALL h5dget_space_f(dset_id, dspace_id, hdferr)
         CALL h5sselect_hyperslab_f (dspace_id, H5S_SELECT_SET_F, offset, count, hdferr, stride, block)
         CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdferr)
         CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdferr)
         CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, x, dimsf, hdferr, &
                        file_space_id=dspace_id, mem_space_id=memspace, xfer_prp=plist_id)
         CALL h5sclose_f(dspace_id, hdferr)
         CALL h5dclose_f(dset_id, hdferr)
         CALL h5pclose_f(plist_id, hdferr)

         ! read y
         CALL h5dopen_f(file_id, dname(2), dset_id, hdferr)
         CALL h5dget_space_f(dset_id, dspace_id, hdferr)
         CALL h5sselect_hyperslab_f (dspace_id, H5S_SELECT_SET_F, offset, count, hdferr, stride, block)
         CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdferr)
         CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdferr)
         CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, y(:,:,:), dimsf, hdferr, &
                        file_space_id=dspace_id, mem_space_id=memspace, xfer_prp=plist_id)
         CALL h5sclose_f(dspace_id, hdferr)
         CALL h5dclose_f(dset_id, hdferr)
         CALL h5pclose_f(plist_id, hdferr)
         ! read z
         CALL h5dopen_f(file_id, dname(3), dset_id, hdferr)
         CALL h5dget_space_f(dset_id, dspace_id, hdferr)
         CALL h5sselect_hyperslab_f (dspace_id, H5S_SELECT_SET_F, offset, count, hdferr, stride, block)
         CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdferr)
         CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdferr)
         CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, z(:,:,:), dimsf, hdferr, &
                        file_space_id=dspace_id, mem_space_id=memspace, xfer_prp=plist_id)
         CALL h5sclose_f(dspace_id, hdferr)
         CALL h5dclose_f(dset_id, hdferr)
         CALL h5pclose_f(plist_id, hdferr)

       CALL h5sclose_f(memspace, hdferr)
       CALL h5fclose_f(file_id, hdferr)

     end subroutine ReadHDF5grid_P

     !
     subroutine WriteHDF51D_p(hslab,buffer)
       implicit none
       type(tp_rdwt_hdf5), intent(in) :: hslab
       real(8), intent(in) :: buffer(hslab%dimsm(1),hslab%dnum)
       integer :: n

       ! Setup file access property list with parallel I/O access.
       CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_file, hdferr)
       CALL h5pset_fapl_mpio_f(plist_file, hslab%comm, hslab%info, hdferr)
       ! Create property list for collective dataset read/write
       CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_cdrw, hdferr)
       CALL h5pset_dxpl_mpio_f(plist_cdrw, H5FD_MPIO_COLLECTIVE_F, hdferr)

       CALL h5fcreate_f(trim(hslab%fname), H5F_ACC_TRUNC_F, file_id, hdferr,access_prp = plist_file)
       CALL h5screate_simple_f(hslab%rank, hslab%dimsm, memspace, hdferr)
       CALL h5screate_simple_f(hslab%rank, hslab%dimsf, dspace_id, hdferr)
       CALL h5gcreate_f(file_id, trim(hslab%gname),group_id, hdferr)
       do n = 1, hslab%dnum
         CALL h5dcreate_f(group_id, trim(hslab%dname(n)), &
                          H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)

         CALL h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, hslab%offset, &
                                    hslab%count, hdferr )
         CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, buffer(:,n),  &
                         hslab%dimsf, hdferr, file_space_id=dspace_id, mem_space_id=memspace, xfer_prp = plist_cdrw)
         CALL h5dclose_f(dset_id, hdferr)
       enddo
       CALL h5sclose_f(dspace_id, hdferr)
       CALL h5sclose_f(memspace, hdferr)
       CALL h5gclose_f(group_id, hdferr)
       CALL h5fclose_f(file_id, hdferr)
       CALL h5pclose_f(plist_file, hdferr)
       CALL h5pclose_f(plist_cdrw, hdferr)


     end subroutine WriteHDF51D_p



     ! parallel version
     ! kdim, idim, jdim for the dimension of vars
     ! kst, ist, jst for the offset
     subroutine ReadHDF5sol_P(fname,kdim,idim,jdim,kst,ist,jst,vars,time)
       character(*), intent(in) :: fname
       integer, intent(in) :: kdim, idim, jdim
       integer, intent(in) :: kst, ist, jst
       real(8), intent(out) :: vars(:,:,:,:)
       real(8), intent(out), optional :: time
       integer(HSIZE_T), dimension(:), allocatable :: dimsf
       integer(HSIZE_T), dimension(3) :: count, offset, stride, block
       character(10), dimension(:), allocatable :: dname
       integer :: rank, num_dset
       integer :: i, j, k

       if(.not.IsHDF5Initialized) then
         print *, 'HDF5 should be Initialized'
         stop
       endif

       rank = 3
       num_dset = 6
       allocate(dname(num_dset))
       allocate(dimsf(rank))

       dname(1) = "u"
       dname(2) = "v"
       dname(3) = "w"
       dname(4) = "p"
       dname(5) = "T"
       dname(6) = "time"
       dimsf(1) = kdim
       dimsf(2) = idim
       dimsf(3) = jdim
       count(1)  = dimsf(1)
       count(2)  = dimsf(2)
       count(3)  = dimsf(3)
       offset(1) = kst - 1
       offset(2) = ist - 1
       offset(3) = jst - 1
       CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdferr)
       CALL h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, hdferr)
       CALL h5fopen_f(trim(fname), H5F_ACC_RDONLY_F, file_id, hdferr, access_prp=plist_id)
       CALL h5pclose_f(plist_id, hdferr)
       CALL h5screate_simple_f(rank, dimsf, memspace, hdferr)
!    print *, 'hdferr = ', hdferr
       do i=1, num_dset - 1
         CALL h5dopen_f(file_id, dname(i), dset_id, hdferr)
         CALL h5dget_space_f(dset_id, dspace_id, hdferr)
         CALL h5sselect_hyperslab_f (dspace_id, H5S_SELECT_SET_F, offset, count, hdferr)
         CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdferr)
         CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdferr)
         CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, vars(:,:,:,i), dimsf, hdferr, &
                        file_space_id=dspace_id, mem_space_id=memspace, xfer_prp=plist_id)
         CALL h5sclose_f(dspace_id, hdferr)
         CALL h5dclose_f(dset_id, hdferr)
         CALL h5pclose_f(plist_id, hdferr)
       enddo

       if(present(time)) then
         CALL h5dopen_f(file_id, trim(dname(6)), dset_id, hdferr)
         CALL h5dget_space_f(dset_id, dspace_id, hdferr)
         CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, time, dimsf, hdferr)
         CALL h5sclose_f(dspace_id, hdferr)
         CALL h5dclose_f(dset_id, hdferr)
       endif
       CALL h5sclose_f(memspace, hdferr)
       CALL h5fclose_f(file_id, hdferr)

     end subroutine ReadHDF5sol_P

     subroutine Read3dHDF5_svariable_P(fname,dname,dim1,dim2,dim3,dim1st,dim2st,dim3st,vars)
       implicit none
       character(*), intent(in) :: fname
       character(*), intent(in) :: dname
       integer, intent(in) :: dim1, dim2, dim3
       integer, intent(in) :: dim1st, dim2st, dim3st
       real(8), intent(out) :: vars(:,:,:)
       integer(HSIZE_T), dimension(3) :: dimsf
       integer(HSIZE_T), dimension(3) :: count, offset, stride, block
       integer :: rank=3

       if(.not.IsHDF5Initialized) then
         print *, 'HDF5 should be Initialized'
         stop
       endif
       if(dim1.ne.size(vars,dim=1).or.dim2.ne.size(vars,dim=2).or.dim3.ne.size(vars,dim=3)) then
         print *, 'The input dimension (dim1,dim2,dim3) should be the same with the dimension of vars'
         stop
       endif

       dimsf = (/dim1,dim2,dim3/)
       block = dimsf
       stride = 1
       count = 1
       offset(1) = dim1st - 1
       offset(2) = dim2st - 1
       offset(3) = dim3st - 1

       CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdferr)
       CALL h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, hdferr)
       CALL h5fopen_f(trim(fname), H5F_ACC_RDONLY_F, file_id, hdferr, access_prp=plist_id)
       CALL h5pclose_f(plist_id, hdferr)
       CALL h5screate_simple_f(rank, dimsf, memspace, hdferr)

       ! read data
       CALL h5dopen_f(file_id, trim(dname), dset_id, hdferr)
       CALL h5dget_space_f(dset_id, dspace_id, hdferr)
       CALL h5sselect_hyperslab_f (dspace_id, H5S_SELECT_SET_F, offset, count, hdferr, stride, block)
       CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdferr)
       CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdferr)
       CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, vars, dimsf, hdferr, &
                      file_space_id=dspace_id, mem_space_id=memspace, xfer_prp=plist_id)
       CALL h5sclose_f(dspace_id, hdferr)
       CALL h5dclose_f(dset_id, hdferr)
       CALL h5pclose_f(plist_id, hdferr)

       CALL h5sclose_f(memspace, hdferr)
       CALL h5fclose_f(file_id, hdferr)

     end subroutine Read3dHDF5_svariable_P


      ! write 3d single variables
      subroutine Write3dHDF5_svariable_P(fname,dname,dim1,dim2,dim3,dim1st,dim2st,dim3st,vars,iexist)
        implicit none
        character(*), intent(in) :: fname
        character(*), intent(in) :: dname
        integer, intent(in) :: dim1, dim2, dim3
        integer, intent(in) :: dim1st, dim2st, dim3st
        real, intent(in) :: vars(:,:,:)
        integer, intent(in) :: iexist
        integer(HSIZE_T), dimension(3) :: dimsf, dimsf_sub
        integer(HSIZE_T), dimension(3) :: count, offset, stride, block
        integer :: rank = 3

        if(.not.IsHDF5Initialized) then
          print *, 'HDF5 should be Initialized'
          stop
        endif

        dimsf = (/dim1,dim2,dim3/)
        dimsf_sub(1) = size(vars,1)
        dimsf_sub(2) = size(vars,2)
        dimsf_sub(3) = size(vars,3)

        block(1) = dimsf_sub(1)
        block(2) = dimsf_sub(2)
        block(3) = dimsf_sub(3)
        count = 1
        stride = 1
        offset(1) = dim1st - 1
        offset(2) = dim2st - 1
        offset(3) = dim3st - 1

        CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdferr)
        CALL h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, hdferr)

        if(iexist.eq.0) then
          CALL h5fcreate_f(trim(fname), H5F_ACC_TRUNC_F, file_id, hdferr, access_prp = plist_id)  
    
        else
          CALL h5fopen_f(trim(fname), H5F_ACC_RDWR_F, file_id, hdferr, access_prp = plist_id )
        endif
        
        CALL h5pclose_f(plist_id, hdferr)
        CALL h5screate_simple_f(rank, dimsf_sub, memspace, hdferr)        

        ! write data
        CALL h5screate_simple_f(rank, dimsf, dspace_id, hdferr)
        CALL h5dcreate_f(file_id, trim(dname), H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)      
        CALL h5sclose_f(dspace_id, hdferr)        
        CALL h5dget_space_f(dset_id, dspace_id, hdferr)        
        CALL h5sselect_hyperslab_f (dspace_id, H5S_SELECT_SET_F, offset, count, hdferr, &
                                    stride, block )                                 

        CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdferr)
        CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdferr)
        CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, vars, dimsf_sub, hdferr, &
                        file_space_id = dspace_id, mem_space_id = memspace, xfer_prp=plist_id)

        CALL h5sclose_f(dspace_id, hdferr)
        CALL h5pclose_f(plist_id, hdferr)
        CALL h5dclose_f(dset_id, hdferr)

        CALL h5sclose_f(memspace, hdferr) 

        CALL h5fclose_f(file_id, hdferr) 

      end subroutine Write3dHDF5_svariable_P


subroutine Write3dHDF5_svariable_S(hslab, buffer,current_num)
  type(tp_rdwt_hdf5), intent(in) :: hslab
  real(8), intent(in) :: buffer(hslab%dimsf(1),hslab%dimsf(2),hslab%dimsf(3))
  integer, intent(in) :: current_num
  integer :: n
  
  if(.not.(IsHDF5Initialized.and.hslab%IsHSInitialized)) then
      print *, 'HDF5 or type hslab NOT intialized !!! STOP'
      stop
  endif

  if(current_num.eq.1) then
    CALL h5fcreate_f(trim(hslab%fname), H5F_ACC_TRUNC_F, file_id, hdferr)
  else
    CALL h5fopen_f(trim(hslab%fname),H5F_ACC_RDWR_F, file_id, hdferr)
  endif

  CALL h5screate_simple_f(hslab%rank, hslab%dimsf, dspace_id, hdferr)
  do n=current_num, current_num !????????????????????????
    CALL h5dcreate_f(file_id, trim(hslab%dname(n)), H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr )
    CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, buffer(:,:,:),hslab%dimsf, hdferr)
    CALL h5dclose_f(dset_id, hdferr)
  enddo
  CALL h5sclose_f(dspace_id, hdferr)
  CALL h5fclose_f(file_id, hdferr)

end subroutine Write3dHDF5_svariable_S


!relocate this later
    subroutine WriteHDF5_3D_S(hslab, buffer)
      type(tp_rdwt_hdf5), intent(in) :: hslab
      real(8), intent(in) :: buffer(hslab%dimsf(1),hslab%dimsf(2),hslab%dimsf(3),hslab%dnum)
      integer :: n

      if(.not.(IsHDF5Initialized.and.hslab%IsHSInitialized)) then
          print *, 'HDF5 or type hslab NOT intialized !!! STOP'
          stop
      endif

      if(.not.hslab%IsMultiGroup) then
        CALL h5fcreate_f(trim(hslab%fname), H5F_ACC_TRUNC_F, file_id, hdferr)
      else
        CALL h5fopen_f(trim(hslab%fname),H5F_ACC_RDWR_F, file_id, hdferr)
      endif

      if(trim(hslab%gname).ne.'/') then
        CALL h5gcreate_f(file_id, trim(hslab%gname),group_id, hdferr)
      endif
      CALL h5screate_simple_f(hslab%rank, hslab%dimsf, dspace_id, hdferr)
      do n=1, hslab%dnum
        if(trim(hslab%gname).ne.'/') then
          CALL h5dcreate_f(group_id, trim(hslab%dname(n)), H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr )
        else
          CALL h5dcreate_f(file_id, trim(hslab%dname(n)), H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr )
        endif
        CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, buffer(:,:,:,n),hslab%dimsf, hdferr)
        CALL h5dclose_f(dset_id, hdferr)
      enddo
      CALL h5sclose_f(dspace_id, hdferr)
      if(trim(hslab%gname).ne.'/') then
        CALL h5gclose_f(group_id, hdferr)
      endif
      CALL h5fclose_f(file_id, hdferr)

    end subroutine WriteHDF5_3D_S


!      ! write 3d single variables
!      subroutine Write3dHDF5_svariable_P(fname,dims,dims_total,dim1st,dim2st,dim3st,dname,vars,iexist)
!        character(*), intent(in) :: fname
!        integer, intent(in) :: dims(3), dims_total(3)
!        integer, intent(in) :: dim1st, dim2st, dim3st
!        character(*), intent(in) :: dname
!        real, intent(in) :: vars(1:dims(1),1:dims(2),1:dims(3))
!        integer, intent(in) :: iexist
!        integer(HSIZE_T), dimension(3) :: dimsf, dimsf_sub
!        integer(HSIZE_T), dimension(3) :: count, offset, stride, block
!        integer :: rank = 3
!
!        dimsf = dims_total
!        dimsf_sub = dims
!
!        ! check dimsf
!       ! print *, 'dimsf_sub = ', dimsf_sub
!       ! print *, 'dimsf = ', dimsf
!
!        count(1) = 1
!        count(2) = 1
!        count(3) = 1
!        block(1) = dimsf_sub(1)
!        block(2) = dimsf_sub(2)
!        block(3) = dimsf_sub(3)
!        stride(1) = 1
!        stride(2) = 1
!        stride(3) = 1
!        offset(1) = dim1st - 1
!        offset(2) = dim2st - 1
!        offset(3) = dim3st - 1
!
!        ! check info
!        !CALL h5fcreate_f(trim(fname), H5F_ACC_TRUNC_F, file_id, hdferr)
!        !CALL h5screate_simple_f(rank, dimsf, dspace_id, hdferr)
!        !CALL h5dcreate_f(file_id, trim(dname), H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)
!        !CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, vars, dimsf, hdferr)
!        !CALL h5sclose_f(dspace_id, hdferr)
!        !CALL h5dclose_f(dset_id, hdferr)
!        !CALL h5fclose_f(file_id, hdferr)
!
!
!        CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdferr)
!        CALL h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, hdferr)
!        if(iexist.eq.0) then
!          CALL h5fcreate_f(trim(fname), H5F_ACC_TRUNC_F, file_id, hdferr, access_prp = plist_id)
!        else
!          CALL h5fopen_f(trim(fname), H5F_ACC_RDWR_F, file_id, hdferr, access_prp = plist_id )
!        endif
!        CALL h5pclose_f(plist_id, hdferr)
!
!        CALL h5screate_simple_f(rank, dimsf_sub, memspace, hdferr)
!
!        CALL h5screate_simple_f(rank, dimsf, dspace_id, hdferr)
!        CALL h5dcreate_f(file_id, trim(dname), H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)
!        CALL h5sclose_f(dspace_id, hdferr)
!        CALL h5dget_space_f(dset_id, dspace_id, hdferr)
!        CALL h5sselect_hyperslab_f (dspace_id, H5S_SELECT_SET_F, offset, count, hdferr, &
!                                    stride, block )
!        CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdferr)
!        CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdferr)
!
!        ! check info
!        ! print *, 'dimsf_sub = ', dimsf_sub
!        ! print *, 'size(vars,1) = ', size(vars,1)
!        ! print *, 'size(vars,2) = ', size(vars,2)
!        ! print *, 'size(vars,3) = ', size(vars,3)
!        ! print *, 'vars(1,2,65) = ', vars(1,2,65)
!
!        CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, vars(1:dims(1),1:dims(2),1:dims(3)), dimsf_sub, hdferr, &
!                        file_space_id = dspace_id, mem_space_id = memspace, xfer_prp=plist_id)
!        CALL h5sclose_f(dspace_id, hdferr)
!        CALL h5pclose_f(plist_id, hdferr)
!        CALL h5dclose_f(dset_id, hdferr)
!        CALL h5sclose_f(memspace, hdferr)
!        CALL h5fclose_f(file_id, hdferr)
!
!print *, 'finish writing'
!
!      end subroutine Write3dHDF5_svariable_P

! -----------------------------------------------------------------------
    ! initialize the variables in digital filter inlet profile
    subroutine InitHDF5_DigFilter(hslab)
        type(tp_rdwt_hdf5),intent(out) :: hslab

        hslab%gname = '/'
        hslab%rank = 2
        hslab%dnum = 11 ! u,v,w,rho,t,uu,vv,ww,uv,uw,vw
        if(allocated(hslab%dimsf)) deallocate(hslab%dimsf)
        if(allocated(hslab%dname)) deallocate(hslab%dname)
        allocate(hslab%dimsf(hslab%rank),hslab%dname(hslab%dnum))
        hslab%dname(1)  = 'uave'
        hslab%dname(2)  = 'vave'
        hslab%dname(3)  = 'wave'
        hslab%dname(4)  = 'rhoave'
        hslab%dname(5)  = 'tave'
        hslab%dname(6)  = 'uu'
        hslab%dname(7)  = 'vv'
        hslab%dname(8)  = 'ww'
        hslab%dname(9)  = 'uv'
        hslab%dname(10) = 'uw'
        hslab%dname(11) = 'vw'
        hslab%IsHSInitialized = .true.
    end subroutine InitHDF5_DigFilter

    ! read 2D data (serial version)
    subroutine ReadHDF5_2D_serial(hslab, buffer)
        type(tp_rdwt_hdf5), intent(in) :: hslab
        real(8), intent(out) :: buffer(hslab%dimsf(1),hslab%dimsf(2),hslab%dnum)
        integer :: n

        if(.not.(IsHDF5Initialized.and.hslab%IsHSInitialized)) then
          print *, 'HDF5 or type hslab NOT intialized !!! STOP'
          stop
        endif

        CALL h5fopen_f(trim(hslab%fname), H5F_ACC_RDONLY_F, file_id, hdferr)
        CALL h5gopen_f(file_id, trim(hslab%gname), group_id, hdferr)

        do n=1, hslab%dnum
        CALL h5dopen_f(group_id, trim(hslab%dname(n)), dset_id, hdferr)
        CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, buffer(:,:,n), hslab%dimsf, hdferr)
        CALL h5dclose_f(dset_id, hdferr)
        enddo

        CALL h5gclose_f(group_id, hdferr)
        CALL h5fclose_f(file_id, hdferr)
    end subroutine ReadHDF5_2D_serial

    ! Write 2D data
    subroutine WriteHDF5_2D_serial(hslab, buffer)
        type(tp_rdwt_hdf5), intent(in) :: hslab
        real(8), intent(in) :: buffer(hslab%dimsf(1),hslab%dimsf(2),hslab%dnum)
        integer :: n

        if(.not.(IsHDF5Initialized.and.hslab%IsHSInitialized)) then
          print *, 'HDF5 or type hslab NOT intialized !!! STOP'
          stop
        endif

        if(.not.hslab%IsMultiGroup) then
            CALL h5fcreate_f(trim(hslab%fname), H5F_ACC_TRUNC_F, file_id, hdferr)
        else
            CALL h5fopen_f(trim(hslab%fname),H5F_ACC_RDWR_F, file_id, hdferr)
        endif

        if(trim(hslab%gname).ne.'/') then
            CALL h5gcreate_f(file_id, trim(hslab%gname),group_id, hdferr)
        endif
            CALL h5screate_simple_f(hslab%rank, hslab%dimsf, dspace_id, hdferr)
        do n=1, hslab%dnum
            if(trim(hslab%gname).ne.'/') then
              CALL h5dcreate_f(group_id, trim(hslab%dname(n)), H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr )
            else
              CALL h5dcreate_f(file_id, trim(hslab%dname(n)), H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr )
            endif
            CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, buffer(:,:,n),hslab%dimsf, hdferr)
            CALL h5dclose_f(dset_id, hdferr)
        enddo
        CALL h5sclose_f(dspace_id, hdferr)
        if(trim(hslab%gname).ne.'/') then
            CALL h5gclose_f(group_id, hdferr)
        endif
        CALL h5fclose_f(file_id, hdferr)
    end subroutine WriteHDF5_2D_serial

    end module MPRWHDF5
