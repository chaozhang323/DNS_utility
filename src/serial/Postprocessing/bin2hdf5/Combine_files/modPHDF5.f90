

    module modPHDF5
      use HDF5
      use modParallel
      implicit none

      integer :: hdferr
      logical :: IsHDF5Initialized = .false.

      type tp_hyperslab
        logical :: IsMyidIncluded
        integer :: rank
        integer(HID_T) :: file_id, group_id
        integer(HID_T) :: memspace
        integer(HID_T)  :: dset_id, dspace_id
        integer(HID_T) :: plist_file, plist_cdrw
        integer(HSIZE_T), dimension(:), pointer :: count
        integer(HSIZE_T), dimension(:), pointer :: offset, stride, block
      end type tp_hyperslab

      type tp_rdwt_hdf5
        integer :: comm, info
        integer :: rank
        integer :: dnum ! number of variables in the given group with group name 'gname'
        character(250) :: fname
        character(100) :: gname
        character(20), dimension(:), allocatable :: dname ! should have a dimension of 'dnum'
        integer(HSIZE_T), dimension(:), allocatable :: dimsf
        integer(HSIZE_T), dimension(:), allocatable :: dimsm
        integer(HSIZE_T), dimension(:), allocatable :: count, offset, block, stride
        logical :: IsHSInitialized = .false.
        logical :: IsMultiGroup = .false.
        logical :: IsSameGroup =  .false.
        logical :: IsSameData  =  .false.

        character(10) :: sname ! scalar name
        ! for attribute
        character(100), dimension(4) :: attr_name
        integer, dimension(3) :: info_1st_array, info_2nd_array
        integer(HSIZE_T), dimension(1) :: attr_dims_info, attr_dims_loc, attr_dims_no = (/1/)
        integer, dimension(:), pointer :: info_loc
        integer, dimension(1) :: attr_time

      end type tp_rdwt_hdf5


contains

      subroutine InitHDF5()
        CALL h5open_f(hdferr)
        IsHDF5Initialized = .True.
      end subroutine InitHDF5

      subroutine FinalizeHDF5()
        CALL h5close_f(hdferr)
      end subroutine FinalizeHDF5

   ! Read 1D array data
    subroutine ReadHDF5_1D_1v(hslab, buffer)
      type(tp_rdwt_hdf5), intent(in) :: hslab
      real(8), intent(out) :: buffer(hslab%dimsf(1))
      integer(HID_T) :: file_id, group_id, memspace
      integer(HID_T) :: dspace_id, dset_id
      integer(HID_T) :: plist_file, plist_cdrw

      if(.not.(IsHDF5Initialized.and.hslab%IsHSInitialized)) then
          print *, 'HDF5 or type hslab NOT intialized !!! STOP'
          stop
      endif

      ! Setup file access property list with parallel I/O access.
      CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_file, hdferr)
      CALL h5pset_fapl_mpio_f(plist_file, hslab%comm, hslab%info, hdferr)

     ! Create property list for collective dataset read/write
      CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_cdrw, hdferr)
      CALL h5pset_dxpl_mpio_f(plist_cdrw, H5FD_MPIO_COLLECTIVE_F, hdferr)

      CALL h5fopen_f(trim(hslab%fname),H5F_ACC_RDONLY_F, file_id, hdferr, access_prp = plist_file)
      CALL h5gopen_f(file_id, trim(hslab%gname), group_id, hdferr)
      CALL h5screate_simple_f(hslab%rank, hslab%dimsm, memspace, hdferr)

      CALL h5dopen_f(group_id, trim(hslab%dname(1)), dset_id, hdferr)

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

    end subroutine ReadHDF5_1D_1v

    ! Write 1D array data
    subroutine WriteHDF5_1D_1v(hslab, buffer)
      type(tp_rdwt_hdf5), intent(in) :: hslab
      real(8), intent(in) :: buffer(hslab%dimsf(1))
      integer(HID_T) :: file_id, group_id, memspace
      integer(HID_T) :: dspace_id, dset_id
      integer(HID_T) :: plist_file, plist_cdrw

      if(.not.(IsHDF5Initialized.and.hslab%IsHSInitialized)) then
        print *, 'HDF5 or type hslab NOT intialized !!! STOP'
        stop
      endif

      ! Setup file access property list with parallel I/O access.
      CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_file, hdferr)
      CALL h5pset_fapl_mpio_f(plist_file, hslab%comm, hslab%info, hdferr)

      ! Create property list for collective dataset read/write
      CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_cdrw, hdferr)
      CALL h5pset_dxpl_mpio_f(plist_cdrw, H5FD_MPIO_COLLECTIVE_F, hdferr)

      CALL h5fopen_f(trim(hslab%fname),H5F_ACC_RDWR_F, file_id, hdferr, access_prp = plist_file)

      CALL h5screate_simple_f(hslab%rank, hslab%dimsm, memspace, hdferr)
      CALL h5screate_simple_f(hslab%rank, hslab%dimsf, dspace_id, hdferr)

      CALL h5dcreate_f(file_id, trim(hslab%dname(1)), &
                        H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)
      CALL h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, hslab%offset, &
                                 hslab%count, hdferr, hslab%stride, hslab%block)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, buffer,  &
                      hslab%dimsf, hdferr, file_space_id=dspace_id, mem_space_id=memspace, xfer_prp = plist_cdrw)
      CALL h5dclose_f(dset_id, hdferr)

      CALL h5sclose_f(dspace_id, hdferr)
      CALL h5sclose_f(memspace, hdferr)
      CALL h5fclose_f(file_id, hdferr)
      CALL h5pclose_f(plist_file, hdferr)
      CALL h5pclose_f(plist_cdrw, hdferr)

    end subroutine WriteHDF5_1D_1v


      subroutine WriteAttribute(hslab)
        implicit none

        type(tp_rdwt_hdf5), intent(in) :: hslab
        integer(HID_T) :: file_id, group_id, memspace
        integer(HID_T) :: plist_file, plist_cdrw
        integer(HID_T) :: attr_id, attr_space, atype_id


        CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_file, hdferr)
        CALL h5pset_fapl_mpio_f(plist_file, hslab%comm, hslab%info, hdferr)

        ! Create property list for collective dataset read/write
        CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_cdrw, hdferr)
        CALL h5pset_dxpl_mpio_f(plist_cdrw, H5FD_MPIO_COLLECTIVE_F, hdferr)

        CALL h5fopen_f(trim(hslab%fname), H5F_ACC_RDWR_F, file_id, hdferr, access_prp = plist_file)
        CALL h5gopen_f(file_id, trim(hslab%gname), group_id, hdferr)

        ! attribute 1
        CALL h5screate_simple_f(1, hslab%attr_dims_no, attr_space, hdferr)
        CALL h5tcopy_f(H5T_NATIVE_INTEGER, atype_id, hdferr)
        CALL h5acreate_f(group_id,hslab%attr_name(1), atype_id,attr_space,attr_id,hdferr)
        CALL h5awrite_f(attr_id,H5T_NATIVE_INTEGER,hslab%attr_time,hslab%attr_dims_no,hdferr)
        CALL h5sclose_f(attr_space, hdferr)
        CALL h5aclose_f(attr_id, hdferr)

        ! attribute 2
        CALL h5screate_simple_f(1, hslab%attr_dims_info, attr_space, hdferr)
        CALL h5tcopy_f(H5T_NATIVE_INTEGER, atype_id, hdferr)
        CALL h5acreate_f(group_id,hslab%attr_name(2), atype_id,attr_space,attr_id,hdferr)
        CALL h5awrite_f(attr_id,H5T_NATIVE_INTEGER,hslab%info_1st_array,hslab%attr_dims_info,hdferr)
        CALL h5sclose_f(attr_space, hdferr)
        CALL h5aclose_f(attr_id, hdferr)

        ! attribute 3
        CALL h5screate_simple_f(1, hslab%attr_dims_info, attr_space, hdferr)
        CALL h5tcopy_f(H5T_NATIVE_INTEGER, atype_id, hdferr)
        CALL h5acreate_f(group_id, hslab%attr_name(3), atype_id,attr_space,attr_id,hdferr)
        CALL h5awrite_f(attr_id,H5T_NATIVE_INTEGER, hslab%info_2nd_array, hslab%attr_dims_info,hdferr)
        CALL h5sclose_f(attr_space, hdferr)
        CALL h5aclose_f(attr_id, hdferr)

        ! attribute 4
        CALL h5screate_simple_f(1, hslab%attr_dims_loc, attr_space, hdferr)
        CALL h5tcopy_f(H5T_NATIVE_INTEGER, atype_id, hdferr)
        CALL h5acreate_f(group_id,hslab%attr_name(4), atype_id,attr_space,attr_id,hdferr)
        CALL h5awrite_f(attr_id,H5T_NATIVE_INTEGER, hslab%info_loc,hslab%attr_dims_loc,hdferr)
        CALL h5sclose_f(attr_space, hdferr)
        CALL h5aclose_f(attr_id, hdferr)

        CALL h5gclose_f(group_id, hdferr)
        CALL h5fclose_f(file_id, hdferr)
        CALL h5pclose_f(plist_file, hdferr)
        CALL h5pclose_f(plist_cdrw, hdferr)

      end subroutine WriteAttribute


    ! Read 4D data. serial read data
    !
    !
    subroutine Serial_ReadHDF5_4D(hslab, buffer)
      type(tp_rdwt_hdf5), intent(in) :: hslab
      real(8), intent(out) :: buffer(hslab%dimsm(1),hslab%dimsm(2),hslab%dimsm(3),hslab%dimsm(4),hslab%dnum)
      integer(HID_T) :: file_id, group_id, memspace
      integer(HID_T) :: dspace_id, dset_id
      integer :: n

      if(.not.(IsHDF5Initialized.and.hslab%IsHSInitialized)) then
          print *, 'HDF5 or type hslab NOT intialized !!! STOP'
          stop
      endif

!   if(myid.eq.0) then
!     print *, 'size(buffer,1) = ', size(buffer,1)
!     print *, 'size(buffer,2) = ', size(buffer,2)
!     print *, 'size(buffer,3) = ', size(buffer,3)
!     print *, 'size(buffer,4) = ', size(buffer,4)
!     print *, 'hslab%offset = ', hslab%offset
!     print *, 'dimsm = ', hslab%dimsm
!     print *, 'fname = ', trim(hslab%fname)
!     print *, 'gname = ', trim(hslab%gname)
!     print *, 'dnum = ', hslab%dnum
!   endif

      CALL h5fopen_f(trim(hslab%fname), H5F_ACC_RDONLY_F, file_id, hdferr)
      CALL h5gopen_f(file_id, trim(hslab%gname), group_id, hdferr)
      CALL h5screate_simple_f(hslab%rank, hslab%dimsm, memspace, hdferr)
      do n=1, hslab%dnum
        CALL h5dopen_f(group_id, trim(hslab%dname(n)), dset_id, hdferr)
        CALL h5dget_space_f(dset_id, dspace_id, hdferr)
        CALL h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, hslab%offset, &
                                   hslab%count, hdferr, hslab%stride, hslab%block)
        CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, buffer(:,:,:,:,n), &
                      hslab%dimsm, hdferr, file_space_id=dspace_id, mem_space_id=memspace)
        CALL h5sclose_f(dspace_id, hdferr)
        CALL h5dclose_f(dset_id, hdferr)
      enddo
      CALL h5sclose_f(memspace, hdferr)
      CALL h5gclose_f(group_id, hdferr)
      CALL h5fclose_f(file_id, hdferr)

    end subroutine Serial_ReadHDF5_4D



      subroutine ReadHDF5_scalar(hslab, scalar)
        type(tp_rdwt_hdf5), intent(in) :: hslab
        real(8), intent(out) :: scalar
        integer(HID_T) :: file_id, group_id
        integer(HID_T) :: dspace_id, dset_id
        integer(HID_T) :: plist_file, plist_cdrw

        ! Setup file access property list with parallel I/O access.
        CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_file, hdferr)
        CALL h5pset_fapl_mpio_f(plist_file, hslab%comm, hslab%info, hdferr)

        CALL h5fopen_f(trim(hslab%fname), H5F_ACC_RDONLY_F, file_id, hdferr, access_prp = plist_file)
        CALL h5gopen_f(file_id, trim(hslab%gname), group_id, hdferr)

        CALL h5dopen_f(group_id, trim(hslab%sname), dset_id, hdferr)
        CALL h5dget_space_f(dset_id, dspace_id, hdferr)
        !   CALL h5dread_f(dset_id, hslab%sd_type, scalar, hslab%dimsm,hdferr)
        CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, scalar, hslab%dimsm,hdferr)

        CALL h5dclose_f(dset_id, hdferr)
        CALL h5sclose_f(dspace_id, hdferr)
        CALL h5gclose_f(group_id, hdferr)
        CALL h5fclose_f(file_id, hdferr)
        CALL h5pclose_f(plist_file, hdferr)

      end subroutine ReadHDF5_scalar

      subroutine ReadHDF5_3D(hslab, buffer)
        type(tp_rdwt_hdf5), intent(in) :: hslab
        real(8), intent(out) :: buffer(hslab%dimsm(1),hslab%dimsm(2),hslab%dimsm(3),hslab%dnum)
        integer :: n
        integer(HID_T) :: file_id, group_id, memspace
        integer(HID_T) :: dspace_id, dset_id
        integer(HID_T) :: plist_file, plist_cdrw

        if(.not.(IsHDF5Initialized.and.hslab%IsHSInitialized)) then
          print *, 'HDF5 or type hslab NOT intialized !!! STOP'
          stop
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
        integer(HID_T) :: file_id, group_id, memspace
        integer(HID_T) :: dspace_id, dset_id
        integer(HID_T) :: plist_file, plist_cdrw

        if(.not.(IsHDF5Initialized.and.hslab%IsHSInitialized)) then
          print *, 'HDF5 or type hslab NOT intialized !!! STOP'
          stop
        endif

        ! Setup file access property list with parallel I/O access.
        CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_file, hdferr)
        CALL h5pset_fapl_mpio_f(plist_file, hslab%comm, hslab%info, hdferr)

        ! Create property list for collective dataset read/write
        CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_cdrw, hdferr)
        CALL h5pset_dxpl_mpio_f(plist_cdrw, H5FD_MPIO_COLLECTIVE_F, hdferr)

        CALL h5fcreate_f(trim(hslab%fname), H5F_ACC_TRUNC_F, file_id, hdferr,access_prp = plist_file)
        ! CALL h5gcreate_f(file_id,trim(hslab%gname), group_id, hdferr)
        CALL h5screate_simple_f(hslab%rank, hslab%dimsm, memspace, hdferr)
        CALL h5screate_simple_f(hslab%rank, hslab%dimsf, dspace_id, hdferr)
        do n = 1, hslab%dnum
         ! CALL h5dcreate_f(group_id, trim(hslab%dname(n)), &
         !                  H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)
         CALL h5dcreate_f(file_id, trim(hslab%dname(n)), &
                          H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)

         CALL h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, hslab%offset, &
                                    hslab%count, hdferr, hslab%stride, hslab%block)
         CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, buffer(:,:,:,n),  &
                         hslab%dimsf, hdferr, file_space_id=dspace_id, mem_space_id=memspace, xfer_prp = plist_cdrw)
         CALL h5dclose_f(dset_id, hdferr)
        enddo
        CALL h5sclose_f(dspace_id, hdferr)
        CALL h5sclose_f(memspace, hdferr)
        ! CALL h5gclose_f(group_id, hdferr)
        CALL h5fclose_f(file_id, hdferr)
        CALL h5pclose_f(plist_file, hdferr)
        CALL h5pclose_f(plist_cdrw, hdferr)

      end subroutine WriteHDF5_3D


      subroutine WriteHDF5_4D(hslab, buffer)
        type(tp_rdwt_hdf5), intent(in) :: hslab
        real(8), intent(in) :: buffer(hslab%dimsm(1),hslab%dimsm(2),hslab%dimsm(3),hslab%dimsm(4),hslab%dnum)
        integer :: n
        integer(HID_T) :: file_id, group_id, memspace
        integer(HID_T) :: dspace_id, dset_id
        integer(HID_T) :: plist_file, plist_cdrw

        if(.not.(IsHDF5Initialized.and.hslab%IsHSInitialized)) then
          print *, 'HDF5 or type hslab NOT intialized !!! STOP'
          stop
        endif

        ! Setup file access property list with parallel I/O access.
        CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_file, hdferr)
        CALL h5pset_fapl_mpio_f(plist_file, hslab%comm, hslab%info, hdferr)

        ! Create property list for collective dataset read/write
        CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_cdrw, hdferr)
        CALL h5pset_dxpl_mpio_f(plist_cdrw, H5FD_MPIO_COLLECTIVE_F, hdferr)

        if((hslab%IsMultiGroup).or.(hslab%IsSameGroup)) then
          CALL h5fopen_f(trim(hslab%fname),H5F_ACC_RDWR_F, file_id, hdferr, access_prp = plist_file)
        else
          CALL h5fcreate_f(trim(hslab%fname), H5F_ACC_TRUNC_F, file_id, hdferr, access_prp = plist_file)
        endif

        if(trim(hslab%gname).ne.'/') then
          if(hslab%IsSameGroup) then
            CALL h5gopen_f(file_id,trim(hslab%gname), group_id, hdferr)
          else
            CALL h5gcreate_f(file_id,trim(hslab%gname), group_id, hdferr)
          endif
        endif
        CALL h5screate_simple_f(hslab%rank, hslab%dimsm, memspace, hdferr)

        do n = 1, hslab%dnum
          if(trim(hslab%gname).ne.'/') then
            if(.not.hslab%IsSameData) then
              CALL h5screate_simple_f(hslab%rank, hslab%dimsf, dspace_id, hdferr)
              CALL h5dcreate_f(group_id, trim(hslab%dname(n)), &
                               H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)
            else
              CALL h5dopen_f(group_id, trim(hslab%dname(n)), dset_id, hdferr)
              CALL h5dget_space_f(dset_id, dspace_id, hdferr)
            endif
          else
            CALL h5screate_simple_f(hslab%rank, hslab%dimsf, dspace_id, hdferr)
            CALL h5dcreate_f(file_id, trim(hslab%dname(n)), &
                             H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)
          endif

          CALL h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, hslab%offset, &
                                    hslab%count, hdferr, hslab%stride, hslab%block)
          CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, buffer(:,:,:,:,n),  &
                         hslab%dimsf, hdferr, file_space_id=dspace_id, mem_space_id=memspace, xfer_prp = plist_cdrw)
          CALL h5dclose_f(dset_id, hdferr)
          CALL h5sclose_f(dspace_id, hdferr)
        enddo ! end n loop
        CALL h5sclose_f(memspace, hdferr)
        if(trim(hslab%gname).ne.'/') then
          CALL h5gclose_f(group_id, hdferr)
        endif
        CALL h5fclose_f(file_id, hdferr)
        CALL h5pclose_f(plist_file, hdferr)
        CALL h5pclose_f(plist_cdrw, hdferr)

      end subroutine WriteHDF5_4D
















































    end module modPHDF5
