!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! InitHDF5()             : open  HDF5
! FinalizeHDF5()         : close HDF5
! InitGridHDF5()         :
! InitFlowHDF5()         :
! InitFlow_PlaneInlet()  :
! InitFlow_tsVol(hslab)  : timeseries volume
! InitFlow_1D()          :
! DetectHDF5()           :
! HDF5_p3d_3D_sol()      : convert HDF5 sol  file to plot3d
! HDF5_p3d_3d_grd()      : convert HDF5 grid file to plot3d
! ReadHDF5_scalar()      :
! WriteHDF5_scalar()     :
! ReadHDF5_1D()          :
! ReadHDF5_2D()          :
! ReadHDF5_3D()          :
! ReadHDF5_4D()          :
! WriteHDF5_1D()         :
! WriteHDF5_2D()         :
! WriteHDF5_3D()         :
! ReadHDF5_3D_1V()       : read 3D one variable data
! WriteHDF5_3D_1V()      : write 3D one variable data
! WriteHDF5_PlaneInlet() :
! WriteHDF5_1D_1v()      : write 1D array data
! ReadtsGridMetrics()    : read timeseries GridMetrics file
! Readtsflow()           : read timeseries flow data
! Write1dHDF5()          : write 1d variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module modRWHDF5
    use HDF5
    implicit none

    integer(HID_T), private :: file_id, group_id, dset_id, dspace_id, memspace_id
    integer, private :: hdferr
    logical :: IsHDF5Initialized = .false.

    type tp_rdwt_hdf5
      character(400) :: fname
      character(100) :: gname
      character(50)  :: sname  ! scalar name
      character(50), dimension(:), allocatable :: dname
      character(50), dimension(:), allocatable :: attr_name
      integer(HSIZE_T), dimension(:), allocatable :: dimsf
      integer :: rank, dnum, anum
      integer(HSIZE_T), dimension(:), allocatable :: dimsm
      integer(HSIZE_T), dimension(:), allocatable :: count, offset, block, stride
      logical :: IsHSInitialized = .false.
      logical :: IsMultiGroup = .false.
!      logical :: IsNewFile = .false.

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

    type tp_plane
      type(tp_rdwt_hdf5), dimension(:,:), allocatable :: TSplane
      type(tp_rdwt_hdf5), dimension(:), allocatable   :: TSgrid
      integer :: nfile, ntpoint_total, npath, nvar_output
      character(10), dimension(:), allocatable :: fname_range
      integer :: iminloc, imaxloc, iavelen
      character(400) :: filepath
      type(tp_DNSIndex) :: DNSIndex
      type(fprop), dimension(:), allocatable :: fileprop

      character(10), dimension(:), allocatable :: varout
      integer :: nxpoint, nypoint, nzpoint, imints, imaxts
      integer :: jmints, jmaxts, jminloc, jmaxloc, javelen
      integer :: nt_buffer
    end type tp_plane

!    type tp_rdwr_planeinlet
!        character(300) :: fname
!        character(20), dimension(:), allocatable :: gname
!        character(20)  :: sname  ! scalar name
!        character(20), dimension(:), allocatable :: dname
!        integer(HSIZE_T), dimension(:), allocatable :: dimsf
!        integer :: rank, gnum, dnum
!        logical :: IsHSInitialized = .false.
!    end type tp_rdwr_planeinlet


contains

    subroutine InitHDF5()
        CALL h5open_f(hdferr)
        IsHDF5Initialized = .True.
    end subroutine InitHDF5

    subroutine FinalizeHDF5()
        CALL h5close_f(hdferr)
    end subroutine FinalizeHDF5


!******************************************************************************
! Initialization: InitGridHDF5(), InitFlowHDF5(), InitFlow_PlaneInlet(), InitFlow_1D()
!******************************************************************************
    subroutine InitGridHDF5(hslab)
      type(tp_rdwt_hdf5), intent(out) :: hslab

      hslab%gname = '/'
      hslab%dnum = 3
      hslab%rank = 3
      if(allocated(hslab%dname)) deallocate(hslab%dname)
      if(allocated(hslab%dimsf)) deallocate(hslab%dimsf)
      if(allocated(hslab%dimsm)) deallocate(hslab%dimsm)
      if(allocated(hslab%count)) deallocate(hslab%count)
      if(allocated(hslab%offset)) deallocate(hslab%offset)
      if(allocated(hslab%block)) deallocate(hslab%block)
      if(allocated(hslab%stride)) deallocate(hslab%stride)
      allocate(hslab%dname(hslab%dnum), hslab%dimsf(hslab%rank))
      hslab%dname(1) = 'x'
      hslab%dname(2) = 'y'
      hslab%dname(3) = 'z'
      hslab%IsHSInitialized = .true.
      allocate(hslab%dimsm(hslab%rank), hslab%count(hslab%rank))
      allocate(hslab%offset(hslab%rank), hslab%block(hslab%rank),hslab%stride(hslab%rank))
      ! default value
      hslab%count = 1
      hslab%stride = 1
      hslab%offset = 0

    end subroutine InitGridHDF5

    subroutine InitHDF5_1V(hslab)
      type(tp_rdwt_hdf5), intent(out) :: hslab

      hslab%gname = '/'
      hslab%dnum = 1
      hslab%rank = 3
      if(allocated(hslab%dname)) deallocate(hslab%dname)
      if(allocated(hslab%dimsf)) deallocate(hslab%dimsf)
      if(allocated(hslab%dimsm)) deallocate(hslab%dimsm)
      if(allocated(hslab%count)) deallocate(hslab%count)
      if(allocated(hslab%offset)) deallocate(hslab%offset)
      if(allocated(hslab%block)) deallocate(hslab%block)
      if(allocated(hslab%stride)) deallocate(hslab%stride)
      allocate(hslab%dname(hslab%dnum), hslab%dimsf(hslab%rank))
      hslab%sname = 'time'
      hslab%IsHSInitialized = .true.
      allocate(hslab%dimsm(hslab%rank), hslab%count(hslab%rank))
      allocate(hslab%offset(hslab%rank), hslab%block(hslab%rank),hslab%stride(hslab%rank))
      ! default value
      hslab%count = 1
      hslab%stride = 1
      hslab%offset = 0

    end subroutine InitHDF5_1V

    subroutine InitFlowHDF5(hslab)
      type(tp_rdwt_hdf5), intent(out) :: hslab

      hslab%gname = '/'
      hslab%dnum = 5
      hslab%rank = 3
      allocate(hslab%dname(hslab%dnum), hslab%dimsf(hslab%rank))
      hslab%dname(1) = 'u'
      hslab%dname(2) = 'v'
      hslab%dname(3) = 'w'
      hslab%dname(4) = 'p'
      hslab%dname(5) = 'T'
      hslab%sname = 'time'
      hslab%IsHSInitialized = .true.
      allocate(hslab%dimsm(hslab%rank), hslab%count(hslab%rank))
      allocate(hslab%offset(hslab%rank), hslab%block(hslab%rank),hslab%stride(hslab%rank))
      ! default value
      hslab%count = 1
      hslab%stride = 1
      hslab%offset = 0

    end subroutine InitFlowHDF5

    subroutine InitFlowHDF5_4D(hslab)
      type(tp_rdwt_hdf5), intent(out) :: hslab

      hslab%gname = '/'
      hslab%dnum = 5
      hslab%rank = 4
      allocate(hslab%dname(hslab%dnum), hslab%dimsf(hslab%rank))
      hslab%dname(1) = 'u'
      hslab%dname(2) = 'v'
      hslab%dname(3) = 'w'
      hslab%dname(4) = 'p'
      hslab%dname(5) = 'T'
      hslab%sname = 'time'
      hslab%IsHSInitialized = .true.
      allocate(hslab%dimsm(hslab%rank), hslab%count(hslab%rank))
      allocate(hslab%offset(hslab%rank), hslab%block(hslab%rank),hslab%stride(hslab%rank))
      ! default value
      hslab%count = 1
      hslab%stride = 1
      hslab%offset = 0

    end subroutine InitFlowHDF5_4D

    subroutine InitFlowHDF5_tsprobe(hslab,indx_var,nvar)
        type(tp_rdwt_hdf5), intent(out) :: hslab
        integer,intent(in) :: nvar
        integer,dimension(nvar),intent(in) :: indx_var
        integer :: n

        hslab%dnum = nvar
        hslab%rank = 4
        allocate(hslab%dname(hslab%dnum), hslab%dimsf(hslab%rank))
        do n = 1,nvar
            if(indx_var(n).eq.1) then
                hslab%dname(n) = 'u'
            endif
            if(indx_var(n).eq.2) then
                hslab%dname(n) = 'v'
            endif
            if(indx_var(n).eq.3) then
                hslab%dname(n) = 'w'
            endif
            if(indx_var(n).eq.4) then
                hslab%dname(n) = 'p'
            endif
            if(indx_var(n).eq.5) then
                hslab%dname(n) = 'T'
            endif
        enddo
        allocate(hslab%dimsm(hslab%rank), hslab%count(hslab%rank))
        allocate(hslab%offset(hslab%rank), hslab%block(hslab%rank),hslab%stride(hslab%rank))
        ! default value
        hslab%count = 1
        hslab%stride = 1
        hslab%block = 1
        hslab%offset = 0
        hslab%IsHSInitialized = .true.
    end subroutine InitFlowHDF5_tsprobe

    subroutine InitFlow_PlaneInlet(hslab)
      type(tp_rdwt_hdf5), intent(out) :: hslab

      hslab%dnum = 5
      hslab%rank = 3
      allocate(hslab%dname(hslab%dnum), hslab%dimsf(hslab%rank))
      hslab%dname(1) = 'u'
      hslab%dname(2) = 'v'
      hslab%dname(3) = 'w'
      hslab%dname(4) = 'rho'
      hslab%dname(5) = 'T'

      hslab%IsHSInitialized = .true.
      allocate(hslab%dimsm(hslab%rank), hslab%count(hslab%rank))
      allocate(hslab%offset(hslab%rank), hslab%block(hslab%rank),hslab%stride(hslab%rank))
      ! default value
      hslab%count = 1
      hslab%stride = 1
      hslab%offset = 0

    end subroutine InitFlow_PlaneInlet

    subroutine InitFlow_tsVol(hslab)
      type(tp_rdwt_hdf5), intent(out) :: hslab

      hslab%gname = 'vol'
      hslab%sname = 'time'
      hslab%dnum = 5
      hslab%rank = 3
      allocate(hslab%dname(hslab%dnum), hslab%dimsf(hslab%rank))
      hslab%dname(1) = 'u'
      hslab%dname(2) = 'v'
      hslab%dname(3) = 'w'
      hslab%dname(4) = 'p'
      hslab%dname(5) = 'T'

      hslab%anum = 3
      allocate(hslab%attr_name(hslab%anum))
      hslab%attr_name(1) = "istart, iend, iskip"
      hslab%attr_name(2) = "jstart, jend, jskip"
      hslab%attr_name(3) = "kstart, kend, kskip"

      hslab%IsHSInitialized = .true.

    end subroutine InitFlow_tsVol


    subroutine InitFlow_1D(hslab)
      type(tp_rdwt_hdf5), intent(out) :: hslab

      hslab%gname = '/'
      hslab%dnum = 1
      hslab%rank = 1
      allocate(hslab%dname(hslab%dnum), hslab%dimsf(hslab%rank))
      hslab%dname(1) = 'time'

      hslab%IsHSInitialized = .true.
      allocate(hslab%dimsm(hslab%rank), hslab%count(hslab%rank))
      allocate(hslab%offset(hslab%rank), hslab%block(hslab%rank),hslab%stride(hslab%rank))
      ! default value
      hslab%count = 1
      hslab%stride = 1
      hslab%offset = 0

    end subroutine InitFlow_1D


  subroutine InitRescalemean(hslab)
    type(tp_rdwt_hdf5), intent(out) :: hslab

    hslab%gname = '/'
    hslab%dnum = 4
    hslab%rank = 3
    allocate(hslab%dname(hslab%dnum),hslab%dimsf(hslab%rank))
    hslab%dname(1) = 'uave_recycle'
    hslab%dname(2) = 'wave_recycle'
    hslab%dname(3) = 'rhoave_recycle'
    hslab%dname(4) = 'tave_recycle'
    allocate(hslab%dimsm(hslab%rank), hslab%count(hslab%rank))
    allocate(hslab%offset(hslab%rank), hslab%block(hslab%rank),hslab%stride(hslab%rank))
    ! default value
    hslab%count = 1
    hslab%stride = 1
    hslab%offset = 0

    hslab%IsHSInitialized = .true.

  end subroutine InitRescalemean


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
! -----------------------------------------------------------------------

! initialize HDF5 for flow 2d
subroutine InitFlow2dHDF5(hslab)
    type(tp_rdwt_hdf5), intent(out) :: hslab

    hslab%gname = '/'
    hslab%dnum = 12
    hslab%rank = 2
    allocate(hslab%dname(hslab%dnum), hslab%dimsf(hslab%rank))
    hslab%dname(1) = 'u'
    hslab%dname(2) = 'v'
    hslab%dname(3) = 'w'
    hslab%dname(4) = 'p'
    hslab%dname(5) = 'T'
    hslab%dname(6) = 'rho'
    hslab%dname(7) = 'uu'
    hslab%dname(8) = 'vv'
    hslab%dname(9) = 'ww'
    hslab%dname(10) = 'uv'
    hslab%dname(11) = 'uw'
    hslab%dname(12) = 'vw'
    hslab%IsHSInitialized = .true.
    allocate(hslab%dimsm(hslab%rank), hslab%count(hslab%rank))
    allocate(hslab%offset(hslab%rank), hslab%block(hslab%rank),hslab%stride(hslab%rank))
    ! default value
    hslab%count = 1
    hslab%stride = 1
    hslab%offset = 0

end subroutine InitFlow2dHDF5


    ! Detect HDF5 file
    subroutine DetectHDF5(hslab)
      implicit none
      type(tp_rdwt_hdf5), intent(inout) :: hslab
      integer(HSIZE_T), dimension(:), allocatable :: dimsf_tmp
      integer :: i, idx

      if(.not.IsHDF5Initialized) then
        print *, 'HDF5 should be Initialized '
        stop
      endif

      if(allocated(hslab%dname)) deallocate(hslab%dname)
      if(allocated(hslab%dimsf)) deallocate(hslab%dimsf)

      CALL h5fopen_f(trim(hslab%fname),H5F_ACC_RDONLY_F, file_id, hdferr)
      CALL h5gopen_f(file_id, trim(hslab%gname), group_id, hdferr)
      CALL h5gn_members_f(group_id, trim(hslab%gname), hslab%dnum, hdferr)

      allocate(hslab%dname(hslab%dnum))
      do idx=0, hslab%dnum-1
        CALL h5gget_obj_info_idx_f(group_id, trim(hslab%gname), idx, hslab%dname(idx+1), H5G_DATASET_F, hdferr)
      enddo

      CALL h5dopen_f(group_id, trim(hslab%dname(1)), dset_id, hdferr)
      CALL h5dget_space_f(dset_id, dspace_id, hdferr)
      CALL h5sget_simple_extent_ndims_f(dspace_id, hslab%rank, hdferr)
      allocate(hslab%dimsf(hslab%rank),dimsf_tmp(hslab%rank))
      CALL h5sget_simple_extent_dims_f(dspace_id, hslab%dimsf, dimsf_tmp, hdferr)
      CALL h5dclose_f(dset_id, hdferr)
      CALL h5sclose_f(dspace_id, hdferr)
      CALL h5gclose_f(group_id, hdferr)
      CALL h5fclose_f(file_id, hdferr)

      if(allocated(hslab%count)) then
        deallocate(hslab%dimsm,hslab%count,hslab%offset,hslab%block,hslab%stride)
      endif
      allocate(hslab%dimsm(hslab%rank), hslab%count(hslab%rank))
      allocate(hslab%offset(hslab%rank), hslab%block(hslab%rank),hslab%stride(hslab%rank))
      ! default value
      hslab%dimsm = hslab%dimsf
      hslab%count = 1
      hslab%stride = 1
      hslab%offset = 0
      hslab%block = hslab%dimsm
      hslab%IsHSInitialized = .true.

    end subroutine DetectHDF5


!     ! Detect HDF5 file. return # of datsset, dataset name, rank, dimsf
!     subroutine DetectHDF5(fname,gname,num_dset,dname,rank,dims)
!       character(*), intent(in) :: fname
!       character(*), intent(in) :: gname
!       integer, intent(out) :: num_dset
!       character(*), dimension(:), allocatable, intent(out) :: dname
!       integer, intent(out) :: rank
!       integer, dimension(:), allocatable, intent(out) :: dims
!       integer(HSIZE_T), dimension(:), allocatable :: dimsf, dimsf_total
!       integer :: i, idx
!
!       if(.not.IsHDF5Initialized) then
!         print *, 'HDF5 should be Initialized'
!         stop
!       endif
!
!       if(allocated(dname)) deallocate(dname)
!       if(allocated(dimsf)) deallocate(dimsf)
!       if(allocated(dimsf_total)) deallocate(dimsf_total)
!       if(allocated(dims)) deallocate(dims)
!
!       CALL h5fopen_f(trim(fname), H5F_ACC_RDONLY_F, file_id, hdferr)
!
!       CALL h5gopen_f(file_id, trim(adjustl(gname)), group_id, hdferr)
!       CALL h5gn_members_f(group_id, trim(adjustl(gname)), num_dset, hdferr)
!
!       allocate(dname(num_dset))
!       do idx=0, num_dset-1
!         CALL h5gget_obj_info_idx_f(group_id, trim(adjustl(gname)), idx, dname(idx+1), H5G_DATASET_F, hdferr)
!       enddo
!
!       CALL h5dopen_f(group_id, trim(dname(1)), dset_id, hdferr)
!       CALL h5dget_space_f(dset_id, dspace_id, hdferr)
!       CALL h5sget_simple_extent_ndims_f(dspace_id, rank, hdferr)
!       allocate(dimsf(rank), dimsf_total(rank))
!       CALL h5sget_simple_extent_dims_f(dspace_id, dimsf, dimsf_total, hdferr)
!       CALL h5dclose_f(dset_id, hdferr)
!       CALL h5gclose_f(group_id, hdferr)
!
!       CALL h5fclose_f(file_id, hdferr)
!
!       allocate(dims(rank))
!       dims = dimsf
!
!     end subroutine DetectHDF5


!******************************************************************************
!     convert HDF5 file to plot3d
!     in: hslab, fname, kioplane
!******************************************************************************
    subroutine HDF5_p3d_3D_sol(hslab,fn,kioplane,time)
      implicit none
      type(tp_rdwt_hdf5), intent(in) :: hslab
      character(*), intent(in) :: fn
      integer, intent(in) :: kioplane
      real(8), intent(in) :: time
      integer :: i, j, k, n, k2, kp, kk
      integer :: idim, jdim, kdim
      integer(HSIZE_T), dimension(:), allocatable :: dimsm
      integer(HSIZE_T), dimension(:), allocatable :: count, offset, stride, block
      real(8), dimension(:,:,:,:), allocatable :: buffer
      integer :: nvar

      if(.not.(IsHDF5Initialized.and.hslab%IsHSInitialized)) then
          print *, 'HDF5 or type hslab NOT intialized !!! STOP'
          stop
      endif

      idim = hslab%dimsf(2)
      jdim = hslab%dimsf(3)
      kdim = hslab%dimsf(1)
      allocate(dimsm(hslab%rank), count(hslab%rank),offset(hslab%rank))
      allocate(stride(hslab%rank),block(hslab%rank))
      allocate(buffer(kioplane,hslab%dimsf(2),hslab%dimsf(3),hslab%dnum+1))
      dimsm(2) = hslab%dimsf(2)
      dimsm(3) = hslab%dimsf(3)
      block(2) = dimsm(2)
      block(3) = dimsm(3)
      offset = 0
      count  = 1
      stride = 1

!      print *, 'Reading file: ', trim(hslab%fname)
!      print *, 'dnum = ', hslab%dnum
!      print *, 'dname = ', hslab%dname
!      print *, 'dimsf = ', hslab%dimsf

      CALL h5fopen_f(trim(hslab%fname), H5F_ACC_RDONLY_F, file_id, hdferr)
      CALL h5gopen_f(file_id, trim(hslab%gname), group_id, hdferr)
      CALL h5screate_simple_f(hslab%rank, dimsm, memspace_id, hdferr)

      nvar = 6
      open(unit=11,file=trim(fn),form='unformatted',status='unknown')
        write(11) 1
        write(11) idim, jdim, kdim, nvar

      nvar = 5
      open(unit=12,file=trim('inlet.sol'),form='unformatted',status='unknown')
        write(12) 1
        write(12) 1, jdim, kdim, nvar

        do k=1, hslab%dimsf(1), kioplane
          k2 = min(k+kioplane-1,hslab%dimsf(1))
          kp = k2-k+1
          dimsm(1) = kp
          block(1) = dimsm(1)
          offset(1) = k-1
          if(k.eq.1.or.kp.ne.kioplane) then
            CALL h5screate_simple_f(hslab%rank, dimsm, memspace_id, hdferr)
          endif

          do n=1, hslab%dnum
            CALL h5dopen_f(group_id, trim(hslab%dname(n)), dset_id, hdferr)
            CALL h5dget_space_f(dset_id, dspace_id, hdferr)
            CALL h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, &
                                       count, hdferr, stride, block)
            CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, buffer(1:kp,:,:,n), &
                           dimsm, hdferr, file_space_id=dspace_id, mem_space_id=memspace_id)
            CALL h5sclose_f(dspace_id, hdferr)
            CALL h5dclose_f(dset_id, hdferr)
          enddo ! end n loop
          if(kp.ne.kioplane) then
            CALL h5sclose_f(memspace_id, hdferr)
          endif
          buffer(1:kp,:,:,hslab%dnum+1) = buffer(1:kp,:,:,4)/(287.d0*buffer(1:kp,:,:,5))
          do kk=1, kp
            write(11) (((buffer(kk,i,j,n),i=1,hslab%dimsf(2)),j=1,hslab%dimsf(3)),n=1,hslab%dnum+1)
          enddo
          do kk=1, kp
            write(12) (((buffer(kk,i,j,n),i=1,1),j=1,hslab%dimsf(3)),n=1,hslab%dnum)
          enddo


        enddo ! end k loop
        write(11) time
      close(11)
      close(12)
      CALL h5gclose_f(group_id, hdferr)
      CALL h5fclose_f(file_id, hdferr)

    end subroutine HDF5_p3d_3D_sol


!******************************************************************************
!     convert HDF5 file to plot3d
!     in: hslab, fname, kioplane
!******************************************************************************
    subroutine HDF5_p3d_3D_grd(hslab,fn,kioplane)
      implicit none
      type(tp_rdwt_hdf5), intent(in) :: hslab
      character(*), intent(in) :: fn
      integer, intent(in) :: kioplane
      integer :: i, j, k, n, k2, kp, kk
      integer :: idim, jdim, kdim
      integer(HSIZE_T), dimension(:), allocatable :: dimsm
      integer(HSIZE_T), dimension(:), allocatable :: count, offset, stride, block
      real(8), dimension(:,:,:,:), allocatable :: buffer

      if(.not.(IsHDF5Initialized.and.hslab%IsHSInitialized)) then
          print *, 'HDF5 or type hslab NOT intialized !!! STOP'
          stop
      endif

      idim = hslab%dimsf(2)
      jdim = hslab%dimsf(3)
      kdim = hslab%dimsf(1)
      allocate(dimsm(hslab%rank), count(hslab%rank),offset(hslab%rank))
      allocate(stride(hslab%rank),block(hslab%rank))
      allocate(buffer(kioplane,hslab%dimsf(2),hslab%dimsf(3),hslab%dnum))
      dimsm(2) = hslab%dimsf(2)
      dimsm(3) = hslab%dimsf(3)
      block(2) = dimsm(2)
      block(3) = dimsm(3)
      offset = 0
      count  = 1
      stride = 1

      CALL h5fopen_f(trim(hslab%fname), H5F_ACC_RDONLY_F, file_id, hdferr)
      CALL h5gopen_f(file_id, trim(hslab%gname), group_id, hdferr)
      CALL h5screate_simple_f(hslab%rank, dimsm, memspace_id, hdferr)

      open(unit=11,file=trim(fn),form='unformatted',status='unknown')
        write(11) 1
        write(11) idim, jdim, kdim

      open(unit=12,file=trim(fn)//'_inlet',form='unformatted',status='unknown')
        write(12) 1
        write(12) 1, jdim, kdim

        do k=1, hslab%dimsf(1), kioplane
          k2 = min(k+kioplane-1,hslab%dimsf(1))
          kp = k2-k+1
          dimsm(1) = kp
          block(1) = dimsm(1)
          offset(1) = k-1
          if(k.eq.1.or.kp.ne.kioplane) then
            CALL h5screate_simple_f(hslab%rank, dimsm, memspace_id, hdferr)
          endif

          do n=1, hslab%dnum
            CALL h5dopen_f(group_id, trim(hslab%dname(n)), dset_id, hdferr)
            CALL h5dget_space_f(dset_id, dspace_id, hdferr)
            CALL h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, &
                                       count, hdferr, stride, block)
            CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, buffer(1:kp,:,:,n), &
                           dimsm, hdferr, file_space_id=dspace_id, mem_space_id=memspace_id)
            CALL h5sclose_f(dspace_id, hdferr)
            CALL h5dclose_f(dset_id, hdferr)
          enddo ! end n loop
          if(kp.ne.kioplane) then
            CALL h5sclose_f(memspace_id, hdferr)
          endif
          do kk=1, kp
            write(11) (((buffer(kk,i,j,n),i=1,hslab%dimsf(2)),j=1,hslab%dimsf(3)),n=1,hslab%dnum)
          enddo
          do kk=1, kp
            write(12) (((buffer(kk,i,j,n),i=1,1),j=1,hslab%dimsf(3)),n=1,hslab%dnum)
          enddo

        enddo ! end k loop

      close(11)
      close(12)
      CALL h5gclose_f(group_id, hdferr)
      CALL h5fclose_f(file_id, hdferr)

    end subroutine HDF5_p3d_3D_grd

    ! Read scalar in the file
    subroutine ReadHDF5_scalar(hslab, scalar)
      type(tp_rdwt_hdf5), intent(in) :: hslab
      real(8), intent(out) :: scalar

      if(.not.(IsHDF5Initialized.and.hslab%IsHSInitialized)) then
          print *, 'HDF5 or type hslab NOT intialized !!! STOP'
          stop
      endif

      CALL h5fopen_f(trim(hslab%fname), H5F_ACC_RDONLY_F, file_id, hdferr)
      CALL h5gopen_f(file_id, trim(hslab%gname), group_id, hdferr)
      CALL h5dopen_f(group_id, trim(hslab%sname), dset_id, hdferr)
      CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, scalar, hslab%dimsf, hdferr)
      CALL h5dclose_f(dset_id, hdferr)
      CALL h5gclose_f(group_id, hdferr)
      CALL h5fclose_f(file_id, hdferr)

    end subroutine ReadHDF5_scalar

    ! Write scalar to the file
    subroutine WriteHDF5_scalar(hslab, scalar)
      type(tp_rdwt_hdf5), intent(in) :: hslab
      real(8), intent(in) :: scalar

      if(.not.(IsHDF5Initialized.and.hslab%IsHSInitialized)) then
          print *, 'HDF5 or type hslab NOT intialized !!! STOP'
          stop
      endif

      CALL h5fopen_f(trim(hslab%fname),H5F_ACC_RDWR_F, file_id, hdferr)
      if(trim(hslab%gname).ne.'/') then
        CALL h5gopen_f(file_id, trim(hslab%gname), group_id, hdferr)
      endif
      CALL h5screate_f(H5S_SCALAR_F, dspace_id, hdferr)
      if(trim(hslab%gname).ne.'/') then
        CALL h5dcreate_f(group_id, trim(hslab%sname),H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)
      else
        CALL h5dcreate_f(file_id, trim(hslab%sname),H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)
      endif
      CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, scalar, hslab%dimsf, hdferr)
      CALL h5dclose_f(dset_id, hdferr)
      CALL h5sclose_f(dspace_id, hdferr)
      if(trim(hslab%gname).ne.'/') then
        CALL h5gclose_f(group_id,hdferr)
      endif
      CALL h5fclose_f(file_id, hdferr)

    end subroutine WriteHDF5_scalar

!    ! Write scalar to the file
!    subroutine WriteHDF5_scalar_group(hslab, scalar)
!      type(tp_rdwt_hdf5), intent(in) :: hslab
!      real(8), intent(in) :: scalar
!
!      if(.not.(IsHDF5Initialized.and.hslab%IsHSInitialized)) then
!          print *, 'HDF5 or type hslab NOT intialized !!! STOP'
!          stop
!      endif
!
!      CALL h5fopen_f(trim(hslab%fname),H5F_ACC_RDWR_F, file_id, hdferr)
!      CALL h5gopen_f(file_id, trim(hslab%gname), group_id, hdferr)
!      CALL h5screate_f(H5S_SCALAR_F, dspace_id, hdferr)
!      CALL h5dcreate_f(group_id, trim(hslab%sname),H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)
!      CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, scalar, hslab%dimsf, hdferr)
!      CALL h5dclose_f(dset_id, hdferr)
!      CALL h5sclose_f(dspace_id, hdferr)
!      CALL h5gclose_f(group_id, hdferr)
!      CALL h5fclose_f(file_id, hdferr)
!
!    end subroutine WriteHDF5_scalar_group

    subroutine ReadHDF5_Attribute(hslab, buffer)
      implicit none
      type(tp_rdwt_hdf5), intent(in) :: hslab
      integer, intent(out) :: buffer(3,hslab%anum)
      integer(HID_T) :: attr_id, attr_space, atype_id
      integer(HSIZE_T), dimension(1) :: attr_dims_info, attr_dims_tmp
      integer :: n

      CALL h5fopen_f(trim(hslab%fname), H5F_ACC_RDONLY_F, file_id, hdferr)
      CALL h5gopen_f(file_id, trim(hslab%gname), group_id, hdferr)
      do n=1, hslab%anum
        CALL h5aopen_f(group_id, hslab%attr_name(n), attr_id, hdferr)
        CALL H5aget_space_f(attr_id, attr_space, hdferr)
        CALL H5sget_simple_extent_dims_f(attr_space, attr_dims_info, attr_dims_tmp,hdferr)
        CALL H5aread_f(attr_id, H5T_NATIVE_INTEGER, buffer(:,n), attr_dims_info, hdferr)
        CALL H5aclose_f(attr_id, hdferr)
      enddo ! end n loop
      CALL H5gclose_f(group_id, hdferr)
      CALL H5Fclose_f(file_id, hdferr)

    end subroutine ReadHDF5_Attribute

    subroutine WriteHDF5_Attribute(hslab, buffer)
      implicit none
      type(tp_rdwt_hdf5), intent(in) :: hslab
      integer, intent(in) :: buffer(3,hslab%anum)
      integer(HID_T) :: attr_id, attr_space, atype_id
      integer(HSIZE_T), dimension(1) :: attr_dims_info
      integer :: n

      attr_dims_info = 3
      CALL h5fopen_f(trim(hslab%fname), H5F_ACC_RDWR_F, file_id, hdferr)
      CALL h5gopen_f(file_id, trim(hslab%gname), group_id, hdferr)
      do n=1, hslab%anum
        CALL h5screate_simple_f(1, attr_dims_info, attr_space, hdferr)
        CALL h5tcopy_f(H5T_NATIVE_INTEGER, atype_id, hdferr)
        CALL h5acreate_f(group_id,hslab%attr_name(n), atype_id,attr_space,attr_id,hdferr)
        CALL h5awrite_f(attr_id,H5T_NATIVE_INTEGER,buffer(:,n),attr_dims_info,hdferr)
        CALL h5sclose_f(attr_space, hdferr)
        CALL h5aclose_f(attr_id, hdferr)
      enddo ! end n loop
      CALL H5gclose_f(group_id, hdferr)
      CALL H5Fclose_f(file_id, hdferr)

    end subroutine WriteHDF5_Attribute


    ! Read 1D data
    subroutine ReadHDF5_1D(hslab, buffer)
      type(tp_rdwt_hdf5), intent(in) :: hslab
      real(8), intent(out) :: buffer(hslab%dimsf(1),hslab%dnum)
      integer :: n

      if(.not.(IsHDF5Initialized.and.hslab%IsHSInitialized)) then
          print *, 'HDF5 or type hslab NOT intialized !!! STOP'
          stop
      endif

      CALL h5fopen_f(trim(hslab%fname), H5F_ACC_RDONLY_F, file_id, hdferr)
      if(hdferr.eq.-1) then
        print*,'File open failed: ',trim(hslab%fname)
        stop
      endif
      CALL h5gopen_f(file_id, trim(hslab%gname), group_id, hdferr)
      if(hdferr.eq.-1) then
        print*,'Group open failed: ',trim(hslab%gname)
        stop
      endif

      do n=1, hslab%dnum
        CALL h5dopen_f(group_id, trim(hslab%dname(n)), dset_id, hdferr)
        CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, buffer(:,n), hslab%dimsf, hdferr)
        CALL h5dclose_f(dset_id, hdferr)
      enddo

      CALL h5gclose_f(group_id, hdferr)
      CALL h5fclose_f(file_id, hdferr)

    end subroutine ReadHDF5_1D

   ! Read 1D array data
    subroutine ReadHDF5_1D_1v(hslab, buffer)
      type(tp_rdwt_hdf5), intent(in) :: hslab
      real(8), intent(out) :: buffer(hslab%dimsf(1))

      if(.not.(IsHDF5Initialized.and.hslab%IsHSInitialized)) then
          print *, 'HDF5 or type hslab NOT intialized !!! STOP'
          stop
      endif

      CALL h5fopen_f(trim(hslab%fname),H5F_ACC_RDONLY_F, file_id, hdferr)
      CALL h5gopen_f(file_id, trim(hslab%gname), group_id, hdferr)

      CALL h5dopen_f(group_id, trim(hslab%dname(1)), dset_id, hdferr)
      CALL h5dget_space_f(dset_id, dspace_id, hdferr)
      CALL h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, hslab%offset, &
                                 hslab%count, hdferr, hslab%stride, hslab%dimsf)
      CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, buffer, &
                     hslab%dimsf, hdferr, file_space_id=dspace_id, mem_space_id=memspace_id)

      CALL h5sclose_f(dspace_id, hdferr)
      CALL h5dclose_f(dset_id, hdferr)
      CALL h5gclose_f(group_id, hdferr)
      CALL h5fclose_f(file_id, hdferr)

    end subroutine ReadHDF5_1D_1v


    ! Read 2D data
    subroutine ReadHDF5_2D(hslab, buffer)
      type(tp_rdwt_hdf5), intent(in) :: hslab
      real(8), intent(out) :: buffer(hslab%dimsf(1),hslab%dimsf(2),hslab%dnum)
      integer :: n

      if(.not.(IsHDF5Initialized.and.hslab%IsHSInitialized)) then
          print *, 'HDF5 or type hslab NOT intialized !!! STOP'
          stop
      endif

      CALL h5fopen_f(trim(hslab%fname), H5F_ACC_RDONLY_F, file_id, hdferr)
      if(hdferr.eq.-1) then
        print*,'File open failed: ',trim(hslab%fname)
        stop
      endif
      CALL h5gopen_f(file_id, trim(hslab%gname), group_id, hdferr)
      if(hdferr.eq.-1) then
        print*,'Group open failed: ',trim(hslab%gname)
        stop
      endif

      do n=1, hslab%dnum
        CALL h5dopen_f(group_id, trim(hslab%dname(n)), dset_id, hdferr)
        CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, buffer(:,:,n), hslab%dimsf, hdferr)
        CALL h5dclose_f(dset_id, hdferr)
      enddo

      CALL h5gclose_f(group_id, hdferr)
      CALL h5fclose_f(file_id, hdferr)

    end subroutine ReadHDF5_2D

    ! Read 2D data
    subroutine ReadHDF5_2D_1V(hslab, buffer,current_num)
      type(tp_rdwt_hdf5), intent(in) :: hslab
      real(8), intent(out) :: buffer(hslab%dimsf(1),hslab%dimsf(2))
      integer, intent(in) :: current_num
      integer :: n

      if(.not.(IsHDF5Initialized.and.hslab%IsHSInitialized)) then
          print *, 'HDF5 or type hslab NOT intialized !!! STOP'
          stop
      endif

      CALL h5fopen_f(trim(hslab%fname), H5F_ACC_RDONLY_F, file_id, hdferr)
      CALL h5gopen_f(file_id, trim(hslab%gname), group_id, hdferr)

      do n=current_num,current_num
        CALL h5dopen_f(group_id, trim(hslab%dname(n)), dset_id, hdferr)
        CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, buffer, hslab%dimsf, hdferr)
        CALL h5dclose_f(dset_id, hdferr)
      enddo

      CALL h5gclose_f(group_id, hdferr)
      CALL h5fclose_f(file_id, hdferr)

    end subroutine ReadHDF5_2D_1V



    ! Read 3D data
    ! hslab%dimsf and hslab%offset (optional) should be provided
    ! by default offset = 0
    subroutine ReadHDF5_3D_test(hslab, buffer)
      type(tp_rdwt_hdf5), intent(in) :: hslab
      real(8), intent(out) :: buffer(hslab%dimsm(1),hslab%dimsm(2),hslab%dimsm(3),hslab%dnum)
      integer :: n

      if(.not.(IsHDF5Initialized.and.hslab%IsHSInitialized)) then
          print *, 'HDF5 or type hslab NOT intialized !!! STOP'
          stop
      endif

      CALL h5fopen_f(trim(hslab%fname), H5F_ACC_RDONLY_F, file_id, hdferr)
      CALL h5gopen_f(file_id, trim(hslab%gname), group_id, hdferr)
      CALL h5screate_simple_f(hslab%rank, hslab%dimsm, memspace_id, hdferr)
      do n=1, hslab%dnum
        CALL h5dopen_f(group_id, trim(hslab%dname(n)), dset_id, hdferr)
        CALL h5dget_space_f(dset_id, dspace_id, hdferr)
        CALL h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, hslab%offset, &
                                   hslab%count, hdferr, hslab%stride, hslab%block)
        CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, buffer(:,:,:,n), &
                       hslab%dimsm, hdferr, file_space_id=dspace_id, mem_space_id=memspace_id)
        CALL h5sclose_f(dspace_id, hdferr)
        CALL h5dclose_f(dset_id, hdferr)
      enddo
      CALL h5sclose_f(memspace_id, hdferr)
      CALL h5gclose_f(group_id, hdferr)
      CALL h5fclose_f(file_id, hdferr)

    end subroutine ReadHDF5_3D_test



    ! Read 3D data
    ! hslab%dimsf and hslab%offset (optional) should be provided
    ! by default offset = 0
    subroutine ReadHDF5_3D(hslab, buffer)
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
      CALL h5screate_simple_f(hslab%rank, hslab%dimsf, memspace_id, hdferr)
      do n=1, hslab%dnum
        CALL h5dopen_f(group_id, trim(hslab%dname(n)), dset_id, hdferr)
        CALL h5dget_space_f(dset_id, dspace_id, hdferr)
        CALL h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, hslab%offset, &
                                   hslab%count, hdferr, hslab%stride, hslab%dimsf)
        CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, buffer(:,:,:,n), &
                       hslab%dimsf, hdferr, file_space_id=dspace_id, mem_space_id=memspace_id)

        ! CALL h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, hslab%offset, &
        !                             hslab%count, hdferr, hslab%stride, hslab%block)
        ! CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, buffer(:,:,:,n), &
        !                hslab%dimsm, hdferr, file_space_id=dspace_id, mem_space_id=memspace_id)
        CALL h5sclose_f(dspace_id, hdferr)
        CALL h5dclose_f(dset_id, hdferr)
      enddo
      CALL h5sclose_f(memspace_id, hdferr)
      CALL h5gclose_f(group_id, hdferr)
      CALL h5fclose_f(file_id, hdferr)

    end subroutine ReadHDF5_3D


    ! Read 3D data with single precision
    ! hslab%dimsf and hslab%offset (optional) should be provided
    ! by default offset = 0
    subroutine ReadHDF5_3D_SPrecision(hslab, buffer)
      type(tp_rdwt_hdf5), intent(in) :: hslab
      real(4), intent(out) :: buffer(hslab%dimsf(1),hslab%dimsf(2),hslab%dimsf(3),hslab%dnum)
      integer :: n

      if(.not.(IsHDF5Initialized.and.hslab%IsHSInitialized)) then
          print *, 'HDF5 or type hslab NOT intialized !!! STOP'
          stop
      endif
      CALL h5fopen_f(trim(hslab%fname), H5F_ACC_RDONLY_F, file_id, hdferr)
      CALL h5gopen_f(file_id, trim(hslab%gname), group_id, hdferr)
      CALL h5screate_simple_f(hslab%rank, hslab%dimsf, memspace_id, hdferr)
      do n=1, hslab%dnum
        CALL h5dopen_f(group_id, trim(hslab%dname(n)), dset_id, hdferr)
        CALL h5dget_space_f(dset_id, dspace_id, hdferr)
        CALL h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, hslab%offset, &
                                   hslab%count, hdferr, hslab%stride, hslab%dimsf)
        CALL h5dread_f(dset_id, H5T_NATIVE_REAL, buffer(:,:,:,n), &
                       hslab%dimsf, hdferr, file_space_id=dspace_id, mem_space_id=memspace_id)

        CALL h5sclose_f(dspace_id, hdferr)
        CALL h5dclose_f(dset_id, hdferr)
      enddo
      CALL h5sclose_f(memspace_id, hdferr)
      CALL h5gclose_f(group_id, hdferr)
      CALL h5fclose_f(file_id, hdferr)

    end subroutine ReadHDF5_3D_SPrecision


!    ! Read 3D data
!    subroutine ReadHDF5_3D(hslab, buffer)
!      type(tp_rdwt_hdf5), intent(in) :: hslab
!      real(8), intent(out) :: buffer(hslab%dimsf(1),hslab%dimsf(2),hslab%dimsf(3),hslab%dnum)
!      integer :: n
!
!      if(.not.(IsHDF5Initialized.and.hslab%IsHSInitialized)) then
!          print *, 'HDF5 or type hslab NOT intialized !!! STOP'
!          stop
!      endif
!
!      CALL h5fopen_f(trim(hslab%fname), H5F_ACC_RDONLY_F, file_id, hdferr)
!      CALL h5gopen_f(file_id, trim(hslab%gname), group_id, hdferr)
!
!      do n=1, hslab%dnum
!        CALL h5dopen_f(group_id, trim(hslab%dname(n)), dset_id, hdferr)
!        CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, buffer(:,:,:,n), hslab%dimsf, hdferr)
!        CALL h5dclose_f(dset_id, hdferr)
!      enddo
!
!      CALL h5gclose_f(group_id, hdferr)
!      CALL h5fclose_f(file_id, hdferr)
!
!    end subroutine ReadHDF5_3D


    ! Read 3D select data
!    subroutine ReadHDF5_3D_select(hslab, pvol, buffer)
!      type(tp_rdwt_hdf5), intent(in) :: hslab
!      type(tp_DNSIndex), intent(in) :: pvol
!      real(8), intent(out) :: buffer(hslab%dimsf(1),hslab%dimsf(2),hslab%dimsf(3),hslab%dnum)
!      integer(HSIZE_T), dimension(:), allocatable :: dimsm
!      integer(HSIZE_T), dimension(:), allocatable :: count, offset, stride, block
!      integer :: n
!
!      if(.not.(IsHDF5Initialized.and.hslab%IsHSInitialized)) then
!          print *, 'HDF5 or type hslab NOT intialized !!! STOP'
!          stop
!      endif
!      allocate(dimsm(hslab%rank), count(hslab%rank), offset(hslab%rank))
!      allocate(stride(hslab%rank), block(hslab%rank))
!
!      dimsm = hslab%dimsf
!      block = dimsm
!      count = 1
!      stride = 1
!      offset(1) = pvol%kbe - 1
!      offset(2) = pvol%ibe - 1
!      offset(3) = pvol%jbe - 1
!
!      if(offset(1).lt.0.or.offset(2).lt.0.or.offset(3).lt.0) then
!        print *, 'Please the check the input dimension. STOP!'
!        stop
!      endif
!
!      CALL h5fopen_f(trim(hslab%fname), H5F_ACC_RDONLY_F, file_id, hdferr)
!      CALL h5gopen_f(file_id, trim(hslab%gname), group_id, hdferr)
!      CALL h5screate_simple_f(hslab%rank, dimsm, memspace_id, hdferr)
!
!      do n=1, hslab%dnum
!        CALL h5dopen_f(group_id, trim(hslab%dname(n)), dset_id, hdferr)
!        CALL h5dget_space_f(dset_id, dspace_id, hdferr)
!        CALL h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, &
!                                    count, hdferr, stride, block)
!        CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, buffer(:,:,:,n), &
!                       dimsm, hdferr, file_space_id=dspace_id, mem_space_id=memspace_id)
!        CALL h5sclose_f(dspace_id, hdferr)
!        CALL h5dclose_f(dset_id, hdferr)
!      enddo
!      CALL h5sclose_f(memspace_id, hdferr)
!      CALL h5gclose_f(group_id, hdferr)
!      CALL h5fclose_f(file_id, hdferr)
!
!    end subroutine ReadHDF5_3D_select


    ! Read TimeseriesVol data
    ! fvol timeseries volume index
    ! pvol selected range index
    subroutine ReadTimeseriesVol(hslab,pvol,buffer,scalar)
      implicit none
      type(tp_rdwt_hdf5), intent(in) :: hslab
      type(tp_DNSIndex), intent(in) :: pvol
      real(8), intent(out) :: buffer(hslab%dimsf(1),hslab%dimsf(2),hslab%dimsf(3),hslab%dnum)
      real(8), intent(out) :: scalar
      integer(HSIZE_T), dimension(:), allocatable :: dimsm
      integer(HSIZE_T), dimension(:), allocatable :: count, offset, stride, block
      integer :: n

      if(.not.(IsHDF5Initialized.and.hslab%IsHSInitialized)) then
          print *, 'HDF5 or type hslab NOT intialized !!! STOP'
          stop
      endif
      allocate(dimsm(hslab%rank), count(hslab%rank), offset(hslab%rank))
      allocate(stride(hslab%rank), block(hslab%rank))

      dimsm = hslab%dimsf
      block = dimsm
      count = 1
      stride = 1
      offset(1) = pvol%kbe - 1
      offset(2) = pvol%ibe - 1
      offset(3) = pvol%jbe - 1

      if(offset(1).lt.0.or.offset(2).lt.0.or.offset(3).lt.0) then
        print *, 'Please the check the input dimension. STOP!'
        stop
      endif

      CALL h5fopen_f(trim(hslab%fname), H5F_ACC_RDONLY_F, file_id, hdferr)
      ! read scalar time
      CALL h5dopen_f(file_id, trim(hslab%sname), dset_id, hdferr)
      CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, scalar, hslab%dimsf, hdferr)
      CALL h5dclose_f(dset_id, hdferr)

      CALL h5gopen_f(file_id, trim(hslab%gname), group_id, hdferr)
      CALL h5screate_simple_f(hslab%rank, dimsm, memspace_id, hdferr)
      do n=1, hslab%dnum
        CALL h5dopen_f(group_id, trim(hslab%dname(n)), dset_id, hdferr)
        CALL h5dget_space_f(dset_id, dspace_id, hdferr)
        CALL h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, &
                                    count, hdferr, stride, block)
        CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, buffer(:,:,:,n), &
                       dimsm, hdferr, file_space_id=dspace_id, mem_space_id=memspace_id)
        CALL h5sclose_f(dspace_id, hdferr)
        CALL h5dclose_f(dset_id, hdferr)
      enddo ! end n loop
      CALL h5sclose_f(memspace_id, hdferr)
      CALL h5gclose_f(group_id, hdferr)
      CALL h5fclose_f(file_id, hdferr)

    end subroutine ReadTimeseriesVol

    ! Read 4D data
    subroutine ReadHDF5_4D(hslab, buffer)
      type(tp_rdwt_hdf5), intent(in) :: hslab
      real(8), intent(out) :: buffer(hslab%dimsf(1),hslab%dimsf(2),hslab%dimsf(3),hslab%dimsf(4),hslab%dnum)
      integer :: n

      if(.not.(IsHDF5Initialized.and.hslab%IsHSInitialized)) then
          print *, 'HDF5 or type hslab NOT intialized !!! STOP'
          stop
      endif

      CALL h5fopen_f(trim(hslab%fname), H5F_ACC_RDONLY_F, file_id, hdferr)
      CALL h5gopen_f(file_id, trim(hslab%gname), group_id, hdferr)
      CALL h5screate_simple_f(hslab%rank, hslab%dimsf, memspace_id, hdferr)
      do n=1, hslab%dnum
        CALL h5dopen_f(group_id, trim(hslab%dname(n)), dset_id, hdferr)
        CALL h5dget_space_f(dset_id, dspace_id, hdferr)
        CALL h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, hslab%offset, &
                                   hslab%count, hdferr, hslab%stride, hslab%dimsf)
        CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, buffer(:,:,:,:,n), &
                       hslab%dimsf, hdferr, file_space_id=dspace_id, mem_space_id=memspace_id)
        !CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, buffer(:,:,:,:,n), hslab%dimsf, hdferr)
        CALL h5sclose_f(dspace_id, hdferr)
        CALL h5dclose_f(dset_id, hdferr)
      enddo
      CALL h5sclose_f(memspace_id, hdferr)
      CALL h5gclose_f(group_id, hdferr)
      CALL h5fclose_f(file_id, hdferr)

    end subroutine ReadHDF5_4D


    ! Write 1D data
    subroutine WriteHDF5_1D(hslab, buffer)
      type(tp_rdwt_hdf5), intent(in) :: hslab
      real(8), intent(in) :: buffer(hslab%dimsf(1),hslab%dnum)
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
        CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, buffer(:,n),hslab%dimsf, hdferr)
        CALL h5dclose_f(dset_id, hdferr)
      enddo
      CALL h5sclose_f(dspace_id, hdferr)
      if(trim(hslab%gname).ne.'/') then
        CALL h5gclose_f(group_id, hdferr)
      endif
      CALL h5fclose_f(file_id, hdferr)

    end subroutine WriteHDF5_1D

    ! Write 2D data
    subroutine WriteHDF5_2D(hslab, buffer)
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

    end subroutine WriteHDF5_2D

    ! Write 2D data
    subroutine WriteHDF5_2D_1V(hslab, buffer,current_num)
      type(tp_rdwt_hdf5), intent(in) :: hslab
      real(8), intent(in) :: buffer(hslab%dimsf(1),hslab%dimsf(2))
      integer, intent(in) :: current_num
      integer :: n

      if(.not.(IsHDF5Initialized.and.hslab%IsHSInitialized)) then
          print *, 'HDF5 or type hslab NOT intialized !!! STOP'
          stop
      endif

      if(.not.hslab%IsMultiGroup .and. current_num.eq.1) then
        CALL h5fcreate_f(trim(hslab%fname), H5F_ACC_TRUNC_F, file_id, hdferr)
        if(trim(hslab%gname).ne.'/') CALL h5gcreate_f(file_id, trim(hslab%gname),group_id, hdferr)
      else
        CALL h5fopen_f(trim(hslab%fname),H5F_ACC_RDWR_F, file_id, hdferr)
        CALL h5gopen_f(file_id,trim(hslab%gname),group_id,hdferr)
      endif

      CALL h5screate_simple_f(hslab%rank, hslab%dimsf, dspace_id, hdferr)
      do n=current_num, current_num
        if(trim(hslab%gname).ne.'/') then
          CALL h5dcreate_f(group_id, trim(hslab%dname(n)), H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr )
        else
          CALL h5dcreate_f(file_id, trim(hslab%dname(n)), H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr )
        endif
        CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, buffer,hslab%dimsf, hdferr)
        CALL h5dclose_f(dset_id, hdferr)
      enddo
      CALL h5sclose_f(dspace_id, hdferr)
      if(trim(hslab%gname).ne.'/') then
        CALL h5gclose_f(group_id, hdferr)
      endif
      CALL h5fclose_f(file_id, hdferr)

    end subroutine WriteHDF5_2D_1V

    ! Write 3D data
    subroutine WriteHDF5_3D(hslab, buffer)
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

    end subroutine WriteHDF5_3D

    ! Write 3D data
    subroutine WriteHDF5_3D_SPrecision(hslab, buffer)
      type(tp_rdwt_hdf5), intent(in) :: hslab
      real(4), intent(in) :: buffer(hslab%dimsf(1),hslab%dimsf(2),hslab%dimsf(3),hslab%dnum)
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
          CALL h5dcreate_f(group_id, trim(hslab%dname(n)), H5T_NATIVE_REAL, dspace_id, dset_id, hdferr )
        else
          CALL h5dcreate_f(file_id, trim(hslab%dname(n)), H5T_NATIVE_REAL, dspace_id, dset_id, hdferr )
        endif
        CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, buffer(:,:,:,n),hslab%dimsf, hdferr)
        CALL h5dclose_f(dset_id, hdferr)
      enddo
      CALL h5sclose_f(dspace_id, hdferr)
      if(trim(hslab%gname).ne.'/') then
        CALL h5gclose_f(group_id, hdferr)
      endif
      CALL h5fclose_f(file_id, hdferr)

    end subroutine WriteHDF5_3D_SPrecision

    ! Write 4D data
    subroutine WriteHDF5_4D(hslab, buffer)
      type(tp_rdwt_hdf5), intent(in) :: hslab
      real(8), intent(in) :: buffer(hslab%dimsf(1),hslab%dimsf(2),hslab%dimsf(3),hslab%dimsf(4),hslab%dnum)
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
        CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, buffer(:,:,:,:,n),hslab%dimsf, hdferr)
        CALL h5dclose_f(dset_id, hdferr)
      enddo
      CALL h5sclose_f(dspace_id, hdferr)
      if(trim(hslab%gname).ne.'/') then
        CALL h5gclose_f(group_id, hdferr)
      endif
      CALL h5fclose_f(file_id, hdferr)

    end subroutine WriteHDF5_4D

    ! Read 3D one variable data
    subroutine ReadHDF5_3D_1V(hslab, buffer,current_num)
      type(tp_rdwt_hdf5), intent(in) :: hslab
      real(8), intent(out) :: buffer(hslab%dimsm(1),hslab%dimsm(2),hslab%dimsm(3))
      integer, intent(in) :: current_num
      integer :: n

      if(.not.(IsHDF5Initialized.and.hslab%IsHSInitialized)) then
          print *, 'HDF5 or type hslab NOT intialized !!! STOP'
          stop
      endif

      CALL h5fopen_f(trim(hslab%fname),H5F_ACC_RDONLY_F, file_id, hdferr)
      CALL h5gopen_f(file_id, trim(hslab%gname), group_id, hdferr)
      CALL h5screate_simple_f(hslab%rank, hslab%dimsm, memspace_id, hdferr)
      do n=current_num, current_num
        CALL h5dopen_f(group_id, trim(hslab%dname(n)), dset_id, hdferr )
        CALL h5dget_space_f(dset_id, dspace_id, hdferr)
        CALL h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, hslab%offset, &
                                   hslab%count, hdferr, hslab%stride, hslab%block)
        CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, buffer, &
                       hslab%dimsm, hdferr, file_space_id=dspace_id, mem_space_id=memspace_id)
        !CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, buffer(:,:,:),hslab%dimsf, hdferr)
        CALL h5sclose_f(dspace_id, hdferr)
        CALL h5dclose_f(dset_id, hdferr)
      enddo
      CALL h5sclose_f(memspace_id, hdferr)
      CALL h5gclose_f(group_id, hdferr)
      CALL h5fclose_f(file_id, hdferr)

    end subroutine ReadHDF5_3D_1V

    ! Write 3D one variable data
    subroutine WriteHDF5_3D_1V(hslab, buffer,current_num)
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

      if(trim(hslab%gname).ne.'/'.and.current_num.eq.1) then
        CALL h5gcreate_f(file_id, trim(hslab%gname),group_id, hdferr)
      endif
      CALL h5screate_simple_f(hslab%rank, hslab%dimsf, dspace_id, hdferr)
      do n=current_num, current_num
        if(trim(hslab%gname).ne.'/') then
          CALL h5dcreate_f(group_id, trim(hslab%dname(n)), H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr )
        else
          CALL h5dcreate_f(file_id, trim(hslab%dname(n)), H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr )
        endif
        !CALL h5dcreate_f(file_id, trim(hslab%dname(n)), H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr )
        CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, buffer(:,:,:),hslab%dimsf, hdferr)
        CALL h5dclose_f(dset_id, hdferr)
      enddo
      CALL h5sclose_f(dspace_id, hdferr)
      CALL h5fclose_f(file_id, hdferr)

    end subroutine WriteHDF5_3D_1V


    ! write planeinlet data
    subroutine WriteHDF5_planeinlet(hslab, buffer, inew)
      type(tp_rdwt_hdf5), intent(in) :: hslab
      real(8), intent(in) :: buffer(hslab%dimsf(1),hslab%dimsf(2),hslab%dnum,hslab%dimsf(3))
      integer, intent(in) :: inew
      integer :: n

      if(.not.(IsHDF5Initialized.and.hslab%IsHSInitialized)) then
          print *, 'HDF5 or type hslab NOT intialized !!! STOP'
          stop
      endif

      if(inew.eq.1) then
        CALL h5fcreate_f(trim(hslab%fname), H5F_ACC_TRUNC_F, file_id, hdferr)
      else
        CALL h5fopen_f(trim(hslab%fname), H5F_ACC_RDWR_F, file_id, hdferr)
      endif

      CALL h5gcreate_f(file_id, trim(hslab%gname),group_id, hdferr)
      CALL h5screate_simple_f(hslab%rank, hslab%dimsf, dspace_id, hdferr)
      do n=1, hslab%dnum
        CALL h5dcreate_f(group_id, trim(hslab%dname(n)), H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr )
        CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, buffer(:,:,n,:),hslab%dimsf, hdferr)
        CALL h5dclose_f(dset_id, hdferr)
      enddo
      CALL h5sclose_f(dspace_id, hdferr)
      CALL h5gclose_f(group_id,hdferr)
      CALL h5fclose_f(file_id, hdferr)

    end subroutine WriteHDF5_planeinlet


   ! Write 1D array data
    subroutine WriteHDF5_1D_1v(hslab, buffer)
      type(tp_rdwt_hdf5), intent(in) :: hslab
      real(8), intent(in) :: buffer(hslab%dimsf(1))
      integer :: n

      if(.not.(IsHDF5Initialized.and.hslab%IsHSInitialized)) then
          print *, 'HDF5 or type hslab NOT intialized !!! STOP'
          stop
      endif

      CALL h5fopen_f(trim(hslab%fname),H5F_ACC_RDWR_F, file_id, hdferr)

      CALL h5screate_simple_f(hslab%rank, hslab%dimsf, dspace_id, hdferr)

      CALL h5dcreate_f(file_id, trim(hslab%dname(1)), H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr )
      CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, buffer,hslab%dimsf, hdferr)
      CALL h5dclose_f(dset_id, hdferr)

      CALL h5sclose_f(dspace_id, hdferr)
      CALL h5fclose_f(file_id, hdferr)

    end subroutine WriteHDF5_1D_1v



!     ! Detect HDF5 file. return # of datsset, dataset name, rank, dimsf
!     subroutine DetectHDF5(fname,gname,num_dset,dname,rank,dims)
!       character(*), intent(in) :: fname
!       character(*), intent(in) :: gname
!       integer, intent(out) :: num_dset
!       character(*), dimension(:), allocatable, intent(out) :: dname
!       integer, intent(out) :: rank
!       integer, dimension(:), allocatable, intent(out) :: dims
!       integer(HSIZE_T), dimension(:), allocatable :: dimsf, dimsf_total
!       integer :: i, idx
!
!       if(.not.IsHDF5Initialized) then
!         print *, 'HDF5 should be Initialized'
!         stop
!       endif
!
!       if(allocated(dname)) deallocate(dname)
!       if(allocated(dimsf)) deallocate(dimsf)
!       if(allocated(dimsf_total)) deallocate(dimsf_total)
!       if(allocated(dims)) deallocate(dims)
!
!       CALL h5fopen_f(trim(fname), H5F_ACC_RDONLY_F, file_id, hdferr)
!
!       CALL h5gopen_f(file_id, trim(adjustl(gname)), group_id, hdferr)
!       CALL h5gn_members_f(group_id, trim(adjustl(gname)), num_dset, hdferr)
!
!       allocate(dname(num_dset))
!       do idx=0, num_dset-1
!         CALL h5gget_obj_info_idx_f(group_id, trim(adjustl(gname)), idx, dname(idx+1), H5G_DATASET_F, hdferr)
!       enddo
!
!       CALL h5dopen_f(group_id, trim(dname(1)), dset_id, hdferr)
!       CALL h5dget_space_f(dset_id, dspace_id, hdferr)
!       CALL h5sget_simple_extent_ndims_f(dspace_id, rank, hdferr)
!       allocate(dimsf(rank), dimsf_total(rank))
!       CALL h5sget_simple_extent_dims_f(dspace_id, dimsf, dimsf_total, hdferr)
!       CALL h5dclose_f(dset_id, hdferr)
!       CALL h5gclose_f(group_id, hdferr)
!
!       CALL h5fclose_f(file_id, hdferr)
!
!       allocate(dims(rank))
!       dims = dimsf
!
!     end subroutine DetectHDF5
!
!     ! read 1d variables with the dataset name
!     subroutine Read1dHDF5(fname,gname,num_dset,dname,rank,dims,buffer1d)
!       character(*), intent(in) :: fname
!       character(*), intent(in) :: gname
!       integer, intent(in) :: num_dset
!       character(*), intent(in) :: dname(num_dset)
!       integer, intent(in) :: rank
!       integer, intent(in) :: dims(rank)
!       real(8), intent(out) :: buffer1d(:,:)
!       integer(HSIZE_T), dimension(:), allocatable :: dimsf
!       integer :: i
!
!       if(allocated(dimsf)) deallocate(dimsf)
!       allocate(dimsf(rank))
!       dimsf = dims
!
!       CALL h5fopen_f(trim(fname), H5F_ACC_RDONLY_F, file_id, hdferr)
!       CALL h5gopen_f(file_id, trim(adjustl(gname)), group_id, hdferr)
!       do i=1, num_dset
!         CALL h5dopen_f(group_id, trim(dname(i)), dset_id, hdferr)
!         CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, buffer1d(:,i), dimsf, hdferr)
!         CALL h5dclose_f(dset_id, hdferr)
!       enddo
!       CALL h5gclose_f(group_id, hdferr)
!       CALL h5fclose_f(file_id, hdferr)
!
!     end subroutine Read1dHDF5
!
!     ! read 2d selected variables with the dataset name
!     subroutine Read2dHDF5(fname,gname,num_dset,dname,rank,dims,buffer2d)
!       character(*), intent(in) :: fname
!       character(*), intent(in) :: gname
!       integer, intent(in) :: num_dset
!       character(*), intent(in) :: dname(num_dset)
!       integer, intent(in) :: rank
!       integer, intent(in) :: dims(rank)
!       real(8), intent(out) :: buffer2d(:,:,:)
!       integer(HSIZE_T), dimension(:), allocatable :: dimsf
!       integer :: i
!
!       if(allocated(dimsf)) deallocate(dimsf)
!       allocate(dimsf(rank))
!       dimsf = dims
!
!       CALL h5fopen_f(trim(fname), H5F_ACC_RDONLY_F, file_id, hdferr)
!       CALL h5gopen_f(file_id, trim(adjustl(gname)), group_id, hdferr)
!       do i=1, num_dset
!         CALL h5dopen_f(group_id, trim(dname(i)), dset_id, hdferr)
!         CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, buffer2d(:,:,i), dimsf, hdferr)
!         CALL h5dclose_f(dset_id, hdferr)
!       enddo
!       CALL h5gclose_f(group_id, hdferr)
!       CALL h5fclose_f(file_id, hdferr)
!
!     end subroutine Read2dHDF5
!
!     ! read 3d selected variables with the dataset name
!     subroutine Read3dHDF5(fname,gname,num_dset,dname,rank,dims,buffer3d)
!       character(*), intent(in) :: fname
!       character(*), intent(in) :: gname
!       integer, intent(in) :: num_dset
!       character(*), intent(in) :: dname(num_dset)
!       integer, intent(in) :: rank
!       integer, intent(in) :: dims(rank)
!       real(8), intent(out) :: buffer3d(:,:,:,:)
!       integer(HSIZE_T), dimension(:), allocatable :: dimsf
!       integer :: i
!
!       if(allocated(dimsf)) deallocate(dimsf)
!       allocate(dimsf(rank))
!       dimsf = dims
!
!       CALL h5fopen_f(trim(fname), H5F_ACC_RDONLY_F, file_id, hdferr)
!       CALL h5gopen_f(file_id, trim(adjustl(gname)), group_id, hdferr)
!       do i=1, num_dset
!         CALL h5dopen_f(group_id, trim(dname(i)), dset_id, hdferr)
!         CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, buffer3d(:,:,:,i), dimsf, hdferr)
!         CALL h5dclose_f(dset_id, hdferr)
!       enddo
!       CALL h5gclose_f(group_id, hdferr)
!       CALL h5fclose_f(file_id, hdferr)
!
!     end subroutine Read3dHDF5
!
!     subroutine ReadAveAcousticGridMetrics(fname,gname,dims,buffer)
!       character(*), intent(in) :: fname
!       character(*), intent(in) :: gname
!       integer, dimension(:), allocatable, intent(out) :: dims
!       real(8), dimension(:,:,:), allocatable, intent(out) :: buffer
!       integer(HSIZE_T), dimension(:), allocatable :: dimsf, dimsf_total
!       character(10), dimension(:), allocatable :: dname
!       integer :: num_dset
!       integer :: rank
!       integer :: i, j, k, idx
!
!       CALL h5fopen_f(trim(fname), H5F_ACC_RDONLY_F, file_id, hdferr)
!       CALL h5gopen_f(file_id, trim(gname), group_id, hdferr)
!
!       CALL h5gn_members_f(group_id, trim(gname), num_dset, hdferr)
!       allocate(dname(num_dset))
!       do idx=0, num_dset-1
!         CALL h5gget_obj_info_idx_f(group_id, trim(gname), idx, dname(idx+1), H5G_DATASET_F, hdferr)
!       enddo
!
!       do i=1, num_dset
!         CALL h5dopen_f(group_id, trim(dname(i)), dset_id, hdferr)
!         CALL h5dget_space_f(dset_id, dspace_id, hdferr)
!         CALL h5sget_simple_extent_ndims_f(dspace_id, rank, hdferr)
!         if(i.eq.1) allocate(dimsf(rank), dimsf_total(rank))
!         CALL h5sget_simple_extent_dims_f(dspace_id, dimsf, dimsf_total, hdferr)
!         if(i.eq.1) allocate(buffer(dimsf(1),dimsf(2),num_dset))
!         CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, buffer(:,:,i), dimsf, hdferr)
!         CALL h5dclose_f(dset_id, hdferr)
!       enddo
!
!       CALL h5gclose_f(group_id, hdferr)
!       CALL h5fclose_f(file_id, hdferr)
!
!       allocate(dims(rank))
!       dims = dimsf
!
!       print *, 'dimsf_grid(1) = ', dimsf(1)
!       print *, 'dimsf_grid(2) = ', dimsf(2)
!
!     end subroutine ReadAveAcousticGridMetrics
!
!
!     subroutine ReadHDF5sol(fname,idim,jdim,kdim,vars,time)
!       character(*), intent(in) :: fname
!       integer, intent(in) :: idim, jdim, kdim
!       real, intent(out) :: vars(:,:,:,:)
!       real(8), intent(out), optional :: time
!       integer(HSIZE_T), dimension(:), allocatable :: dimsf
!       character(10), dimension(:), allocatable :: dname
!       integer :: rank, num_dset
!       integer :: i, j, k
!
!       if(.not.IsHDF5Initialized) then
!         print *, 'HDF5 should be Initialized'
!         stop
!       endif
!
!       rank = 3
!       num_dset = 6
!       allocate(dname(num_dset))
!       allocate(dimsf(rank))
!
!       dname(1) = "u"
!       dname(2) = "v"
!       dname(3) = "w"
!       dname(4) = "p"
!       dname(5) = "T"
!       dname(6) = "time"
!       dimsf(1) = kdim
!       dimsf(2) = idim
!       dimsf(3) = jdim
!
!       CALL h5fopen_f(trim(fname), H5F_ACC_RDONLY_F, file_id, hdferr)
!       do i=1, num_dset - 1
!         CALL h5dopen_f(file_id, dname(i), dset_id, hdferr)
!         CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, vars(:,:,:,i), dimsf, hdferr)
!         CALL h5dclose_f(dset_id, hdferr)
!       enddo
!
!       if(present(time)) then
!         CALL h5dopen_f(file_id, trim(dname(6)), dset_id, hdferr)
!         CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, time, dimsf, hdferr)
!         CALL h5dclose_f(dset_id, hdferr)
!       endif
!
!       CALL h5fclose_f(file_id, hdferr)
!
!     end subroutine ReadHDF5sol
!
!     subroutine WriteHDF5sol(fname,idim,jdim,kdim,vars,time)
!       character(*), intent(in) :: fname
!       integer, intent(in) :: idim, jdim, kdim
!       real, intent(in) :: vars(:,:,:,:)
!       real(8), intent(in), optional :: time
!       integer(HSIZE_T), dimension(:), allocatable :: dimsf
!       character(10), dimension(:), allocatable :: dname
!       integer :: rank, num_dset
!       integer :: i, j, k
!
!       if(.not.IsHDF5Initialized) then
!         print *, 'HDF5 should be Initialized'
!         stop
!       endif
!
!       rank = 3
!       num_dset = 6
!       allocate(dname(num_dset))
!       allocate(dimsf(rank))
!
!       dname(1) = "u"
!       dname(2) = "v"
!       dname(3) = "w"
!       dname(4) = "p"
!       dname(5) = "T"
!       dname(6) = "time"
!       dimsf(1) = kdim
!       dimsf(2) = idim
!       dimsf(3) = jdim
!
!       CALL h5fcreate_f(trim(fname), H5F_ACC_TRUNC_F, file_id, hdferr)
!       CALL h5screate_simple_f(rank, dimsf, dspace_id, hdferr)
!       do i=1, num_dset - 1
!         CALL h5dcreate_f(file_id, dname(i), H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)
!         CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, vars(:,:,:,i), dimsf, hdferr)
!         CALL h5dclose_f(dset_id, hdferr)
!       enddo
!       CALL h5sclose_f(dspace_id, hdferr)
!
!       if(present(time)) then
!         CALL h5screate_f(H5S_SCALAR_F, dspace_id, hdferr)
!         CALL h5dcreate_f(file_id, dname(num_dset), H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)
!         CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, time, dimsf, hdferr)
!         CALL h5sclose_f(dspace_id, hdferr)
!         CALL h5dclose_f(dset_id, hdferr)
!       endif
!
!       CALL h5fclose_f(file_id, hdferr)
!
!     end subroutine WriteHDF5sol
!
!
!
!!     subroutine ReadtsGrid_jplane(fname,imax,kmax,num_dset,vars)
!!       character(*), intent(in) :: fname
!!       integer, intent(in) :: imax, kmax, num_dset
!!       real, intent(out) :: vars(:,:,:,:)         !!!!
!!       integer :: i, j, k, rank
!!       character(10), dimension(:), allocatable :: dname
!!       integer(HSIZE_T), dimension(:), allocatable :: dimsf
!
!!       rank = 3
!!       allocate(dname(num_dset))
!!       allocate(dimsf(rank))
!
!!       dimsf(1) = imax
!!       dimsf(2) = kmax
!!       dimsf(3) = 1      !!!!1
!
!!       dname(1) = "x"
!!       dname(2) = "y"
!!       dname(3) = "z"
!!       dname(4) = "didx"
!!       dname(5) = "djdx"
!!       dname(6) = "dkdx"
!!       dname(7) = "didy"
!!       dname(8) = "djdy"
!!       dname(9) = "dkdy"
!!       dname(10) = "didz"
!!       dname(11) = "djdz"
!!       dname(12) = "dkdz"
!
!!       CALL h5fopen_f(trim(fname), H5F_ACC_RDONLY_F, file_id, hdferr)
!!       CALL h5gopen_f(file_id,'/jplane',group_id,hdferr)
!!       do i=1, num_dset
!!         CALL h5dopen_f(group_id, dname(i), dset_id, hdferr)
!!         CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, vars(:,:,:,i), dimsf, hdferr)
!!         CALL h5dclose_f(dset_id, hdferr)
!!       enddo
!!       CALL h5gclose_f(group_id, hdferr)
!!       CALL h5fclose_f(file_id, hdferr)
!
!!     end subroutine ReadtsGrid_jplane
!
     ! read timeseries GridMetrics file
     ! file name: fname. group name: gname.
     ! dim1, dim2, dim3 corresponding to the dimsf in the timeseries_GridMetrics.h5
     subroutine ReadtsGridMetrics(fname,gname,dim1,dim2,dim3,vars)
       character(*), intent(in) :: fname
       character(*), intent(in) :: gname
       integer, intent(in) :: dim1, dim2, dim3
       real, intent(out) :: vars(:,:,:,:)         !!!!
       character(10), dimension(:), allocatable :: dname
       integer(HSIZE_T), dimension(:), allocatable :: dimsf
       integer :: num_dset, rank
       integer :: i, j, k

       if(dim1.ne.size(vars,dim=1).or.dim2.ne.size(vars,dim=2).or.dim3.ne.size(vars,dim=3)) then
         print *, 'The input dimension (dim1,dim2,dim3) should be the same with the dimension of vars'
         stop
       endif

       if(.not.IsHDF5Initialized) then
         print *, 'HDF5 should be Initialized'
         stop
       endif

       rank = 3
       num_dset = 12
       allocate(dname(num_dset))
       allocate(dimsf(rank))

       dimsf(1) = dim1
       dimsf(2) = dim2
       dimsf(3) = dim3   !!!!1

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

       CALL h5fopen_f(trim(fname), H5F_ACC_RDONLY_F, file_id, hdferr)
       CALL h5gopen_f(file_id,trim(gname),group_id,hdferr)
       do i=1, num_dset
         CALL h5dopen_f(group_id, dname(i), dset_id, hdferr)
         CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, vars(:,:,:,i), dimsf, hdferr)
         CALL h5dclose_f(dset_id, hdferr)
       enddo
       CALL h5gclose_f(group_id, hdferr)
       CALL h5fclose_f(file_id, hdferr)

     end subroutine ReadtsGridMetrics

     ! read timeseries GridMetrics file
     ! file name: fname. group name: gname.
     ! dim1, dim2, dim3 corresponding to the dimsf in the timeseries_GridMetrics.h5
     subroutine ReadtsGridMetrics_Subset(hslab,vars)
       type(tp_rdwt_hdf5), intent(in) :: hslab
       real, intent(out) :: vars(:,:,:,:)         !!!!
       integer :: num_dset, rank
       integer :: i, j, k

       if(hslab%dimsm(1).ne.size(vars,dim=1).or.hslab%dimsm(2).ne.size(vars,dim=2).or.hslab%dimsm(3).ne.size(vars,dim=3)) then
         print *, 'The input dimension (dim1,dim2,dim3) should be the same with the dimension of vars'
         print*,'Input: ',(hslab%dimsm(i),i=1,3)
         print*,'vars: ',(size(vars,dim=i),i=1,3)
         stop
       endif

       if(hslab%dnum.ne.size(vars,dim=4)) then
           print*,'# of variables is not consistent in vars'
           stop
       endif

       if(.not.IsHDF5Initialized) then
         print *, 'HDF5 should be Initialized'
         stop
       endif

       num_dset = hslab%dnum

       CALL h5fopen_f(trim(hslab%fname), H5F_ACC_RDONLY_F, file_id, hdferr)
       CALL h5gopen_f(file_id,trim(hslab%gname),group_id,hdferr)
       CALL h5screate_simple_f(hslab%rank, hslab%dimsm, memspace_id, hdferr)
       do i=1, num_dset
         CALL h5dopen_f(group_id, hslab%dname(i), dset_id, hdferr)
         CALL h5dget_space_f(dset_id, dspace_id, hdferr)
         CALL h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, hslab%offset, &
                                    hslab%count, hdferr, hslab%stride, hslab%block)
         CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, vars(:,:,:,i), &
                        hslab%dimsm,hdferr,file_space_id=dspace_id,mem_space_id=memspace_id)
         CALL h5sclose_f(dspace_id, hdferr)
         CALL h5dclose_f(dset_id, hdferr)
       enddo
       CALL h5sclose_f(memspace_id, hdferr)
       CALL h5gclose_f(group_id, hdferr)
       CALL h5fclose_f(file_id, hdferr)

     end subroutine ReadtsGridMetrics_Subset

     ! read timeseries flow data.
     ! file name: fname. group name: gname
     ! nt for # of time steps
     ! n: 1 for kplane, 2 for iplane, 3 for jplane
     !subroutine Readtsflow(fname,gname,nt,dim1,dim2,dim3,vars,n)
     subroutine Readtsflow(fname,gname,nt,dim1,dim2,dim3,ioffset,vars,n)
       character(*), intent(in) :: fname
       character(*), intent(in) :: gname
       integer, intent(in) :: nt, dim1, dim2, dim3, ioffset, n
       real(8), intent(out) :: vars(:,:,:,:,:)
       character(10), dimension(:), allocatable :: dname
       integer(HSIZE_T), dimension(:), allocatable :: dimsf
       integer(HSIZE_T), dimension(:), allocatable :: dimsm
       integer(HSIZE_T), dimension(:), allocatable :: count, offset, block, stride
       integer :: num_dset, rank
       integer :: i, j, k

       if(nt.ne.size(vars,dim=1).or.dim1.ne.size(vars,dim=2).or.dim2.ne.size(vars,dim=3).or.dim3.ne.size(vars,dim=4)) then
         print *, 'The input dimension (nt,dim1,dim2,dim3) should be the same with the dimension of vars'
         stop
       endif

       rank = 4
       num_dset = 10
       allocate(dname(num_dset))
       allocate(dimsf(rank))

       dimsf(1) = nt
       dimsf(2) = dim1
       dimsf(3) = dim2
       dimsf(4) = dim3

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
         print *, ' n should be 1: kplane, 2: iplane, 3: jplane'
         stop
       endif

       allocate( dimsm(rank),count(rank))
       allocate(offset(rank),block(rank),stride(rank))
       dimsm = (/nt,dim1,dim2,1/)
       block = dimsm
       count = 1
       stride = 1
       offset(1) = 0
       offset(2) = 0
       offset(3) = 0
       offset(4) = ioffset - 1

       CALL h5fopen_f(trim(fname), H5F_ACC_RDONLY_F, file_id, hdferr)
       CALL h5gopen_f(file_id,trim(gname),group_id,hdferr)
       CALL h5screate_simple_f(rank, dimsm, memspace_id, hdferr)
       do i=1, num_dset
         CALL h5dopen_f(group_id, dname(i), dset_id, hdferr)
         CALL h5dget_space_f(dset_id, dspace_id, hdferr)
         CALL h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, &
                                    count, hdferr, stride, block)
         CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, vars(:,:,:,:,i), dimsm, hdferr,&
                            file_space_id=dspace_id, mem_space_id=memspace_id)
         CALL h5sclose_f(dspace_id, hdferr)
         CALL h5dclose_f(dset_id, hdferr)
       enddo
       CALL h5sclose_f(memspace_id, hdferr)
       CALL h5gclose_f(group_id, hdferr)
       CALL h5fclose_f(file_id, hdferr)

       deallocate(dname)
       deallocate(dimsf)
!       print *, 'min(rho) = ', minval(vars(:,:,:,:,4)/(287.0*vars(:,:,:,:,5)))
     end subroutine Readtsflow

     ! read timeseries flow data.
     ! file name: fname. group name: gname
     ! nt for # of time steps
     ! n: 1 for kplane, 2 for iplane, 3 for jplane
     !subroutine Readtsflow(fname,gname,nt,dim1,dim2,dim3,vars,n)
     subroutine Readtsflow_Subset(hslab,vars)
       type(tp_rdwt_hdf5), intent(in) :: hslab
       real(8), intent(out) :: vars(:,:,:,:,:)
       integer :: num_dset, rank
       integer :: i, j, k

       if(hslab%dimsm(1).ne.size(vars,dim=1).or.hslab%dimsm(2).ne.size(vars,dim=2).or.hslab%dimsm(3).ne.size(vars,dim=3)) then
         print *, 'The input dimension (dim1,dim2,dim3) should be the same with the dimension of vars'
         print*,'Input: ',(hslab%dimsm(i),i=1,3)
         print*,'vars: ',(size(vars,dim=i),i=1,3)
         stop
       endif

       if(hslab%dnum.ne.size(vars,dim=5)) then
           print*,'# of variables is not consistent in vars'
           stop
       endif

       num_dset = hslab%dnum
!
!       dimsf(1) = nt
!       dimsf(2) = dim1
!       dimsf(3) = dim2
!       dimsf(4) = dim3
!
!       dname(1) = "u"
!       dname(2) = "v"
!       dname(3) = "w"
!       dname(4) = "p"
!       dname(5) = "T"
!       if(n.eq.1) then
!         dname(6) = "uk"
!         dname(7) = "vk"
!         dname(8) = "wk"
!         dname(9) = "pk"
!         dname(10) = "Tk"
!       elseif(n.eq.2) then
!         dname(6) = "ui"
!         dname(7) = "vi"
!         dname(8) = "wi"
!         dname(9) = "pi"
!         dname(10) = "Ti"
!       elseif(n.eq.3) then
!         dname(6) = "uj"
!         dname(7) = "vj"
!         dname(8) = "wj"
!         dname(9) = "pj"
!         dname(10) = "Tj"
!       else
!         print *, ' n should be 1: kplane, 2: iplane, 3: jplane'
!         stop
!       endif
!
!       allocate( dimsm(rank),count(rank))
!       allocate(offset(rank),block(rank),stride(rank))
!       dimsm = (/nt,dim1,dim2,1/)
!       block = dimsm
!       count = 1
!       stride = 1
!       offset(1) = 0
!       offset(2) = 0
!       offset(3) = 0
!       offset(4) = ioffset - 1

       CALL h5fopen_f(trim(hslab%fname), H5F_ACC_RDONLY_F, file_id, hdferr)
       CALL h5gopen_f(file_id,trim(hslab%gname),group_id,hdferr)
       CALL h5screate_simple_f(hslab%rank, hslab%dimsm, memspace_id, hdferr)
       do i=1, num_dset
         CALL h5dopen_f(group_id, hslab%dname(i), dset_id, hdferr)
         CALL h5dget_space_f(dset_id, dspace_id, hdferr)
         CALL h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, hslab%offset, &
                                    hslab%count, hdferr, hslab%stride, hslab%block)
         CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, vars(:,:,:,:,i), hslab%dimsm, hdferr,&
                            file_space_id=dspace_id, mem_space_id=memspace_id)
         CALL h5sclose_f(dspace_id, hdferr)
         CALL h5dclose_f(dset_id, hdferr)
       enddo
       CALL h5sclose_f(memspace_id, hdferr)
       CALL h5gclose_f(group_id, hdferr)
       CALL h5fclose_f(file_id, hdferr)

!       print *, 'min(rho) = ', minval(vars(:,:,:,:,4)/(287.0*vars(:,:,:,:,5)))
     end subroutine Readtsflow_Subset

      ! write 1d variables
      subroutine Write1dHDF5(fname,dims,num_dset,dname,buffer1d,iexist)
        character(*), intent(in) :: fname
        integer, intent(in) :: dims(1)    !!!!!!!!!!!!!!!!
        integer, intent(in) :: num_dset
        character(*), intent(in) :: dname(num_dset)
!        real(8), intent(in) :: buffer1d(:,:)        !!!!!
        real(8), intent(in) :: buffer1d(:)
        integer(HSIZE_T), dimension(:), allocatable :: dimsf
        integer :: i, j, k, rank
        integer, intent(in) :: iexist

        rank = 1
        if(allocated(dimsf)) deallocate(dimsf)
        allocate(dimsf(rank))
        dimsf(1) = dims(1)


        if(iexist.eq.0) then
          CALL h5fcreate_f(trim(fname), H5F_ACC_TRUNC_F, file_id, hdferr)
        else
          CALL h5fopen_f(trim(fname), H5F_ACC_RDWR_F, file_id, hdferr )
        endif

        CALL h5screate_simple_f(rank, dimsf, dspace_id, hdferr)
        do i=1, num_dset
          CALL h5dcreate_f(file_id, trim(dname(i)), H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)
          CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, buffer1d, dimsf, hdferr)
          CALL h5dclose_f(dset_id, hdferr)
        enddo

        CALL h5sclose_f(dspace_id, hdferr)
        CALL h5fclose_f(file_id, hdferr)

      end subroutine Write1dHDF5


     subroutine ReadTSiplane(TSiplane,iloc, buffer)
       implicit none
       type(tp_plane), intent(inout) :: TSiplane
       integer, intent(in) :: iloc
       real(8), intent(out) :: buffer(TSiplane%ntpoint_total, TSiplane%nypoint, TSiplane%nzpoint,1,TSiplane%nvar_output)
       !real(8), intent(out), optional :: buffer_grid(nypoint_i, nzpoint_i,1,12)
       integer :: n, nnn, itmp(1), i
       character(4) :: fnum1
       integer :: nn
       integer :: ntpoint_tmp

       ! ********************************
       ! print *, 'iloc = ', iloc
       ! print *, 'npath = ', npath
       ! print *, 'jmints_i = ', jmints_i
       ! print *, 'nfile = ', nfile
       ! print *, 'trim(fileprop(1)%filepath) = ', trim(fileprop(1)%filepath)
       ! print *, 'TSiplane(1,1)%dimsf = ', TSiplane(1,1)%dimsf
       ! print *, 'TSiplane(1,1)%offset = ', TSiplane(1,1)%offset
       ! ********************************

       ntpoint_tmp = 0
       do nn=1, TSiplane%npath
         nnn = TSiplane%jmints
         do n=1, TSiplane%nfile
           write(unit=fnum1,fmt='(I04.4)') iloc
           TSiplane%TSplane(n,nn)%fname = trim(TSiplane%fileprop(nn)%filepath)//'timeseries_iplane'//fnum1//'_'//TSiplane%fname_range(n)//'.h5'
           print *, 'Reading file: ', trim(TSiplane%TSplane(n,nn)%fname)
           call ReadHDF5_4D(TSiplane%TSplane(n,nn), buffer(ntpoint_tmp+1:ntpoint_tmp+TSiplane%TSplane(n,nn)%dimsf(1), nnn:(nnn+TSiplane%TSplane(n,nn)%dimsf(2)-1), &
                              1:TSiplane%TSplane(n,nn)%dimsf(3),1:TSiplane%TSplane(n,nn)%dimsf(4),1:TSiplane%TSplane(n,nn)%dnum) )
           nnn = nnn + TSiplane%TSplane(n,nn)%dimsf(2)
         enddo
         ntpoint_tmp = ntpoint_tmp + TSiplane%TSplane(1,nn)%dimsf(1)
       enddo

!       if(present(buffer_grid)) then
!         TSigrid(1)%fname = trim(fileprop(1)%filepath)//'timeseries_GridMetrics.h5'
!         call ReadTSHDF5_3D(TSigrid(1), buffer_grid)
!       endif

     end subroutine ReadTSiplane

     subroutine ReadTSiplane_grid(TSiplane,buffer_grid,igrid_offset)
       implicit none
       type(tp_plane), intent(inout) :: TSiplane
       real(8), intent(out) :: buffer_grid(TSiplane%nypoint, TSiplane%nzpoint,1,12)
       integer, intent(in) :: igrid_offset
       integer :: n, nnn, itmp(1), i

       TSiplane%TSgrid(1)%fname = trim(TSiplane%fileprop(1)%filepath)//'timeseries_GridMetrics.h5'
       print *, 'reading file: ', trim(TSiplane%TSgrid(1)%fname)
       call ReadHDF5_3D(TSiplane%TSgrid(1), buffer_grid)

     end subroutine ReadTSiplane_grid

     subroutine ReadTSiplane_buffer(TSiplane,iloc, buffer, nt_offset)
       implicit none
       type(tp_plane), intent(inout) :: TSiplane
       integer, intent(in) :: iloc
       real(8), intent(out) :: buffer(TSiplane%fileprop(1)%ntpoint, TSiplane%nypoint, TSiplane%nzpoint,1,TSiplane%nvar_output)
       integer, intent(in) :: nt_offset
       integer :: n, nnn, itmp(1), i
       character(4) :: fnum1
       integer :: nn
       integer :: ntpoint_tmp

       nnn = TSiplane%jmints
       do n=1, TSiplane%nfile
         write(unit=fnum1,fmt='(I04.4)') iloc
         TSiplane%TSplane(n,1)%fname = trim(TSiplane%fileprop(1)%filepath)//'timeseries_iplane'//fnum1//'_'//TSiplane%fname_range(n)//'.h5'
         print *, 'Reading file: ', trim(TSiplane%TSplane(n,1)%fname)
         TSiplane%TSplane(n,1)%offset(1) = nt_offset
         call ReadHDF5_4D(TSiplane%TSplane(n,1), buffer(1:TSiplane%TSplane(n,1)%dimsf(1), nnn:(nnn+TSiplane%TSplane(n,1)%dimsf(2)-1), &
                            1:TSiplane%TSplane(n,1)%dimsf(3),1:TSiplane%TSplane(n,1)%dimsf(4),1:TSiplane%TSplane(n,1)%dnum) )
         nnn = nnn + TSiplane%TSplane(n,1)%dimsf(2)
       enddo

     end subroutine ReadTSiplane_buffer

     subroutine InitTSiplane(nt, np, DNSIndex, jbe_ave_inp, jend_ave_inp, kbe_ave_inp, kend_ave_inp, &
                             nvarout, varindex_output, fp, TSiplane )
       implicit none
       integer, intent(in) :: nt, np
       type(tp_DNSIndex), intent(in) :: DNSIndex
       integer, intent(in) :: jbe_ave_inp, jend_ave_inp, kbe_ave_inp, kend_ave_inp
       integer, intent(in) :: nvarout, varindex_output(nvarout)
       type(tp_plane), intent(inout) :: TSiplane
       type(fprop), intent(in) :: fp(:)
       !integer, intent(in), optional :: ireadgrid !!!
       type(tp_DNSIndex) :: DNSIndex_i

       character(10), dimension(10) :: varname
       integer :: jbe_DNS, jend_DNS, jskip_DNS
       integer :: kbe_DNS, kend_DNS, kskip_DNS
       integer :: jbe_ave, jend_ave, kbe_ave, kend_ave
       integer :: jbuffer !!!
       integer :: n_total, joffset_lower, jnum_lower, jnum_upper
       integer :: jbe_pfile, jend_pfile, n_lower, n_upper
       integer :: n, nn, nnn
       character(4) :: fnum_jbe, fnum_jend
       integer :: nzloc !!
       real(8) :: dys !!!!!!!!!!!

       DNSIndex_i = DNSIndex
       TSiplane%ntpoint_total = nt
       TSiplane%nvar_output = nvarout
       TSiplane%npath = np
       allocate(TSiplane%fileprop(TSiplane%npath))
       TSiplane%fileprop = fp

       varname = (/'u','v','w','p','T','ui','vi','wi','pi','Ti'/)

       jbuffer    = DNSIndex_i%jbuffer
       jbe_DNS    = DNSIndex_i%jbe
       jend_DNS   = DNSIndex_i%jend
       jskip_DNS  = DNSIndex_i%jskip
       kbe_DNS    = DNSIndex_i%kbe
       kend_DNS   = DNSIndex_i%kend
       kskip_DNS  = DNSIndex_i%kskip

       if(jbe_ave_inp.gt.jend_ave_inp.or.kbe_ave_inp.gt.kend_ave_inp.or. &
         jbe_ave_inp.lt.jbe_DNS.or.jend_ave_inp.gt.jend_DNS.or.kbe_ave_inp.lt.kbe_DNS.or.kend_ave_inp.gt.kend_DNS) then
         print *, 'Ranges for averaging is out of bound ... STOP'
         print *, 'jbe_ave_inp = ', jbe_ave_inp, 'jend_ave_inp = ', jend_ave_inp
         print *, 'kbe_ave_inp = ', kbe_ave_inp, 'kend_ave_inp = ', kend_ave_inp
         if(jbe_ave_inp.lt.jbe_DNS) print *, 'jbe_ave_inp < jbe_DNS'
         if(jend_ave_inp.gt.jend_DNS) print *, 'jend_ave_inp > jend_DNS'
         if(kbe_ave_inp.lt.kbe_DNS) print *, 'kbe_ave_inp < kbe_DNS'
         if(kend_ave_inp.gt.kend_DNS) print *, 'kend_ave_inp > kend_DNS'
         stop
       endif

       if(jbe_DNS.gt.jend_DNS.or.kbe_DNS.gt.kend_DNS.or.jskip_DNS.le.0.or.kskip_DNS.le.0.or.jbuffer.le.0) then
         print *, 'Input for DNS timeseries files out of bound ... STOP'
         print *, 'jbe_DNS = ', jbe_DNS, 'jend_DNS = ', jend_DNS, 'jskip_DNS = ', jskip_DNS
         print *, 'kbe_DNS = ', kbe_DNS, 'kend_DNS = ', kend_DNS, 'kskip_DNS = ', kskip_DNS
         print *, 'jbuffer = ', jbuffer
         stop
       endif

       ! begin, end ,and number of spatial points for Averaging

       if(mod(jbe_ave_inp-jbe_DNS,jskip_DNS).ne.0) then
         jbe_ave = jbe_ave_inp + jskip_DNS - mod(jbe_ave_inp-jbe_DNS,jskip_DNS)
       else
         jbe_ave = jbe_ave_inp
       endif
       TSiplane%nypoint = (jend_ave_inp - jbe_ave)/jskip_DNS + 1
       jend_ave = jbe_ave + (TSiplane%nypoint-1)*jskip_DNS

       if(mod(kbe_ave_inp-kbe_DNS,kskip_DNS).ne.0) then
         kbe_ave = kbe_ave_inp + kskip_DNS - mod(kbe_ave_inp-kbe_DNS,kskip_DNS)
       else
         kbe_ave = kbe_ave_inp
       endif
       TSiplane%nzpoint = (kend_ave_inp-kbe_ave)/kskip_DNS + 1
       kend_ave = kbe_ave + (TSiplane%nzpoint-1)*kskip_DNS

       print *, 'DNS Output timeseries spatial index range'
       print *, 'jbe = ', jbe_DNS, 'jend = ', jend_DNS
       print *, 'kbe = ', kbe_DNS, 'kend = ', kend_DNS
       print *, 'Actual Spatial Average range'
       print *, 'jbe_ave = ', jbe_ave, 'jend_ave = ', jend_ave
       print *, 'kbe_ave = ', kbe_ave, 'kend_ave = ', kend_ave
       print *, 'Number of spatial points for Averaging'
       print *, 'nypoint_i = ', TSiplane%nypoint, 'nzpoint_i = ', TSiplane%nzpoint

       ! Dimension range for tsdata
       TSiplane%jmints = 1 - min(0,(jbe_ave-jbe_DNS)/jskip_DNS)
       TSiplane%jmaxts = TSiplane%nypoint + min(0,(jend_DNS-jend_ave)/jskip_DNS)
       TSiplane%javelen = TSiplane%jmaxts - TSiplane%jmints + 1
       print *, 'spanwise Index range for iplane buffer'
       print *, 'jmints_i = ', TSiplane%jmints, 'jmaxts_i = ', TSiplane%jmaxts, 'javelen = ', TSiplane%javelen

       TSiplane%jminloc = max(jbe_ave,jbe_DNS)
       TSiplane%jmaxloc = min(jend_ave,jend_DNS)
       print *, 'Physical j-index range read TS'
       print *, 'jminloc = ', TSiplane%jminloc, 'jmaxloc = ', TSiplane%jmaxloc

       ! Total number of files that cover the full range of DNS output in j dir
       n_total = ceiling( ( (jend_DNS - jbe_DNS)/jskip_DNS + 1 )/real(jbuffer) )

       ! Determine the range of files that need to be read based on jbe_ave & jend_ave
       n_lower = 0; n_upper = 0 ! File indexes that are partially read
       joffset_lower = 0; jnum_lower = 0; jnum_upper = 0
       do n=1, n_total
         jbe_pfile = jbe_DNS + ((n-1)*jbuffer)*jskip_DNS
         jend_pfile = jbe_pfile + (jbuffer-1)*jskip_DNS
         if(jbe_pfile.le.TSiplane%jminloc.and.jend_pfile.ge.TSiplane%jminloc) then
           n_lower = n
           joffset_lower = (TSiplane%jminloc - jbe_pfile)/jskip_DNS
           jnum_lower = min((jend_pfile - TSiplane%jminloc)/jskip_DNS+1, TSiplane%javelen )
         endif
         if(jbe_pfile.le.TSiplane%jmaxloc.and.jend_pfile.ge.TSiplane%jmaxloc) then
           n_upper = n
           if(n_upper.gt.n_lower) jnum_upper = (TSiplane%jmaxloc - jbe_pfile)/jskip_DNS + 1
           exit
         endif
       enddo ! end n loop

       TSiplane%nfile = n_upper - n_lower + 1 ! Total number of files that need to be read

       allocate(TSiplane%fname_range(TSiplane%nfile))

       nn = 0
       do n = n_lower, n_upper
         nn = nn + 1
         jbe_pfile = jbe_DNS + ((n-1)*jbuffer)*jskip_DNS
         jend_pfile = min(jbe_pfile + (jbuffer-1)*jskip_DNS, jend_DNS)
         write(unit=fnum_jbe,fmt='(I04.4)') jbe_pfile
         write(unit=fnum_jend,fmt='(I04.4)') jend_pfile
         TSiplane%fname_range(nn) = 'j'//fnum_jbe//'-'//fnum_jend
         print *, 'fname_jrange = ', TSiplane%fname_range(nn)
       enddo ! end n loop

       allocate(TSiplane%TSplane(TSiplane%nfile,TSiplane%npath))
       TSiplane%TSplane%gname = '/iplane'
       TSiplane%TSplane%rank = 4

       do nnn=1, TSiplane%npath
         do nn=1, TSiplane%nfile
           allocate(TSiplane%TSplane(nn,nnn)%dname(TSiplane%nvar_output))
           allocate(TSiplane%TSplane(nn,nnn)%dimsf(TSiplane%TSplane(nn,nnn)%rank))
           allocate(TSiplane%TSplane(nn,nnn)%dimsm(TSiplane%TSplane(nn,nnn)%rank))
           allocate(TSiplane%TSplane(nn,nnn)%count(TSiplane%TSplane(nn,nnn)%rank))
           allocate(TSiplane%TSplane(nn,nnn)%offset(TSiplane%TSplane(nn,nnn)%rank))
           allocate(TSiplane%TSplane(nn,nnn)%block(TSiplane%TSplane(nn,nnn)%rank))
           allocate(TSiplane%TSplane(nn,nnn)%stride(TSiplane%TSplane(nn,nnn)%rank))
         enddo
       enddo

       ! Files that are fully read
       do nnn=1, TSiplane%npath
         do nn=1, TSiplane%nfile
           TSiplane%TSplane(nn,nnn)%dimsf(1) = TSiplane%fileprop(nnn)%ntpoint
           TSiplane%TSplane(nn,nnn)%dimsf(2) = jbuffer
           TSiplane%TSplane(nn,nnn)%dimsf(3) = (kend_ave - kbe_ave)/kskip_DNS + 1
           TSiplane%TSplane(nn,nnn)%dimsf(4) = 1
           TSiplane%TSplane(nn,nnn)%offset(1) = 0
           TSiplane%TSplane(nn,nnn)%offset(2) = 0
           TSiplane%TSplane(nn,nnn)%offset(3) = kbe_ave - 1
           TSiplane%TSplane(nn,nnn)%offset(4) = 0
         enddo
       enddo

       ! Files that are partially read in j-dir
       do nnn=1, TSiplane%npath
         TSiplane%TSplane(1,nnn)%dimsf(2) = jnum_lower
         TSiplane%TSplane(1,nnn)%offset(2) = joffset_lower
       enddo
       if(TSiplane%nfile.gt.1) then
         do nnn=1, TSiplane%npath
           TSiplane%TSplane(TSiplane%nfile,nnn)%dimsf(2) = jnum_upper
         enddo
       endif

       do nnn=1, TSiplane%npath
         do nn=1, TSiplane%nfile
           TSiplane%TSplane(nn,nnn)%dimsm(1:4) = TSiplane%TSplane(nn,nnn)%dimsf(1:4)
           TSiplane%TSplane(nn,nnn)%block(1:4) = TSiplane%TSplane(nn,nnn)%dimsf(1:4)
           TSiplane%TSplane(nn,nnn)%count = 1
           TSiplane%TSplane(nn,nnn)%stride = 1
         enddo
       enddo

       allocate(TSiplane%varout(TSiplane%nvar_output))
       do nnn=1, TSiplane%npath
         do nn=1, TSiplane%nfile
           do n=1, TSiplane%nvar_output
             TSiplane%TSplane(nn,nnn)%dname(n) = varname(varindex_output(n))
           enddo
           TSiplane%TSplane(nn,nnn)%dnum = TSiplane%nvar_output
         enddo
       enddo
       TSiplane%varout = TSiplane%TSplane(1,1)%dname
       print *, 'Variable name read in: ', (trim(TSiplane%TSplane(1,1)%dname(n))//',', n=1, TSiplane%nvar_output)
       TSiplane%TSplane%IsHSInitialized = .true.

       allocate(TSiplane%TSgrid(1))
       TSiplane%TSgrid%gname = "/iplane"
       TSiplane%TSgrid%rank = 3

       allocate( TSiplane%TSgrid(1)%dname(12))
       allocate( TSiplane%TSgrid(1)%dimsf(TSiplane%TSgrid(1)%rank))
       allocate( TSiplane%TSgrid(1)%dimsm(TSiplane%TSgrid(1)%rank))
       allocate( TSiplane%TSgrid(1)%count(TSiplane%TSgrid(1)%rank))
       allocate( TSiplane%TSgrid(1)%offset(TSiplane%TSgrid(1)%rank))
       allocate( TSiplane%TSgrid(1)%block(TSiplane%TSgrid(1)%rank))
       allocate( TSiplane%TSgrid(1)%stride(TSiplane%TSgrid(1)%rank))
       TSiplane%TSgrid(1)%dnum = 12
       TSiplane%TSgrid(1)%dname(1) = "x"
       TSiplane%TSgrid(1)%dname(2) = "y"
       TSiplane%TSgrid(1)%dname(3) = "z"
       TSiplane%TSgrid(1)%dname(4) = "didx"
       TSiplane%TSgrid(1)%dname(5) = "djdx"
       TSiplane%TSgrid(1)%dname(6) = "dkdx"
       TSiplane%TSgrid(1)%dname(7) = "didy"
       TSiplane%TSgrid(1)%dname(8) = "djdy"
       TSiplane%TSgrid(1)%dname(9) = "dkdy"
       TSiplane%TSgrid(1)%dname(10) = "didz"
       TSiplane%TSgrid(1)%dname(11) = "djdz"
       TSiplane%TSgrid(1)%dname(12) = "dkdz"

       TSiplane%TSgrid(1)%dimsf(1) = TSiplane%nypoint ! jend_ave - jbe_ave + 1
       TSiplane%TSgrid(1)%dimsf(2) = TSiplane%nzpoint ! kend_ave - kbe_ave + 1
       TSiplane%TSgrid(1)%dimsf(3) = 1                       !! should be changed
       TSiplane%TSgrid(1)%offset(1) = jbe_ave - jbe_DNS
       TSiplane%TSgrid(1)%offset(2) = kbe_ave - kbe_DNS
       TSiplane%TSgrid(1)%offset(3) = 0                      !! should be changed

       TSiplane%TSgrid(1)%dimsm(1:3) = TSiplane%TSgrid(1)%dimsf(1:3)
       TSiplane%TSgrid(1)%block(1:3) = TSiplane%TSgrid(1)%dimsf(1:3)
       TSiplane%TSgrid(1)%count = 1
       TSiplane%TSgrid(1)%stride = 1
       TSiplane%TSgrid(1)%IsHSInitialized = .true.

     end subroutine InitTSiplane

















!
!
!
!
!!    end module modReadHDF5
!
!!    module modWriteHDF5
!!      use modVolume  !!!!!!!
!!      use HDF5
!!      implicit none
!
!!      integer(HID_T), private :: file_id, group_id, dset_id, dspace_id
!!      integer, private :: hdferr
!
!!   contains
!
!      ! write 1d variables
!      subroutine Write1dHDF5(fname,dims,num_dset,dname,buffer1d,iexist)
!        character(*), intent(in) :: fname
!        integer, intent(in) :: dims(1)    !!!!!!!!!!!!!!!!
!        integer, intent(in) :: num_dset
!        character(*), intent(in) :: dname(num_dset)
!        real(8), intent(in) :: buffer1d(:,:)        !!!!!
!        integer(HSIZE_T), dimension(:), allocatable :: dimsf
!        integer :: i, j, k, rank
!        integer, intent(in) :: iexist
!
!        rank = 1
!        if(allocated(dimsf)) deallocate(dimsf)
!        allocate(dimsf(rank))
!        dimsf(1) = dims(1)
!
!        if(iexist.eq.0) then
!          CALL h5fcreate_f(trim(fname), H5F_ACC_TRUNC_F, file_id, hdferr)
!        else
!          CALL h5fopen_f(trim(fname), H5F_ACC_RDWR_F, file_id, hdferr )
!        endif
!
!        CALL h5screate_simple_f(rank, dimsf, dspace_id, hdferr)
!        do i=1, num_dset
!          CALL h5dcreate_f(file_id, trim(dname(i)), H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)
!          CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, buffer1d(:,i), dimsf, hdferr)
!          CALL h5dclose_f(dset_id, hdferr)
!        enddo
!
!        CALL h5sclose_f(dspace_id, hdferr)
!        CALL h5fclose_f(file_id, hdferr)
!
!      end subroutine Write1dHDF5
!
!      ! write 2d variables
!      subroutine Write2dHDF5(fname,dims,num_dset,dname,buffer2d,iexist)
!        character(*), intent(in) :: fname
!        integer, intent(in) :: dims(2)
!        integer, intent(in) :: num_dset
!        character(*), intent(in) :: dname(num_dset)
!        real(8), intent(in) :: buffer2d(:,:,:)        !!!!!
!        integer(HSIZE_T), dimension(:), allocatable :: dimsf
!        integer :: i, j, k, rank
!        integer, intent(in) :: iexist
!
!        rank = 2
!        if(allocated(dimsf)) deallocate(dimsf)
!        allocate(dimsf(rank))
!        dimsf = dims
!
!        if(iexist.eq.0) then
!          CALL h5fcreate_f(trim(fname), H5F_ACC_TRUNC_F, file_id, hdferr)
!        else
!          CALL h5fopen_f(trim(fname), H5F_ACC_RDWR_F, file_id, hdferr )
!        endif
!
!        CALL h5screate_simple_f(rank, dimsf, dspace_id, hdferr)
!        do i=1, num_dset
!          CALL h5dcreate_f(file_id, trim(dname(i)), H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)
!          CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, buffer2d(:,:,i), dimsf, hdferr)
!          CALL h5dclose_f(dset_id, hdferr)
!        enddo
!
!        CALL h5sclose_f(dspace_id, hdferr)
!        CALL h5fclose_f(file_id, hdferr)
!
!      end subroutine Write2dHDF5
!
!
!      ! write 3d variables
!      subroutine Write3dHDF5(fname,dims,num_dset,dname,buffer3d,iexist)
!        character(*), intent(in) :: fname
!        integer, intent(in) :: dims(3)
!        integer, intent(in) :: num_dset
!        character(*), intent(in) :: dname(num_dset)
!        real(8), intent(in) :: buffer3d(:,:,:,:)        !!!!!
!        integer(HSIZE_T), dimension(:), allocatable :: dimsf
!        integer :: i, rank
!        integer, intent(in) :: iexist
!
!        rank = 3
!        if(allocated(dimsf)) deallocate(dimsf)
!        allocate(dimsf(rank))
!        dimsf = dims
!
!        if(iexist.eq.0) then
!          CALL h5fcreate_f(trim(fname), H5F_ACC_TRUNC_F, file_id, hdferr)
!        else
!          CALL h5fopen_f(trim(fname), H5F_ACC_RDWR_F, file_id, hdferr )
!        endif
!
!        CALL h5screate_simple_f(rank, dimsf, dspace_id, hdferr)
!        do i=1, num_dset
!          CALL h5dcreate_f(file_id, trim(dname(i)), H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)
!          CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, buffer3d(:,:,:,i), dimsf, hdferr)
!          CALL h5dclose_f(dset_id, hdferr)
!        enddo
!
!        CALL h5sclose_f(dspace_id, hdferr)
!        CALL h5fclose_f(file_id, hdferr)
!
!      end subroutine Write3dHDF5
!
!
!
!      ! write 3d single variables
!      subroutine Write3dHDF5_svariable(fname,dims,dname,vars,iexist)
!        character(*), intent(in) :: fname
!        integer, intent(in) :: dims(:)
!        character(*), intent(in) :: dname
!        real, intent(in) :: vars(:,:,:)
!        integer, intent(in) :: iexist
!        integer(HSIZE_T) :: dimsf(3)
!        integer :: rank = 3
!
!        dimsf = dims
!
!        if(iexist.eq.0) then
!          CALL h5fcreate_f(trim(fname), H5F_ACC_TRUNC_F, file_id, hdferr)
!        else
!          CALL h5fopen_f(trim(fname), H5F_ACC_RDWR_F, file_id, hdferr )
!        endif
!
!        CALL h5screate_simple_f(rank, dimsf, dspace_id, hdferr)
!        CALL h5dcreate_f(file_id, trim(dname), H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)
!        CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, vars, dimsf, hdferr)
!        CALL h5sclose_f(dspace_id, hdferr)
!        CALL h5dclose_f(dset_id, hdferr)
!
!        CALL h5fclose_f(file_id, hdferr)
!
!      end subroutine Write3dHDF5_svariable


!    end module modWriteHDF5






 ! read all the 1d variables without the dataset name.
!     subroutine Read1dHDF5(fname,gname,num_dset,dname,buffer1d)
!       character(*), intent(in) :: fname
!       character(*), intent(in) :: gname
!       integer, intent(out) :: num_dset
!       character(*), dimension(:), allocatable, intent(out) :: dname
!       real(8), dimension(:,:), allocatable, intent(out) :: buffer1d
!       integer(HSIZE_T), dimension(:), allocatable :: dimsf, dimsf_total
!       integer :: rank
!       integer :: i, j, k, idx

!       if(allocated(dname)) deallocate(dname)
!       if(allocated(dimsf)) deallocate(dimsf)
!       if(allocated(dimsf_total)) deallocate(dimsf_total)
!       if(allocated(buffer1d)) deallocate(buffer1d)

!       CALL h5fopen_f(trim(fname), H5F_ACC_RDONLY_F, file_id, hdferr)
!       CALL h5gopen_f(file_id, trim(gname), group_id, hdferr)

!       CALL h5gn_members_f(group_id, trim(gname), num_dset, hdferr)
!       allocate(dname(num_dset))
!       do idx=0, num_dset-1
!         CALL h5gget_obj_info_idx_f(group_id, trim(gname), idx, dname(idx+1), H5G_DATASET_F, hdferr)
!       enddo

!       do i=1, num_dset
!         CALL h5dopen_f(group_id, trim(dname(i)), dset_id, hdferr)
!         CALL h5dget_space_f(dset_id, dspace_id, hdferr)
!         CALL h5sget_simple_extent_ndims_f(dspace_id, rank, hdferr)
!         if(i.eq.1) allocate(dimsf(rank), dimsf_total(rank))
!         CALL h5sget_simple_extent_dims_f(dspace_id, dimsf, dimsf_total, hdferr)
!         if(i.eq.1) allocate(buffer1d(dimsf(1), num_dset))
!         CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, buffer1d(:,i), dimsf, hdferr)
!         CALL h5dclose_f(dset_id, hdferr)
!       enddo

!       CALL h5gclose_f(group_id, hdferr)
!       CALL h5fclose_f(file_id, hdferr)

!     end subroutine Read1dHDF5

 ! read all the 2d variables without the dataset name.
!     subroutine Read2dHDF5(fname,gname,num_dset,dname,buffer2d)
!       character(*), intent(in) :: fname
!       character(*), intent(in) :: gname
!       integer, intent(out) :: num_dset
!       character(*), dimension(:), allocatable, intent(out) :: dname
!       real(8), dimension(:,:,:), allocatable, intent(out) :: buffer2d
!       integer(HSIZE_T), dimension(:), allocatable :: dimsf, dimsf_total
!       integer :: rank
!       integer :: i, j, k, idx

!       if(allocated(dname)) deallocate(dname)
!       if(allocated(dimsf)) deallocate(dimsf)
!       if(allocated(dimsf_total)) deallocate(dimsf_total)
!       if(allocated(buffer2d)) deallocate(buffer2d)

!       CALL h5fopen_f(trim(fname), H5F_ACC_RDONLY_F, file_id, hdferr)
!       CALL h5gopen_f(file_id, trim(gname), group_id, hdferr)

!       CALL h5gn_members_f(group_id, trim(gname), num_dset, hdferr)
!       allocate(dname(num_dset))
!       do idx=0, num_dset-1
!         CALL h5gget_obj_info_idx_f(group_id, trim(gname), idx, dname(idx+1), H5G_DATASET_F, hdferr)
!       enddo

!       do i=1, num_dset
!         CALL h5dopen_f(group_id, trim(dname(i)), dset_id, hdferr)
!         CALL h5dget_space_f(dset_id, dspace_id, hdferr)
!         CALL h5sget_simple_extent_ndims_f(dspace_id, rank, hdferr)
!         if(i.eq.1) allocate(dimsf(rank), dimsf_total(rank))
!         CALL h5sget_simple_extent_dims_f(dspace_id, dimsf, dimsf_total, hdferr)
!         if(i.eq.1) allocate(buffer2d(dimsf(1),dimsf(2),num_dset))
!         CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, buffer2d(:,:,i), dimsf, hdferr)
!         CALL h5dclose_f(dset_id, hdferr)
!       enddo

!       CALL h5gclose_f(group_id, hdferr)
!       CALL h5fclose_f(file_id, hdferr)

!     end subroutine Read2dHDF5

 ! read all the 3d variables without the dataset name.
!     subroutine Read3dHDF5(fname,gname,num_dset,dname,buffer3d)
!       character(*), intent(in) :: fname
!       character(*), intent(in) :: gname
!       integer, intent(out) :: num_dset
!       character(*), dimension(:), allocatable, intent(out) :: dname
!       real(8), dimension(:,:,:,:), allocatable, intent(out) :: buffer3d
!       integer(HSIZE_T), dimension(:), allocatable :: dimsf, dimsf_total
!       integer :: rank
!       integer :: i, j, k, idx

!       if(allocated(dname)) deallocate(dname)
!       if(allocated(dimsf)) deallocate(dimsf)
!       if(allocated(dimsf_total)) deallocate(dimsf_total)
!       if(allocated(buffer3d)) deallocate(buffer3d)

!       CALL h5fopen_f(trim(fname), H5F_ACC_RDONLY_F, file_id, hdferr)
!       CALL h5gopen_f(file_id, trim(gname), group_id, hdferr)

!       CALL h5gn_members_f(group_id, trim(gname), num_dset, hdferr)
!       allocate(dname(num_dset))
!       do idx=0, num_dset-1
!         CALL h5gget_obj_info_idx_f(group_id, trim(gname), idx, dname(idx+1), H5G_DATASET_F, hdferr)
!       enddo

!       do i=1, num_dset
!         CALL h5dopen_f(group_id, trim(dname(i)), dset_id, hdferr)
!         CALL h5dget_space_f(dset_id, dspace_id, hdferr)
!         CALL h5sget_simple_extent_ndims_f(dspace_id, rank, hdferr)
!         if(i.eq.1) allocate(dimsf(rank), dimsf_total(rank))
!         CALL h5sget_simple_extent_dims_f(dspace_id, dimsf, dimsf_total, hdferr)
!         if(i.eq.1) allocate(buffer3d(dimsf(1),dimsf(2),dimsf(3),num_dset))
!         CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, buffer3d(:,:,:,i), dimsf, hdferr)
!         CALL h5dclose_f(dset_id, hdferr)
!       enddo

!       CALL h5gclose_f(group_id, hdferr)
!       CALL h5fclose_f(file_id, hdferr)

!     end subroutine Read3dHDF5




end module modRWHDF5

