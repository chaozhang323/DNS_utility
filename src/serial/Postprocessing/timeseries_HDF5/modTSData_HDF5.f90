module MTSHDF5
   use HDF5
   implicit none
   integer :: hdferr
   logical :: IsHDF5Initialized = .false.
   type tp_hyperslab
        integer :: rank
        integer :: dnum ! number of variables in the given group with group name 'gname'
        character(400) :: fname
        character(100) :: gname
        character(20), dimension(:), allocatable :: dname ! should have a dimension of 'dnum'
        integer(HSIZE_T), dimension(:), allocatable :: dimsf
        integer(HSIZE_T), dimension(:), allocatable :: dimsm  !!!!!!
        integer(HSIZE_T), dimension(:), allocatable :: count, offset, block, stride
        logical :: IsHSInitialized = .false.
        logical :: IsMultiGroup = .false.
        integer :: anum ! number of attribute variables
        character(100), dimension(:), allocatable :: attr_name
        integer(HSIZE_T), dimension(1) :: attr_dims_no=(/1/)
        integer :: dims_integer, dims_real

   end type tp_hyperslab

   type fprop
     character(400) :: filepath
     integer :: ntpoint, nskip
     integer :: file_be, file_end, file_skip
   end type fprop

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

     real(8) :: uinf, rhoinf, Tinf
     integer :: iFileType

!     type(fprop), dimension(:), allocatable :: fileprop
!     type(tp_DNSIndex) :: DNSIndex_k, DNSIndex_j, DNSIndex_i

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

   subroutine ReadTSHDF5_4D(TSplane, buffer)
      type(tp_hyperslab), intent(in) :: TSplane
      real(8), intent(out) :: buffer(TSplane%dimsf(1),TSplane%dimsf(2),TSplane%dimsf(3),TSplane%dimsf(4),TSplane%dnum)

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
         CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, buffer(:,:,:,:,n),  &
                             TSplane%dimsf, hdferr, file_space_id=dspace_id, mem_space_id=memspace)
         CALL h5sclose_f(dspace_id, hdferr)
         CALL h5dclose_f(dset_id, hdferr)
     enddo
     CALL h5sclose_f(memspace, hdferr)
     CALL h5gclose_f(group_id, hdferr)
     CALL h5fclose_f(file_id, hdferr)

   end subroutine ReadTSHDF5_4D

   subroutine ReadTSHDF5_3D(TSplane, buffer)
      type(tp_hyperslab), intent(in) :: TSplane
      real(8), intent(out) :: buffer(TSplane%dimsf(1),TSplane%dimsf(2),TSplane%dimsf(3),TSplane%dnum)

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
         CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, buffer(:,:,:,n),  &
                             TSplane%dimsf, hdferr, file_space_id=dspace_id, mem_space_id=memspace)
         CALL h5sclose_f(dspace_id, hdferr)
         CALL h5dclose_f(dset_id, hdferr)
     enddo
     CALL h5sclose_f(memspace, hdferr)
     CALL h5gclose_f(group_id, hdferr)
     CALL h5fclose_f(file_id, hdferr)

   end subroutine ReadTSHDF5_3D

   subroutine ReadTSHDF5_1D(TSplane, buffer)
      type(tp_hyperslab), intent(in) :: TSplane
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



   subroutine ReadTSHDF5_1D_integer(TSplane, buffer)
      type(tp_hyperslab), intent(in) :: TSplane
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

   subroutine WriteTSHDF5_2D(TSplane, buffer)
      type(tp_hyperslab), intent(in) :: TSplane
      real(8), intent(in) :: buffer(TSplane%dimsf(1),TSplane%dimsf(2),TSplane%dnum)

      integer :: n
      integer(HID_T) :: file_id, group_id, memspace
      integer(HID_T) :: dspace_id, dset_id

      if(.not.(IsHDF5Initialized.and.TSplane%IsHSInitialized)) then
          print *, 'HDF5 or ty_Hyperslab NOT intialized !!! STOP in WriteTSHDF5_2D'
          stop
      endif

      if(.not.TSplane%IsMultiGroup) then
        CALL h5fcreate_f(trim(TSplane%fname), H5F_ACC_TRUNC_F, file_id, hdferr)
      else
        CALL h5fopen_f(trim(TSplane%fname), H5F_ACC_RDWR_F, file_id, hdferr)
      endif
      if(trim(TSplane%gname).ne.'/') then
        CALL h5gcreate_f(file_id,trim(TSplane%gname), group_id, hdferr)
      endif
      CALL h5screate_simple_f(TSplane%rank, TSplane%dimsm, memspace, hdferr)
      CALL h5screate_simple_f(TSplane%rank, TSplane%dimsf, dspace_id, hdferr)
      do n = 1, TSplane%dnum
        if(trim(TSplane%gname).ne.'/') then
          CALL h5dcreate_f(group_id, trim(TSplane%dname(n)), &
                           H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)
        else
          CALL h5dcreate_f(file_id, trim(TSplane%dname(n)), &
                           H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)
        endif
        CALL h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, TSplane%offset, &
                                   TSplane%count, hdferr, TSplane%stride, TSplane%block)
        CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, buffer(:,:,n),  &
                            TSplane%dimsf, hdferr, file_space_id=dspace_id, mem_space_id=memspace)
        CALL h5dclose_f(dset_id, hdferr)
      enddo
      CALL h5sclose_f(dspace_id, hdferr)
      CALL h5sclose_f(memspace, hdferr)
      if(trim(TSplane%gname).ne.'/') then
        CALL h5gclose_f(group_id, hdferr)
      endif
      CALL h5fclose_f(file_id, hdferr)

   end subroutine WriteTSHDF5_2D


   subroutine WriteTSHDF5_1D(TSplane, buffer)
      type(tp_hyperslab), intent(in) :: TSplane
      real(8), intent(in) :: buffer(TSplane%dimsf(1),TSplane%dnum)

      integer :: n
      integer(HID_T) :: file_id, group_id, memspace
      integer(HID_T) :: dspace_id, dset_id

      if(.not.(IsHDF5Initialized.and.TSplane%IsHSInitialized)) then
          print *, 'HDF5 or ty_Hyperslab NOT intialized !!! STOP in WriteTSHDF5_1D'
          stop
      endif

      if(.not.TSplane%IsMultiGroup) then
        CALL h5fcreate_f(trim(TSplane%fname), H5F_ACC_TRUNC_F, file_id, hdferr)
      else
        CALL h5fopen_f(trim(TSplane%fname), H5F_ACC_RDWR_F, file_id, hdferr)
      endif
      if(trim(TSplane%gname).ne.'/') then
        CALL h5gcreate_f(file_id,trim(TSplane%gname), group_id, hdferr)
      endif

      CALL h5screate_simple_f(TSplane%rank, TSplane%dimsm, memspace, hdferr)
      CALL h5screate_simple_f(TSplane%rank, TSplane%dimsf, dspace_id, hdferr)
      do n = 1, TSplane%dnum
         if(trim(TSplane%gname).ne.'/') then
           CALL h5dcreate_f(group_id, trim(TSplane%dname(n)), &
                            H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)
         else
           CALL h5dcreate_f(file_id, trim(TSplane%dname(n)), &
                            H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)
         endif
         CALL h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, TSplane%offset, &
                                    TSplane%count, hdferr, TSplane%stride, TSplane%block)
         CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, buffer(:,n),  &
                             TSplane%dimsf, hdferr, file_space_id=dspace_id, mem_space_id=memspace)
         CALL h5dclose_f(dset_id, hdferr)
     enddo
     CALL h5sclose_f(dspace_id, hdferr)
     CALL h5sclose_f(memspace, hdferr)
     if(trim(TSplane%gname).ne.'/') then
       CALL h5gclose_f(group_id, hdferr)
     endif
     CALL h5fclose_f(file_id, hdferr)

   end subroutine WriteTSHDF5_1D


   subroutine WriteTSHDF5_Attribute(TSplane,buffer_integer,buffer_real)
     type(tp_hyperslab), intent(in) :: TSplane
     integer, intent(in) :: buffer_integer(TSplane%dims_integer)
     real(8), intent(in) :: buffer_real(TSplane%dims_real)
     integer(HID_T) :: file_id, group_id
     integer(HID_T) :: attr_id, attr_space, atype_id
     integer :: i

     CALL h5fopen_f(trim(TSplane%fname), H5F_ACC_RDWR_F, file_id, hdferr)
     CALL h5gopen_f(file_id, trim(TSplane%gname), group_id, hdferr)
     ! write integer attr
     do i=1, TSplane%dims_integer
       CALL h5screate_simple_f(1, TSplane%attr_dims_no, attr_space, hdferr)
       CALL h5tcopy_f(H5T_NATIVE_INTEGER, atype_id, hdferr)
       CALL h5acreate_f(group_id,TSplane%attr_name(i), atype_id,attr_space,attr_id,hdferr)
       CALL h5awrite_f(attr_id,H5T_NATIVE_INTEGER,buffer_integer(i),TSplane%attr_dims_no,hdferr)
       CALL h5sclose_f(attr_space, hdferr)
       CALL h5aclose_f(attr_id, hdferr)
     enddo

     ! write real attr
     do i=1, TSplane%dims_real
        CALL h5screate_simple_f(1, TSplane%attr_dims_no, attr_space, hdferr)
       CALL h5tcopy_f(H5T_NATIVE_DOUBLE, atype_id, hdferr)
       CALL h5acreate_f(group_id,TSplane%attr_name(i+TSplane%dims_integer), atype_id,attr_space,attr_id,hdferr)
       CALL h5awrite_f(attr_id,H5T_NATIVE_DOUBLE,buffer_real(i),TSplane%attr_dims_no,hdferr)
       CALL h5sclose_f(attr_space, hdferr)
       CALL h5aclose_f(attr_id, hdferr)
     enddo
     CALL h5gclose_f(group_id, hdferr)
     CALL h5fclose_f(file_id, hdferr)

   end subroutine WriteTSHDF5_Attribute

    subroutine ReadHDF5_Attribute(hslab, buffer)
      implicit none
      type(tp_hyperslab), intent(in) :: hslab
      integer, intent(out) :: buffer(3,hslab%anum)
      integer(HID_T) :: attr_id, attr_space, atype_id
      integer(HID_T) :: file_id, group_id
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


end module MTSHDF5


module MTSVOL
     use MTSHDF5
     implicit none

     type(tp_hyperslab), private :: TSVol, TSVolgrid
     integer, private :: nfile ! Total number of files to read data from
     integer, private :: ntpoint_total, nvar_output, iminloc, imaxloc, iavelen, kminloc, kmaxloc, kavelen
     character(200), private :: filepath
     type(tp_DNSIndex), private :: DNSIndex_Vol
     type(fprop), private :: fileprop

     ! public variables
     character(10), dimension(5) :: varname_Vol
     character(10), dimension(:), allocatable :: varout_Vol
     integer :: nxpoint_Vol, nypoint_Vol, nzpoint_Vol, imints_Vol, imaxts_Vol, kmints_Vol, kmaxts_Vol
contains

!     subroutine ReadDNS_index_Vol(nt,dt,DNSIndex,filepath)
!       integer, intent(out) :: nt
!       real(8), intent(out) :: dt
!       type(tp_DNSIndex), intent(out) :: DNSIndex
!       character(*), intent(in) :: filepath
!       type(tp_hyperslab) :: Index_Vol
!       real(8) :: time1, time2
!
!     end subroutine ReadDNS_index_Vol

     subroutine ReadTSVol(buffer)
       implicit none
       real(8), intent(out) :: buffer(ntpoint_total,kmints_Vol:kmaxts_Vol,imints_Vol:imaxts_Vol,nypoint_Vol,nvar_output)
       character(8) :: fnum1
       integer :: n, nn, ntpoint_tmp, num_file

       ntpoint_tmp = 1
       num_file = (fileprop%file_end - fileprop%file_be)/fileprop%file_skip + 1

       do n=1, num_file
         write(unit=fnum1,fmt='(I08.8)') fileprop%file_be + (n-1)*fileprop%file_skip
         TSVol%fname = trim(fileprop%filepath)//'timeseriesVol_'//fnum1//'.h5'
         print *, 'reading file: ', trim(TSVol%fname)
         call ReadTSHDF5_3D(TSVol,buffer(ntpoint_tmp,kmints_Vol:kmaxts_Vol,imints_Vol:imaxts_Vol,1:nypoint_Vol,1:TSVol%dnum)  )
         ntpoint_tmp = ntpoint_tmp + 1
       enddo

     end subroutine ReadTSVol


     subroutine InitTSVol(nt, np, DNSIndex, ibe_ave_inp, iend_ave_inp, jbe_ave_inp, jend_ave_inp, kbe_ave_inp, kend_ave_inp, &
                          iwinl, iwinr, kwind, kwinu, nvarout, varindex_output, fp, ireadgrid )
       integer, intent(in) :: nt, np
       type(tp_DNSIndex), intent(in) :: DNSIndex
       integer, intent(in) :: ibe_ave_inp, iend_ave_inp, jbe_ave_inp, jend_ave_inp, kbe_ave_inp, kend_ave_inp, iwinl, iwinr, kwind, kwinu
       integer, intent(in) :: nvarout, varindex_output(nvarout)
       type(fprop), intent(in) :: fp
       integer, intent(in), optional :: ireadgrid

       integer :: ibe_DNS, iend_DNS, iskip_DNS
       integer :: jbe_DNS, jend_DNS, jskip_DNS
       integer :: kbe_DNS, kend_DNS, kskip_DNS
       integer :: ibe_ave, iend_ave, jbe_ave, jend_ave, kbe_ave, kend_ave
       integer :: n_total
       integer :: n

       DNSIndex_Vol = DNSIndex
       ntpoint_total = nt
       nvar_output = nvarout

       fileprop = fp

       varname_Vol(1) = 'u'
       varname_Vol(2) = 'v'
       varname_Vol(3) = 'w'
       varname_Vol(4) = 'p'
       varname_Vol(5) = 'T'

       ibe_DNS    = DNSIndex_Vol%ibe
       iend_DNS   = DNSIndex_Vol%iend
       iskip_DNS  = 1 !DNSIndex_Vol%iskip
       jbe_DNS    = DNSIndex_Vol%jbe
       jend_DNS   = DNSIndex_Vol%jend
       jskip_DNS  = 1 !DNSIndex_Vol%jskip
       kbe_DNS    = DNSIndex_Vol%kbe
       kend_DNS   = DNSIndex_Vol%kend
       kskip_DNS  = 1 !DNSIndex_Vol%kskip

       if(ibe_ave_inp.gt.iend_ave_inp.or.jbe_ave_inp.gt.jend_ave_inp.or. &
          kbe_ave_inp.gt.kend_ave_inp.or. &
          ibe_ave_inp.lt.ibe_DNS.or.iend_ave_inp.gt.iend_DNS.or. &
          jbe_ave_inp.lt.jbe_DNS.or.jend_ave_inp.gt.jend_DNS.or. &
          kbe_ave_inp.lt.kbe_DNS.or.kend_ave_inp.gt.kend_DNS ) then
         print *, 'Ranges for averaging is out of bound ... STOP'
         print *, 'ibe_ave_inp =', ibe_ave_inp, 'iend_ave_inp = ', iend_ave_inp
         print *, 'jbe_ave_inp =', jbe_ave_inp, 'jend_ave_inp = ', jend_ave_inp
         print *, 'kbe_ave_inp =', kbe_ave_inp, 'kend_ave_inp = ', kend_ave_inp
         if(ibe_ave_inp.lt.ibe_DNS) print *, 'ibe_ave_inp < ibe_DNS'
         if(iend_ave_inp.gt.iend_DNS) print *, 'iend_ave_inp > iend_DNS'
         if(jbe_ave_inp.lt.jbe_DNS) print *, 'jbe_ave_inp < jbe_DNS'
         if(jend_ave_inp.gt.jend_DNS) print *, 'jend_ave_inp > jend_DNS'
         if(kbe_ave_inp.lt.kbe_DNS) print *, 'kbe_ave_inp < kbe_DNS'
         if(kend_ave_inp.gt.kend_DNS) print *, 'kend_ave_inp > kend_DNS'
         stop
       endif

       if(ibe_DNS.gt.iend_DNS.or.jbe_DNS.gt.jend_DNS.or.kbe_DNS.gt.kend_DNS.or.iskip_DNS.le.0.or.jskip_DNS.le.0.or.kskip_DNS.le.0) then
         print *, 'Input for DNS timeseries files out of bound ... STOP'
         print *, 'ibe_DNS =', ibe_DNS, 'iend_DNS = ', iend_DNS, 'iskip_DNS =', iskip_DNS
         print *, 'jbe_DNS =', jbe_DNS, 'jend_DNS = ', jend_DNS, 'jskip_DNS =', jskip_DNS
         print *, 'kbe_DNS =', kbe_DNS, 'kend_DNS = ', kend_DNS, 'kskip_DNS =', kskip_DNS
         stop
       endif

       ! begin, end and number of spatial points for Averaging
       if(mod(ibe_ave_inp-ibe_DNS,iskip_DNS).ne.0) then
         ibe_ave = ibe_ave_inp + iskip_DNS  - mod(ibe_ave_inp-ibe_DNS,iskip_DNS)
       else
         ibe_ave = ibe_ave_inp
       endif
       nxpoint_Vol = (iend_ave_inp-ibe_ave)/iskip_DNS +1
       iend_ave = ibe_ave + (nxpoint_Vol-1)*iskip_DNS

       if(mod(jbe_ave_inp-jbe_DNS,jskip_DNS).ne.0) then
         jbe_ave = jbe_ave_inp + jskip_DNS - mod(jbe_ave_inp-jbe_DNS,jskip_DNS)
       else
         jbe_ave = jbe_ave_inp
       endif
       nypoint_Vol = (jend_ave_inp-jbe_ave)/jskip_DNS + 1
       jend_ave = jbe_ave + (nypoint_Vol-1)*jskip_DNS

       if(mod(kbe_ave_inp-kbe_DNS,kskip_DNS).ne.0) then
         kbe_ave = kbe_ave_inp + kskip_DNS - mod(kbe_ave_inp-kbe_DNS,kskip_DNS)
       else
         kbe_ave = kbe_ave_inp
       endif
       nzpoint_Vol = (kend_ave_inp-kbe_ave)/kskip_DNS + 1
       kend_ave = kbe_ave + (nzpoint_Vol-1)*kskip_DNS

       print *, '######################################################'
       print *, 'DNS Output timeseries spatial index range'
       print *, 'ibe = ', ibe_DNS, 'iend = ', iend_DNS
       print *, 'jbe = ', jbe_DNS, 'jend = ', jend_DNS
       print *, 'kbe = ', kbe_DNS, 'kend = ', kend_DNS
       print *, 'Actual Spatial Average range'
       print *, 'ibe_ave = ', ibe_ave, 'iend_ave = ', iend_ave
       print *, 'jbe_ave = ', jbe_ave, 'jend_ave = ', jend_ave
       print *, 'kbe_ave = ', kbe_ave, 'kend_ave = ', kend_ave
       print *, 'Number of spatial points for Averaging'
       print *, ' nxpoint_Vol = ', nxpoint_Vol, 'nypoint_Vol = ', nypoint_Vol, 'nzpoint_Vol = ', nzpoint_Vol


       ! Dimension range for tsdata
       imints_Vol = 1 - min(iwinl, (ibe_ave-ibe_DNS)/iskip_DNS)
       imaxts_Vol = nxpoint_Vol + min(iwinr, (iend_DNS-iend_ave)/iskip_DNS)
       iavelen = imaxts_Vol - imints_Vol + 1
       print *, 'streamwise Index range for TimeSeries Volume buffer'
       print *, ' imints_Vol = ', imints_Vol, ' imaxts_Vol = ', imaxts_Vol, 'iavelen (buffer length) = ', iavelen

       iminloc = max(ibe_ave  - iwinl*iskip_DNS,ibe_DNS)
       imaxloc = min(iend_ave + iwinr*iskip_DNS,iend_DNS)
       print *, 'Physical i-index range read TS (including xwindow)'
       print *, '   iminloc=', iminloc, 'imaxloc=', imaxloc

       kmints_Vol = 1 - min(kwind, (kbe_ave-kbe_DNS)/kskip_DNS)
       kmaxts_Vol = nzpoint_Vol + min(kwinu, (kend_DNS-kend_ave)/kskip_DNS)
       kavelen = kmaxts_Vol - kmints_Vol + 1
       print *, 'wall-normal Index range for TimeSeries Volume buffer'
       print *, ' kmints_Vol = ', kmints_Vol, ' kmaxts_Vol = ', kmaxts_Vol, 'kavelen (buffer length) = ', kavelen

       kminloc = max(kbe_ave  - kwind*kskip_DNS,kbe_DNS)
       kmaxloc = min(kend_ave + kwinu*kskip_DNS,kend_DNS)
       print *, 'Physical k-index range read TS (including kwindow)'
       print *, '   kminloc=', kminloc, 'kmaxloc=', kmaxloc

       TSVol%gname = '/vol'
       TSVol%rank = 3

       allocate(TSVol%dname(nvar_output))
       allocate(TSVol%dimsf(TSVol%rank),TSVol%dimsm(TSVol%rank))
       allocate(TSVol%count(TSVol%rank),TSVol%offset(TSVol%rank))
       allocate(TSVol%block(TSVol%rank),TSVol%stride(TSVol%rank))

       TSVol%dimsf(1) = kmaxts_Vol - kmints_Vol + 1 !kend_ave - kbe_ave + 1
       TSVol%dimsf(2) = imaxts_Vol - imints_Vol + 1
       TSVol%dimsf(3) = jend_ave - jbe_ave + 1
       TSVol%offset(1) = kbe_ave - kbe_DNS - kwind !kbe_ave - kbe_DNS
       TSVol%offset(2) = ibe_ave - ibe_DNS - iwinl
       TSVol%offset(3) = jbe_ave - jbe_DNS
       TSVol%dimsm = TSVol%dimsf
       TSVol%block = TSVol%dimsm
       TSVol%count = 1
       TSVol%stride = 1

       allocate(varout_Vol(nvar_output))
       do n=1, nvar_output
         TSVol%dname(n) = varname_Vol(varindex_output(n))
       enddo
       TSVol%dnum = nvar_output

       varout_Vol = TSVol%dname
       print *, 'Variable name read in:  ', (trim(TSVol%dname(n))//',', n=1, nvar_output) !!!!!


!print *, 'TSVol%dimsf = ', TSVol%dimsf
!print *, 'TSVol%dimsm = ', TSVol%dimsm
!print *, 'TSVol%offset = ', TSVol%offset
!print *, 'gname = ', TSVol%gname


       TSVol%anum = 3
       allocate(TSVol%attr_name(TSVol%anum))
       TSVol%attr_name(1) = "istart, iend, iskip"
       TSVol%attr_name(2) = "jstart, jend, jskip"
       TSVol%attr_name(3) = "kstart, kend, kskip"

       TSVol%IsHSInitialized = .true.

     end subroutine InitTSVol





end module MTSVOL


module MTSiplane
     use MTSHDF5

     ! private variables
!     character(10), dimension(10), private :: varname_iplane
     type(tp_hyperslab), dimension(:,:), allocatable, private :: TSiplane
     type(tp_hyperslab), dimension(:), allocatable, private :: TSigrid
     integer, private :: nfile ! Total number of files to read data from
     character(10), dimension(:), allocatable, private :: fname_jrange
     integer, private :: ntpoint_total, npath, nvar_output, jminloc, jmaxloc, javelen
     character(200), private :: filepath
     type(tp_DNSIndex), private :: DNSIndex_i
     type(fprop), dimension(:), allocatable, private :: fileprop

     ! public variables
     character(10), dimension(10) :: varname_iplane
     character(10), dimension(:), allocatable :: varout_i
     integer :: nypoint_i, nzpoint_i, jmints_i, jmaxts_i
contains

    subroutine ReadDNS_index_iplane(nt, dy, dz, dt, DNSIndex, filepath)
       integer, intent(out) :: nt
       real(8), intent(out) :: dy, dz, dt
       type(tp_DNSIndex), intent(out) :: DNSIndex
       character(*), intent(in) :: filepath
       type(tp_hyperslab) :: Index_iplane
       real(8), dimension(:,:), allocatable :: buffer_real
       integer, dimension(:,:), allocatable :: buffer_integer
       integer, dimension(:,:), allocatable :: buffer_ilocs
       real(8) :: time1, time2

       open(11,file=trim(filepath)//'series_time_ascii.dat',status='old')
         read(11,*) time1
         read(11,*) time2
       close(11)
       dt = time2 - time1

       Index_iplane%gname = '/iplane'
       Index_iplane%rank  = 1

       !Index_iplane%dnum  = 1
       Index_iplane%dnum  = 2
       allocate(Index_iplane%dname(Index_iplane%dnum),Index_iplane%dimsf(Index_iplane%rank))
       allocate(Index_iplane%dimsm(Index_iplane%rank),Index_iplane%count(Index_iplane%rank))
       allocate(Index_iplane%offset(Index_iplane%rank),Index_iplane%block(Index_iplane%rank))
       allocate(Index_iplane%stride(Index_iplane%rank))
       Index_iplane%dname(1)  = 'dy'
       Index_iplane%dname(2)  = 'dz'
       Index_iplane%dimsf = (/1/)
       Index_iplane%dimsm = Index_iplane%dimsf
       Index_iplane%block = Index_iplane%dimsm
       Index_iplane%offset = 0
       Index_iplane%count  = 1
       Index_iplane%stride = 1
       Index_iplane%IsHSInitialized = .true.
       allocate(buffer_real(1,2))
       !allocate(buffer_real(1,1))

       Index_iplane%fname = trim(filepath)//'DNS_index.h5'

       call ReadTSHDF5_1D(Index_iplane,buffer_real)

       dy = buffer_real(1,1)
       dz = buffer_real(1,2)
       !dz = 0.

       deallocate(Index_iplane%dname)
       Index_iplane%dnum  = 8
       allocate(Index_iplane%dname(Index_iplane%dnum))

       Index_iplane%dname(1)  = 'jstart'
       Index_iplane%dname(2)  = 'jend'
       Index_iplane%dname(3)  = 'jskip'
       Index_iplane%dname(4)  = 'kstart'
       Index_iplane%dname(5)  = 'kend'
       Index_iplane%dname(6)  = 'kskip'
       Index_iplane%dname(7)  = 'niplane'
       Index_iplane%dname(8)  = 'ntpoint'

       allocate(buffer_integer(1,8))

       call ReadTSHDF5_1D_integer(Index_iplane,buffer_integer)

       DNSIndex%jbe     = buffer_integer(1,1)
       DNSIndex%jend    = buffer_integer(1,2)
       DNSIndex%jskip   = buffer_integer(1,3)
       DNSIndex%kbe     = buffer_integer(1,4)
       DNSIndex%kend    = buffer_integer(1,5)
       DNSIndex%kskip   = buffer_integer(1,6)
       DNSIndex%niplane = buffer_integer(1,7)
       nt = buffer_integer(1,8)

       deallocate(Index_iplane%dname)
       Index_iplane%dnum  = 1
       allocate(Index_iplane%dname(Index_iplane%dnum))
       Index_iplane%dname(1) = 'ilocs'
       Index_iplane%dimsf = (/DNSIndex%niplane/)
       Index_iplane%dimsm = Index_iplane%dimsf
       Index_iplane%block = Index_iplane%dimsm
       allocate(buffer_ilocs(DNSIndex%niplane,1))

       call ReadTSHDF5_1D_integer(Index_iplane,buffer_ilocs)

       allocate(DNSIndex%iplane(DNSIndex%niplane))
       DNSIndex%iplane = buffer_ilocs(:,1)

     end subroutine ReadDNS_index_iplane

     subroutine ReadVOLTSIplane_files(i, iloc, buffer, buffer_grid)
        implicit none
        integer, intent(in) :: i, iloc
        real(8), intent(out) :: buffer(ntpoint_total, nypoint_i, nzpoint_i, 1, nvar_output) ! ntpoint
        real(8), intent(out), optional :: buffer_grid(nypoint_i, nzpoint_i,  1, 12)
        integer :: n, nnn, ktmp(1)
        character(8) :: fnum1
        integer :: nn
        integer :: ntpoint_tmp, num_file

         ntpoint_tmp = 0
         do nn=1, npath
           num_file = (fileprop(nn)%file_end - fileprop(nn)%file_be)/fileprop(nn)%file_skip + 1

           TSiplane(1,nn)%offset(4) = i - 1
           do n=1, num_file
             write(unit=fnum1,fmt='(I08.8)') fileprop(nn)%file_be + (n-1)*fileprop(nn)%file_skip
             TSiplane(1,nn)%fname = trim(fileprop(nn)%filepath)//'timeseries_'//fnum1//'.h5'

!             print *, 'Reading file: ', trim(TSkplane(1,nn)%fname)

             call ReadTSHDF5_4D(TSIplane(1,nn), buffer((ntpoint_tmp+1):(ntpoint_tmp+TSIplane(1,nn)%dimsf(1)), 1:TSIplane(1,nn)%dimsf(2), &
                                1:TSIplane(1,nn)%dimsf(3), 1:TSIplane(1,nn)%dimsf(4), 1:TSIplane(1,nn)%dnum) )

           ntpoint_tmp = ntpoint_tmp + TSIplane(1,nn)%dimsf(1)

           enddo ! end n loop
         enddo ! end nn loop

     end subroutine ReadVOLTSIplane_files


     subroutine ReadVOLTSiplane(fn, iloc, buffer, buffer_grid)
       implicit none
       character(*), intent(in) :: fn
       integer, intent(in) :: iloc
       real(8), intent(out) :: buffer(ntpoint_total, nypoint_i, nzpoint_i, 1, nvar_output) ! ntpoint
       real(8), intent(out), optional :: buffer_grid(nypoint_i, nzpoint_i, 1, 12)

       TSiplane(1,1)%fname = fn
       TSiplane(1,1)%offset(4) = iloc - 1
       TSigrid(1)%offset(3) = iloc - 1

       call ReadTSHDF5_4D(TSiplane(1,1), buffer)

       if(present(buffer_grid)) then
         if(.not.allocated(TSigrid)) allocate(TSigrid(1))
         TSigrid(1)%fname = trim(fileprop(1)%filepath)//'timeseries_GridMetrics.h5'
         call ReadTSHDF5_3D(TSigrid(1), buffer_grid)  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       endif

     end subroutine ReadVOLTSiplane


     subroutine ReadTSiplane(iloc, buffer, buffer_grid)
       implicit none
       integer, intent(in) :: iloc
       real(8), intent(out) :: buffer(ntpoint_total, nypoint_i, nzpoint_i,1,nvar_output)
       real(8), intent(out), optional :: buffer_grid(nypoint_i, nzpoint_i,1,12)
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
       do nn=1, npath
         nnn = jmints_i
         do n=1, nfile
           write(unit=fnum1,fmt='(I04.4)') iloc
           TSiplane(n,nn)%fname = trim(fileprop(nn)%filepath)//'timeseries_iplane'//fnum1//'_'//fname_jrange(n)//'.h5'
           print *, 'Reading file: ', trim(TSiplane(n,nn)%fname)
           call ReadTSHDF5_4D(TSiplane(n,nn), buffer(ntpoint_tmp+1:ntpoint_tmp+TSiplane(n,nn)%dimsf(1), nnn:(nnn+TSiplane(n,nn)%dimsf(2)-1), &
                              1:TSiplane(n,nn)%dimsf(3),1:TSiplane(n,nn)%dimsf(4),1:TSiplane(n,nn)%dnum) )
           nnn = nnn + TSiplane(n,nn)%dimsf(2)
         enddo
         ntpoint_tmp = ntpoint_tmp + TSiplane(1,nn)%dimsf(1)
       enddo

       if(present(buffer_grid)) then
         TSigrid(1)%fname = trim(fileprop(1)%filepath)//'timeseries_GridMetrics.h5'
         call ReadTSHDF5_3D(TSigrid(1), buffer_grid)
       endif

     end subroutine ReadTSiplane


     subroutine InitTSiplane(nt, np, DNSIndex, jbe_ave_inp, jend_ave_inp, kbe_ave_inp, kend_ave_inp, &
                             nvarout, varindex_output, fp, ireadgrid )
       implicit none
       integer, intent(in) :: nt, np
       type(tp_DNSIndex), intent(in) :: DNSIndex
       integer, intent(in) :: jbe_ave_inp, jend_ave_inp, kbe_ave_inp, kend_ave_inp
       integer, intent(in) :: nvarout, varindex_output(nvarout)
       type(fprop), intent(in) :: fp(:)
       integer, intent(in), optional :: ireadgrid !!!

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
       ntpoint_total = nt
       nvar_output = nvarout
       npath = np
       allocate(fileprop(npath))
       fileprop = fp

       varname_iplane(1) = 'u'
       varname_iplane(2) = 'v'
       varname_iplane(3) = 'w'
       varname_iplane(4) = 'p'
       varname_iplane(5) = 'T'
       varname_iplane(6) = 'ui'
       varname_iplane(7) = 'vi'
       varname_iplane(8) = 'wi'
       varname_iplane(9) = 'pi'
       varname_iplane(10) = 'Ti'

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
       nypoint_i = (jend_ave_inp - jbe_ave)/jskip_DNS + 1
       jend_ave = jbe_ave + (nypoint_i-1)*jskip_DNS

       if(mod(kbe_ave_inp-kbe_DNS,kskip_DNS).ne.0) then
         kbe_ave = kbe_ave_inp + kskip_DNS - mod(kbe_ave_inp-kbe_DNS,kskip_DNS)
       else
         kbe_ave = kbe_ave_inp
       endif
       nzpoint_i = (kend_ave_inp-kbe_ave)/kskip_DNS + 1
       kend_ave = kbe_ave + (nzpoint_i-1)*kskip_DNS

       print *, 'DNS Output timeseries spatial index range'
       print *, 'jbe = ', jbe_DNS, 'jend = ', jend_DNS
       print *, 'kbe = ', kbe_DNS, 'kend = ', kend_DNS
       print *, 'Actual Spatial Average range'
       print *, 'jbe_ave = ', jbe_ave, 'jend_ave = ', jend_ave
       print *, 'kbe_ave = ', kbe_ave, 'kend_ave = ', kend_ave
       print *, 'Number of spatial points for Averaging'
       print *, 'nypoint_i = ', nypoint_i, 'nzpoint_i = ', nzpoint_i

       ! Dimension range for tsdata
       jmints_i = 1 - min(0,(jbe_ave-jbe_DNS)/jskip_DNS)
       jmaxts_i = nypoint_i + min(0,(jend_DNS-jend_ave)/jskip_DNS)
       javelen = jmaxts_i - jmints_i + 1
       print *, 'spanwise Index range for iplane buffer'
       print *, 'jmints_i = ', jmints_i, 'jmaxts_i = ', jmaxts_i, 'javelen = ', javelen

       jminloc = max(jbe_ave,jbe_DNS)
       jmaxloc = min(jend_ave,jend_DNS)
       print *, 'Physical j-index range read TS'
       print *, 'jminloc = ', jminloc, 'jmaxloc = ', jmaxloc

       ! Total number of files that cover the full range of DNS output in j dir
       n_total = ceiling( ( (jend_DNS - jbe_DNS)/jskip_DNS + 1 )/real(jbuffer) )

       ! Determine the range of files that need to be read based on jbe_ave & jend_ave
       n_lower = 0; n_upper = 0 ! File indexes that are partially read
       joffset_lower = 0; jnum_lower = 0; jnum_upper = 0
       do n=1, n_total
         jbe_pfile = jbe_DNS + ((n-1)*jbuffer)*jskip_DNS
         jend_pfile = jbe_pfile + (jbuffer-1)*jskip_DNS
         if(jbe_pfile.le.jminloc.and.jend_pfile.ge.jminloc) then
           n_lower = n
           joffset_lower = (jminloc - jbe_pfile)/jskip_DNS
           jnum_lower = min((jend_pfile - jminloc)/jskip_DNS+1, javelen )
         endif
         if(jbe_pfile.le.jmaxloc.and.jend_pfile.ge.jmaxloc) then
           n_upper = n
           if(n_upper.gt.n_lower) jnum_upper = (jmaxloc - jbe_pfile)/jskip_DNS + 1
           exit
         endif
       enddo ! end n loop

       nfile = n_upper - n_lower + 1 ! Total number of files that need to be read

       allocate(fname_jrange(nfile))

       nn = 0
       do n = n_lower, n_upper
         nn = nn + 1
         jbe_pfile = jbe_DNS + ((n-1)*jbuffer)*jskip_DNS
         jend_pfile = min(jbe_pfile + (jbuffer-1)*jskip_DNS, jend_DNS)
         write(unit=fnum_jbe,fmt='(I04.4)') jbe_pfile
         write(unit=fnum_jend,fmt='(I04.4)') jend_pfile
         fname_jrange(nn) = 'j'//fnum_jbe//'-'//fnum_jend
         print *, 'fname_jrange = ', fname_jrange(nn)
       enddo ! end n loop

       allocate(TSiplane(nfile,npath))
       TSiplane%gname = '/iplane'
       TSiplane%rank = 4

       do nnn=1, npath
         do nn=1, nfile
           allocate(TSiplane(nn,nnn)%dname(nvar_output))
           allocate(TSiplane(nn,nnn)%dimsf(TSiplane(nn,nnn)%rank))
           allocate(TSiplane(nn,nnn)%dimsm(TSiplane(nn,nnn)%rank))
           allocate(TSiplane(nn,nnn)%count(TSiplane(nn,nnn)%rank))
           allocate(TSiplane(nn,nnn)%offset(TSiplane(nn,nnn)%rank))
           allocate(TSiplane(nn,nnn)%block(TSiplane(nn,nnn)%rank))
           allocate(TSiplane(nn,nnn)%stride(TSiplane(nn,nnn)%rank))
         enddo
       enddo

       ! Files that are fully read
       do nnn=1, npath
         do nn=1, nfile
           TSiplane(nn,nnn)%dimsf(1) = fileprop(nnn)%ntpoint
           TSiplane(nn,nnn)%dimsf(2) = jbuffer
           TSiplane(nn,nnn)%dimsf(3) = (kend_ave - kbe_ave)/kskip_DNS + 1
           TSiplane(nn,nnn)%dimsf(4) = 1
           TSiplane(nn,nnn)%offset(1) = 0
           TSiplane(nn,nnn)%offset(2) = 0
           TSiplane(nn,nnn)%offset(3) = kbe_ave/kskip_DNS - 1       !!!!!!!!!!
           TSiplane(nn,nnn)%offset(4) = 0
         enddo
       enddo

       ! Files that are partially read in j-dir
       do nnn=1, npath
         TSiplane(1,nnn)%dimsf(2) = jnum_lower
         TSiplane(1,nnn)%offset(2) = joffset_lower
       enddo
       if(nfile.gt.1) then
         do nnn=1, npath
           TSiplane(nfile,nnn)%dimsf(2) = jnum_upper
         enddo
       endif

       do nnn=1, npath
         do nn=1, nfile
           TSiplane(nn,nnn)%dimsm = TSiplane(nn,nnn)%dimsf
           TSiplane(nn,nnn)%block = TSiplane(nn,nnn)%dimsf
           TSiplane(nn,nnn)%count = 1
           TSiplane(nn,nnn)%stride = 1
         enddo
       enddo

       allocate(varout_i(nvar_output))
       do nnn=1, npath
         do nn=1, nfile
           do n=1, nvar_output
             TSiplane(nn,nnn)%dname(n) = varname_iplane(varindex_output(n))
           enddo
           TSiplane(nn,nnn)%dnum = nvar_output
         enddo
       enddo
       varout_i = TSiplane(1,1)%dname
       print *, 'Variable name read in: ', (trim(TSiplane(1,1)%dname(n))//',', n=1, nvar_output)
       TSiplane%IsHSInitialized = .true.

       if(present(ireadgrid)) then
         if(ireadgrid.eq.1) then
           allocate(TSigrid(1))
           TSigrid%gname = '/iplane'
           TSigrid%rank = 3

           allocate( TSigrid(1)%dname(12) )
           allocate( TSigrid(1)%dimsf(TSigrid(1)%rank) )
           allocate( TSigrid(1)%dimsm(TSigrid(1)%rank) )
           allocate( TSigrid(1)%count(TSigrid(1)%rank) )
           allocate( TSigrid(1)%offset(TSigrid(1)%rank) )
           allocate( TSigrid(1)%block(TSigrid(1)%rank) )
           allocate( TSigrid(1)%stride(TSigrid(1)%rank) )
           TSigrid(1)%dnum = 12
           TSigrid(1)%dname(1) = "x"
           TSigrid(1)%dname(2) = "y"
           TSigrid(1)%dname(3) = "z"
           TSigrid(1)%dname(4) = "didx"
           TSigrid(1)%dname(5) = "djdx"
           TSigrid(1)%dname(6) = "dkdx"
           TSigrid(1)%dname(7) = "didy"
           TSigrid(1)%dname(8) = "djdy"
           TSigrid(1)%dname(9) = "dkdy"
           TSigrid(1)%dname(10) = "didz"
           TSigrid(1)%dname(11) = "djdz"
           TSigrid(1)%dname(12) = "dkdz"

           TSigrid(1)%dimsf(1) = jend_ave - jbe_ave + 1
           TSigrid(1)%dimsf(2) = kend_ave - kbe_ave + 1
           TSigrid(1)%dimsf(3) = 1                       !! should be changed
           TSigrid(1)%offset(1) = jbe_ave - jbe_DNS
           TSigrid(1)%offset(2) = kbe_ave - kbe_DNS
           TSigrid(1)%offset(3) = 0                      !! should be changed

           TSigrid(1)%dimsm(1) = TSigrid(1)%dimsf(1)
           TSigrid(1)%dimsm(2) = TSigrid(1)%dimsf(2)
           TSigrid(1)%dimsm(3) = TSigrid(1)%dimsf(3)
           TSigrid(1)%block(1) = TSigrid(1)%dimsf(1)
           TSigrid(1)%block(2) = TSigrid(1)%dimsf(2)
           TSigrid(1)%block(3) = TSigrid(1)%dimsf(3)
           TSigrid(1)%count = 1
           TSigrid(1)%stride = 1
           TSigrid(1)%IsHSInitialized = .true.

         endif
       endif


     end subroutine InitTSiplane


     subroutine InitVOLTSiplane(nt, np, DNSIndex, jbe_ave_inp, jend_ave_inp, kbe_ave_inp, kend_ave_inp, &
                  iwinl, iwinr, nvarout, varindex_output, fp, ireadgrid)
        implicit none

        integer, intent(in) :: nt, np
        type(tp_DNSIndex), intent(in) :: DNSIndex
        integer, intent(in) :: jbe_ave_inp, jend_ave_inp, kbe_ave_inp, kend_ave_inp, iwinl, iwinr
        integer, intent(in) :: nvarout, varindex_output(nvarout)
        type(fprop), intent(in) :: fp(:)
        integer, intent(in), optional :: ireadgrid

        integer :: jbe_DNS, jend_DNS, jskip_DNS
        integer :: kbe_DNS, kend_DNS, kskip_DNS
        integer :: jbe_ave, jend_ave, kbe_ave, kend_ave

        integer :: n, nn, nnn
!        character(4) :: fnum_ibe, fnum_iend
        integer :: njloc

        DNSIndex_i = DNSIndex
        ntpoint_total = nt
        nvar_output = nvarout
        npath = np
        allocate(fileprop(npath))
        fileprop = fp

        varname_iplane(1) = 'u'
        varname_iplane(2) = 'v'
        varname_iplane(3) = 'w'
        varname_iplane(4) = 'p'
        varname_iplane(5) = 'T'
        varname_iplane(6) = 'ui'
        varname_iplane(7) = 'vi'
        varname_iplane(8) = 'wi'
        varname_iplane(9) = 'pi'
        varname_iplane(10) = 'Ti'

        jbe_DNS   = DNSIndex_i%jbe
        jend_DNS  = DNSIndex_i%jend
        jskip_DNS = DNSIndex_i%jskip
        kbe_DNS   = DNSIndex_i%kbe
        kend_DNS  = DNSIndex_i%kend
        kskip_DNS = DNSIndex_i%kskip

        if(jbe_ave_inp.gt.jend_ave_inp.or.kbe_ave_inp.gt.kend_ave_inp.or. &
           jbe_ave_inp.lt.jbe_DNS.or.jend_ave_inp.gt.jend_DNS.or.kbe_ave_inp.lt.kbe_DNS.or.kend_ave_inp.gt.kend_DNS) then
            print *, 'Ranges for averaging is out of bound ... STOP'
            print *, 'jbe_ave_inp =', jbe_ave_inp, 'jend_ave_inp = ', jend_ave_inp
            print *, 'kbe_ave_inp =', kbe_ave_inp, 'kend_ave_inp = ', kend_ave_inp
            if(jbe_ave_inp.lt.jbe_DNS) print *, 'jbe_ave_inp < jbe_DNS'
            if(jend_ave_inp.gt.jend_DNS) print *, 'jend_ave_inp > jend_DNS'
            if(kbe_ave_inp.lt.kbe_DNS) print *, 'kbe_ave_inp < kbe_DNS'
            if(kend_ave_inp.gt.kend_DNS) print *, 'kend_ave_inp > kend_DNS'
            stop
        endif

        if(jbe_DNS.gt.jend_DNS.or.kbe_DNS.gt.kend_DNS.or.jskip_DNS.le.0.or.kskip_DNS.le.0) then
            print *, 'Input for DNS timeseries files out of bound ... STOP'
            print *, 'jbe_DNS =', jbe_DNS, 'jend_DNS = ', jend_DNS, 'jskip_DNS =', jskip_DNS
            print *, 'kbe_DNS =', kbe_DNS, 'kend_DNS = ', kend_DNS, 'kskip_DNS =', kskip_DNS
            stop
        endif

        ! begin, end, and number of spatial points for Averaging

        if(mod(jbe_ave_inp-jbe_DNS,jskip_DNS).ne.0) then
          jbe_ave = jbe_ave_inp + jskip_DNS  - mod(jbe_ave_inp-jbe_DNS,jskip_DNS)
        else
          jbe_ave = jbe_ave_inp
        endif
        nypoint_i = (jend_ave_inp-jbe_ave)/jskip_DNS +1
        jend_ave = jbe_ave + (nypoint_i-1)*jskip_DNS


        if(mod(kbe_ave_inp-kbe_DNS,kskip_DNS).ne.0) then
          kbe_ave = kbe_ave_inp + kskip_DNS - mod(kbe_ave_inp-kbe_DNS,kskip_DNS)
        else
          kbe_ave = kbe_ave_inp
        endif
        nzpoint_i = (kend_ave_inp-kbe_ave)/kskip_DNS + 1
        kend_ave = kbe_ave + (nzpoint_i-1)*kskip_DNS

        print *, 'DNS Output timeseries spatial index range'
        print *, 'jbe = ', jbe_DNS, 'jend = ', jend_DNS
        print *, 'kbe = ', kbe_DNS, 'kend = ', kend_DNS
        print *, 'Actual Spatial Average range'
        print *, 'jbe_ave = ', jbe_ave, 'jend_ave = ', jend_ave
        print *, 'kbe_ave = ', kbe_ave, 'kend_ave = ', kend_ave
        print *, 'Number of spatial points for Averaging'
        print *, ' nypoint_i = ', nypoint_i, 'nzpoint_i = ', nzpoint_i

        ! Dimension range for tsdata

!        imints_j = 1 - min(iwinl, (ibe_ave-ibe_DNS)/iskip_DNS)
!        imaxts_j = nxpoint_j + min(iwinr, (iend_DNS-iend_ave)/iskip_DNS)
!        iavelen = imaxts_j - imints_j + 1
!        print *, 'streamwise Index range for jplane buffer'
!        print *, ' imints_j = ', imints_j, ' imaxts_j = ', imaxts_j, 'iavelen (buffer length) = ', iavelen

!        iminloc = max(ibe_ave  - iwinl*iskip_DNS,ibe_DNS)
!        imaxloc = min(iend_ave + iwinr*iskip_DNS,iend_DNS)
!        print *, 'Physical i-index range read TS (including xwindow)'
!        print *, '   iminloc=', iminloc, 'imaxloc=', imaxloc


       allocate(TSiplane(1,npath))
        TSiplane%gname = '/iplane'
        TSiplane%rank = 4

        do nnn=1, npath
          allocate(TSiplane(1,nnn)%dname(nvar_output))
          allocate(TSiplane(1,nnn)%dimsf(TSiplane(1,nnn)%rank))
          allocate(TSiplane(1,nnn)%dimsm(TSiplane(1,nnn)%rank))
          allocate(TSiplane(1,nnn)%count(TSiplane(1,nnn)%rank))
          allocate(TSiplane(1,nnn)%offset(TSiplane(1,nnn)%rank))
          allocate(TSiplane(1,nnn)%block(TSiplane(1,nnn)%rank))
          allocate(TSiplane(1,nnn)%stride(TSiplane(1,nnn)%rank))
       enddo

       do nnn=1, npath
            TSiplane(1,nnn)%dimsf(1)  = fileprop(nnn)%ntpoint     ! ntpoint
            TSiplane(1,nnn)%dimsf(2)  = nypoint_i ! imaxts_j - imints_j + 1
!            TSkplane(1,nnn)%dimsf(3)  = iend_ave - ibe_ave + 1
            TSiplane(1,nnn)%dimsf(3)  = nzpoint_i                 ! check
            TSiplane(1,nnn)%dimsf(4)  = 1                         !!!!!!!!!    should be changed
            TSiplane(1,nnn)%offset(1) = 0
            TSiplane(1,nnn)%offset(2) = jbe_ave - jbe_DNS         !!!!
            TSiplane(1,nnn)%offset(3) = kbe_ave - kbe_DNS         !!!!
            print *, 'offset(2) = ', jbe_ave - jbe_DNS
            TSiplane(1,nnn)%offset(4) = 0                         !!!!!!!!!!!! should be changed
            TSiplane(1,nnn)%dimsm(1)  = TSiplane(1,nnn)%dimsf(1)
            TSiplane(1,nnn)%dimsm(2)  = TSiplane(1,nnn)%dimsf(2)
            TSiplane(1,nnn)%dimsm(3)  = TSiplane(1,nnn)%dimsf(3)
            TSiplane(1,nnn)%dimsm(4)  = TSiplane(1,nnn)%dimsf(4)
            TSiplane(1,nnn)%block(1)  = TSiplane(1,nnn)%dimsf(1)
            TSiplane(1,nnn)%block(2)  = TSiplane(1,nnn)%dimsf(2)
            TSiplane(1,nnn)%block(3)  = TSiplane(1,nnn)%dimsf(3)
            TSiplane(1,nnn)%block(4)  = TSiplane(1,nnn)%dimsf(4)
            TSiplane(1,nnn)%count  = 1
            TSiplane(1,nnn)%stride = 1

       enddo ! end nnn loop


       allocate(varout_i(nvar_output))
       do nnn=1, npath
         do n=1, nvar_output
           TSiplane(1,nnn)%dname(n) = varname_iplane(varindex_output(n))
         enddo
         TSiplane(1,nnn)%dnum       = nvar_output
       enddo

       varout_i = TSiplane(1,1)%dname
       print *, 'Variable name read in:  ', (trim(TSiplane(1,1)%dname(n))//',', n=1, nvar_output) !!!!!
       TSiplane%IsHSInitialized = .true.

       if(present(ireadgrid)) then
          if(ireadgrid.eq.1) then
             allocate(TSigrid(1))
             TSigrid%gname = "/iplane"
             TSigrid%rank = 3

             allocate( TSigrid(1)%dname(12))
             allocate( TSigrid(1)%dimsf(TSigrid(1)%rank))
             allocate( TSigrid(1)%dimsm(TSigrid(1)%rank))
             allocate( TSigrid(1)%count(TSigrid(1)%rank))
             allocate( TSigrid(1)%offset(TSigrid(1)%rank))
             allocate( TSigrid(1)%block(TSigrid(1)%rank))
             allocate( TSigrid(1)%stride(TSigrid(1)%rank))
             TSigrid(1)%dnum = 12
             TSigrid(1)%dname(1) = "x"
             TSigrid(1)%dname(2) = "y"
             TSigrid(1)%dname(3) = "z"
             TSigrid(1)%dname(4) = "didx"
             TSigrid(1)%dname(5) = "djdx"
             TSigrid(1)%dname(6) = "dkdx"
             TSigrid(1)%dname(7) = "didy"
             TSigrid(1)%dname(8) = "djdy"
             TSigrid(1)%dname(9) = "dkdy"
             TSigrid(1)%dname(10) = "didz"
             TSigrid(1)%dname(11) = "djdz"
             TSigrid(1)%dname(12) = "dkdz"

             TSigrid(1)%dimsf(1) = nypoint_i ! jend_ave - jbe_ave + 1
             TSigrid(1)%dimsf(2) = nzpoint_i ! kend_ave - kbe_ave + 1
             TSigrid(1)%dimsf(3) = 1                       !! should be changed
             TSigrid(1)%offset(1) = jbe_ave - jbe_DNS
             TSigrid(1)%offset(2) = kbe_ave - kbe_DNS
             TSigrid(1)%offset(3) = 0                      !! should be changed

             TSigrid(1)%dimsm(1) = TSigrid(1)%dimsf(1)
             TSigrid(1)%dimsm(2) = TSigrid(1)%dimsf(2)
             TSigrid(1)%dimsm(3) = TSigrid(1)%dimsf(3)
             TSigrid(1)%block(1) = TSigrid(1)%dimsf(1)
             TSigrid(1)%block(2) = TSigrid(1)%dimsf(2)
             TSigrid(1)%block(3) = TSigrid(1)%dimsf(3)
             TSigrid(1)%count = 1
             TSigrid(1)%stride = 1
             TSigrid(1)%IsHSInitialized = .true.
           endif
         endif


     end subroutine InitVOLTSiplane

end module MTSiplane



module MTSjplane
     use MTSHDF5

     character(10), dimension(10), private :: varname_jplane
     type(tp_hyperslab), dimension(:,:), allocatable, private :: TSJplane
     type(tp_hyperslab), dimension(:), allocatable, private :: TSJgrid
     integer, private :: nfile ! Total number of files to read data from
     integer, private :: ntpoint_total, npath, nvar_output, iminloc, imaxloc, iavelen
     integer, private :: kminloc, kmaxloc, kavelen
!     character(200), private :: filepath
     type(tp_DNSIndex), private :: DNSIndex_j
     type(fprop), dimension(:), allocatable, private :: fileprop

     character(10), dimension(:), allocatable :: varout_j
     integer :: nxpoint_j, nzpoint_j, imints_j, imaxts_j, kmints_j, kmaxts_j
contains


    subroutine ReadDNS_index_jplane(nt, dx, dz, dt, DNSIndex, filepath)
       integer, intent(out) :: nt
       real(8), intent(out) :: dx, dz, dt
       type(tp_DNSIndex), intent(out) :: DNSIndex
       character(*), intent(in) :: filepath
       type(tp_hyperslab) :: Index_jplane
       real(8), dimension(:,:), allocatable :: buffer_real
       integer, dimension(:,:), allocatable :: buffer_integer
       integer, dimension(:,:), allocatable :: buffer_jlocs
       real(8) :: time1, time2

       dz = 0.d0

       open(11,file=trim(filepath)//'series_time_ascii.dat',status='old')
         read(11,*) time1
         read(11,*) time2
       close(11)
       dt = time2 - time1

       Index_jplane%gname = '/jplane'
       Index_jplane%rank  = 1

       !Index_jplane%dnum  = 2
       Index_jplane%dnum  = 1
       allocate(Index_jplane%dname(Index_jplane%dnum),Index_jplane%dimsf(Index_jplane%rank))
       allocate(Index_jplane%dimsm(Index_jplane%rank),Index_jplane%count(Index_jplane%rank))
       allocate(Index_jplane%offset(Index_jplane%rank),Index_jplane%block(Index_jplane%rank))
       allocate(Index_jplane%stride(Index_jplane%rank))
       Index_jplane%dname(1)  = 'dx'
       !Index_jplane%dname(2)  = 'dz'
       Index_jplane%dimsf = (/1/)
       Index_jplane%dimsm = Index_jplane%dimsf
       Index_jplane%block = Index_jplane%dimsm
       Index_jplane%offset = 0
       Index_jplane%count  = 1
       Index_jplane%stride = 1
       Index_jplane%IsHSInitialized = .true.
       !allocate(buffer_real(1,2))
       allocate(buffer_real(1,1))

!       Index_jplane%fname = 'DNS_index.h5'
       Index_jplane%fname = trim(filepath)//'DNS_index.h5'
!       Index_kplane%fname = trim(filepath)//'DNS_index.h5'

       call ReadTSHDF5_1D(Index_jplane,buffer_real)

       dx = buffer_real(1,1)
       !dz = buffer_real(1,2)

       deallocate(Index_jplane%dname)
       Index_jplane%dnum  = 8
       allocate(Index_jplane%dname(Index_jplane%dnum))

       Index_jplane%dname(1)  = 'istart'
       Index_jplane%dname(2)  = 'iend'
       Index_jplane%dname(3)  = 'iskip'
       Index_jplane%dname(4)  = 'kstart'
       Index_jplane%dname(5)  = 'kend'
       Index_jplane%dname(6)  = 'kskip'
       Index_jplane%dname(7)  = 'njplane'
       Index_jplane%dname(8)  = 'ntpoint'

       allocate(buffer_integer(1,8))

       call ReadTSHDF5_1D_integer(Index_jplane,buffer_integer)

       DNSIndex%ibe     = buffer_integer(1,1)
       DNSIndex%iend    = buffer_integer(1,2)
       DNSIndex%iskip   = buffer_integer(1,3)
       DNSIndex%kbe     = buffer_integer(1,4)
       DNSIndex%kend    = buffer_integer(1,5)
       DNSIndex%kskip   = buffer_integer(1,6)
       DNSIndex%njplane = buffer_integer(1,7)
       nt = buffer_integer(1,8)


!print *, '', DNSIndex%

       deallocate(Index_jplane%dname)
       Index_jplane%dnum  = 1
       allocate(Index_jplane%dname(Index_jplane%dnum))
       Index_jplane%dname(1) = 'jlocs'
       Index_jplane%dimsf = (/DNSIndex%njplane/)
       Index_jplane%dimsm = Index_jplane%dimsf
       Index_jplane%block = Index_jplane%dimsm
       allocate(buffer_jlocs(DNSIndex%njplane,1))

       call ReadTSHDF5_1D_integer(Index_jplane,buffer_jlocs)

       allocate(DNSIndex%jplane(DNSIndex%njplane))
       DNSIndex%jplane = buffer_jlocs(:,1)

     end subroutine ReadDNS_index_jplane

     subroutine ReadVOLTSJplane_files(j, jloc, buffer, buffer_grid)
        implicit none
        integer, intent(in) :: j, jloc
        real(8), intent(out) :: buffer(ntpoint_total, imints_j:imaxts_j, kmints_j:kmaxts_j, 1, nvar_output) ! ntpoint
        real(8), intent(out), optional :: buffer_grid(imints_j:imaxts_j, kmints_j:kmaxts_j, 1, 12)
        integer :: n, nnn, ktmp(1), i
        real(8) :: dxi, dyi, dzi
        character(8) :: fnum1
        integer :: nn
        integer :: ntpoint_tmp, num_file

         ntpoint_tmp = 0
         do nn=1, npath
           nnn = imints_j
           num_file = (fileprop(nn)%file_end - fileprop(nn)%file_be)/fileprop(nn)%file_skip + 1

           TSjplane(1,nn)%offset(4) = j - 1
           do n=1, num_file
             write(unit=fnum1,fmt='(I08.8)') fileprop(nn)%file_be + (n-1)*fileprop(nn)%file_skip
             TSjplane(1,nn)%fname = trim(fileprop(nn)%filepath)//'timeseries_'//fnum1//'.h5'

!             print *, 'Reading file: ', trim(TSkplane(1,nn)%fname)

             call ReadTSHDF5_4D(TSJplane(1,nn), buffer((ntpoint_tmp+1):(ntpoint_tmp+TSJplane(1,nn)%dimsf(1)), nnn:(nnn+TSjplane(1,nn)%dimsf(2)-1), &
                                kmints_j:kmints_j+TSjplane(1,nn)%dimsf(3)-1, 1:TSjplane(1,nn)%dimsf(4), 1:TSjplane(1,nn)%dnum) )
 ! should be changed
 !            if(present(buffer_grid).and.TSkgrid(1)%IsHSInitialized.and.n.eq.1.and.nn.eq.1) then
 !              TSkgrid(1)%fname = trim(fileprop(1)%filepath)//'timeseries_GridMetrics.h5'
 !              call ReadTSHDF5_3D(TSkgrid(1), buffer_grid(1:TSkgrid(1)%dimsf(1), &
 !                                 nnn:(nnn+TSkgrid(1)%dimsf(2)-1), 1:TSkgrid(1)%dimsf(3), 1:TSkgrid(1)%dnum) )
 !            endif
           ntpoint_tmp = ntpoint_tmp + TSJplane(1,nn)%dimsf(1)

           enddo ! end n loop
         enddo ! end nn loop

     end subroutine ReadVOLTSJplane_files


     subroutine ReadVOLTSJplane(fn, jloc, buffer, buffer_grid)
       implicit none
!       integer, intent(in) :: jloc
       character(*), intent(in) :: fn
       integer :: jloc
       real(8), intent(out) :: buffer(ntpoint_total, imints_j:imaxts_j, nzpoint_j, 1, nvar_output) ! ntpoint
       real(8), intent(out), optional :: buffer_grid(imints_j:imaxts_j, nzpoint_j, 1, 12)
!       integer :: n, nnn, jtmp(1), i
!       real(8) :: dxi, dyi, dzi !!!!!!!!!!
!       character(8) :: fnum
!       integer :: nn
!       integer :: num_file
       TSjplane(1,1)%fname = fn
       TSjplane(1,1)%offset(4) = jloc - 1
       TSjgrid(1)%offset(3) = jloc - 1

       call ReadTSHDF5_4D(TSjplane(1,1), buffer)

       if(present(buffer_grid)) then
         if(.not.allocated(TSJgrid)) allocate(TSJgrid(1))
         TSJgrid(1)%fname = trim(fileprop(1)%filepath)//'timeseries_GridMetrics.h5'
         call ReadTSHDF5_3D(TSJgrid(1), buffer_grid)  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       endif

     end subroutine ReadVOLTSJplane

     subroutine InitVOLTSJplane(nt, np, DNSIndex, ibe_ave_inp, iend_ave_inp, kbe_ave_inp, kend_ave_inp, &
                  iwinl, iwinr, kwind, kwinu, nvarout, varindex_output, fp, ireadgrid)
        implicit none

        integer, intent(in) :: nt, np
        type(tp_DNSIndex), intent(in) :: DNSIndex
        integer, intent(in) :: ibe_ave_inp, iend_ave_inp, kbe_ave_inp, kend_ave_inp, iwinl, iwinr, kwind, kwinu
        integer, intent(in) :: nvarout, varindex_output(nvarout)
        type(fprop), intent(in) :: fp(:)
        integer, intent(in), optional :: ireadgrid

        integer :: ibe_DNS, iend_DNS, iskip_DNS
        integer :: kbe_DNS, kend_DNS, kskip_DNS
        integer :: ibe_ave, iend_ave, kbe_ave, kend_ave

        integer :: n, nn, nnn
        character(4) :: fnum_ibe, fnum_iend
        integer :: njloc

        DNSIndex_j = DNSIndex
        ntpoint_total = nt
        nvar_output = nvarout
        npath = np
        allocate(fileprop(npath))
        fileprop = fp

        varname_jplane(1) = 'u'
        varname_jplane(2) = 'v'
        varname_jplane(3) = 'w'
        varname_jplane(4) = 'p'
        varname_jplane(5) = 'T'
        varname_jplane(6) = 'uj'
        varname_jplane(7) = 'vj'
        varname_jplane(8) = 'wj'
        varname_jplane(9) = 'pj'
        varname_jplane(10) = 'Tj'

        ibe_DNS   = DNSIndex_j%ibe
        iend_DNS  = DNSIndex_j%iend
        iskip_DNS = DNSIndex_j%iskip
        kbe_DNS   = DNSIndex_j%kbe
        kend_DNS  = DNSIndex_j%kend
        kskip_DNS = DNSIndex_j%kskip

        if(ibe_ave_inp.gt.iend_ave_inp.or.kbe_ave_inp.gt.kend_ave_inp.or. &
           ibe_ave_inp.lt.ibe_DNS.or.iend_ave_inp.gt.iend_DNS.or.kbe_ave_inp.lt.kbe_DNS.or.kend_ave_inp.gt.kend_DNS) then
            print *, 'Ranges for averaging is out of bound ... STOP'
            print *, 'ibe_ave_inp =', ibe_ave_inp, 'iend_ave_inp = ', iend_ave_inp
            print *, 'kbe_ave_inp =', kbe_ave_inp, 'kend_ave_inp = ', kend_ave_inp
            if(ibe_ave_inp.lt.ibe_DNS) print *, 'ibe_ave_inp < ibe_DNS'
            if(iend_ave_inp.gt.iend_DNS) print *, 'iend_ave_inp > iend_DNS'
            if(kbe_ave_inp.lt.kbe_DNS) print *, 'kbe_ave_inp < kbe_DNS'
            if(kend_ave_inp.gt.kend_DNS) print *, 'kend_ave_inp > kend_DNS'
            stop
        endif

        if(ibe_DNS.gt.iend_DNS.or.kbe_DNS.gt.kend_DNS.or.iskip_DNS.le.0.or.kskip_DNS.le.0) then
            print *, 'Input for DNS timeseries files out of bound ... STOP'
            print *, 'ibe_DNS =', ibe_DNS, 'iend_DNS = ', iend_DNS, 'iskip_DNS =', iskip_DNS
            print *, 'kbe_DNS =', kbe_DNS, 'kend_DNS = ', kend_DNS, 'kskip_DNS =', kskip_DNS
            stop
        endif


        ! begin, end, and number of spatial points for Averaging

        if(mod(ibe_ave_inp-ibe_DNS,iskip_DNS).ne.0) then
          ibe_ave = ibe_ave_inp + iskip_DNS  - mod(ibe_ave_inp-ibe_DNS,iskip_DNS)
        else
          ibe_ave = ibe_ave_inp
        endif
        nxpoint_j = (iend_ave_inp-ibe_ave)/iskip_DNS +1
        iend_ave = ibe_ave + (nxpoint_j-1)*iskip_DNS

        if(mod(kbe_ave_inp-kbe_DNS,kskip_DNS).ne.0) then
          kbe_ave = kbe_ave_inp + kskip_DNS - mod(kbe_ave_inp-kbe_DNS,kskip_DNS)
        else
          kbe_ave = kbe_ave_inp
        endif
        nzpoint_j = (kend_ave_inp-kbe_ave)/kskip_DNS + 1
        kend_ave = kbe_ave + (nzpoint_j-1)*kskip_DNS

        print *, 'DNS Output timeseries spatial index range'
        print *, 'ibe = ', ibe_DNS, 'iend = ', iend_DNS
        print *, 'kbe = ', kbe_DNS, 'kend = ', kend_DNS
        print *, 'Actual Spatial Average range'
        print *, 'ibe_ave = ', ibe_ave, 'iend_ave = ', iend_ave
        print *, 'kbe_ave = ', kbe_ave, 'kend_ave = ', kend_ave
        print *, 'Number of spatial points for Averaging'
        print *, ' nxpoint_j = ', nxpoint_j, 'nzpoint_j = ', nzpoint_j

        ! Dimension range for tsdata

        imints_j = 1 - min(iwinl, (ibe_ave-ibe_DNS)/iskip_DNS)
        imaxts_j = nxpoint_j + min(iwinr, (iend_DNS-iend_ave)/iskip_DNS)
        iavelen = imaxts_j - imints_j + 1
        print *, 'streamwise Index range for jplane buffer'
        print *, ' imints_j = ', imints_j, ' imaxts_j = ', imaxts_j, 'iavelen (buffer length) = ', iavelen

        iminloc = max(ibe_ave  - iwinl*iskip_DNS,ibe_DNS)
        imaxloc = min(iend_ave + iwinr*iskip_DNS,iend_DNS)
        print *, 'Physical i-index range read TS (including xwindow)'
        print *, '   iminloc=', iminloc, 'imaxloc=', imaxloc

        kmints_j = 1 -  min(kwind, (kbe_ave-kbe_DNS)/kskip_DNS)
        kmaxts_j = nzpoint_j + min(kwinu, (kend_DNS-kend_ave)/kskip_DNS)
        kavelen = kmaxts_j - kmints_j + 1
        print *, 'wall-normal Index range for jplane buffer'
        print *, ' kmints_j = ', kmints_j, ' kmaxts_j = ', kmaxts_j, 'kavelen (buffer length) = ', kavelen

       allocate(TSjplane(1,npath))
        TSjplane%gname = '/jplane'
        TSjplane%rank = 4

        do nnn=1, npath
          allocate(TSjplane(1,nnn)%dname(nvar_output))
          allocate(TSjplane(1,nnn)%dimsf(TSjplane(1,nnn)%rank))
          allocate(TSjplane(1,nnn)%dimsm(TSjplane(1,nnn)%rank))
          allocate(TSjplane(1,nnn)%count(TSjplane(1,nnn)%rank))
          allocate(TSjplane(1,nnn)%offset(TSjplane(1,nnn)%rank))
          allocate(TSjplane(1,nnn)%block(TSjplane(1,nnn)%rank))
          allocate(TSjplane(1,nnn)%stride(TSjplane(1,nnn)%rank))
       enddo

       do nnn=1, npath
            TSjplane(1,nnn)%dimsf(1)  = fileprop(nnn)%ntpoint     ! ntpoint
            TSjplane(1,nnn)%dimsf(2)  = imaxts_j - imints_j + 1
!            TSkplane(1,nnn)%dimsf(3)  = iend_ave - ibe_ave + 1
            TSjplane(1,nnn)%dimsf(3)  = kmaxts_j - kmints_j + 1   ! check
            TSjplane(1,nnn)%dimsf(4)  = 1                         !!!!!!!!!    should be changed
            TSjplane(1,nnn)%offset(1) = 0
            TSjplane(1,nnn)%offset(2) = ibe_ave - ibe_DNS - iwinl
            !TSjplane(1,nnn)%offset(3) = kbe_ave - kbe_DNS
            TSjplane(1,nnn)%offset(3) = kbe_ave - kbe_DNS - kwind
            print *, 'offset(2) (i) = ', ibe_ave - ibe_DNS - iwinl
            print *, 'offset(3) (k) = ', kbe_ave - kbe_DNS - kwind
            TSjplane(1,nnn)%offset(4) = 0                         !!!!!!!!!!!! should be changed
            if(TSjplane(1,nnn)%offset(2).lt.0) then
              print *, 'iwinl is bigger ... Stop! '
            endif
            if(TSjplane(1,nnn)%offset(4).lt.0) then
              print *, 'kwind is bigger ... Stop! '
            endif
            TSjplane(1,nnn)%dimsm(1)  = TSjplane(1,nnn)%dimsf(1)
            TSjplane(1,nnn)%dimsm(2)  = TSjplane(1,nnn)%dimsf(2)
            TSjplane(1,nnn)%dimsm(3)  = TSjplane(1,nnn)%dimsf(3)
            TSjplane(1,nnn)%dimsm(4)  = TSjplane(1,nnn)%dimsf(4)
            TSjplane(1,nnn)%block(1)  = TSjplane(1,nnn)%dimsf(1)
            TSjplane(1,nnn)%block(2)  = TSjplane(1,nnn)%dimsf(2)
            TSjplane(1,nnn)%block(3)  = TSjplane(1,nnn)%dimsf(3)
            TSjplane(1,nnn)%block(4)  = TSjplane(1,nnn)%dimsf(4)
            TSjplane(1,nnn)%count  = 1
            TSjplane(1,nnn)%stride = 1

       enddo ! end nnn loop


       allocate(varout_j(nvar_output))
       do nnn=1, npath
         do n=1, nvar_output
           TSJplane(1,nnn)%dname(n) = varname_jplane(varindex_output(n))
         enddo
         TSJplane(1,nnn)%dnum       = nvar_output
       enddo

       varout_j = TSJplane(1,1)%dname
       print *, 'Variable name read in:  ', (trim(TSJplane(1,1)%dname(n))//',', n=1, nvar_output) !!!!!
       TSJplane%IsHSInitialized = .true.

       if(present(ireadgrid)) then
          if(ireadgrid.eq.1) then
             allocate(TSJgrid(1))
             TSJgrid%gname = "/jplane"
             TSJgrid%rank = 3

             allocate( TSjgrid(1)%dname(12))
             allocate( TSjgrid(1)%dimsf(TSjgrid(1)%rank))
             allocate( TSjgrid(1)%dimsm(TSjgrid(1)%rank))
             allocate( TSjgrid(1)%count(TSjgrid(1)%rank))
             allocate( TSjgrid(1)%offset(TSjgrid(1)%rank))
             allocate( TSjgrid(1)%block(TSjgrid(1)%rank))
             allocate( TSjgrid(1)%stride(TSjgrid(1)%rank))
             TSjgrid(1)%dnum = 12
             TSjgrid(1)%dname(1) = "x"
             TSjgrid(1)%dname(2) = "y"
             TSjgrid(1)%dname(3) = "z"
             TSjgrid(1)%dname(4) = "didx"
             TSjgrid(1)%dname(5) = "djdx"
             TSjgrid(1)%dname(6) = "dkdx"
             TSjgrid(1)%dname(7) = "didy"
             TSjgrid(1)%dname(8) = "djdy"
             TSjgrid(1)%dname(9) = "dkdy"
             TSjgrid(1)%dname(10) = "didz"
             TSjgrid(1)%dname(11) = "djdz"
             TSjgrid(1)%dname(12) = "dkdz"

             TSjgrid(1)%dimsf(1) = imaxts_j - imints_j + 1
             !TSjgrid(1)%dimsf(2) = kend_ave - kbe_ave + 1
             TSjgrid(1)%dimsf(2) = kmaxts_j - kmints_j + 1
             TSjgrid(1)%dimsf(3) = 1                       !! should be changed
             TSjgrid(1)%offset(1) = ibe_ave - ibe_DNS - iwinl
             !TSjgrid(1)%offset(2) = kbe_ave - kbe_DNS
             TSjgrid(1)%offset(2) = kbe_ave - kbe_DNS - kwind
             TSjgrid(1)%offset(3) = 0                      !! should be changed

             TSjgrid(1)%dimsm(1) = TSjgrid(1)%dimsf(1)
             TSjgrid(1)%dimsm(2) = TSjgrid(1)%dimsf(2)
             TSjgrid(1)%dimsm(3) = TSjgrid(1)%dimsf(3)
             TSjgrid(1)%block(1) = TSjgrid(1)%dimsf(1)
             TSjgrid(1)%block(2) = TSjgrid(1)%dimsf(2)
             TSjgrid(1)%block(3) = TSjgrid(1)%dimsf(3)
             TSjgrid(1)%count = 1
             TSjgrid(1)%stride = 1
             TSjgrid(1)%IsHSInitialized = .true.
           endif
         endif


     end subroutine InitVOLTSjplane

end module MTSjplane


module MTSKplane
     use MTSHDF5
!     type tp_DNSIndex
!        integer :: ibe, iend, iskip
!        integer :: jbe, jend, jskip
!        integer :: ibuffer
!        integer :: nzloc
!        integer, dimension(:), allocatable :: kplane
!        real(8), dimension(:), allocatable :: zloc, dzdk
!        real(8) :: dx, dy
!     end type tp_DNSIndex
     ! private variables (accessible only by subroutines within this module
!     character(10), dimension(10), private :: varname_kplane
     type(tp_hyperslab), dimension(:,:), allocatable, private :: TSkplane
     type(tp_hyperslab), dimension(:), allocatable, private :: TSkgrid
     integer, private :: nfile ! Total number of files to read data from
     character(10), dimension(:), allocatable, private :: fname_irange
     integer, private :: ntpoint_total, npath, nvar_output, iminloc, imaxloc, iavelen
     character(200), private :: filepath
     type(tp_DNSIndex), private :: DNSIndex_k
     type(fprop), dimension(:), allocatable, private :: fileprop  !!!!!!!!

     ! public variables
     character(10), dimension(10) :: varname_kplane
     character(10), dimension(:), allocatable :: varout_k
     integer :: nxpoint_k, nypoint_k, imints_k, imaxts_k
contains


    subroutine ReadDNS_index_kplane(nt, dx, dy, dt, DNSIndex, filepath)
       integer, intent(out) :: nt
       real(8), intent(out) :: dx, dy, dt
       type(tp_DNSIndex), intent(out) :: DNSIndex
       character(*), intent(in) :: filepath
       type(tp_hyperslab) :: Index_kplane
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

       if(allocated(Index_kplane%dname)) deallocate(Index_kplane%dname)
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

       if(allocated(Index_kplane%dname)) deallocate(Index_kplane%dname)
       Index_kplane%dnum  = 1
       allocate(Index_kplane%dname(Index_kplane%dnum))
       Index_kplane%dname(1)  = 'zloc'
       Index_kplane%dimsf = (/DNSIndex%nkplane/)
       Index_kplane%dimsm = Index_kplane%dimsf
       Index_kplane%block = Index_kplane%dimsm
       if(allocated(buffer_real)) deallocate(buffer_real)
       if(allocated(DNSIndex%zloc)) deallocate(DNSIndex%zloc)
       allocate(buffer_real(DNSIndex%nkplane,1))
       allocate(DNSIndex%zloc(DNSIndex%nkplane))
       call ReadTSHDF5_1D(Index_kplane,buffer_real)
       DNSIndex%zloc = buffer_real(:,1)

     end subroutine ReadDNS_index_kplane



     subroutine ReadVOLTSkplane_perFile(fn, kloc, buffer, buffer_grid)
       implicit none

       character(*), intent(in) :: fn
       integer :: kloc !!!
       real(8), intent(out) :: buffer(ntpoint_total, nypoint_k,imints_k:imaxts_k,1,nvar_output)
       real(8), intent(out), optional :: buffer_grid( nypoint_k,imints_k:imaxts_k,1,12)

!  print *, 'dnum = ', TSkplane(1,1)%dnum
!  print *, 'dname = ', TSkplane(1,1)%dname
!  print *, 'dimsf = ', TSkplane(1,1)%dimsf

       TSkplane(1,1)%fname = trim(fn)
       TSkplane(1,1)%offset(4) = kloc - 1

       if(present(buffer_grid)) TSkgrid(1)%offset(3) = kloc - 1
!       print *, 'offset(4) = ', TSkplane(1,1)%offset(4)



       call ReadTSHDF5_4D(TSkplane(1,1), buffer)

       if(present(buffer_grid)) then
         if(.not.allocated(TSkgrid)) allocate(TSkgrid(1))
         TSkgrid(1)%fname = trim(fileprop(1)%filepath)//'timeseries_GridMetrics.h5'
         call ReadTSHDF5_3D(TSkgrid(1), buffer_grid)  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       endif

     end subroutine ReadVOLTSkplane_perFile

     subroutine ReadVOLTSkplane(k, kloc, buffer, buffer_grid)
        implicit none
        integer, intent(in) :: k, kloc
        real(8), intent(out) :: buffer(ntpoint_total, nypoint_k, imints_k:imaxts_k, 1, nvar_output) ! ntpoint
        real(8), intent(out), optional :: buffer_grid(nypoint_k, imints_k:imaxts_k, 1, 12)
        integer :: n, nnn, ktmp(1), i
        real(8) :: dxi, dyi, dzi
        character(8) :: fnum1
        integer :: nn
        integer :: ntpoint_tmp, num_file

         ntpoint_tmp = 0
         do nn=1, npath
           nnn = imints_k
           num_file = (fileprop(nn)%file_end - fileprop(nn)%file_be)/fileprop(nn)%file_skip + 1

           TSkplane(1,nn)%offset(4) = k - 1
           do n=1, num_file
             write(unit=fnum1,fmt='(I08.8)') fileprop(nn)%file_be + (n-1)*fileprop(nn)%file_skip
             if(iFileType.eq.1) then
               TSkplane(1,nn)%fname = trim(fileprop(nn)%filepath)//'timeseries_'//fnum1//'.h5'
             elseif(iFileType.eq.2) then
               TSkplane(1,nn)%fname = trim(fileprop(nn)%filepath)//'timeseries_kplane_'//fnum1//'.h5'
             endif
!             print *, 'Reading file: ', trim(TSkplane(1,nn)%fname)

             call ReadTSHDF5_4D(TSkplane(1,nn), buffer((ntpoint_tmp+1):(ntpoint_tmp+TSkplane(1,nn)%dimsf(1)), 1:TSkplane(1,nn)%dimsf(2), &
                                nnn:(nnn+TSkplane(1,nn)%dimsf(3)-1), 1:TSkplane(1,nn)%dimsf(4), 1:TSkplane(1,nn)%dnum) )
 ! should be changed
 !            if(present(buffer_grid).and.TSkgrid(1)%IsHSInitialized.and.n.eq.1.and.nn.eq.1) then
 !              TSkgrid(1)%fname = trim(fileprop(1)%filepath)//'timeseries_GridMetrics.h5'
 !              call ReadTSHDF5_3D(TSkgrid(1), buffer_grid(1:TSkgrid(1)%dimsf(1), &
 !                                 nnn:(nnn+TSkgrid(1)%dimsf(2)-1), 1:TSkgrid(1)%dimsf(3), 1:TSkgrid(1)%dnum) )
 !            endif
           ntpoint_tmp = ntpoint_tmp + TSkplane(1,nn)%dimsf(1)

           enddo ! end n loop
         enddo ! end nn loop

     end subroutine ReadVOLTSkplane

     subroutine ReadTSkplane(kloc, buffer, buffer_grid, kgrid_offset)
        implicit none
        integer, intent(in) :: kloc
        real(8), intent(out) :: buffer(ntpoint_total, nypoint_k, imints_k:imaxts_k, 1, nvar_output) ! ntpoint
        real(8), intent(out), optional :: buffer_grid(nypoint_k, imints_k:imaxts_k, 1, 12)
        integer, intent(in), optional :: kgrid_offset
        integer :: n, nnn, ktmp(1), i
        real(8) :: dxi, dyi, dzi
        character(4) :: fnum1
        integer :: nn
        integer :: ntpoint_tmp

        ! test
        ! print *, 'dimsf =  ', TSkplane(1,1)%dimsf
        ! print *, 'dimsm = ', TSkplane(1,1)%dimsm
        ! print *, 'offset = ', TSkplane(1,1)%offset

        ! print *, 'dimsf =  ', TSkplane(2,1)%dimsf
        ! print *, 'dimsm = ', TSkplane(2,1)%dimsm
        ! print *, 'offset = ', TSkplane(2,1)%offset

!         nnn = imints_k
         ntpoint_tmp = 0
         do nn=1, npath
           nnn = imints_k
           do n=1, nfile
             write(unit=fnum1,fmt='(I04.4)') kloc
             TSkplane(n,nn)%fname = trim(fileprop(nn)%filepath)//'timeseries_kplane'//fnum1//'_'//fname_irange(n)//'.h5'

             call ReadTSHDF5_4D(TSkplane(n,nn), buffer(ntpoint_tmp+1:ntpoint_tmp+TSkplane(n,nn)%dimsf(1),1:TSkplane(n,nn)%dimsf(2), &
                                        nnn:(nnn+TSkplane(n,nn)%dimsf(3)-1),1:TSkplane(n,nn)%dimsf(4),1:TSkplane(n,nn)%dnum) )
   ! should be changed !!!!!!!!!!!!
   !          if(present(buffer_grid).and.TSkgrid(n)%IsHSInitialized.and.nn.eq.1) then
   !             TSkgrid(n)%fname = trim(fileprop(nn)%filepath)//'timeseries_kgrid'//fnum1//'_'//fname_irange(n)//'.h5'
   !             call ReadTSHDF5_3D(TSkgrid(n), buffer_grid(1:TSkgrid(n)%dimsf(1),&
   !                                 nnn:(nnn+TSkgrid(n)%dimsf(2)-1),1:TSkgrid(n)%dimsf(3),1:TSkgrid(n)%dnum) )
   !          endif
             nnn = nnn+TSkplane(n,nn)%dimsf(3)
           enddo ! End loop n = 1, nfile
           ntpoint_tmp = ntpoint_tmp + TSkplane(1,nn)%dimsf(1)
         enddo ! end nn loop

         if(present(buffer_grid)) then
           TSkgrid(1)%offset(3) =  kgrid_offset
           TSkgrid(1)%fname = trim(fileprop(1)%filepath)//'timeseries_GridMetrics.h5'
           call ReadTSHDF5_3D(TSkgrid(1), buffer_grid)
         endif

   !      if(present(buffer_grid).and.(.not.TSkgrid(1)%IsHSInitialized)) then
   !         buffer_grid(:,:, 1, 3) = 0.d0
   !         do i = imints_k, imaxts_k
   !            buffer_grid(:,i, 1, 1) = (i-1)*DNSIndex_k%dx       ! x
   !         enddo
!             TSkgrid(nn)%dname(1) = "x"
!             TSkgrid(nn)%dname(2) = "y"
!             TSkgrid(nn)%dname(3) = "z"
!             TSkgrid(nn)%dname(4) = "didx"
!             TSkgrid(nn)%dname(5) = "djdx"
!             TSkgrid(nn)%dname(6) = "dkdx"
!             TSkgrid(nn)%dname(7) = "didy"
!             TSkgrid(nn)%dname(8) = "djdy"
!             TSkgrid(nn)%dname(9) = "dkdy"
!             TSkgrid(nn)%dname(10) = "didz"
!             TSkgrid(nn)%dname(11) = "djdz"
!             TSkgrid(nn)%dname(12) = "dkdz"

   !         ktmp = minloc(DNSIndex_k%kplane-kloc)
   !         buffer_grid(:,:, 1, 3)  =  DNSIndex_k%zloc(ktmp(1)) ! z

   !         dxi = 1.d0/DNSIndex_k%dx
   !         dyi = 1.d0/DNSIndex_k%dy
   !         dzi = 1.d0/DNSIndex_k%dzdk(ktmp(1))
   !         buffer_grid(:,:, 1, 4)  =  dxi  ! didx
   !         buffer_grid(:,:, 1, 8)  =  dyi  ! djdy
   !         buffer_grid(:,:, 1, 12) =  dzi  ! dkdz
   !      endif
      end subroutine ReadTSkplane

      subroutine InitVOLTSKplane(nt, np, DNSIndex, ibe_ave_inp, iend_ave_inp, jbe_ave_inp, jend_ave_inp, &
                  iwinl, iwinr, nvarout, varindex_output, fp, ireadgrid)
        implicit none
        integer, intent(in) :: nt, np
        type(tp_DNSIndex), intent(in) :: DNSIndex
        integer, intent(in) :: ibe_ave_inp, iend_ave_inp, jbe_ave_inp, jend_ave_inp, iwinl, iwinr
        integer, intent(in) :: nvarout, varindex_output(nvarout)
!        character(*), intent(in) :: fpath
        type(fprop), intent(in) :: fp(:)
        integer, intent(in), optional :: ireadgrid

        integer :: ibe_DNS, iend_DNS, iskip_DNS
        integer :: jbe_DNS, jend_DNS, jskip_DNS
        integer :: ibe_ave, iend_ave, jbe_ave, jend_ave
!        integer :: ibuffer
        integer :: n_total, ioffset_lower, inum_lower, inum_upper !!!!!!!!!!!!!!
        integer :: ibe_pfile, iend_pfile, n_lower, n_upper        !!!!!!!!!!!!!!
        integer :: n, nn, nnn
        character(4) :: fnum_ibe, fnum_iend
        integer :: nzloc
        real(8) :: dxs, dys

        DNSIndex_k = DNSIndex
        ntpoint_total = nt                 ! ntpoint
        nvar_output = nvarout
        npath = np
        if(allocated(fileprop)) deallocate(fileprop)
        allocate(fileprop(npath))
        fileprop = fp

        varname_kplane(1) = 'u'
        varname_kplane(2) = 'v'
        varname_kplane(3) = 'w'
        varname_kplane(4) = 'p'
        varname_kplane(5) = 'T'
        varname_kplane(6) = 'uk'
        varname_kplane(7) = 'vk'
        varname_kplane(8) = 'wk'
        varname_kplane(9) = 'pk'
        varname_kplane(10) = 'Tk'

!        ibuffer   = DNSIndex_k%ibuffer
        ibe_DNS   = DNSIndex_k%ibe
        iend_DNS  = DNSIndex_k%iend
        iskip_DNS = DNSIndex_k%iskip
        jbe_DNS   = DNSIndex_k%jbe
        jend_DNS  = DNSIndex_k%jend
        jskip_DNS = DNSIndex_k%jskip

        if(ibe_ave_inp.gt.iend_ave_inp.or.jbe_ave_inp.gt.jend_ave_inp.or. &
           ibe_ave_inp.lt.ibe_DNS.or.iend_ave_inp.gt.iend_DNS.or.jbe_ave_inp.lt.jbe_DNS.or.jend_ave_inp.gt.jend_DNS) then
            print *, 'Ranges for averaging is out of bound ... STOP'
            print *, 'ibe_ave_inp =', ibe_ave_inp, 'iend_ave_inp = ', iend_ave_inp
            print *, 'jbe_ave_inp =', jbe_ave_inp, 'jend_ave_inp = ', jend_ave_inp
            if(ibe_ave_inp.lt.ibe_DNS) print *, 'ibe_ave_inp < ibe_DNS'
            if(iend_ave_inp.gt.iend_DNS) print *, 'iend_ave_inp > iend_DNS'
            if(jbe_ave_inp.lt.jbe_DNS) print *, 'jbe_ave_inp < jbe_DNS'
            if(jend_ave_inp.gt.jend_DNS) print *, 'jend_ave_inp > jend_DNS'
            stop
        endif

        if(ibe_DNS.gt.iend_DNS.or.jbe_DNS.gt.jend_DNS.or.iskip_DNS.le.0.or.jskip_DNS.le.0) then
            print *, 'Input for DNS timeseries files out of bound ... STOP'
            print *, 'ibe_DNS =', ibe_DNS, 'iend_DNS = ', iend_DNS, 'iskip_DNS =', iskip_DNS
            print *, 'jbe_DNS =', jbe_DNS, 'jend_DNS = ', jend_DNS, 'jskip_DNS =', jskip_DNS
!            print *, 'ibuffer =', ibuffer
            stop
        endif

        ! begin, end, and number of spatial points for Averaging

        if(mod(ibe_ave_inp-ibe_DNS,iskip_DNS).ne.0) then
          ibe_ave = ibe_ave_inp + iskip_DNS  - mod(ibe_ave_inp-ibe_DNS,iskip_DNS)
        else
          ibe_ave = ibe_ave_inp
        endif
        nxpoint_k = (iend_ave_inp-ibe_ave)/iskip_DNS +1
        iend_ave = ibe_ave + (nxpoint_k-1)*iskip_DNS

        if(mod(jbe_ave_inp-jbe_DNS,jskip_DNS).ne.0) then
          jbe_ave = jbe_ave_inp + jskip_DNS - mod(jbe_ave_inp-jbe_DNS,jskip_DNS)
        else
          jbe_ave = jbe_ave_inp
        endif
        nypoint_k = (jend_ave_inp-jbe_ave)/jskip_DNS + 1
        jend_ave = jbe_ave + (nypoint_k-1)*jskip_DNS

        print *, '######################################################'
        print *, 'DNS Output timeseries spatial index range'
        print *, 'ibe = ', ibe_DNS, 'iend = ', iend_DNS
        print *, 'jbe = ', jbe_DNS, 'jend = ', jend_DNS
        print *, 'Actual Spatial Average range'
        print *, 'ibe_ave = ', ibe_ave, 'iend_ave = ', iend_ave
        print *, 'jbe_ave = ', jbe_ave, 'jend_ave = ', jend_ave
        print *, 'Number of spatial points for Averaging'
        print *, ' nxpoint_k = ', nxpoint_k, 'nypoint_k = ', nypoint_k

        ! Dimension range for tsdata

        imints_k = 1 - min(iwinl, (ibe_ave-ibe_DNS)/iskip_DNS)
        imaxts_k = nxpoint_k + min(iwinr, (iend_DNS-iend_ave)/iskip_DNS)
        iavelen = imaxts_k - imints_k + 1
        print *, 'streamwise Index range for kplane buffer'
        print *, ' imints_k = ', imints_k, ' imaxts_k = ', imaxts_k, 'iavelen (buffer length) = ', iavelen


        iminloc = max(ibe_ave  - iwinl*iskip_DNS,ibe_DNS)
        imaxloc = min(iend_ave + iwinr*iskip_DNS,iend_DNS)
        print *, 'Physical i-index range read TS (including xwindow)'
        print *, '   iminloc=', iminloc, 'imaxloc=', imaxloc

        if(allocated(TSkplane)) deallocate(TSkplane)
        allocate(TSkplane(1,npath))
        TSkplane%gname = '/kplane'
        TSkplane%rank = 4



        do nnn=1, npath
          allocate(TSkplane(1,nnn)%dname(nvar_output))
          allocate(TSkplane(1,nnn)%dimsf(TSkplane(1,nnn)%rank))
          allocate(TSkplane(1,nnn)%dimsm(TSkplane(1,nnn)%rank))
          allocate(TSkplane(1,nnn)%count(TSkplane(1,nnn)%rank))
          allocate(TSkplane(1,nnn)%offset(TSkplane(1,nnn)%rank))
          allocate(TSkplane(1,nnn)%block(TSkplane(1,nnn)%rank))
          allocate(TSkplane(1,nnn)%stride(TSkplane(1,nnn)%rank))
       enddo

       do nnn=1, npath
            TSkplane(1,nnn)%dimsf(1)  = fileprop(nnn)%ntpoint                      ! ntpoint
            TSkplane(1,nnn)%dimsf(2)  = jend_ave - jbe_ave + 1
!            TSkplane(1,nnn)%dimsf(3)  = iend_ave - ibe_ave + 1
            TSkplane(1,nnn)%dimsf(3)  = imaxts_k - imints_k + 1
            TSkplane(1,nnn)%dimsf(4)  = 1                    !!!!!!!!!    should be changed
            TSkplane(1,nnn)%offset(1) = 0
            TSkplane(1,nnn)%offset(2) = jbe_ave - jbe_DNS          !!!!
            TSkplane(1,nnn)%offset(3) = ibe_ave - ibe_DNS - iwinl        !!!!
    print *, 'offset(3) = ', ibe_ave - ibe_DNS - iwinl
            TSkplane(1,nnn)%offset(4) = 0                   !!!!!!!!!!!! should be changed
            TSkplane(1,nnn)%dimsm(1)  = TSkplane(1,nnn)%dimsf(1)
            TSkplane(1,nnn)%dimsm(2)  = TSkplane(1,nnn)%dimsf(2)
            TSkplane(1,nnn)%dimsm(3)  = TSkplane(1,nnn)%dimsf(3)
            TSkplane(1,nnn)%dimsm(4)  = TSkplane(1,nnn)%dimsf(4)
            TSkplane(1,nnn)%block(1)  = TSkplane(1,nnn)%dimsf(1)
            TSkplane(1,nnn)%block(2)  = TSkplane(1,nnn)%dimsf(2)
            TSkplane(1,nnn)%block(3)  = TSkplane(1,nnn)%dimsf(3)
            TSkplane(1,nnn)%block(4)  = TSkplane(1,nnn)%dimsf(4)
            TSkplane(1,nnn)%count(1)  = 1
            TSkplane(1,nnn)%count(2)  = 1
            TSkplane(1,nnn)%count(3)  = 1
            TSkplane(1,nnn)%count(4)  = 1
            TSkplane(1,nnn)%stride(1) = 1
            TSkplane(1,nnn)%stride(2) = 1
            TSkplane(1,nnn)%stride(3) = 1
            TSkplane(1,nnn)%stride(4) = 1
       enddo

       if(allocated(varout_k)) deallocate(varout_k)
       allocate(varout_k(nvar_output))
       do nnn=1, npath
         do n=1, nvar_output
           TSkplane(1,nnn)%dname(n) = varname_kplane(varindex_output(n))
         enddo
         TSkplane(1,nnn)%dnum       = nvar_output
       enddo

       varout_k = TSkplane(1,1)%dname
       print *, 'Variable name read in:  ', (trim(TSkplane(1,1)%dname(n))//',', n=1, nvar_output) !!!!!
       TSkplane%IsHSInitialized = .true.

       if(present(ireadgrid)) then
          if(ireadgrid.eq.1) then
             allocate(TSkgrid(1))
             TSkgrid%gname = "/kplane"
             TSkgrid%rank = 3

             allocate( TSkgrid(1)%dname(12))
             allocate( TSkgrid(1)%dimsf(TSkgrid(1)%rank))
             allocate( TSkgrid(1)%dimsm(TSkgrid(1)%rank))
             allocate( TSkgrid(1)%count(TSkgrid(1)%rank))
             allocate( TSkgrid(1)%offset(TSkgrid(1)%rank))
             allocate( TSkgrid(1)%block(TSkgrid(1)%rank))
             allocate( TSkgrid(1)%stride(TSkgrid(1)%rank))
             TSkgrid(1)%dnum = 12
             TSkgrid(1)%dname(1) = "x"
             TSkgrid(1)%dname(2) = "y"
             TSkgrid(1)%dname(3) = "z"
             TSkgrid(1)%dname(4) = "didx"
             TSkgrid(1)%dname(5) = "djdx"
             TSkgrid(1)%dname(6) = "dkdx"
             TSkgrid(1)%dname(7) = "didy"
             TSkgrid(1)%dname(8) = "djdy"
             TSkgrid(1)%dname(9) = "dkdy"
             TSkgrid(1)%dname(10) = "didz"
             TSkgrid(1)%dname(11) = "djdz"
             TSkgrid(1)%dname(12) = "dkdz"

             TSkgrid(1)%dimsf(1) = jend_ave - jbe_ave + 1
!             TSkgrid(1)%dimsf(2) = iend_ave - ibe_ave + 1
             TSkgrid(1)%dimsf(2) = imaxts_k - imints_k + 1
             TSkgrid(1)%dimsf(3) = 1                       !! should be changed
             TSkgrid(1)%offset(1) = jbe_ave - jbe_DNS
             TSkgrid(1)%offset(2) = ibe_ave - ibe_DNS - iwinl
             TSkgrid(1)%offset(3) = 0                      !! should be changed

             TSkgrid(1)%dimsm(1) = TSkgrid(1)%dimsf(1)
             TSkgrid(1)%dimsm(2) = TSkgrid(1)%dimsf(2)
             TSkgrid(1)%dimsm(3) = TSkgrid(1)%dimsf(3)
             TSkgrid(1)%block(1) = TSkgrid(1)%dimsf(1)
             TSkgrid(1)%block(2) = TSkgrid(1)%dimsf(2)
             TSkgrid(1)%block(3) = TSkgrid(1)%dimsf(3)
             TSkgrid(1)%count(1) = 1
             TSkgrid(1)%count(2) = 1
             TSkgrid(1)%count(3) = 1
             TSkgrid(1)%stride(1) = 1
             TSkgrid(1)%stride(2) = 1
             TSkgrid(1)%stride(3) = 1
             TSkgrid(1)%IsHSInitialized = .true.
           endif
         endif


      end subroutine InitVOLTSKplane


      subroutine InitVOLTSKplane_NoPrintInfo(nt, np, DNSIndex, ibe_ave_inp, iend_ave_inp, jbe_ave_inp, jend_ave_inp, &
                  iwinl, iwinr, nvarout, varindex_output, fp, ireadgrid)
        implicit none
        integer, intent(in) :: nt, np
        type(tp_DNSIndex), intent(in) :: DNSIndex
        integer, intent(in) :: ibe_ave_inp, iend_ave_inp, jbe_ave_inp, jend_ave_inp, iwinl, iwinr
        integer, intent(in) :: nvarout, varindex_output(nvarout)
        type(fprop), intent(in) :: fp(:)
        integer, intent(in), optional :: ireadgrid

        integer :: ibe_DNS, iend_DNS, iskip_DNS
        integer :: jbe_DNS, jend_DNS, jskip_DNS
        integer :: ibe_ave, iend_ave, jbe_ave, jend_ave
        integer :: n_total, ioffset_lower, inum_lower, inum_upper !!!!!!!!!!!!!!
        integer :: ibe_pfile, iend_pfile, n_lower, n_upper        !!!!!!!!!!!!!!
        integer :: n, nn, nnn
        character(4) :: fnum_ibe, fnum_iend
        integer :: nzloc
        real(8) :: dxs, dys

        DNSIndex_k = DNSIndex
        ntpoint_total = nt                 ! ntpoint
        nvar_output = nvarout
        npath = np
        if(allocated(fileprop)) deallocate(fileprop)
        allocate(fileprop(npath))
        fileprop = fp

        varname_kplane(1) = 'u'
        varname_kplane(2) = 'v'
        varname_kplane(3) = 'w'
        varname_kplane(4) = 'p'
        varname_kplane(5) = 'T'
        varname_kplane(6) = 'uk'
        varname_kplane(7) = 'vk'
        varname_kplane(8) = 'wk'
        varname_kplane(9) = 'pk'
        varname_kplane(10) = 'Tk'

        ibe_DNS   = DNSIndex_k%ibe
        iend_DNS  = DNSIndex_k%iend
        iskip_DNS = DNSIndex_k%iskip
        jbe_DNS   = DNSIndex_k%jbe
        jend_DNS  = DNSIndex_k%jend
        jskip_DNS = DNSIndex_k%jskip

        if(ibe_ave_inp.gt.iend_ave_inp.or.jbe_ave_inp.gt.jend_ave_inp.or. &
           ibe_ave_inp.lt.ibe_DNS.or.iend_ave_inp.gt.iend_DNS.or.jbe_ave_inp.lt.jbe_DNS.or.jend_ave_inp.gt.jend_DNS) then
            print *, 'Ranges for averaging is out of bound ... STOP'
            print *, 'ibe_ave_inp =', ibe_ave_inp, 'iend_ave_inp = ', iend_ave_inp
            print *, 'jbe_ave_inp =', jbe_ave_inp, 'jend_ave_inp = ', jend_ave_inp
            if(ibe_ave_inp.lt.ibe_DNS) print *, 'ibe_ave_inp < ibe_DNS'
            if(iend_ave_inp.gt.iend_DNS) print *, 'iend_ave_inp > iend_DNS'
            if(jbe_ave_inp.lt.jbe_DNS) print *, 'jbe_ave_inp < jbe_DNS'
            if(jend_ave_inp.gt.jend_DNS) print *, 'jend_ave_inp > jend_DNS'
            stop
        endif

        if(ibe_DNS.gt.iend_DNS.or.jbe_DNS.gt.jend_DNS.or.iskip_DNS.le.0.or.jskip_DNS.le.0) then
            print *, 'Input for DNS timeseries files out of bound ... STOP'
            print *, 'ibe_DNS =', ibe_DNS, 'iend_DNS = ', iend_DNS, 'iskip_DNS =', iskip_DNS
            print *, 'jbe_DNS =', jbe_DNS, 'jend_DNS = ', jend_DNS, 'jskip_DNS =', jskip_DNS
            stop
        endif

        ! begin, end, and number of spatial points for Averaging

        if(mod(ibe_ave_inp-ibe_DNS,iskip_DNS).ne.0) then
          ibe_ave = ibe_ave_inp + iskip_DNS  - mod(ibe_ave_inp-ibe_DNS,iskip_DNS)
        else
          ibe_ave = ibe_ave_inp
        endif
        nxpoint_k = (iend_ave_inp-ibe_ave)/iskip_DNS +1
        iend_ave = ibe_ave + (nxpoint_k-1)*iskip_DNS

        if(mod(jbe_ave_inp-jbe_DNS,jskip_DNS).ne.0) then
          jbe_ave = jbe_ave_inp + jskip_DNS - mod(jbe_ave_inp-jbe_DNS,jskip_DNS)
        else
          jbe_ave = jbe_ave_inp
        endif
        nypoint_k = (jend_ave_inp-jbe_ave)/jskip_DNS + 1
        jend_ave = jbe_ave + (nypoint_k-1)*jskip_DNS

        print *, '######################################################'
        print *, 'Actual Spatial Average range'
        print *, 'ibe_ave = ', ibe_ave, 'iend_ave = ', iend_ave
        print *, 'jbe_ave = ', jbe_ave, 'jend_ave = ', jend_ave

        ! Dimension range for tsdata

        imints_k = 1 - min(iwinl, (ibe_ave-ibe_DNS)/iskip_DNS)
        imaxts_k = nxpoint_k + min(iwinr, (iend_DNS-iend_ave)/iskip_DNS)
        iavelen = imaxts_k - imints_k + 1

        iminloc = max(ibe_ave  - iwinl*iskip_DNS,ibe_DNS)
        imaxloc = min(iend_ave + iwinr*iskip_DNS,iend_DNS)

        if(allocated(TSkplane)) deallocate(TSkplane)
        allocate(TSkplane(1,npath))
        TSkplane%gname = '/kplane'
        TSkplane%rank = 4

        do nnn=1, npath
          allocate(TSkplane(1,nnn)%dname(nvar_output))
          allocate(TSkplane(1,nnn)%dimsf(TSkplane(1,nnn)%rank))
          allocate(TSkplane(1,nnn)%dimsm(TSkplane(1,nnn)%rank))
          allocate(TSkplane(1,nnn)%count(TSkplane(1,nnn)%rank))
          allocate(TSkplane(1,nnn)%offset(TSkplane(1,nnn)%rank))
          allocate(TSkplane(1,nnn)%block(TSkplane(1,nnn)%rank))
          allocate(TSkplane(1,nnn)%stride(TSkplane(1,nnn)%rank))
       enddo

       do nnn=1, npath
            TSkplane(1,nnn)%dimsf(1)  = fileprop(nnn)%ntpoint                      ! ntpoint
            TSkplane(1,nnn)%dimsf(2)  = jend_ave - jbe_ave + 1

            TSkplane(1,nnn)%dimsf(3)  = imaxts_k - imints_k + 1
            TSkplane(1,nnn)%dimsf(4)  = 1                    !!!!!!!!!    should be changed
            TSkplane(1,nnn)%offset(1) = 0
            TSkplane(1,nnn)%offset(2) = jbe_ave - jbe_DNS
            TSkplane(1,nnn)%offset(3) = ibe_ave - ibe_DNS - iwinl
            TSkplane(1,nnn)%offset(4) = 0                   !!!!!!!!!!!! should be changed
            TSkplane(1,nnn)%dimsm(1)  = TSkplane(1,nnn)%dimsf(1)
            TSkplane(1,nnn)%dimsm(2)  = TSkplane(1,nnn)%dimsf(2)
            TSkplane(1,nnn)%dimsm(3)  = TSkplane(1,nnn)%dimsf(3)
            TSkplane(1,nnn)%dimsm(4)  = TSkplane(1,nnn)%dimsf(4)
            TSkplane(1,nnn)%block(1)  = TSkplane(1,nnn)%dimsf(1)
            TSkplane(1,nnn)%block(2)  = TSkplane(1,nnn)%dimsf(2)
            TSkplane(1,nnn)%block(3)  = TSkplane(1,nnn)%dimsf(3)
            TSkplane(1,nnn)%block(4)  = TSkplane(1,nnn)%dimsf(4)
            TSkplane(1,nnn)%count(1)  = 1
            TSkplane(1,nnn)%count(2)  = 1
            TSkplane(1,nnn)%count(3)  = 1
            TSkplane(1,nnn)%count(4)  = 1
            TSkplane(1,nnn)%stride(1) = 1
            TSkplane(1,nnn)%stride(2) = 1
            TSkplane(1,nnn)%stride(3) = 1
            TSkplane(1,nnn)%stride(4) = 1
       enddo

       if(allocated(varout_k)) deallocate(varout_k)
       allocate(varout_k(nvar_output))
       do nnn=1, npath
         do n=1, nvar_output
           TSkplane(1,nnn)%dname(n) = varname_kplane(varindex_output(n))
         enddo
         TSkplane(1,nnn)%dnum       = nvar_output
       enddo

       varout_k = TSkplane(1,1)%dname
       !  print *, 'Variable name read in:  ', (trim(TSkplane(1,1)%dname(n))//',', n=1, nvar_output) !!!!!
       TSkplane%IsHSInitialized = .true.

      end subroutine InitVOLTSKplane_NoPrintInfo


      subroutine InitTSkplane(nt, np, DNSIndex, ibe_ave_inp, iend_ave_inp, jbe_ave_inp, jend_ave_inp, &
                  iwinl, iwinr, nvarout, varindex_output, fp, ireadgrid)
        implicit none
        integer, intent(in) :: nt, np
        type(tp_DNSIndex), intent(in) :: DNSIndex
        integer, intent(in) :: ibe_ave_inp, iend_ave_inp, jbe_ave_inp, jend_ave_inp, iwinl, iwinr
        integer, intent(in) :: nvarout, varindex_output(nvarout)
!        character(*), intent(in) :: fpath
        type(fprop), intent(in) :: fp(:)
        integer, intent(in), optional :: ireadgrid

        integer :: ibe_DNS, iend_DNS, iskip_DNS
        integer :: jbe_DNS, jend_DNS, jskip_DNS
        integer :: ibe_ave, iend_ave, jbe_ave, jend_ave
        integer :: ibuffer
        integer :: n_total, ioffset_lower, inum_lower, inum_upper
        integer :: ibe_pfile, iend_pfile, n_lower, n_upper
        integer :: n, nn, nnn
        character(4) :: fnum_ibe, fnum_iend
        integer :: nzloc
        real(8) :: dxs, dys

        DNSIndex_k = DNSIndex
        ntpoint_total = nt                 ! ntpoint
        nvar_output = nvarout
        npath = np
        allocate(fileprop(npath))
        fileprop = fp
!        filepath = trim(fpath)

        varname_kplane(1) = 'u'
        varname_kplane(2) = 'v'
        varname_kplane(3) = 'w'
        varname_kplane(4) = 'p'
        varname_kplane(5) = 'T'
        varname_kplane(6) = 'uk'
        varname_kplane(7) = 'vk'
        varname_kplane(8) = 'wk'
        varname_kplane(9) = 'pk'
        varname_kplane(10) = 'Tk'

        ibuffer   = DNSIndex_k%ibuffer
        ibe_DNS   = DNSIndex_k%ibe
        iend_DNS  = DNSIndex_k%iend
        iskip_DNS = DNSIndex_k%iskip
        jbe_DNS   = DNSIndex_k%jbe
        jend_DNS  = DNSIndex_k%jend
        jskip_DNS = DNSIndex_k%jskip

        if(ibe_ave_inp.gt.iend_ave_inp.or.jbe_ave_inp.gt.jend_ave_inp.or. &
           ibe_ave_inp.lt.ibe_DNS.or.iend_ave_inp.gt.iend_DNS.or.jbe_ave_inp.lt.jbe_DNS.or.jend_ave_inp.gt.jend_DNS) then
            print *, 'Ranges for averaging is out of bound ... STOP'
            print *, 'ibe_ave_inp =', ibe_ave_inp, 'iend_ave_inp = ', iend_ave_inp
            print *, 'jbe_ave_inp =', jbe_ave_inp, 'jend_ave_inp = ', jend_ave_inp
            if(ibe_ave_inp.lt.ibe_DNS) print *, 'ibe_ave_inp < ibe_DNS'
            if(iend_ave_inp.gt.iend_DNS) print *, 'iend_ave_inp > iend_DNS'
            if(jbe_ave_inp.lt.jbe_DNS) print *, 'jbe_ave_inp < jbe_DNS'
            if(jend_ave_inp.gt.jend_DNS) print *, 'jend_ave_inp > jend_DNS'
            stop
        endif

        if(ibe_DNS.gt.iend_DNS.or.jbe_DNS.gt.jend_DNS.or.iskip_DNS.le.0.or.jskip_DNS.le.0.or.ibuffer.le.0) then
            print *, 'Input for DNS timeseries files out of bound ... STOP'
            print *, 'ibe_DNS =', ibe_DNS, 'iend_DNS = ', iend_DNS, 'iskip_DNS =', iskip_DNS
            print *, 'jbe_DNS =', jbe_DNS, 'jend_DNS = ', jend_DNS, 'jskip_DNS =', jskip_DNS
            print *, 'ibuffer =', ibuffer
            stop
        endif

        ! begin, end, and number of spatial points for Averaging

        if(mod(ibe_ave_inp-ibe_DNS,iskip_DNS).ne.0) then
          ibe_ave = ibe_ave_inp + iskip_DNS  - mod(ibe_ave_inp-ibe_DNS,iskip_DNS)
        else
          ibe_ave = ibe_ave_inp
        endif
        nxpoint_k = (iend_ave_inp-ibe_ave)/iskip_DNS +1
        iend_ave = ibe_ave + (nxpoint_k-1)*iskip_DNS

        if(mod(jbe_ave_inp-jbe_DNS,jskip_DNS).ne.0) then
          jbe_ave = jbe_ave_inp + jskip_DNS - mod(jbe_ave_inp-jbe_DNS,jskip_DNS)
        else
          jbe_ave = jbe_ave_inp
        endif
        nypoint_k = (jend_ave_inp-jbe_ave)/jskip_DNS + 1
        jend_ave = jbe_ave + (nypoint_k-1)*jskip_DNS

!  print *, 'jskip_DNS = ', jskip_DNS

        print *, 'DNS Output timeseries spatial index range'
        print *, 'ibe = ', ibe_DNS, 'iend = ', iend_DNS
        print *, 'jbe = ', jbe_DNS, 'jend = ', jend_DNS
        print *, 'Actual Spatial Average range'
        print *, 'ibe_ave = ', ibe_ave, 'iend_ave = ', iend_ave
        print *, 'jbe_ave = ', jbe_ave, 'jend_ave = ', jend_ave
        print *, 'Number of spatial points for Averaging'
        print *, ' nxpoint_k = ', nxpoint_k, 'nypoint_k = ', nypoint_k

        ! Dimension range for tsdata

        imints_k = 1 - min(iwinl, (ibe_ave-ibe_DNS)/iskip_DNS)
        imaxts_k = nxpoint_k + min(iwinr, (iend_DNS-iend_ave)/iskip_DNS)
        iavelen = imaxts_k - imints_k + 1
        print *, 'streamwise Index range for kplane buffer'
        print *, ' imints_k = ', imints_k, ' imaxts_k = ', imaxts_k, 'iavelen (buffer length) = ', iavelen


        iminloc = max(ibe_ave  - iwinl*iskip_DNS,ibe_DNS)
        imaxloc = min(iend_ave + iwinr*iskip_DNS,iend_DNS)
        print *, 'Physical i-index range read TS (including xwindow)'
        print *, '   iminloc=', iminloc, 'imaxloc=', imaxloc


        ! Total number of files that cover the full range of DNS output in i dir
        n_total = ceiling( ((iend_DNS - ibe_DNS)/iskip_DNS + 1)/real(ibuffer)  )

        ! Determine the range of files that need to be read based on ibe_ave & iend_ave
        n_lower = 0; n_upper = 0 ! File indexes that are partially read
        ioffset_lower = 0; inum_lower = 0; inum_upper = 0
        do n = 1, n_total
           ibe_pfile = ibe_DNS + ((n-1)*ibuffer)*iskip_DNS ! Could also be obtained by reading from the file name
           iend_pfile = ibe_pfile + (ibuffer-1)*iskip_DNS  ! Could also be obtained by reading from the file name
           if(ibe_pfile.le.iminloc.and.iend_pfile.ge.iminloc) then
              n_lower = n
              ioffset_lower = (iminloc - ibe_pfile)/iskip_DNS
              inum_lower = min( (iend_pfile - iminloc)/iskip_DNS + 1, iavelen)
           endif
           if(ibe_pfile.le.imaxloc.and.iend_pfile.ge.imaxloc) then
              n_upper = n
              if(n_upper.gt.n_lower) inum_upper = (imaxloc - ibe_pfile)/iskip_DNS + 1
              exit
           endif
        enddo

       nfile = n_upper - n_lower + 1 ! Total number of files that need to be read

       allocate(fname_irange(nfile))
       nn = 0
       do n = n_lower, n_upper
          nn = nn + 1
          ibe_pfile  = ibe_DNS + ((n-1)*ibuffer)*iskip_DNS
          iend_pfile = min(ibe_pfile + (ibuffer-1)*iskip_DNS, iend_DNS)
          write(unit=fnum_ibe,fmt='(I04.4)') ibe_pfile
          write(unit=fnum_iend,fmt='(I04.4)') iend_pfile
          fname_irange(nn) =   'i'//fnum_ibe//'-'//fnum_iend
          print *, 'fname_irang =', fname_irange(nn)
       enddo


       allocate(TSkplane(nfile,npath))
       TSkplane%gname = '/kplane'
       TSkplane%rank = 4

       do nnn=1, npath
         do nn = 1, nfile
            allocate(TSkplane(nn,nnn)%dname(nvar_output))
            allocate(TSkplane(nn,nnn)%dimsf(TSkplane(nn,nnn)%rank))
            allocate(TSkplane(nn,nnn)%dimsm(TSkplane(nn,nnn)%rank))
            allocate(TSkplane(nn,nnn)%count(TSkplane(nn,nnn)%rank))
            allocate(TSkplane(nn,nnn)%offset(TSkplane(nn,nnn)%rank))
            allocate(TSkplane(nn,nnn)%block(TSkplane(nn,nnn)%rank))
            allocate(TSkplane(nn,nnn)%stride(TSkplane(nn,nnn)%rank))
        enddo
      enddo

       ! Files that are fully read
       do nnn=1, npath
         do nn = 1, nfile
            TSkplane(nn,nnn)%dimsf(1) = fileprop(nnn)%ntpoint                      ! ntpoint
            TSkplane(nn,nnn)%dimsf(2) = (jend_ave - jbe_ave)/jskip_DNS + 1
            TSkplane(nn,nnn)%dimsf(3) = ibuffer
            TSkplane(nn,nnn)%dimsf(4) = 1
            TSkplane(nn,nnn)%offset(1) = 0
            TSkplane(nn,nnn)%offset(2) = jbe_ave - 1
            TSkplane(nn,nnn)%offset(3) = 0
            TSkplane(nn,nnn)%offset(4) = 0
         enddo
       enddo

       ! Files that are partially read in i-dir
       do nnn=1, npath
         TSkplane(1,nnn)%dimsf(3) = inum_lower
         TSkplane(1,nnn)%offset(3) = ioffset_lower
       enddo
       if(nfile.gt.1) then
         do nnn=1, npath
           TSkplane(nfile,nnn)%dimsf(3) = inum_upper
         enddo
       endif

       do nnn=1, npath
         do nn = 1, nfile
            TSkplane(nn,nnn)%dimsm(1)  = TSkplane(nn,nnn)%dimsf(1)
            TSkplane(nn,nnn)%dimsm(2)  = TSkplane(nn,nnn)%dimsf(2)
            TSkplane(nn,nnn)%dimsm(3)  = TSkplane(nn,nnn)%dimsf(3)
            TSkplane(nn,nnn)%dimsm(4)  = TSkplane(nn,nnn)%dimsf(4)
            TSkplane(nn,nnn)%block(1)  = TSkplane(nn,nnn)%dimsf(1)
            TSkplane(nn,nnn)%block(2)  = TSkplane(nn,nnn)%dimsf(2)
            TSkplane(nn,nnn)%block(3)  = TSkplane(nn,nnn)%dimsf(3)
            TSkplane(nn,nnn)%block(4)  = TSkplane(nn,nnn)%dimsf(4)
            TSkplane(nn,nnn)%count(1)  = 1
            TSkplane(nn,nnn)%count(2)  = 1
            TSkplane(nn,nnn)%count(3)  = 1
            TSkplane(nn,nnn)%count(4)  = 1
            TSkplane(nn,nnn)%stride(1) = 1
            TSkplane(nn,nnn)%stride(2) = 1
            TSkplane(nn,nnn)%stride(3) = 1
            TSkplane(nn,nnn)%stride(4) = 1
         enddo
       enddo

       allocate(varout_k(nvar_output))
       do nnn=1, npath
         do nn = 1, nfile
            do n=1, nvar_output
               TSkplane(nn,nnn)%dname(n) = varname_kplane(varindex_output(n))
            enddo
            TSkplane(nn,nnn)%dnum      = nvar_output
         enddo
       enddo

!       varout_k = TSkplane(1)%dname  !
!       print *, 'Variable name read in:  ', (trim(TSkplane(1)%dname(n))//',', n=1, nvar_output) !!!!!
       varout_k = TSkplane(1,1)%dname
       print *, 'Variable name read in:  ', (trim(TSkplane(1,1)%dname(n))//',', n=1, nvar_output) !!!!!
       TSkplane%IsHSInitialized = .true.

       if(present(ireadgrid)) then
          if(ireadgrid.eq.1) then
            allocate(TSkgrid(1))
            TSkgrid(1)%gname = "/kplane"
            TSkgrid(1)%rank = 3
            allocate( TSkgrid(1)%dname(12))
            allocate( TSkgrid(1)%dimsf(TSkgrid(1)%rank))
            allocate( TSkgrid(1)%dimsm(TSkgrid(1)%rank))
            allocate( TSkgrid(1)%count(TSkgrid(1)%rank))
            allocate( TSkgrid(1)%offset(TSkgrid(1)%rank))
            allocate( TSkgrid(1)%block(TSkgrid(1)%rank))
            allocate( TSkgrid(1)%stride(TSkgrid(1)%rank))
            TSkgrid(1)%dnum = 12
            TSkgrid(1)%dname(1) = "x"
            TSkgrid(1)%dname(2) = "y"
            TSkgrid(1)%dname(3) = "z"
            TSkgrid(1)%dname(4) = "didx"
            TSkgrid(1)%dname(5) = "djdx"
            TSkgrid(1)%dname(6) = "dkdx"
            TSkgrid(1)%dname(7) = "didy"
            TSkgrid(1)%dname(8) = "djdy"
            TSkgrid(1)%dname(9) = "dkdy"
            TSkgrid(1)%dname(10) = "didz"
            TSkgrid(1)%dname(11) = "djdz"
            TSkgrid(1)%dname(12) = "dkdz"

            TSkgrid(1)%dimsf(1) = jend_ave - jbe_ave + 1
            TSkgrid(1)%dimsf(2) = imaxts_k - imints_k + 1
            TSkgrid(1)%dimsf(3) = 1
            TSkgrid(1)%offset(1) = jbe_ave - jbe_DNS
            TSkgrid(1)%offset(2) = ibe_ave - ibe_DNS - 2
            TSkgrid(1)%offset(3) = 0

            TSkgrid(1)%dimsm(1) = TSkgrid(1)%dimsf(1)
            TSkgrid(1)%dimsm(2) = TSkgrid(1)%dimsf(2)
            TSkgrid(1)%dimsm(3) = TSkgrid(1)%dimsf(3)
            TSkgrid(1)%block(1) = TSkgrid(1)%dimsf(1)
            TSkgrid(1)%block(2) = TSkgrid(1)%dimsf(2)
            TSkgrid(1)%block(3) = TSkgrid(1)%dimsf(3)
            TSkgrid(1)%count = 1
            TSkgrid(1)%stride = 1
            TSkgrid(1)%IsHSInitialized = .true.


!             allocate(TSkgrid(nfile))
!             TSkgrid%gname = "/kplane"
!             TSkgrid%rank = 3
!
!             do nn = 1, nfile
!                allocate( TSkgrid(nn)%dname(12))
!                allocate( TSkgrid(nn)%dimsf(TSkgrid(nn)%rank))
!                allocate( TSkgrid(nn)%dimsm(TSkgrid(nn)%rank))
!                allocate( TSkgrid(nn)%count(TSkgrid(nn)%rank))
!                allocate( TSkgrid(nn)%offset(TSkgrid(nn)%rank))
!                allocate( TSkgrid(nn)%block(TSkgrid(nn)%rank))
!                allocate( TSkgrid(nn)%stride(TSkgrid(nn)%rank))
!                TSkgrid(nn)%dnum = 12
!                TSkgrid(nn)%dname(1) = "x"
!                TSkgrid(nn)%dname(2) = "y"
!                TSkgrid(nn)%dname(3) = "z"
!                TSkgrid(nn)%dname(4) = "didx"
!                TSkgrid(nn)%dname(5) = "djdx"
!                TSkgrid(nn)%dname(6) = "dkdx"
!                TSkgrid(nn)%dname(7) = "didy"
!                TSkgrid(nn)%dname(8) = "djdy"
!                TSkgrid(nn)%dname(9) = "dkdy"
!                TSkgrid(nn)%dname(10) = "didz"
!                TSkgrid(nn)%dname(11) = "djdz"
!                TSkgrid(nn)%dname(12) = "dkdz"
!             enddo
!
!             do nn = 1, nfile
!                TSkgrid(nn)%dimsf(1) = jend_ave - jbe_ave + 1
!                TSkgrid(nn)%dimsf(2) = ibuffer
!                TSkgrid(nn)%dimsf(3) = 1
!                TSkgrid(nn)%offset(1) = jbe_ave - 1
!                TSkgrid(nn)%offset(2) = 0
!                TSkgrid(nn)%offset(3) = 0
!             enddo
!
!             TSkgrid(1)%dimsf(2) = inum_lower
!             TSkgrid(1)%offset(2) = ioffset_lower
!             if(nfile.gt.1) TSkgrid(nfile)%dimsf(2) = inum_upper
!
!             ! grid file info
!             do nn = 1, nfile
!                TSkgrid(nn)%dimsm(1) = TSkgrid(nn)%dimsf(1)
!                TSkgrid(nn)%dimsm(2) = TSkgrid(nn)%dimsf(2)
!                TSkgrid(nn)%dimsm(3) = TSkgrid(nn)%dimsf(3)
!                TSkgrid(nn)%block(1) = TSkgrid(nn)%dimsf(1)
!                TSkgrid(nn)%block(2) = TSkgrid(nn)%dimsf(2)
!                TSkgrid(nn)%block(3) = TSkgrid(nn)%dimsf(3)
!                TSkgrid(nn)%count(1) = 1
!                TSkgrid(nn)%count(2) = 1
!                TSkgrid(nn)%count(3) = 1
!                TSkgrid(nn)%stride(1) = 1
!                TSkgrid(nn)%stride(2) = 1
!                TSkgrid(nn)%stride(3) = 1
!                TSkgrid(nn)%IsHSInitialized = .true.
!             enddo
!          else  ! Assume Cartesian grid
!             DNSIndex_k%kplane = nzloc
!             open(101,file=trim(filepath)//'kplane_index.dat')  !!!!!!!!!!!!!!!1
!             rewind(101)
!             read(101,*)
!             read(101,*) dxs, dys, nzloc
!             DNSIndex_k%dx = dxs
!             DNSIndex_k%dy = dys
!             DNSIndex_k%nzloc = nzloc
!             allocate(DNSIndex_k%kplane(nzloc), DNSIndex_k%zloc(nzloc), DNSIndex_k%dzdk(nzloc))
!             read(101,*)
!             do nn = 1, nzloc
!                read(101,*) DNSIndex_k%kplane(nn), DNSIndex_k%zloc(nn), DNSIndex_k%dzdk(nn)
!             enddo
!             close(101)
          endif
       endif
    end subroutine InitTSkplane

end module MTSKplane

