!> convert timeseries volume to separated data files


    Program convert_timeseries

      use HDF5
      implicit none

      type hyperslab
        character(100) :: gname
        character(100), dimension(:), pointer :: dname
        character(100), dimension(5) :: attr_name
        integer(HID_T) :: group1_id, memspace1
        integer(HID_T) :: dspace1_id, dset1_id
        integer(HID_T) :: group2_id, dspace2_id, dset2_id, memspace2
        integer(HSIZE_T), dimension(:), allocatable :: dimsf_total
        integer(HSIZE_T), dimension(:), allocatable :: dimsf
        integer(HSIZE_T), dimension(:), allocatable :: dimsf_sub
        integer(HSIZE_T), dimension(:), allocatable :: dimsm
        integer(HSIZE_T), dimension(:), allocatable :: count
        integer(HSIZE_T), dimension(:), allocatable :: offset_rd
        integer(HSIZE_T), dimension(:), allocatable :: offset_wt
        integer(HSIZE_T), dimension(:), allocatable :: block
        integer(HSIZE_T), dimension(:), allocatable :: stride
        integer :: rank
        integer :: num_dset
        integer, dimension(3) :: info_1st_array, info_2nd_array
        integer(HSIZE_T), dimension(1) :: attr_dims_info, attr_dims_loc, attr_dims_no
        integer, dimension(:), pointer :: info_loc
        integer :: ibe, iend, iskip, jbe, jend, jskip, kbe, kend, kskip
      end type hyperslab
      type(hyperslab) :: TSkplane, TSiplane, TSjplane    !!!!!!!
      type(hyperslab) :: kgrid, igrid, jgrid             !!!!!!!!!!!!!!!
      integer(HID_T) :: file1_id, file2_id
      integer(HID_T) :: attr_id, attr_space, atype_id
      integer(HSIZE_T), dimension(1) :: dimsf_time, dimsf_sub_time
      integer(HID_T) :: dspace_id_time, dset_id_time
      character(100) :: dname_time = 'time'
      integer, dimension(1) :: attr_time
      integer, dimension(1) :: attr_ntpoint
      integer, dimension(1) :: attr_numkp, attr_numip, attr_numjp
      real(8), dimension(:,:,:,:,:), allocatable :: buffer_kplane, buffer_iplane, buffer_jplane
      real(8), dimension(:,:,:,:), allocatable :: buffer_kgrid, buffer_igrid, buffer_jgrid      !!!!!!!1
      real(8), dimension(:,:), allocatable :: buffer_time
      real(8), dimension(:), allocatable :: buffer_1D_time
      integer :: hdferr
      character(400) :: fname1, fname2
      character(4) :: fnum1, fnum_ibe, fnum_iend
      character(8) :: fnum2
      integer :: i, j, k, n, nn, num_file
      !!!!!!
      integer :: ireadkp, ireadip, ireadjp
      integer :: nv = 10 
      character(200) :: filepath
      integer :: file_be, file_end, file_skip
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer :: idx, ibuffer_k, jbuffer_i, ibuffer_j
      integer :: ibe_k, iend_k, iskip_k
      integer :: ibe_k_inp, iend_k_inp, iskip_k_inp
      integer :: ibe_k_file, iend_k_file, iskip_k_file
      integer :: jbe_i, jend_i, jskip_i
      integer :: jbe_i_inp, jend_i_inp, jskip_i_inp
      integer :: jbe_i_file, jend_i_file, jskip_i_file
      integer :: ibe_j, iend_j, iskip_j
      integer :: ibe_j_inp, iend_j_inp, iskip_j_inp
      integer :: ibe_j_file, iend_j_file, iskip_j_file
      integer :: num_kplane_inp, num_iplane_inp, num_jplane_inp
      integer :: nn_total
      integer :: num_kplane_be, num_kplane_end
      integer :: num_iplane_be, num_iplane_end
      integer :: num_jplane_be, num_jplane_end
      integer, dimension(1) :: ibuffer_k_tmp         !!!!!!!
   

      call Init()

      num_file = (file_end - file_be)/file_skip + 1

      call InitHDF5()

      if(ireadkp.eq.1) then
        call ReadkplaneAttr()
        call GetkplaneInfo()
        call GetkgridInfo()
        dimsf_time(1) = TSkplane%dimsf(1)*num_file
        dimsf_sub_time(1) = TSkplane%dimsf(1)
        allocate(buffer_time(TSkplane%dimsf(1), num_file), buffer_1D_time(dimsf_time(1)))
      endif


      if(ireadip.eq.1) then
        call ReadiplaneAttr()
        call GetiplaneInfo()
        if(.not.allocated(buffer_time)) then
          dimsf_time(1) = TSiplane%dimsf(1)*num_file
          dimsf_sub_time(1) = TSiplane%dimsf(1)
          allocate(buffer_time(TSiplane%dimsf(1), num_file), buffer_1D_time(dimsf_time(1)))
        endif
      endif


      if(ireadjp.eq.1) then
        call ReadjplaneAttr()
        call GetjplaneInfo()
        if(.not.allocated(buffer_time)) then
          dimsf_time(1) = TSjplane%dimsf(1)*num_file
          dimsf_sub_time(1) = TSjplane%dimsf(1)
          allocate(buffer_time(TSjplane%dimsf(1), num_file), buffer_1D_time(dimsf_time(1)))
        endif
      endif



      if(ireadkp.eq.1) then
 
        print *, 'number of files: ', num_file

        call InitTSkplane_HDF5()

        ! do k=1, attr_numkp(1)
        do k=num_kplane_be, num_kplane_end

          if(ibe_k_inp.lt.ibe_k_file.or.iend_k_inp.gt.iend_k_file) then
            print *, 'The input i index is not in the domain'
            stop 
          endif
          if(ibe_k_inp.gt.iend_k_file.or.iend_k_inp.lt.ibe_k_file) then
            print *, 'The input i index is not in the domain'
            stop 
          endif

          ! adjust the i index
          if(mod(ibe_k_inp - ibe_k_file, iskip_k_file).ne.0) then
            ibe_k_inp = ibe_k_inp + iskip_k_file - mod(ibe_k_inp - ibe_k_file, iskip_k_file)
          endif
          iend_k_inp = ibe_k_inp + ibuffer_k*(((iend_k_inp - ibe_k_inp)/iskip_k_file+1)/ibuffer_k) - 1

          if(mod((iend_k_inp - ibe_k_inp)/iskip_k_file+1, ibuffer_k).eq.0) then
            nn_total = ((iend_k_inp-ibe_k_inp)/iskip_k_file + 1)/ibuffer_k
          else
            nn_total = ((iend_k_inp-ibe_k_inp)/iskip_k_file + 1)/ibuffer_k + 1
          endif

 
          print *, '******************************************'
          print *, 'Number of sub files for one kplane = ', nn_total
          print *, 'DNS Output timeseries spatial index range'
          print *, 'ibe_k_file = ', ibe_k_file, 'iend_k_file = ', iend_k_file, 'iskip_k = ', iskip_k_file
          print *, 'Actual Spatial Convert range'
          print *, '     ibe_k = ', ibe_k_inp, '     iend_k = ', iend_k_inp, 'iskip_k = ', iskip_k_file
          print *, 'Converting kplane location: ', TSkplane%info_loc(k)  
          print *, '******************************************'


          do nn=1, nn_total
            ibe_k = ((nn-1)*ibuffer_k)*iskip_k_file + ibe_k_inp 
            iend_k = ibe_k + (ibuffer_k-1)*iskip_k_file

            if(nn.eq.nn_total.and.iend_k.gt.iend_k_inp) then
              iend_k = iend_k_inp
            endif

            print *, 'ibe_k = ', ibe_k
            print *, 'iend_k = ', iend_k

            TSkplane%offset_rd(1) = 0
            TSkplane%offset_rd(2) = 0
            TSkplane%offset_rd(3) = (nn-1)*ibuffer_k + ibe_k_inp - ibe_k_file
            TSkplane%offset_rd(4) = (k-1)

            write(unit=fnum_ibe,fmt='(I04.4)') ibe_k
            write(unit=fnum_iend,fmt='(I04.4)') iend_k

            write(unit=fnum1,fmt='(I04.4)')  TSkplane%info_loc(k)      !kplane(k)
            fname1 = 'timeseries_kplane'//fnum1//'_i'//fnum_ibe//'-'//fnum_iend//'.h5'

            CALL h5fcreate_f(trim(fname1), H5F_ACC_TRUNC_F, file1_id, hdferr)

            ! create 'time'
            CALL h5screate_simple_f(1, dimsf_time, dspace_id_time, hdferr)
            CALL h5dcreate_f(file1_id, trim(dname_time), H5T_NATIVE_DOUBLE, dspace_id_time, dset_id_time, hdferr)
            CALL h5dclose_f(dset_id_time, hdferr)
            CALL h5sclose_f(dspace_id_time, hdferr)
            ! create group '/kplane'
            CALL h5gcreate_f(file1_id, trim(TSkplane%gname), TSkplane%group1_id, hdferr)
            CALL h5screate_simple_f(TSkplane%rank, TSkplane%dimsm, TSkplane%memspace1, hdferr)
            ! create the dataset 'u, v, w, p, T, uk, vk, wk, pk, Tk'
            do i=1, TSkplane%num_dset
              CALL h5screate_simple_f(TSkplane%rank, TSkplane%dimsf_total, TSkplane%dspace1_id, hdferr)
              CALL h5dcreate_f(TSkplane%group1_id, trim(TSkplane%dname(i)), H5T_NATIVE_DOUBLE, &
                               TSkplane%dspace1_id, TSkplane%dset1_id, hdferr)
              CALL h5dclose_f(TSkplane%dset1_id, hdferr)
              CALL h5sclose_f(TSkplane%dspace1_id, hdferr)
            enddo
            ! read time and dataset from files
            do n=1, num_file
              write(unit=fnum2,fmt='(I08.8)') (file_be + (n-1)*file_skip)
              fname2 = trim(filepath)//'timeseries_'//fnum2//'.h5'
              CALL h5fopen_f(trim(fname2), H5F_ACC_RDONLY_F, file2_id, hdferr)
            
              ! read time
              CALL h5dopen_f(file2_id, dname_time, dset_id_time, hdferr)
              CALL h5dread_f(dset_id_time, H5T_NATIVE_DOUBLE, buffer_time(:,n), dimsf_sub_time, hdferr)
              CALL h5dclose_f(dset_id_time, hdferr)
              ! if(nn.eq.nn_total) print *, 'buffer_time = ', buffer_time

              CALL h5gopen_f(file2_id, TSkplane%gname, TSkplane%group2_id, hdferr)
              CALL h5screate_simple_f(TSkplane%rank, TSkplane%dimsm, TSkplane%memspace2, hdferr)
              do i=1, TSkplane%num_dset 
                CALL h5dopen_f(TSkplane%group2_id, trim(TSkplane%dname(i)), TSkplane%dset2_id, hdferr)
                CALL h5dget_space_f(TSkplane%dset2_id, TSkplane%dspace2_id, hdferr)
                CALL h5sselect_hyperslab_f(TSkplane%dspace2_id, H5S_SELECT_SET_F, TSkplane%offset_rd, & 
                                           TSkplane%count, hdferr, TSkplane%stride, TSkplane%block)
                CALL h5dread_f(TSkplane%dset2_id, H5T_NATIVE_DOUBLE, buffer_kplane(:,:,:,1,i), TSkplane%dimsf_sub, hdferr, &
                               file_space_id=TSkplane%dspace2_id, mem_space_id=TSkplane%memspace2) !!
                CALL h5sclose_f(TSkplane%dspace2_id, hdferr)
                CALL h5dclose_f(TSkplane%dset2_id, hdferr)
                ! write the data to the new file
                CALL h5dopen_f(TSkplane%group1_id, trim(TSkplane%dname(i)), TSkplane%dset1_id, hdferr)
                CALL h5dget_space_f(TSkplane%dset1_id, TSkplane%dspace1_id, hdferr)
                TSkplane%offset_wt(1) = TSkplane%dimsf(1)*(n-1)
                TSkplane%offset_wt(2) = 0
                TSkplane%offset_wt(3) = 0
                TSkplane%offset_wt(4) = 0
                CALL h5sselect_hyperslab_f(TSkplane%dspace1_id, H5S_SELECT_SET_F, &
                                           TSkplane%offset_wt, TSkplane%count, hdferr, TSkplane%stride, TSkplane%block)
                CALL h5dwrite_f(TSkplane%dset1_id, H5T_NATIVE_DOUBLE, buffer_kplane(:,:,:,1,i), TSkplane%dimsf_total, hdferr, &
                                file_space_id=TSkplane%dspace1_id, mem_space_id=TSkplane%memspace1)
                CALL h5dclose_f(TSkplane%dset1_id, hdferr)
                CALL h5sclose_f(TSkplane%dspace1_id, hdferr)
           
              enddo ! end i loop

              CALL h5gclose_f(TSkplane%group2_id, hdferr)
              CALL h5fclose_f(file2_id, hdferr)
            enddo ! end n loop

            ! write kplane attribute
            ! attribute 1
            attr_time(1) = TSkplane%dimsf(1)*num_file
            CALL h5screate_simple_f(1, TSkplane%attr_dims_no, attr_space, hdferr)
            CALL h5tcopy_f(H5T_NATIVE_INTEGER, atype_id, hdferr)
            CALL h5acreate_f(TSkplane%group1_id,TSkplane%attr_name(1), atype_id,attr_space,attr_id,hdferr) 
            CALL h5awrite_f(attr_id,H5T_NATIVE_INTEGER,attr_time,TSkplane%attr_dims_no,hdferr)
            CALL h5sclose_f(attr_space, hdferr)
            CALL h5aclose_f(attr_id, hdferr)  

            TSkplane%info_1st_array(1) =  ibe_k  
            TSkplane%info_1st_array(2) =  iend_k
            TSkplane%info_1st_array(3) =  iskip_k_file
 
            ! attribute 2
            CALL h5screate_simple_f(1, TSkplane%attr_dims_info, attr_space, hdferr)
            CALL h5tcopy_f(H5T_NATIVE_INTEGER, atype_id, hdferr)
            CALL h5acreate_f(TSkplane%group1_id,TSkplane%attr_name(2), atype_id,attr_space,attr_id,hdferr)
            CALL h5awrite_f(attr_id,H5T_NATIVE_INTEGER,TSkplane%info_1st_array,TSkplane%attr_dims_info,hdferr)
            CALL h5sclose_f(attr_space, hdferr)
            CALL h5aclose_f(attr_id, hdferr)  

            ! attribute 3
            CALL h5screate_simple_f(1, TSkplane%attr_dims_info, attr_space, hdferr)
            CALL h5tcopy_f(H5T_NATIVE_INTEGER, atype_id, hdferr)
            CALL h5acreate_f(TSkplane%group1_id,TSkplane%attr_name(3), atype_id,attr_space,attr_id,hdferr)
            CALL h5awrite_f(attr_id,H5T_NATIVE_INTEGER,TSkplane%info_2nd_array,TSkplane%attr_dims_info,hdferr)
            CALL h5sclose_f(attr_space, hdferr)
            CALL h5aclose_f(attr_id, hdferr)  

            ! attribute 4
            CALL h5screate_simple_f(1, TSkplane%attr_dims_no, attr_space, hdferr)
            CALL h5tcopy_f(H5T_NATIVE_INTEGER, atype_id, hdferr)
            CALL h5acreate_f(TSkplane%group1_id,TSkplane%attr_name(4), atype_id,attr_space,attr_id,hdferr)
            CALL h5awrite_f(attr_id,H5T_NATIVE_INTEGER,TSkplane%info_loc(k),TSkplane%attr_dims_no,hdferr)
            CALL h5sclose_f(attr_space, hdferr)
            CALL h5aclose_f(attr_id, hdferr) 

            CALL h5sclose_f(TSkplane%memspace1, hdferr)
            CALL h5gclose_f(TSkplane%group1_id, hdferr)

            do n=1, num_file
              do i=1, TSkplane%dimsf(1)
                buffer_1D_time(i+(n-1)*TSkplane%dimsf(1)) = buffer_time(i,n)
              enddo
            enddo

            ! write time
            CALL h5dopen_f(file1_id, trim(dname_time), dset_id_time, hdferr)
            CALL h5dget_space_f(dset_id_time, dspace_id_time, hdferr)
            CALL h5dwrite_f(dset_id_time, H5T_NATIVE_DOUBLE, buffer_1D_time, dimsf_sub_time, hdferr)
            CALL h5sclose_f(dspace_id_time, hdferr)
            CALL h5dclose_f(dset_id_time, hdferr)

            CALL h5fclose_f(file1_id, hdferr)

            ! read grid Info

!            kgrid%offset_rd(1) = 0
!            kgrid%offset_rd(2) = (nn-1)*ibuffer_k + ibe_k_inp - ibe_k_file
!            kgrid%offset_rd(3) = (k-1)
!
!            fname1 = 'timeseries_kgrid'//fnum1//'_i'//fnum_ibe//'-'//fnum_iend//'.h5'
!            fname2 = trim(filepath)//'timeseries_GridMetrics.h5'                         !!!!!!!!!
!
!            ! read grid from the file timeseries_GridMetrics.h5
!            CALL h5fopen_f(trim(fname2), H5F_ACC_RDONLY_F, file2_id, hdferr)
!            CALL h5gopen_f(file2_id, kgrid%gname, kgrid%group2_id, hdferr)
!            CALL h5screate_simple_f(kgrid%rank, kgrid%dimsm, kgrid%memspace2, hdferr)
!            do i=1, kgrid%num_dset
!              CALL h5dopen_f(kgrid%group2_id, trim(kgrid%dname(i)), kgrid%dset2_id, hdferr)
!              CALL h5dget_space_f(kgrid%dset2_id, kgrid%dspace2_id, hdferr)
!              CALL h5sselect_hyperslab_f(kgrid%dspace2_id, H5S_SELECT_SET_F, kgrid%offset_rd, &
!                                         kgrid%count, hdferr, kgrid%stride, kgrid%block)
!              CALL h5dread_f(kgrid%dset2_id, H5T_NATIVE_DOUBLE, buffer_kgrid(:,:,1,i), kgrid%dimsf_sub, hdferr, &
!                             file_space_id=kgrid%dspace2_id, mem_space_id=kgrid%memspace2) !!
!              CALL h5sclose_f(kgrid%dspace2_id, hdferr)
!              CALL h5dclose_f(kgrid%dset2_id, hdferr)
!            enddo
!            CALL h5gclose_f(kgrid%group2_id, hdferr)
!            CALL h5fclose_f(file2_id, hdferr)
!
!            ! write the grid to the new file
!            CALL h5fcreate_f(trim(fname1), H5F_ACC_TRUNC_F, file1_id, hdferr)
!            CALL h5gcreate_f(file1_id, trim(kgrid%gname), kgrid%group1_id, hdferr)
!            ! CALL h5screate_simple_f(kgrid%rank, kgird%dimsm, kgrid%memspace1, hdferr)
!            do i=1, kgrid%num_dset
!              CALL h5screate_simple_f(kgrid%rank, kgrid%dimsf_total, kgrid%dspace1_id, hdferr)
!              CALL h5dcreate_f(kgrid%group1_id, trim(kgrid%dname(i)), H5T_NATIVE_DOUBLE, &
!                               kgrid%dspace1_id, kgrid%dset1_id, hdferr)
!              CALL h5dwrite_f(kgrid%dset1_id, H5T_NATIVE_DOUBLE, buffer_kgrid(:,:,1,i), kgrid%dimsf, hdferr)
!              CALL h5dclose_f(kgrid%dset1_id, hdferr)
!              CALL h5sclose_f(kgrid%dspace1_id, hdferr)
!            enddo
!
!            ! write kgrid attribute
!            ! attribute 1: istart, iend, iskip for the global
!
!            kgrid%info_1st_array(1) = ibe_k_inp
!            kgrid%info_1st_array(2) = iend_k_inp
!            kgrid%info_1st_array(3) = iskip_k_file
!
!            CALL h5screate_simple_f(1, TSkplane%attr_dims_info, attr_space, hdferr)
!            CALL h5tcopy_f(H5T_NATIVE_INTEGER, atype_id, hdferr)
!            CALL h5acreate_f(kgrid%group1_id,kgrid%attr_name(1), atype_id,attr_space,attr_id,hdferr)
!            CALL h5awrite_f(attr_id,H5T_NATIVE_INTEGER,kgrid%info_1st_array,TSkplane%attr_dims_info,hdferr)
!            CALL h5sclose_f(attr_space, hdferr)
!            CALL h5aclose_f(attr_id, hdferr)
!
!            ! attribute 2: jstart, jend, jskip
!            CALL h5screate_simple_f(1, kgrid%attr_dims_info, attr_space, hdferr)
!            CALL h5tcopy_f(H5T_NATIVE_INTEGER, atype_id, hdferr)
!            CALL h5acreate_f(kgrid%group1_id,kgrid%attr_name(2), atype_id,attr_space,attr_id,hdferr)
!            CALL h5awrite_f(attr_id,H5T_NATIVE_INTEGER,TSkplane%info_2nd_array,kgrid%attr_dims_info,hdferr)
!            CALL h5sclose_f(attr_space, hdferr)
!            CALL h5aclose_f(attr_id, hdferr)
!
!            ! attribute 3: k locations
!            CALL h5screate_simple_f(1, kgrid%attr_dims_no, attr_space, hdferr)
!            CALL h5tcopy_f(H5T_NATIVE_INTEGER, atype_id, hdferr)
!            CALL h5acreate_f(kgrid%group1_id,kgrid%attr_name(3), atype_id,attr_space,attr_id,hdferr)
!            CALL h5awrite_f(attr_id,H5T_NATIVE_INTEGER,ibuffer_k_tmp,kgrid%attr_dims_no,hdferr)
!            CALL h5sclose_f(attr_space, hdferr)
!            CALL h5aclose_f(attr_id, hdferr)
!
!
!            CALL h5fclose_f(file1_id, hdferr)

          enddo ! end nn loop 

        enddo ! end k loop
      endif



     call FinalizeHDF5()

 contains


    subroutine GetkgridInfo()
      implicit none


 !     TSkplane%attr_name(2) = "istart, iend, iskip"
 !     TSkplane%attr_name(3) = "jstart, jend, jskip"
 !     TSkplane%attr_name(4) = "K locations"


      kgrid%gname = '/kplane'
      kgrid%attr_name(1) = "istart, iend, iskip for this conversion"
      kgrid%attr_name(2) = "jstart, jend, jskip"
      kgrid%attr_name(3) = "ibuffer_k"
      kgrid%attr_dims_no = (/1/)
      kgrid%attr_dims_info = 3
      ibuffer_k_tmp(1) = ibuffer_k

      fname2 = trim(filepath)//'timeseries_GridMetrics.h5'

      CALL h5fopen_f(trim(fname2), H5F_ACC_RDONLY_F, file2_id, hdferr)
      CALL h5gopen_f(file2_id, kgrid%gname, kgrid%group2_id, hdferr)
      CALL h5gn_members_f(kgrid%group2_id, trim(kgrid%gname), kgrid%num_dset, hdferr)
      allocate(kgrid%dname(kgrid%num_dset))
      do idx =0, kgrid%num_dset-1
        CALL h5gget_obj_info_idx_f(kgrid%group2_id, trim(kgrid%gname), idx, kgrid%dname(idx+1), H5G_DATASET_F, hdferr)
      enddo

      CALL h5dopen_f(kgrid%group2_id, trim(kgrid%dname(1)), kgrid%dset2_id, hdferr)
      CALL h5dget_space_f(kgrid%dset2_id, kgrid%dspace2_id, hdferr)
      CALL h5sget_simple_extent_ndims_f(kgrid%dspace2_id, kgrid%rank, hdferr)
      allocate(kgrid%dimsf(kgrid%rank), kgrid%dimsf_total(kgrid%rank))
      allocate(kgrid%dimsf_sub(kgrid%rank))
      CALL h5sget_simple_extent_dims_f(kgrid%dspace2_id, kgrid%dimsf, kgrid%dimsf_total, hdferr)

      ! allocate(buffer_kplane(TSkplane%dimsf(1), TSkplane%dimsf(2), ibuffer_k, 1, TSkplane%num_dset)) !!!!!!!!!!!!!!!!
      allocate(buffer_kgrid(kgrid%dimsf(1),ibuffer_k,1,kgrid%num_dset))

      CALL h5dclose_f(kgrid%dset2_id, hdferr)
      CALL h5gclose_f(kgrid%group2_id, hdferr)
      CALL h5fclose_f(file2_id, hdferr)


    end subroutine GetkgridInfo


    subroutine GetkplaneInfo()
      implicit none

      TSkplane%gname = '/kplane'

      write(unit=fnum2,fmt='(I08.8)') file_be 
      fname2 = trim(filepath)//'timeseries_'//fnum2//'.h5'

      CALL h5fopen_f(trim(fname2), H5F_ACC_RDONLY_F, file2_id, hdferr)
      CALL h5gopen_f(file2_id, TSkplane%gname, TSkplane%group2_id, hdferr)

      CALL h5gn_members_f(TSkplane%group2_id, trim(TSkplane%gname), TSkplane%num_dset, hdferr)
      allocate(TSkplane%dname(TSkplane%num_dset))
      do idx =0, TSkplane%num_dset-1
        CALL h5gget_obj_info_idx_f(TSkplane%group2_id, trim(TSkplane%gname), idx, TSkplane%dname(idx+1), H5G_DATASET_F, hdferr)
      enddo
      CALL h5dopen_f(TSkplane%group2_id, trim(TSkplane%dname(1)), TSkplane%dset2_id, hdferr)
      CALL h5dget_space_f(TSkplane%dset2_id, TSkplane%dspace2_id, hdferr)
      CALL h5sget_simple_extent_ndims_f(TSkplane%dspace2_id, TSkplane%rank, hdferr)
      allocate(TSkplane%dimsf(TSkplane%rank), TSkplane%dimsf_total(TSkplane%rank))
      allocate(TSkplane%dimsf_sub(TSkplane%rank))
      CALL h5sget_simple_extent_dims_f(TSkplane%dspace2_id, TSkplane%dimsf, TSkplane%dimsf_total, hdferr)

      ! get the # of kplane
      ! num_kplane = TSkplane%dimsf(4)
      allocate(buffer_kplane(TSkplane%dimsf(1), TSkplane%dimsf(2), ibuffer_k, 1, TSkplane%num_dset)) !!!!!!!!!!!!!!!!
   
      CALL h5dclose_f(TSkplane%dset2_id, hdferr)
      CALL h5gclose_f(TSkplane%group2_id, hdferr)
      CALL h5fclose_f(file2_id, hdferr)

    end subroutine GetkplaneInfo

    subroutine GetiplaneInfo()
      implicit none

      TSiplane%gname = '/iplane'  !!!!!

      write(unit=fnum2,fmt='(I08.8)') file_be 
      fname2 = trim(filepath)//'timeseries_'//fnum2//'.h5'

      CALL h5fopen_f(trim(fname2), H5F_ACC_RDONLY_F, file2_id, hdferr)
      CALL h5gopen_f(file2_id, TSiplane%gname, TSiplane%group2_id, hdferr)

      CALL h5gn_members_f(TSiplane%group2_id, trim(TSiplane%gname), TSiplane%num_dset, hdferr)
      allocate(TSiplane%dname(TSiplane%num_dset))
      do idx =0, TSiplane%num_dset-1
        CALL h5gget_obj_info_idx_f(TSiplane%group2_id, trim(TSiplane%gname), idx, TSiplane%dname(idx+1), H5G_DATASET_F, hdferr)
      enddo

      CALL h5dopen_f(TSiplane%group2_id, trim(TSiplane%dname(1)), TSiplane%dset2_id, hdferr)
      CALL h5dget_space_f(TSiplane%dset2_id, TSiplane%dspace2_id, hdferr)
      CALL h5sget_simple_extent_ndims_f(TSiplane%dspace2_id, TSiplane%rank, hdferr)

      allocate(TSiplane%dimsf(TSiplane%rank), TSiplane%dimsf_total(TSiplane%rank))
      allocate(TSiplane%dimsf_sub(TSiplane%rank))
      CALL h5sget_simple_extent_dims_f(TSiplane%dspace2_id, TSiplane%dimsf, TSiplane%dimsf_total, hdferr)

      allocate(buffer_iplane(TSiplane%dimsf(1), jbuffer_i, TSiplane%dimsf(3), 1, TSiplane%num_dset)) !!!!!!!!!!!!!!!!

      CALL h5dclose_f(TSiplane%dset2_id, hdferr)
      CALL h5gclose_f(TSiplane%group2_id, hdferr)
      CALL h5fclose_f(file2_id, hdferr)
      ! CALL h5sclose_f(TSiplane%dspace2_id, hdferr)

    end subroutine GetiplaneInfo



    subroutine GetjplaneInfo()
      implicit none

      TSjplane%gname = '/jplane'  !!!!!

      write(unit=fnum2,fmt='(I08.8)') file_be 
      fname2 = trim(filepath)//'timeseries_'//fnum2//'.h5'

      CALL h5fopen_f(trim(fname2), H5F_ACC_RDONLY_F, file2_id, hdferr)
      CALL h5gopen_f(file2_id, TSjplane%gname, TSjplane%group2_id, hdferr)

      CALL h5gn_members_f(TSjplane%group2_id, trim(TSjplane%gname), TSjplane%num_dset, hdferr)
      allocate(TSjplane%dname(TSjplane%num_dset))
      do idx =0, TSjplane%num_dset-1
        CALL h5gget_obj_info_idx_f(TSjplane%group2_id, trim(TSjplane%gname), idx, TSjplane%dname(idx+1), H5G_DATASET_F, hdferr)
      enddo

      CALL h5dopen_f(TSjplane%group2_id, trim(TSjplane%dname(1)), TSjplane%dset2_id, hdferr)
      CALL h5dget_space_f(TSjplane%dset2_id, TSjplane%dspace2_id, hdferr)
      CALL h5sget_simple_extent_ndims_f(TSjplane%dspace2_id, TSjplane%rank, hdferr)

      allocate(TSjplane%dimsf(TSjplane%rank), TSjplane%dimsf_total(TSjplane%rank))
      allocate(TSjplane%dimsf_sub(TSjplane%rank))
      CALL h5sget_simple_extent_dims_f(TSjplane%dspace2_id, TSjplane%dimsf, TSjplane%dimsf_total, hdferr)

      allocate(buffer_jplane(TSjplane%dimsf(1), ibuffer_j, TSjplane%dimsf(3), 1, TSjplane%num_dset)) !!!!!!!!!!!!!!!!

      CALL h5dclose_f(TSjplane%dset2_id, hdferr)
      CALL h5gclose_f(TSjplane%group2_id, hdferr)
      CALL h5fclose_f(file2_id, hdferr)
      ! CALL h5sclose_f(TSiplane%dspace2_id, hdferr)

    end subroutine GetjplaneInfo

    subroutine ReadKplaneAttr()
      implicit none

      TSkplane%gname = '/kplane'
      TSkplane%attr_name(1) = "number of time steps"
      TSkplane%attr_name(2) = "istart, iend, iskip"
      TSkplane%attr_name(3) = "jstart, jend, jskip"
      TSkplane%attr_name(4) = "K locations"
      ! TSkplane%attr_name(5) = "number of K planes"
      TSkplane%attr_dims_no = (/1/)
      TSkplane%attr_dims_info = 3

      write(unit=fnum2,fmt='(I08.8)') file_be 
      fname2 = trim(filepath)//'timeseries_'//fnum2//'.h5'
   print *, 'fname2 = ', trim(fname2)

      CALL h5fopen_f(trim(fname2), H5F_ACC_RDONLY_F, file2_id, hdferr)
      CALL h5gopen_f(file2_id, TSkplane%gname, TSkplane%group2_id, hdferr)

      ! read attribute 5
!      CALL h5aopen_f(TSkplane%group2_id, TSkplane%attr_name(5), attr_id, hdferr)
!      CALL h5aget_type_f(attr_id, atype_id, hdferr)
!      CALL h5aread_f(attr_id, atype_id, attr_numkp, TSkplane%attr_dims_no, hdferr)
!      CALL h5aclose_f(attr_id, hdferr)   

!      TSkplane%attr_dims_loc = attr_numkp(1)
!      attr_numkp(1) = 8

      attr_numkp = num_kplane_inp
      TSkplane%attr_dims_loc = attr_numkp
      allocate(TSkplane%info_loc(attr_numkp(1)))

      ! read attribute 1
!      CALL h5aopen_f(TSkplane%group2_id, TSkplane%attr_name(1), attr_id, hdferr)
!      CALL h5aget_type_f(attr_id, atype_id, hdferr)
!      CALL h5aread_f(attr_id, atype_id, attr_ntpoint, TSkplane%attr_dims_no, hdferr)
!      CALL h5aclose_f(attr_id, hdferr)     
!      print *, 'number of time steps per file = ', attr_ntpoint

      ! read attribute 2
      CALL h5aopen_f(TSkplane%group2_id, TSkplane%attr_name(2), attr_id, hdferr)
      CALL h5aget_type_f(attr_id, atype_id, hdferr)
      CALL h5aread_f(attr_id, atype_id, TSkplane%info_1st_array, TSkplane%attr_dims_info, hdferr)
      CALL h5aclose_f(attr_id, hdferr)  
      print *, 'istart, iend, iskip = ', TSkplane%info_1st_array

      ibe_k_file = TSkplane%info_1st_array(1)
      iend_k_file = TSkplane%info_1st_array(2)
      iskip_k_file = TSkplane%info_1st_array(3)

      ! read attribute 3
      CALL h5aopen_f(TSkplane%group2_id, TSkplane%attr_name(3), attr_id, hdferr)
      CALL h5aget_type_f(attr_id, atype_id, hdferr)
      CALL h5aread_f(attr_id, atype_id, TSkplane%info_2nd_array, TSkplane%attr_dims_info, hdferr)
      CALL h5aclose_f(attr_id, hdferr)
      print *, 'jstart, jend, jskip = ', TSkplane%info_2nd_array

      ! read attribute 4
      CALL h5aopen_f(TSkplane%group2_id, TSkplane%attr_name(4), attr_id, hdferr)
      CALL h5aget_type_f(attr_id, atype_id, hdferr)
      CALL h5aread_f(attr_id, atype_id, TSkplane%info_loc, TSkplane%attr_dims_loc, hdferr)
      CALL h5aclose_f(attr_id, hdferr)

      CALL h5gclose_f(TSkplane%group2_id, hdferr)
      CALL h5fclose_f(file2_id, hdferr)

      print *, 'number of K planes = ', attr_numkp
      print *, 'K locations = ', TSkplane%info_loc

    end subroutine ReadKplaneAttr



    subroutine ReadiplaneAttr()

      implicit none

      TSiplane%gname = '/iplane'
      TSiplane%attr_name(1) = "number of time steps"
      TSiplane%attr_name(2) = "jstart, jend, jskip"
      TSiplane%attr_name(3) = "kstart, kend, kskip"
      TSiplane%attr_name(4) = "I locations"
!      TSiplane%attr_name(5) = "number of I planes"
      TSiplane%attr_dims_no = (/1/)
      TSiplane%attr_dims_info = 3

      write(unit=fnum2,fmt='(I08.8)') file_be 
      fname2 = trim(filepath)//'timeseries_'//fnum2//'.h5'
      CALL h5fopen_f(trim(fname2), H5F_ACC_RDONLY_F, file2_id, hdferr)
      CALL h5gopen_f(file2_id, TSiplane%gname, TSiplane%group2_id, hdferr)

      ! read attribute 5
!      CALL h5aopen_f(TSiplane%group2_id, TSiplane%attr_name(5), attr_id, hdferr)
!      CALL h5aget_type_f(attr_id, atype_id, hdferr)
!      CALL h5aread_f(attr_id, atype_id, attr_numip, TSiplane%attr_dims_no, hdferr)
!      CALL h5aclose_f(attr_id, hdferr)   

      attr_numip = num_iplane_inp
      TSiplane%attr_dims_loc = attr_numip(1)
      allocate(TSiplane%info_loc(attr_numip(1)))

      ! read attribute 1
!      CALL h5aopen_f(TSiplane%group2_id, TSiplane%attr_name(1), attr_id, hdferr)
!      CALL h5aget_type_f(attr_id, atype_id, hdferr)
!      CALL h5aread_f(attr_id, atype_id, attr_ntpoint, TSiplane%attr_dims_no, hdferr)
!      CALL h5aclose_f(attr_id, hdferr)     
!      print *, 'number of time steps per file = ', attr_ntpoint

      ! read attribute 2
      CALL h5aopen_f(TSiplane%group2_id, TSiplane%attr_name(2), attr_id, hdferr)
      CALL h5aget_type_f(attr_id, atype_id, hdferr)
      CALL h5aread_f(attr_id, atype_id, TSiplane%info_1st_array, TSiplane%attr_dims_info, hdferr)
      CALL h5aclose_f(attr_id, hdferr)  
      print *, 'jstart, jend, jskip = ', TSiplane%info_1st_array

      jbe_i_file = TSiplane%info_1st_array(1)
      jend_i_file = TSiplane%info_1st_array(2)
      jskip_i_file = TSiplane%info_1st_array(3) 

      ! read attribute 3
      CALL h5aopen_f(TSiplane%group2_id, TSiplane%attr_name(3), attr_id, hdferr)
      CALL h5aget_type_f(attr_id, atype_id, hdferr)
      CALL h5aread_f(attr_id, atype_id, TSiplane%info_2nd_array, TSiplane%attr_dims_info, hdferr)
      CALL h5aclose_f(attr_id, hdferr)
      print *, 'kstart, kend, kskip = ', TSiplane%info_2nd_array

      ! read attribute 4
      CALL h5aopen_f(TSiplane%group2_id, TSiplane%attr_name(4), attr_id, hdferr)
      CALL h5aget_type_f(attr_id, atype_id, hdferr)
      CALL h5aread_f(attr_id, atype_id, TSiplane%info_loc, TSiplane%attr_dims_loc, hdferr)
      CALL h5aclose_f(attr_id, hdferr)

      CALL h5gclose_f(TSiplane%group2_id, hdferr)
      CALL h5fclose_f(file2_id, hdferr)

      print *, 'number of I planes = ', attr_numip
      print *, 'I locations = ', TSiplane%info_loc

    end subroutine ReadiplaneAttr




    subroutine ReadjplaneAttr()
      implicit none

      TSjplane%gname = '/jplane'
      TSjplane%attr_name(1) = "number of time steps"
      TSjplane%attr_name(2) = "istart, iend, iskip"
      TSjplane%attr_name(3) = "kstart, kend, kskip"
      TSjplane%attr_name(4) = "J locations"
!      TSjplane%attr_name(5) = "number of J planes"
      TSjplane%attr_dims_no = (/1/)
      TSjplane%attr_dims_info = 3

      write(unit=fnum2,fmt='(I08.8)') file_be 
      fname2 = trim(filepath)//'timeseries_'//fnum2//'.h5'
      CALL h5fopen_f(trim(fname2), H5F_ACC_RDONLY_F, file2_id, hdferr)
      CALL h5gopen_f(file2_id, TSjplane%gname, TSjplane%group2_id, hdferr)

      ! read attribute 5
!      CALL h5aopen_f(TSjplane%group2_id, TSjplane%attr_name(5), attr_id, hdferr)
!      CALL h5aget_type_f(attr_id, atype_id, hdferr)
!      CALL h5aread_f(attr_id, atype_id, attr_numjp, TSjplane%attr_dims_no, hdferr)
!      CALL h5aclose_f(attr_id, hdferr)   

      attr_numjp = num_jplane_inp
      TSjplane%attr_dims_loc = attr_numjp(1)
      allocate(TSjplane%info_loc(attr_numjp(1)))

      ! read attribute 1
!      CALL h5aopen_f(TSjplane%group2_id, TSjplane%attr_name(1), attr_id, hdferr)
!      CALL h5aget_type_f(attr_id, atype_id, hdferr)
!      CALL h5aread_f(attr_id, atype_id, attr_ntpoint, TSjplane%attr_dims_no, hdferr)
!      CALL h5aclose_f(attr_id, hdferr)     
!      print *, 'number of time steps per file = ', attr_ntpoint

      ! read attribute 2
      CALL h5aopen_f(TSjplane%group2_id, TSjplane%attr_name(2), attr_id, hdferr)
      CALL h5aget_type_f(attr_id, atype_id, hdferr)
      CALL h5aread_f(attr_id, atype_id, TSjplane%info_1st_array, TSjplane%attr_dims_info, hdferr)
      CALL h5aclose_f(attr_id, hdferr)  
      print *, 'istart, iend, iskip = ', TSjplane%info_1st_array

      ibe_j_file = TSjplane%info_1st_array(1)
      iend_j_file = TSjplane%info_1st_array(2)
      iskip_j_file = TSjplane%info_1st_array(3)

      ! read attribute 3
      CALL h5aopen_f(TSjplane%group2_id, TSjplane%attr_name(3), attr_id, hdferr)
      CALL h5aget_type_f(attr_id, atype_id, hdferr)
      CALL h5aread_f(attr_id, atype_id, TSjplane%info_2nd_array, TSjplane%attr_dims_info, hdferr)
      CALL h5aclose_f(attr_id, hdferr)
      print *, 'kstart, kend, kskip = ', TSjplane%info_2nd_array


      ! read attribute 4
      CALL h5aopen_f(TSjplane%group2_id, TSjplane%attr_name(4), attr_id, hdferr)
      CALL h5aget_type_f(attr_id, atype_id, hdferr)
      CALL h5aread_f(attr_id, atype_id, TSjplane%info_loc, TSjplane%attr_dims_loc, hdferr)
      CALL h5aclose_f(attr_id, hdferr)

      CALL h5gclose_f(TSjplane%group2_id, hdferr)
      CALL h5fclose_f(file2_id, hdferr)

      print *, 'number of J planes = ', attr_numjp
      print *, 'J locations = ', TSjplane%info_loc

    end subroutine ReadjplaneAttr




    subroutine InitTSkplane_HDF5()
      implicit none

!      allocate(TSkplane%dimsf_sub(TSkplane%rank))
!      allocate(TSkplane%dimsf_total(TSkplane%rank))
      allocate(TSkplane%dimsm(TSkplane%rank))
      allocate(TSkplane%stride(TSkplane%rank))
      allocate(TSkplane%count(TSkplane%rank))
      allocate(TSkplane%block(TSkplane%rank))
      allocate(TSkplane%offset_rd(TSkplane%rank))
      allocate(TSkplane%offset_wt(TSkplane%rank))

      TSkplane%dimsf_sub(1) = TSkplane%dimsf(1)
      TSkplane%dimsf_sub(2) = TSkplane%dimsf(2)
      TSkplane%dimsf_sub(3) = ibuffer_k
      TSkplane%dimsf_sub(4) = 1

      TSkplane%dimsf_total(1) = TSkplane%dimsf(1)*num_file
      TSkplane%dimsf_total(2) = TSkplane%dimsf(2)
      TSkplane%dimsf_total(3) = ibuffer_k
      TSkplane%dimsf_total(4) = 1
      TSkplane%dimsm(1) = TSkplane%dimsf(1)
      TSkplane%dimsm(2) = TSkplane%dimsf(2)
      TSkplane%dimsm(3) = ibuffer_k
      TSkplane%dimsm(4) = 1
      TSkplane%stride(1) = 1
      TSkplane%stride(2) = 1
      TSkplane%stride(3) = 1
      TSkplane%stride(4) = 1
      TSkplane%count(1) = 1
      TSkplane%count(2) = 1
      TSkplane%count(3) = 1
      TSkplane%count(4) = 1
      TSkplane%block(1) = TSkplane%dimsf(1)
      TSkplane%block(2) = TSkplane%dimsf(2)
      TSkplane%block(3) = ibuffer_k
      TSkplane%block(4) = 1

      ! kgrid info
!      allocate(kgrid%dimsf_sub(kgrid%rank))
!      allocate(kgrid%dimsf_total(kgrid%rank))
      allocate(kgrid%dimsm(kgrid%rank))
      allocate(kgrid%stride(kgrid%rank))
      allocate(kgrid%count(kgrid%rank))
      allocate(kgrid%block(kgrid%rank))
      allocate(kgrid%offset_rd(kgrid%rank))
      allocate(kgrid%offset_wt(kgrid%rank))

      kgrid%dimsf_sub(1) = kgrid%dimsf(1)
      kgrid%dimsf_sub(2) = ibuffer_k
      kgrid%dimsf_sub(3) = 1

      kgrid%dimsf_total(1) = kgrid%dimsf(1)
      kgrid%dimsf_total(2) = ibuffer_k
      kgrid%dimsf_total(3) = 1
      kgrid%dimsm(1) = kgrid%dimsf(1)
      kgrid%dimsm(2) = ibuffer_k
      kgrid%dimsm(3) = 1
      kgrid%stride(1) = 1
      kgrid%stride(2) = 1
      kgrid%stride(3) = 1
      kgrid%count(1) = 1
      kgrid%count(2) = 1
      kgrid%count(3) = 1
      kgrid%block(1) = kgrid%dimsf(1)
      kgrid%block(2) = ibuffer_k
      kgrid%block(3) = 1

    end subroutine InitTSkplane_HDF5



    subroutine InitTSiplane_HDF5()
      implicit none
      integer :: nv = 10

      TSiplane%dimsf_sub(1) = TSiplane%dimsf(1)
      TSiplane%dimsf_sub(2) = jbuffer_i
      TSiplane%dimsf_sub(3) = TSiplane%dimsf(3)
      TSiplane%dimsf_sub(4) = 1
 
      TSiplane%dimsf_total(1) = TSiplane%dimsf(1)*num_file
      TSiplane%dimsf_total(2) = jbuffer_i
      TSiplane%dimsf_total(3) = TSiplane%dimsf(3)
      TSiplane%dimsf_total(4) = 1
      TSiplane%dimsm(1) = TSiplane%dimsf(1)
      TSiplane%dimsm(2) = jbuffer_i
      TSiplane%dimsm(3) = TSiplane%dimsf(3)
      TSiplane%dimsm(4) = 1
      TSiplane%stride(1) = 1
      TSiplane%stride(2) = 1
      TSiplane%stride(3) = 1
      TSiplane%stride(4) = 1
      TSiplane%count(1) = 1
      TSiplane%count(2) = 1
      TSiplane%count(3) = 1
      TSiplane%count(4) = 1
      TSiplane%block(1) = TSiplane%dimsf(1)
      TSiplane%block(2) = jbuffer_i
      TSiplane%block(3) = TSiplane%dimsf(3)
      TSiplane%block(4) = 1

    end subroutine InitTSiplane_HDF5



    subroutine InitTSjplane_HDF5()
      implicit none

      TSjplane%dimsf_sub(1) = TSjplane%dimsf(1)
      TSjplane%dimsf_sub(2) = ibuffer_j
      TSjplane%dimsf_sub(3) = TSjplane%dimsf(3)
      TSjplane%dimsf_sub(4) = 1

      TSjplane%dimsf_total(1) = TSjplane%dimsf(1)*num_file
      TSjplane%dimsf_total(2) = ibuffer_j
      TSjplane%dimsf_total(3) = TSjplane%dimsf(3)
      TSjplane%dimsf_total(4) = 1
      TSjplane%dimsm(1) = TSjplane%dimsf(1)
      TSjplane%dimsm(2) = ibuffer_j
      TSjplane%dimsm(3) = TSjplane%dimsf(3)
      TSjplane%dimsm(4) = 1
      TSjplane%stride(1) = 1
      TSjplane%stride(2) = 1
      TSjplane%stride(3) = 1
      TSjplane%stride(4) = 1
      TSjplane%count(1) = 1
      TSjplane%count(2) = 1
      TSjplane%count(3) = 1
      TSjplane%count(4) = 1
      TSjplane%block(1) = TSjplane%dimsf(1)
      TSjplane%block(2) = ibuffer_j
      TSjplane%block(3) = TSjplane%dimsf(3)
      TSjplane%block(4) = 1

    end subroutine InitTSjplane_HDF5




    subroutine InitHDF5()
      ! Initialize HDF5 library and Fortran interfaces.
      CALL h5open_f(hdferr)
    end subroutine InitHDF5

    subroutine FinalizeHDF5()
      ! Close FORTRAN interfaces and HDF5 library.
      CALL h5close_f(hdferr)
    end subroutine FinalizeHDF5


    subroutine Init()
      implicit none
      integer, parameter :: nid=5
!      integer, parameter :: nid=10
      integer :: i,j
!      if (nid.eq.5) then
!        print*,'please input parameters or redirect from a file'
!      else
!        open(nid,file='convert.inp',status='old')
!      end if

      read(nid,*)
      read(nid,'(a)')filepath
      read(nid,*)
      read(nid,*)
      read(nid,*) file_be, file_end, file_skip
      read(nid,*)
      read(nid,*) ireadkp, ibe_k_inp, iend_k_inp, ibuffer_k, num_kplane_inp, num_kplane_be, num_kplane_end
      read(nid,*)
      read(nid,*) ireadip, jbe_i_inp, jend_i_inp, jbuffer_i, num_iplane_inp, num_iplane_be, num_iplane_end
      read(nid,*)
      read(nid,*) ireadjp, ibe_j_inp, iend_j_inp, ibuffer_j, num_jplane_inp, num_jplane_be, num_jplane_end

      if(ireadkp.eq.1.and.(iend_k_inp-ibe_k_inp + 1 )/ibuffer_k.eq.0) then
        print *, '(iend_k_inp - ibe_k_inp) should be larger or equal ibuffer_k'
        stop
      endif
      if(ireadip.eq.1.and.(jend_i_inp-jbe_i_inp + 1 )/jbuffer_i.eq.0) then
        print *, '(jend_i_inp - jbe_i_inp) should be larger or equal jbuffer_i'
        stop
      endif
      if(ireadjp.eq.1.and.(iend_j_inp-ibe_j_inp + 1 )/ibuffer_j.eq.0) then
        print *, '(iend_j_inp - ibe_j_inp) should be larger or equal ibuffer_j'
        stop
      endif
      if(ireadkp.eq.1.and.(num_kplane_end - num_kplane_be).lt.0) then
        print *, 'num_kplane_end should be larger or equal num_kplane_be'
        stop
      endif
      if(ireadkp.eq.1.and.num_kplane_end.gt.num_kplane_inp) then
        print *, 'kplane_be and kplane_end should be smaller than num_kplane_inp'
        stop
      endif
      if(ireadip.eq.1.and.(num_iplane_end - num_iplane_be).lt.0) then
        print *, 'num_iplane_end should be larger or equal num_iplane_be'
        stop
      endif
      if(ireadip.eq.1.and.num_iplane_end.gt.num_iplane_inp) then
        print *, 'iplane_be and iplane_end should be smaller than num_iplane_inp'
        stop
      endif
      if(ireadjp.eq.1.and.(num_jplane_end - num_jplane_be).lt.0) then
        print *, 'num_jplane_end should be larger or equal num_jplane_be'
        stop
      endif
      if(ireadjp.eq.1.and.num_jplane_end.gt.num_jplane_inp) then
        print *, 'jplane_be and jplane_end should be smaller than num_jplane_inp'
        stop
      endif

    end subroutine Init


  End Program convert_timeseries
