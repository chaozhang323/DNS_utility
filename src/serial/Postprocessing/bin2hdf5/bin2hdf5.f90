!> convert old format to the HDF5 format with buffer size  (timeseries data)


program convert
  use MTSTimeData
  use MTSKplane
  use MTSIplane
  use MTSJplane
!  use MDerivative
  use HDF5
  implicit none
  integer, parameter :: nvar_output = 10
  integer :: iascii
  integer :: ireadkp, ireadip, ireadjp, nodes_k, nodes_i, nodes_j
  integer :: nskip, ixwindowl, ixwindowr
  integer :: num_kplane, iwindow
  integer :: num_iplane, num_jplane
  integer :: nxpoint, nypoint, nzpoint, imints, imaxts
  integer, dimension(:), allocatable :: kplane, iplane, jplane
  real(8), dimension(:,:,:,:), allocatable :: tsdata
  integer :: varindex_output(nvar_output)
  character(3) :: varname_output(nvar_output)
  integer :: ibe_k, iend_k, jbe_k, jend_k
  integer :: jbe_i, jend_i, kbe_i, kend_i
  integer :: ibe_j, iend_j, kbe_j, kend_j
  character(200) :: filepath
  character(200) :: fileinput, fileoutput
  character(4) :: fnum, fnum1, fnum2
  integer :: i, j, k, n, numpt
  real(8) :: dx_kp, dy_kp, dx_ip, dy_ip, dx_jp, dy_jp
  real(8) :: zloc_kp, dzdk_kp
  real(8),pointer :: zloc_ip(:), dzdk_ip(:), zloc_jp(:), dzdk_jp(:)
  logical :: IsZgridRead = .false.
  integer :: klen, idim
  real(8) :: dt_sample
  integer :: numpt1, numpt2
  real(8) :: uvave1, uvave2
  real(8) :: tmp(-2:2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  type hyperslab
    character(100) :: gname
    character(100), dimension(:), pointer :: dname
    character(100), dimension(4) :: attr_name
    integer(HID_T) :: group_id, memspace
    integer(HID_T) :: dspace_id, dset_id
    integer(HSIZE_T), dimension(:), pointer :: dimsf
    integer :: rank
    integer, dimension(3) :: info_1st_array, info_2nd_array
    integer(HSIZE_T), dimension(1) :: attr_dims_info, attr_dims_loc, attr_dims_no = (/1/)
    integer, dimension(:), pointer :: info_loc
  end type hyperslab
  type(hyperslab) :: TSkplane, TSiplane, TSjplane 
  integer(HID_T) :: file_id, attr_id, attr_space, atype_id
  integer(HSIZE_T), dimension(1) :: dimsf_time
  integer(HID_T) :: dspace_id_time, dset_id_time
  character(100) :: dname_time
  integer, dimension(1) :: attr_time
  real(8), dimension(:,:,:,:), allocatable :: buffer_kplane, buffer_iplane, buffer_jplane
  integer :: hdferr
  character(100) :: fname
  integer :: nn
  integer :: ibe_k_inp, iend_k_inp, iskip_k_inp, jbe_k_inp, jend_k_inp, jskip_k_inp
  integer :: jbe_i_inp, jend_i_inp, jskip_i_inp, kbe_i_inp, kend_i_inp, kskip_i_inp
  integer :: ibe_j_inp, iend_j_inp, iskip_j_inp, kbe_j_inp, kend_j_inp, kskip_j_inp
  integer :: nnn, nnn_total, ibuffer_k, jbuffer_i, ibuffer_j
  integer :: convert_file_be_k, convert_file_be_i, convert_file_be_j

!  (/'u','v','w','p','t','ux','uy','uz','vx','vy','vz','wx','wy','wz',&
!    'uxx','uxy','uxz','vxy','vyy','vyz','wxz','wyz','wzz'/))
  
  call Input()
  call ReadTSTime(filepath,ntpoint,nskip)
  dt_sample = Getdts()
  print *, 'ntpoint =', ntpoint
  varindex_output = (/( i,i=1,nvar_output) /)
  print *, 'varindex_output =', varindex_output
  do n = 1, nvar_output
     varname_output(n) = varname_Kplane(varindex_output(n)) ! Need to be updated
  enddo
  print *, 'Variable name read in:   ', (varname_output(n),n=1,nvar_output)

  attr_time = ntpoint
  
  call InitHDF5()

  if(ireadkp.eq.1.and.num_kplane.gt.0) then  ! Using constant-k plane time series

    call GetDNSInfo_kplane(nodes_k,iascii, nskip, dx_kp, dy_kp, ibe_k_inp, iend_k_inp, iskip_k_inp, &
                            jbe_k_inp, jend_k_inp, jskip_k_inp )

    do k = 1, num_kplane

      if(mod((iend_k_inp-ibe_k_inp)/iskip_k_inp + 1, ibuffer_k).eq.0) then
        nnn_total = ((iend_k_inp-ibe_k_inp)/iskip_k_inp + 1)/ibuffer_k
      else
        nnn_total = ((iend_k_inp-ibe_k_inp)/iskip_k_inp + 1)/ibuffer_k + 1
      endif

      do nnn=convert_file_be_k, nnn_total

        ibe_k = ((nnn-1)*ibuffer_k)*iskip_k_inp + ibe_k_inp
        iend_k = ibe_k + (ibuffer_k - 1)*iskip_k_inp
        if(nnn.eq.nnn_total.and.iend_k.gt.((iend_k_inp-ibe_k_inp)/iskip_k_inp + 1)) then  !!!!!!!!!!!!!
          iend_k = iend_k_inp
        endif


        call InitTSkplane(nxpoint,nypoint, ibe_k,iend_k)
        allocate(tsdata(nvar_output,nxpoint,nypoint,ntpoint))

        call InitTSkplane_HDF5()


        call ReadTSKplane(kplane(k),nvar_output,varindex_output,tsdata, zloc_kp, dzdk_kp)

        write(unit=fnum,fmt='(I04.4)') kplane(k)
        write(unit=fnum1,fmt='(I04.4)') ibe_k
        write(unit=fnum2,fmt='(I04.4)') iend_k
        fname = 'timeseries_kplane'//fnum//'_i'//fnum1//'-'//fnum2//'.h5'
        CALL h5fcreate_f(trim(fname), H5F_ACC_TRUNC_F, file_id, hdferr)
       ! write time
        CALL h5screate_simple_f(1,dimsf_time,dspace_id_time, hdferr)
        CALL h5dcreate_f(file_id,trim(dname_time),H5T_NATIVE_DOUBLE,dspace_id_time,dset_id_time,hdferr)
        CALL h5dwrite_f(dset_id_time,H5T_NATIVE_DOUBLE,tstime,dimsf_time,hdferr)
        CALL h5sclose_f(dspace_id_time,hdferr)
        CALL h5dclose_f(dset_id_time,hdferr)
        CALL h5gcreate_f(file_id,trim(TSkplane%gname),TSkplane%group_id, hdferr)
       ! write (u,v,w,p,t,uk,vk,wk,pk,tk)
        do nn=1, 10
          do n=1, ntpoint
            do j=1, nypoint
              do i=1, nxpoint
              buffer_kplane(n,j,i,1) = tsdata(nn,i,j,n)
              enddo
            enddo
          enddo

          CALL h5screate_simple_f(TSkplane%rank, TSkplane%dimsf, TSkplane%dspace_id, hdferr)
          CALL h5dcreate_f(TSkplane%group_id, trim(TSkplane%dname(nn)), &
                           H5T_NATIVE_DOUBLE, TSkplane%dspace_id, &
                           TSkplane%dset_id, hdferr)
          CALL h5dwrite_f(TSkplane%dset_id, H5T_NATIVE_DOUBLE, buffer_kplane, TSkplane%dimsf, hdferr)
          CALL h5sclose_f(TSkplane%dspace_id, hdferr)
          CALL h5dclose_f(TSkplane%dset_id, hdferr)
        enddo

       ! write kplane attribute
       ! attribute 1
        CALL h5screate_simple_f(1, TSkplane%attr_dims_no, attr_space, hdferr)
        CALL h5tcopy_f(H5T_NATIVE_INTEGER, atype_id, hdferr)
        CALL h5acreate_f(TSkplane%group_id,TSkplane%attr_name(1), atype_id,attr_space,attr_id,hdferr)
      
        CALL h5awrite_f(attr_id,H5T_NATIVE_INTEGER,attr_time,TSkplane%attr_dims_no,hdferr)
        CALL h5sclose_f(attr_space, hdferr)
        CALL h5aclose_f(attr_id, hdferr)

        TSkplane%info_1st_array(1) =  ibe_k
        TSkplane%info_1st_array(2) =  iend_k
        TSkplane%info_1st_array(3) =  iskip_k_inp
 
       ! attribute 2
        CALL h5screate_simple_f(1, TSkplane%attr_dims_info, attr_space, hdferr)
        CALL h5tcopy_f(H5T_NATIVE_INTEGER, atype_id, hdferr)
        CALL h5acreate_f(TSkplane%group_id,TSkplane%attr_name(2), atype_id,attr_space,attr_id,hdferr)
        CALL h5awrite_f(attr_id,H5T_NATIVE_INTEGER,TSkplane%info_1st_array,TSkplane%attr_dims_info,hdferr)
        CALL h5sclose_f(attr_space, hdferr)
        CALL h5aclose_f(attr_id, hdferr)

       ! attribute 3
        CALL h5screate_simple_f(1, TSkplane%attr_dims_info, attr_space, hdferr)
        CALL h5tcopy_f(H5T_NATIVE_INTEGER, atype_id, hdferr)
        CALL h5acreate_f(TSkplane%group_id,TSkplane%attr_name(3), atype_id,attr_space,attr_id,hdferr)
        CALL h5awrite_f(attr_id,H5T_NATIVE_INTEGER,TSkplane%info_2nd_array,TSkplane%attr_dims_info,hdferr)
        CALL h5sclose_f(attr_space, hdferr)
        CALL h5aclose_f(attr_id, hdferr)

       ! attribute 4
        CALL h5screate_simple_f(1, TSkplane%attr_dims_no, attr_space, hdferr)
        CALL h5tcopy_f(H5T_NATIVE_INTEGER, atype_id, hdferr)
        CALL h5acreate_f(TSkplane%group_id,TSkplane%attr_name(4), atype_id,attr_space,attr_id,hdferr)
        CALL h5awrite_f(attr_id,H5T_NATIVE_INTEGER,TSkplane%info_loc(k),TSkplane%attr_dims_no,hdferr)
        CALL h5sclose_f(attr_space, hdferr)
        CALL h5aclose_f(attr_id, hdferr)

        CALL h5gclose_f(TSkplane%group_id, hdferr)
        CALL h5fclose_f(file_id, hdferr)
        deallocate(tsdata)
      enddo ! end nnn loop
    enddo  ! end looping over k planes


  endif  ! end calculating const-k plane

  if(ireadip.eq.1.and.num_iplane.gt.0) then  ! Using constant-i plane time series
    ! Reading i plane data
    if(allocated(tsdata)) deallocate(tsdata)



    do i = 1, num_iplane
      call GetDNSInfo_iplane(iplane(i),nodes_i,iascii,nskip,dx_ip,dy_ip,jbe_i_inp,jend_i_inp,jskip_i_inp, &
                             kbe_i_inp, kend_i_inp, kskip_i_inp)
     
       if(mod((jend_i_inp-jbe_i_inp)/jskip_i_inp + 1, jbuffer_i).eq.0) then
         nnn_total = ((jend_i_inp-jbe_i_inp)/jskip_i_inp + 1)/jbuffer_i
       else
         nnn_total = ((jend_i_inp-jbe_i_inp)/jskip_i_inp + 1)/jbuffer_i + 1
       endif

    do nnn=convert_file_be_i, nnn_total
      jbe_i = ((nnn-1)*jbuffer_i)*jskip_i_inp + jbe_i_inp
      jend_i = jbe_i + (jbuffer_i - 1)*jskip_i_inp
      if(nnn.eq.nnn_total.and.jend_i.gt.((jend_i_inp-jbe_i_inp)/jskip_i_inp + 1)) then
        jend_i = jend_i_inp
      endif

!      call GetDNSInfo_iplane(iplane(i),nodes_i,iascii,nskip,dx_ip,dy_ip,jbe_i,jend_i, &
!                            kbe_i_inp, kend_i_inp)


       call InitTSIplane(nypoint,nzpoint, jbe_i, jend_i)
       allocate(tsdata(nvar_output,nypoint,nzpoint,ntpoint))
       call ReadTSIplane(iplane(i), nvar_output, varindex_output,tsdata)
       write(unit=fnum,fmt='(I04.4)') iplane(i)
       write(unit=fnum1,fmt='(I04.4)') jbe_i
       write(unit=fnum2,fmt='(I04.4)') jend_i
       fname = 'timeseries_iplane'//fnum//'_j'//fnum1//'-'//fnum2//'.h5'
       call InitTSiplane_HDF5()

       CALL h5fcreate_f(trim(fname), H5F_ACC_TRUNC_F, file_id, hdferr)
       ! write time
       CALL h5screate_simple_f(1,dimsf_time,dspace_id_time, hdferr)
       CALL h5dcreate_f(file_id,trim(dname_time),H5T_NATIVE_DOUBLE,dspace_id_time,dset_id_time,hdferr)
       CALL h5dwrite_f(dset_id_time,H5T_NATIVE_DOUBLE,tstime,dimsf_time,hdferr)
       CALL h5sclose_f(dspace_id_time,hdferr)
       CALL h5dclose_f(dset_id_time,hdferr)

       CALL h5gcreate_f(file_id,trim(TSiplane%gname),TSiplane%group_id, hdferr)

       ! write (u,v,w,p,t,uk,vk,wk,pk,tk)

       do nn=1, 10
         do n=1, ntpoint
           do j=1, nypoint
             do k=1, nzpoint
             buffer_iplane(n,j,k,1) = tsdata(nn,j,k,n)
             enddo
           enddo
         enddo
         CALL h5screate_simple_f(TSiplane%rank, TSiplane%dimsf, TSiplane%dspace_id, hdferr)
         CALL h5dcreate_f(TSiplane%group_id, trim(TSiplane%dname(nn)), &
                          H5T_NATIVE_DOUBLE, TSiplane%dspace_id, &
                          TSiplane%dset_id, hdferr)
         CALL h5dwrite_f(TSiplane%dset_id, H5T_NATIVE_DOUBLE, buffer_iplane, TSiplane%dimsf, hdferr)
         CALL h5sclose_f(TSiplane%dspace_id, hdferr)
         CALL h5dclose_f(TSiplane%dset_id, hdferr)
       enddo

       ! write iplane attribute
       ! attribute 1
       CALL h5screate_simple_f(1, TSiplane%attr_dims_no, attr_space, hdferr)
       CALL h5tcopy_f(H5T_NATIVE_INTEGER, atype_id, hdferr)
       CALL h5acreate_f(TSiplane%group_id,TSiplane%attr_name(1), atype_id,attr_space,attr_id,hdferr)
       CALL h5awrite_f(attr_id,H5T_NATIVE_INTEGER,attr_time,TSiplane%attr_dims_no,hdferr)
       CALL h5sclose_f(attr_space, hdferr)
       CALL h5aclose_f(attr_id, hdferr)  

      TSiplane%info_1st_array(1) =  jbe_i
      TSiplane%info_1st_array(2) =  jend_i
      TSiplane%info_1st_array(3) =  jskip_i_inp

       ! attribute 2
       CALL h5screate_simple_f(1, TSiplane%attr_dims_info, attr_space, hdferr)
       CALL h5tcopy_f(H5T_NATIVE_INTEGER, atype_id, hdferr)
       CALL h5acreate_f(TSiplane%group_id,TSiplane%attr_name(2), atype_id,attr_space,attr_id,hdferr)
       CALL h5awrite_f(attr_id,H5T_NATIVE_INTEGER,TSiplane%info_1st_array,TSiplane%attr_dims_info,hdferr)
       CALL h5sclose_f(attr_space, hdferr)
       CALL h5aclose_f(attr_id, hdferr)  

       ! attribute 3
       CALL h5screate_simple_f(1, TSiplane%attr_dims_info, attr_space, hdferr)
       CALL h5tcopy_f(H5T_NATIVE_INTEGER, atype_id, hdferr)
       CALL h5acreate_f(TSiplane%group_id,TSiplane%attr_name(3), atype_id,attr_space,attr_id,hdferr)
       CALL h5awrite_f(attr_id,H5T_NATIVE_INTEGER,TSiplane%info_2nd_array,TSiplane%attr_dims_info,hdferr)
       CALL h5sclose_f(attr_space, hdferr)
       CALL h5aclose_f(attr_id, hdferr)  

       ! attribute 4
       CALL h5screate_simple_f(1, TSiplane%attr_dims_no, attr_space, hdferr)
       CALL h5tcopy_f(H5T_NATIVE_INTEGER, atype_id, hdferr)
       CALL h5acreate_f(TSiplane%group_id,TSiplane%attr_name(4), atype_id,attr_space,attr_id,hdferr)
       CALL h5awrite_f(attr_id,H5T_NATIVE_INTEGER,TSiplane%info_loc(i),TSiplane%attr_dims_no,hdferr)
       CALL h5sclose_f(attr_space, hdferr)
       CALL h5aclose_f(attr_id, hdferr) 

       CALL h5gclose_f(TSiplane%group_id, hdferr)
       CALL h5fclose_f(file_id, hdferr)
     
      deallocate(buffer_iplane) 
      deallocate(TSiplane%dimsf)
      deallocate(TSiplane%dname)  
      deallocate(tsdata) 
    enddo ! end nnn loop
    enddo  ! end looping over i planes

  endif  ! end calculating const-i plane


  if(ireadjp.eq.1.and.num_jplane.gt.0) then  ! Using constant-j plane time series
    if(allocated(tsdata)) deallocate(tsdata)
    ! Reading j plane data
    do j = 1, num_jplane
       call GetDNSInfo_jplane(jplane(j),nodes_j,iascii,nskip,dx_jp,dy_jp,ibe_j_inp,iend_j_inp,iskip_j_inp, &
                              kbe_j_inp, kend_j_inp, kskip_j_inp)

       if(mod((iend_j_inp-ibe_j_inp)/iskip_j_inp + 1, ibuffer_j).eq.0) then
         nnn_total = ((iend_j_inp-ibe_j_inp)/iskip_j_inp + 1)/ibuffer_j
       else
         nnn_total = ((iend_j_inp-ibe_j_inp)/iskip_j_inp + 1)/ibuffer_j + 1
       endif


    do nnn=convert_file_be_j, nnn_total
       ibe_j = ((nnn-1)*ibuffer_j)*iskip_j_inp + ibe_j_inp
       iend_j = ibe_j + (ibuffer_j - 1 )*iskip_j_inp
         if(nnn.eq.nnn_total.and.iend_j.gt.((iend_j_inp-ibe_j_inp)/iskip_j_inp + 1)) then
           iend_j = iend_j_inp
         endif

!       call GetDNSInfo_jplane(jplane(j),nodes_j,iascii,nskip,dx_jp,dy_jp,ibe_j,iend_j, &
!                              kbe_j_inp, kend_j_inp)

       call InitTSJplane(nxpoint,nzpoint, ibe_j, iend_j)
       allocate(tsdata(nvar_output,nxpoint,nzpoint,ntpoint))
       idim = nxpoint
       call ReadTSJplane(jplane(j), nvar_output, varindex_output,tsdata)

       write(unit=fnum,fmt='(I04.4)') jplane(j)
       write(unit=fnum1,fmt='(I04.4)') ibe_j
       write(unit=fnum2,fmt='(I04.4)') iend_j

       fname = 'timeseries_jplane'//fnum//'_i'//fnum1//'-'//fnum2//'.h5'

       call InitTSjplane_HDF5()

       CALL h5fcreate_f(trim(fname), H5F_ACC_TRUNC_F, file_id, hdferr)
       ! write time
       CALL h5screate_simple_f(1,dimsf_time,dspace_id_time, hdferr)
       CALL h5dcreate_f(file_id,trim(dname_time),H5T_NATIVE_DOUBLE,dspace_id_time,dset_id_time,hdferr)
       CALL h5dwrite_f(dset_id_time,H5T_NATIVE_DOUBLE,tstime,dimsf_time,hdferr)
       CALL h5sclose_f(dspace_id_time,hdferr)
       CALL h5dclose_f(dset_id_time,hdferr)

       CALL h5gcreate_f(file_id,trim(TSjplane%gname),TSjplane%group_id, hdferr)

       ! write (u,v,w,p,t,uk,vk,wk,pk,tk)

       do nn=1, 10
         do n=1, ntpoint
           do i=1, nxpoint
             do k=1, nzpoint
             buffer_jplane(n,i,k,1) = tsdata(nn,i,k,n)
             enddo
           enddo
         enddo
         CALL h5screate_simple_f(TSjplane%rank, TSjplane%dimsf, TSjplane%dspace_id, hdferr)
         CALL h5dcreate_f(TSjplane%group_id, trim(TSjplane%dname(nn)), &
                          H5T_NATIVE_DOUBLE, TSjplane%dspace_id, &
                          TSjplane%dset_id, hdferr)
         CALL h5dwrite_f(TSjplane%dset_id, H5T_NATIVE_DOUBLE, buffer_jplane, TSjplane%dimsf, hdferr)
         CALL h5sclose_f(TSjplane%dspace_id, hdferr)
         CALL h5dclose_f(TSjplane%dset_id, hdferr)
       enddo

       ! write iplane attribute
       ! attribute 1
       CALL h5screate_simple_f(1, TSjplane%attr_dims_no, attr_space, hdferr)
       CALL h5tcopy_f(H5T_NATIVE_INTEGER, atype_id, hdferr)
       CALL h5acreate_f(TSjplane%group_id,TSjplane%attr_name(1), atype_id,attr_space,attr_id,hdferr)
       CALL h5awrite_f(attr_id,H5T_NATIVE_INTEGER,attr_time,TSjplane%attr_dims_no,hdferr)
       CALL h5sclose_f(attr_space, hdferr)
       CALL h5aclose_f(attr_id, hdferr)  

      TSjplane%info_1st_array(1) =  ibe_j  
      TSjplane%info_1st_array(2) =  iend_j
      TSjplane%info_1st_array(3) =  iskip_j_inp

       ! attribute 2
       CALL h5screate_simple_f(1, TSjplane%attr_dims_info, attr_space, hdferr)
       CALL h5tcopy_f(H5T_NATIVE_INTEGER, atype_id, hdferr)
       CALL h5acreate_f(TSjplane%group_id,TSjplane%attr_name(2), atype_id,attr_space,attr_id,hdferr)
       CALL h5awrite_f(attr_id,H5T_NATIVE_INTEGER,TSjplane%info_1st_array,TSjplane%attr_dims_info,hdferr)
       CALL h5sclose_f(attr_space, hdferr)
       CALL h5aclose_f(attr_id, hdferr)  

       ! attribute 3
       CALL h5screate_simple_f(1, TSjplane%attr_dims_info, attr_space, hdferr)
       CALL h5tcopy_f(H5T_NATIVE_INTEGER, atype_id, hdferr)
       CALL h5acreate_f(TSjplane%group_id,TSjplane%attr_name(3), atype_id,attr_space,attr_id,hdferr)
       CALL h5awrite_f(attr_id,H5T_NATIVE_INTEGER,TSjplane%info_2nd_array,TSjplane%attr_dims_info,hdferr)
       CALL h5sclose_f(attr_space, hdferr)
       CALL h5aclose_f(attr_id, hdferr)  

       ! attribute 4
       CALL h5screate_simple_f(1, TSjplane%attr_dims_no, attr_space, hdferr)
       CALL h5tcopy_f(H5T_NATIVE_INTEGER, atype_id, hdferr)
       CALL h5acreate_f(TSjplane%group_id,TSjplane%attr_name(4), atype_id,attr_space,attr_id,hdferr)
       CALL h5awrite_f(attr_id,H5T_NATIVE_INTEGER,TSjplane%info_loc(j),TSjplane%attr_dims_no,hdferr)
       CALL h5sclose_f(attr_space, hdferr)
       CALL h5aclose_f(attr_id, hdferr) 

       CALL h5gclose_f(TSjplane%group_id, hdferr)
       CALL h5fclose_f(file_id, hdferr)
     
      deallocate(buffer_jplane)  
      deallocate(TSjplane%dimsf)
      deallocate(TSjplane%dname)  
      deallocate(tsdata)
    enddo ! end nnn loop
    enddo  ! end looping over j planes
  
  endif  ! end calculating const-j plane


  call FinalizeHDF5() 



  contains

    subroutine InitTSkplane_HDF5()
      implicit none
      integer :: nv = 10
 
      TSkplane%gname = '/kplane'
 !     if(allocated(TSkplane%dname)) deallocate(TSkplane%dname)
      allocate(TSkplane%dname(nv))
      TSkplane%dname(1) = "u"
      TSkplane%dname(2) = "v"
      TSkplane%dname(3) = "w"
      TSkplane%dname(4) = "p"
      TSkplane%dname(5) = "T"
      TSkplane%dname(6) = "uk"
      TSkplane%dname(7) = "vk"
      TSkplane%dname(8) = "wk"
      TSkplane%dname(9) = "pk"
      TSkplane%dname(10) = "Tk"

      TSkplane%attr_name(1) = "number of time steps"
      TSkplane%attr_name(2) = "istart, iend, iskip"
      TSkplane%attr_name(3) = "jstart, jend, jskip"
      TSkplane%attr_name(4) = "K locations" 

      TSkplane%rank = 4
!      if(allocated(TSkplane%dimsf)) deallocate(TSkplane%dimsf)
      allocate(TSkplane%dimsf(TSkplane%rank))
      TSkplane%dimsf(1) = ntpoint
      TSkplane%dimsf(2) = nypoint
      TSkplane%dimsf(3) = nxpoint
      TSkplane%dimsf(4) = 1

      if(allocated(buffer_kplane)) deallocate(buffer_kplane)
      allocate(buffer_kplane(TSkplane%dimsf(1), TSkplane%dimsf(2), TSkplane%dimsf(3), TSkplane%dimsf(4)))
 
      dname_time = 'time'
      dimsf_time(1) = ntpoint

!      TSkplane%info_1st_array(1) =  ibe_k_inp  
!      TSkplane%info_1st_array(2) =  iend_k_inp
!      TSkplane%info_1st_array(3) =  iskip_k_inp
      TSkplane%info_2nd_array(1) =  jbe_k_inp
      TSkplane%info_2nd_array(2) =  jend_k_inp  
      TSkplane%info_2nd_array(3) =  jskip_k_inp
      TSkplane%attr_dims_info = 3
      TSkplane%attr_dims_loc = num_kplane
      allocate(TSkplane%info_loc(num_kplane))
      TSkplane%info_loc = kplane
 
    end subroutine InitTSkplane_HDF5


    subroutine InitTSiplane_HDF5()
      implicit none
      integer :: nv = 10
 
      TSiplane%gname = '/iplane'
      
      allocate(TSiplane%dname(nv))
      TSiplane%dname(1) = "u"
      TSiplane%dname(2) = "v"
      TSiplane%dname(3) = "w"
      TSiplane%dname(4) = "p"
      TSiplane%dname(5) = "T"
      TSiplane%dname(6) = "ui"
      TSiplane%dname(7) = "vi"
      TSiplane%dname(8) = "wi"
      TSiplane%dname(9) = "pi"
      TSiplane%dname(10) = "Ti"

      TSiplane%attr_name(1) = "number of time steps"
      TSiplane%attr_name(2) = "jstart, jend, jskip"
      TSiplane%attr_name(3) = "kstart, kend, kskip"
      TSiplane%attr_name(4) = "I locations" 

      TSiplane%rank = 4
      allocate(TSiplane%dimsf(TSiplane%rank))
      TSiplane%dimsf(1) = ntpoint
      TSiplane%dimsf(2) = nypoint
      TSiplane%dimsf(3) = nzpoint
      TSiplane%dimsf(4) = 1

      allocate(buffer_iplane(TSiplane%dimsf(1), TSiplane%dimsf(2), TSiplane%dimsf(3), TSiplane%dimsf(4)))

      dname_time = 'time'
      dimsf_time(1) = ntpoint

!      TSiplane%info_1st_array(1) =  jbe_i_inp
!      TSiplane%info_1st_array(2) =  jend_i_inp
!      TSiplane%info_1st_array(3) =  jskip_i_inp
      TSiplane%info_2nd_array(1) =  kbe_i_inp  
      TSiplane%info_2nd_array(2) =  kend_i_inp  
      TSiplane%info_2nd_array(3) =  kskip_i_inp
      TSiplane%attr_dims_info = 3
      TSiplane%attr_dims_loc = num_iplane
      allocate(TSiplane%info_loc(num_iplane))
      TSiplane%info_loc = iplane  !!!!!!!!!!!!!!!
 
    end subroutine InitTSiplane_HDF5


    subroutine InitTSjplane_HDF5()
      implicit none
      integer :: nv = 10
 
      TSjplane%gname = '/jplane'
      allocate(TSjplane%dname(nv))
      TSjplane%dname(1) = "u"
      TSjplane%dname(2) = "v"
      TSjplane%dname(3) = "w"
      TSjplane%dname(4) = "p"
      TSjplane%dname(5) = "T"
      TSjplane%dname(6) = "uj"
      TSjplane%dname(7) = "vj"
      TSjplane%dname(8) = "wj"
      TSjplane%dname(9) = "pj"
      TSjplane%dname(10) = "Tj"

      TSjplane%attr_name(1) = "number of time steps"
      TSjplane%attr_name(2) = "istart, iend, iskip"
      TSjplane%attr_name(3) = "kstart, kend, kskip"
      TSjplane%attr_name(4) = "J locations" 

      TSjplane%rank = 4
      allocate(TSjplane%dimsf(TSjplane%rank))
      TSjplane%dimsf(1) = ntpoint
      TSjplane%dimsf(2) = nxpoint
      TSjplane%dimsf(3) = nzpoint
      TSjplane%dimsf(4) = 1

      allocate(buffer_jplane(TSjplane%dimsf(1), TSjplane%dimsf(2), TSjplane%dimsf(3), TSjplane%dimsf(4)))

      dname_time = 'time'
      dimsf_time(1) = ntpoint

!      TSjplane%info_1st_array(1) =  ibe_j_inp  
!      TSjplane%info_1st_array(2) =  iend_j_inp
!      TSjplane%info_1st_array(3) =  iskip_j_inp
      TSjplane%info_2nd_array(1) =  kbe_j_inp  
      TSjplane%info_2nd_array(2) =  kend_j_inp  
      TSjplane%info_2nd_array(3) =  kskip_j_inp
      TSjplane%attr_dims_info = 3
      TSjplane%attr_dims_loc = num_jplane
      allocate(TSjplane%info_loc(num_jplane))
      TSjplane%info_loc = jplane  
 
    end subroutine InitTSjplane_HDF5


    subroutine InitHDF5()
      ! Initialize HDF5 library and Fortran interfaces.
      CALL h5open_f(hdferr)
    end subroutine InitHDF5

    subroutine FinalizeHDF5()
      ! Close FORTRAN interfaces and HDF5 library.
      CALL h5close_f(hdferr)
    end subroutine FinalizeHDF5


    subroutine Input()
      integer, parameter :: nid=5
!      integer, parameter :: nid=10
      integer :: i,j
      if (nid.eq.5) then
        print*,'please input parameters or redirect from a file'
      else
        open(nid,file='convert.inp',status='old')
      end if
      read(nid,*)
      read(nid,*)ntpoint, nskip, iascii
      read(nid,*)
      read(nid,'(a)')filepath
      read(nid,*)
      read(nid,*)
      read(nid,*) ireadkp, convert_file_be_k
      read(nid,*)
      read(nid,*) ibe_k_inp, iend_k_inp, jbe_k_inp, jend_k_inp, nodes_k, ibuffer_k
      read(nid,*)
      read(nid,*) num_kplane
      allocate(kplane(num_kplane))
      read(nid,*) (kplane(i), i = 1, num_kplane)
      read(nid,*)
      read(nid,*)
      read(nid,*) ireadip, convert_file_be_i
      read(nid,*)
      read(nid,*) jbe_i_inp, jend_i_inp, kbe_i_inp, kend_i_inp, nodes_i, jbuffer_i
      read(nid,*)
      read(nid,*) num_iplane
      allocate(iplane(num_iplane))
      read(nid,*) (iplane(i), i = 1, num_iplane)
      read(nid,*)
      read(nid,*)
      read(nid,*) ireadjp, convert_file_be_j
      read(nid,*)
      read(nid,*) ibe_j_inp, iend_j_inp, kbe_j_inp, kend_j_inp, nodes_j, ibuffer_j
      read(nid,*)
      read(nid,*) num_jplane
      allocate(jplane(num_jplane))
      read(nid,*) (jplane(j), j = 1, num_jplane)
      
       if(convert_file_be_k.eq.0.or.convert_file_be_i.eq.0.or.convert_file_be_j.eq.0) then
        print *, 'convert_file_be_k or convert_file_be_i convert_file_be_j should not be 0'
        stop
      endif
!      jbe_k  = jbe_k_inp
!      jend_k = jend_k_inp
!      kbe_i  = kbe_i_inp
!      kend_i = kend_i_inp
!      kbe_j  = kbe_j_inp
!      kend_j = kend_j_inp

      if (nid.ne.5) close(nid)
      ntpoint = ntpoint - nskip
      if(ireadkp.eq.1) then
         print *, 'number of wall-normal locations =', num_kplane
         print *, 'Wall-normal plane indexes k =', kplane
      endif
      if(ireadip.eq.1) then
         print *, 'number of const-i locations =', num_iplane
         print *, 'const-i plane indexes i =', iplane
      endif
      if(ireadjp.eq.1) then
         print *, 'number of const-j locations =', num_jplane
         print *, 'const-j plane indexes j =', jplane
      endif
      write(*,*)
    end subroutine input

end program convert
