
! Read1dHDF5: read all the 1d variables without the dataset name.
! Read2dHDF5: read all the 2d variables without the dataset name.
! Read3dHDF5: read all the 3d variables without the dataset name.

!    module modReadHDF5
    module modRWHDF5
      use HDF5
      implicit none

      integer(HID_T), private :: file_id, group_id, dset_id, dspace_id
      integer, private :: hdferr
      logical :: IsHDF5Initialized = .false.

  contains


     ! Detect HDF5 file. return # of datsset, dataset name, rank, dimsf
     subroutine DetectHDF5(fname,gname,num_dset,dname,rank,dims)
       character(*), intent(in) :: fname
       character(*), intent(in) :: gname
       integer, intent(out) :: num_dset
       character(*), dimension(:), allocatable, intent(out) :: dname
       integer, intent(out) :: rank
       integer, dimension(:), allocatable, intent(out) :: dims
       integer(HSIZE_T), dimension(:), allocatable :: dimsf, dimsf_total
       integer :: i, idx

       if(allocated(dname)) deallocate(dname)
       if(allocated(dimsf)) deallocate(dimsf)
       if(allocated(dimsf_total)) deallocate(dimsf_total)
       if(allocated(dims)) deallocate(dims)

       CALL h5fopen_f(trim(fname), H5F_ACC_RDONLY_F, file_id, hdferr)

       CALL h5gopen_f(file_id, trim(adjustl(gname)), group_id, hdferr)
       CALL h5gn_members_f(group_id, trim(adjustl(gname)), num_dset, hdferr)

       allocate(dname(num_dset))
       do idx=0, num_dset-1
         CALL h5gget_obj_info_idx_f(group_id, trim(adjustl(gname)), idx, dname(idx+1), H5G_DATASET_F, hdferr)
       enddo

       CALL h5dopen_f(group_id, trim(dname(1)), dset_id, hdferr)
       CALL h5dget_space_f(dset_id, dspace_id, hdferr)
       CALL h5sget_simple_extent_ndims_f(dspace_id, rank, hdferr)
       allocate(dimsf(rank), dimsf_total(rank))
       CALL h5sget_simple_extent_dims_f(dspace_id, dimsf, dimsf_total, hdferr)
       CALL h5dclose_f(dset_id, hdferr)
       CALL h5gclose_f(group_id, hdferr)

       CALL h5fclose_f(file_id, hdferr)

       allocate(dims(rank))
       dims = dimsf

  !     print *, 'rank = ', rank
 !      print *, 'dims = ', dims

     end subroutine DetectHDF5





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
     subroutine Read2dHDF5(fname,gname,num_dset,dname,buffer2d)
       character(*), intent(in) :: fname
       character(*), intent(in) :: gname
       integer, intent(out) :: num_dset
       character(*), dimension(:), allocatable, intent(out) :: dname
       real(8), dimension(:,:,:), allocatable, intent(out) :: buffer2d
       integer(HSIZE_T), dimension(:), allocatable :: dimsf, dimsf_total
       integer :: rank
       integer :: i, j, k, idx

       if(allocated(dname)) deallocate(dname)
       if(allocated(dimsf)) deallocate(dimsf)
       if(allocated(dimsf_total)) deallocate(dimsf_total)
       if(allocated(buffer2d)) deallocate(buffer2d)

       CALL h5fopen_f(trim(fname), H5F_ACC_RDONLY_F, file_id, hdferr)
       CALL h5gopen_f(file_id, trim(gname), group_id, hdferr)

       CALL h5gn_members_f(group_id, trim(gname), num_dset, hdferr)
       allocate(dname(num_dset))
       do idx=0, num_dset-1
         CALL h5gget_obj_info_idx_f(group_id, trim(gname), idx, dname(idx+1), H5G_DATASET_F, hdferr)
       enddo

       do i=1, num_dset
         CALL h5dopen_f(group_id, trim(dname(i)), dset_id, hdferr)
         CALL h5dget_space_f(dset_id, dspace_id, hdferr)
         CALL h5sget_simple_extent_ndims_f(dspace_id, rank, hdferr)
         if(i.eq.1) allocate(dimsf(rank), dimsf_total(rank))
         CALL h5sget_simple_extent_dims_f(dspace_id, dimsf, dimsf_total, hdferr)
         if(i.eq.1) allocate(buffer2d(dimsf(1),dimsf(2),num_dset))
         CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, buffer2d(:,:,i), dimsf, hdferr)
         CALL h5dclose_f(dset_id, hdferr)
       enddo

       CALL h5gclose_f(group_id, hdferr)
       CALL h5fclose_f(file_id, hdferr)

     end subroutine Read2dHDF5

     ! read all the 3d variables without the dataset name.
     subroutine Read3dHDF5(fname,gname,num_dset,dname,buffer3d)
       character(*), intent(in) :: fname
       character(*), intent(in) :: gname
       integer, intent(out) :: num_dset
       character(*), dimension(:), allocatable, intent(out) :: dname
       real(8), dimension(:,:,:,:), allocatable, intent(out) :: buffer3d
       integer(HSIZE_T), dimension(:), allocatable :: dimsf, dimsf_total
!       character(10), dimension(:), allocatable :: dname
       integer :: rank
       integer :: i, j, k, idx

       if(allocated(dname)) deallocate(dname)
       if(allocated(dimsf)) deallocate(dimsf)
       if(allocated(dimsf_total)) deallocate(dimsf_total)
       if(allocated(buffer3d)) deallocate(buffer3d)

       CALL h5fopen_f(trim(fname), H5F_ACC_RDONLY_F, file_id, hdferr)
       CALL h5gopen_f(file_id, trim(gname), group_id, hdferr)

       CALL h5gn_members_f(group_id, trim(gname), num_dset, hdferr)
       allocate(dname(num_dset))
       do idx=0, num_dset-1
         CALL h5gget_obj_info_idx_f(group_id, trim(gname), idx, dname(idx+1), H5G_DATASET_F, hdferr)
       enddo

       do i=1, num_dset
         CALL h5dopen_f(group_id, trim(dname(i)), dset_id, hdferr)
         CALL h5dget_space_f(dset_id, dspace_id, hdferr)
         CALL h5sget_simple_extent_ndims_f(dspace_id, rank, hdferr)
         if(i.eq.1) allocate(dimsf(rank), dimsf_total(rank))
         CALL h5sget_simple_extent_dims_f(dspace_id, dimsf, dimsf_total, hdferr)
         if(i.eq.1) allocate(buffer3d(dimsf(1),dimsf(2),dimsf(3),num_dset))
         CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, buffer3d(:,:,:,i), dimsf, hdferr)
         CALL h5dclose_f(dset_id, hdferr)
       enddo

       CALL h5gclose_f(group_id, hdferr)
       CALL h5fclose_f(file_id, hdferr)

     end subroutine Read3dHDF5

     ! read 1d variables with the dataset name
     subroutine Read1dHDF5_partial(fname,gname,num_dset,dname,buffer1d)
       character(*), intent(in) :: fname
       character(*), intent(in) :: gname
       integer, intent(in) :: num_dset
       character(*), intent(in) :: dname(num_dset)
       real(8), dimension(:,:), allocatable, intent(out) :: buffer1d
       integer(HSIZE_T), dimension(:), allocatable :: dimsf, dimsf_total
       integer :: rank
       integer :: i

       if(allocated(dimsf)) deallocate(dimsf)
       if(allocated(dimsf_total)) deallocate(dimsf_total)
       if(allocated(buffer1d)) deallocate(buffer1d)

       CALL h5fopen_f(trim(fname), H5F_ACC_RDONLY_F, file_id, hdferr)
       CALL h5gopen_f(file_id, trim(gname), group_id, hdferr)
       do i=1, num_dset
         CALL h5dopen_f(group_id, trim(dname(i)), dset_id, hdferr)
         CALL h5dget_space_f(dset_id, dspace_id, hdferr)
         CALL h5sget_simple_extent_ndims_f(dspace_id, rank, hdferr)
         if(i.eq.1) allocate(dimsf(rank), dimsf_total(rank))
         CALL h5sget_simple_extent_dims_f(dspace_id, dimsf, dimsf_total, hdferr)
         if(i.eq.1) allocate(buffer1d(dimsf(1),num_dset))
         CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, buffer1d(:,i), dimsf, hdferr)
         CALL h5dclose_f(dset_id, hdferr)
       enddo
       CALL h5gclose_f(group_id, hdferr)
       CALL h5fclose_f(file_id, hdferr)

     end subroutine Read1dHDF5_partial

     ! read 2d selected variables with the dataset name
     subroutine Read2dHDF5_partial(fname,gname,num_dset,dname,buffer2d)
       character(*), intent(in) :: fname
       character(*), intent(in) :: gname
       integer, intent(in) :: num_dset
       character(*), intent(in) :: dname(num_dset)
       real(8), dimension(:,:,:), allocatable, intent(out) :: buffer2d
       integer(HSIZE_T), dimension(:), allocatable :: dimsf, dimsf_total
       integer :: rank
       integer :: i

       if(allocated(dimsf)) deallocate(dimsf)
       if(allocated(dimsf_total)) deallocate(dimsf_total)
       if(allocated(buffer2d)) deallocate(buffer2d)

       CALL h5fopen_f(trim(fname), H5F_ACC_RDONLY_F, file_id, hdferr)
       CALL h5gopen_f(file_id, trim(gname), group_id, hdferr)
       do i=1, num_dset
         CALL h5dopen_f(group_id, trim(dname(i)), dset_id, hdferr)
         CALL h5dget_space_f(dset_id, dspace_id, hdferr)
         CALL h5sget_simple_extent_ndims_f(dspace_id, rank, hdferr)
         if(i.eq.1) allocate(dimsf(rank), dimsf_total(rank))
         CALL h5sget_simple_extent_dims_f(dspace_id, dimsf, dimsf_total, hdferr)
         if(i.eq.1) allocate(buffer2d(dimsf(1),dimsf(2),num_dset))
         CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, buffer2d(:,:,i), dimsf, hdferr)
         CALL h5dclose_f(dset_id, hdferr)
       enddo
       CALL h5gclose_f(group_id, hdferr)
       CALL h5fclose_f(file_id, hdferr)

     end subroutine Read2dHDF5_partial

     ! read 3d selected variables with the dataset name
     subroutine Read3dHDF5_partial(fname,gname,num_dset,dname,buffer3d)
       character(*), intent(in) :: fname
       character(*), intent(in) :: gname
       integer, intent(in) :: num_dset
       character(*), intent(in) :: dname(num_dset)
       real(8), dimension(:,:,:,:), allocatable, intent(out) :: buffer3d
       integer(HSIZE_T), dimension(:), allocatable :: dimsf, dimsf_total
       integer :: rank
       integer :: i

       if(allocated(dimsf)) deallocate(dimsf)
       if(allocated(dimsf_total)) deallocate(dimsf_total)
       if(allocated(buffer3d)) deallocate(buffer3d)

       CALL h5fopen_f(trim(fname), H5F_ACC_RDONLY_F, file_id, hdferr)
       CALL h5gopen_f(file_id, trim(gname), group_id, hdferr)
       do i=1, num_dset
         CALL h5dopen_f(group_id, trim(dname(i)), dset_id, hdferr)
         CALL h5dget_space_f(dset_id, dspace_id, hdferr)
         CALL h5sget_simple_extent_ndims_f(dspace_id, rank, hdferr)
         if(i.eq.1) allocate(dimsf(rank), dimsf_total(rank))
         CALL h5sget_simple_extent_dims_f(dspace_id, dimsf, dimsf_total, hdferr)
         if(i.eq.1) allocate(buffer3d(dimsf(1),dimsf(2),dimsf(3),num_dset))
         CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, buffer3d(:,:,:,i), dimsf, hdferr)
         CALL h5dclose_f(dset_id, hdferr)
       enddo
       CALL h5gclose_f(group_id, hdferr)
       CALL h5fclose_f(file_id, hdferr)

     end subroutine Read3dHDF5_partial

     subroutine ReadAveAcousticGridMetrics(fname,gname,buffer)
       character(*), intent(in) :: fname
       character(*), intent(in) :: gname
       real(8), dimension(:,:,:), allocatable, intent(out) :: buffer
       integer(HSIZE_T), dimension(:), allocatable :: dimsf, dimsf_total
       character(10), dimension(:), allocatable :: dname
       integer :: num_dset
       integer :: rank
       integer :: i, j, k, idx

       CALL h5fopen_f(trim(fname), H5F_ACC_RDONLY_F, file_id, hdferr)
       CALL h5gopen_f(file_id, trim(gname), group_id, hdferr)

       CALL h5gn_members_f(group_id, trim(gname), num_dset, hdferr)
       allocate(dname(num_dset))
       do idx=0, num_dset-1
         CALL h5gget_obj_info_idx_f(group_id, trim(gname), idx, dname(idx+1), H5G_DATASET_F, hdferr)
       enddo

       do i=1, num_dset
         CALL h5dopen_f(group_id, trim(dname(i)), dset_id, hdferr)
         CALL h5dget_space_f(dset_id, dspace_id, hdferr)
         CALL h5sget_simple_extent_ndims_f(dspace_id, rank, hdferr)
         if(i.eq.1) allocate(dimsf(rank), dimsf_total(rank))
         CALL h5sget_simple_extent_dims_f(dspace_id, dimsf, dimsf_total, hdferr)
         if(i.eq.1) allocate(buffer(dimsf(1),dimsf(2),num_dset))
         CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, buffer(:,:,i), dimsf, hdferr)
         CALL h5dclose_f(dset_id, hdferr)
       enddo

       CALL h5gclose_f(group_id, hdferr)
       CALL h5fclose_f(file_id, hdferr)

       print *, 'dimsf_grid(1) = ', dimsf(1)
       print *, 'dimsf_grid(2) = ', dimsf(2)

     end subroutine ReadAveAcousticGridMetrics








     subroutine ReadHDF5sol(fname,idim,jdim,kdim,num_dset,vars)
       character(*), intent(in) :: fname
       integer, intent(in) :: idim, jdim, kdim, num_dset
       real, intent(out) :: vars(kdim,idim,jdim,num_dset)
       integer(HSIZE_T), dimension(:), allocatable :: dimsf
       character(10), dimension(:), allocatable :: dname
       integer :: rank
       integer :: i, j, k

       if(.not.IsHDF5Initialized) then
         print *, 'HDF5 should be Initialized'
         stop
       endif

       rank = 3
       allocate(dname(num_dset))
       allocate(dimsf(rank))

       dname(1) = "u"
       dname(2) = "v"
       dname(3) = "w"
       dname(4) = "p"
       dname(5) = "T"
       dimsf(1) = kdim
       dimsf(2) = idim
       dimsf(3) = jdim

       CALL h5fopen_f(trim(fname), H5F_ACC_RDONLY_F, file_id, hdferr)
       do i=1, num_dset
         CALL h5dopen_f(file_id, dname(i), dset_id, hdferr)
         CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, vars(:,:,:,i), dimsf, hdferr)
         CALL h5dclose_f(dset_id, hdferr)
       enddo
       CALL h5fclose_f(file_id, hdferr)

     end subroutine ReadHDF5sol

!     subroutine ReadtsGrid_jplane(fname,imax,kmax,num_dset,vars)
!       character(*), intent(in) :: fname
!       integer, intent(in) :: imax, kmax, num_dset
!       real, intent(out) :: vars(:,:,:,:)         !!!!
!       integer :: i, j, k, rank
!       character(10), dimension(:), allocatable :: dname
!       integer(HSIZE_T), dimension(:), allocatable :: dimsf

!       rank = 3
!       allocate(dname(num_dset))
!       allocate(dimsf(rank))

!       dimsf(1) = imax
!       dimsf(2) = kmax
!       dimsf(3) = 1      !!!!1

!       dname(1) = "x"
!       dname(2) = "y"
!       dname(3) = "z"
!       dname(4) = "didx"
!       dname(5) = "djdx"
!       dname(6) = "dkdx"
!       dname(7) = "didy"
!       dname(8) = "djdy"
!       dname(9) = "dkdy"
!       dname(10) = "didz"
!       dname(11) = "djdz"
!       dname(12) = "dkdz"

!       CALL h5fopen_f(trim(fname), H5F_ACC_RDONLY_F, file_id, hdferr)
!       CALL h5gopen_f(file_id,'/jplane',group_id,hdferr)
!       do i=1, num_dset
!         CALL h5dopen_f(group_id, dname(i), dset_id, hdferr)
!         CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, vars(:,:,:,i), dimsf, hdferr)
!         CALL h5dclose_f(dset_id, hdferr)
!       enddo
!       CALL h5gclose_f(group_id, hdferr)
!       CALL h5fclose_f(file_id, hdferr)

!     end subroutine ReadtsGrid_jplane

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



     ! read selected variables
     subroutine Readtsflow_partial(fname,gname,nt,dim1,dim2,dim3,num_dset,dname,vars)
       character(*), intent(in) :: fname
       character(*), intent(in) :: gname
       integer, intent(in) :: nt, dim1, dim2, dim3, num_dset
       real(8), intent(out) :: vars(:,:,:,:,:)
       character(*), intent(in) :: dname(:) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       integer(HSIZE_T), dimension(:), allocatable :: dimsf
       integer :: rank
       integer :: i

       rank = 4
       allocate(dimsf(rank))

       dimsf(1) = nt
       dimsf(2) = dim1
       dimsf(3) = dim2
       dimsf(4) = dim3

       CALL h5fopen_f(trim(fname), H5F_ACC_RDONLY_F, file_id, hdferr)
       CALL h5gopen_f(file_id,trim(gname),group_id,hdferr)
       do i=1, num_dset
         CALL h5dopen_f(group_id, dname(i), dset_id, hdferr)
         CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, vars(:,:,:,:,i), dimsf, hdferr)
         CALL h5dclose_f(dset_id, hdferr)
       enddo
       CALL h5gclose_f(group_id, hdferr)
       CALL h5fclose_f(file_id, hdferr)

       deallocate(dimsf)

     end subroutine Readtsflow_partial






















     ! read timeseries flow data.
     ! file name: fname. group name: gname
     ! nt for # of time steps
     ! n: 1 for kplane, 2 for iplane, 3 for jplane
     subroutine Readtsflow(fname,gname,nt,dim1,dim2,dim3,vars,n)
       character(*), intent(in) :: fname
       character(*), intent(in) :: gname
       integer, intent(in) :: nt, dim1, dim2, dim3, n
       real(8), intent(out) :: vars(:,:,:,:,:)
       character(10), dimension(:), allocatable :: dname
       integer(HSIZE_T), dimension(:), allocatable :: dimsf
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
         print *, ' n should be 1: kplane, 2: iplane, 3: iplane'
         stop
       endif

       CALL h5fopen_f(trim(fname), H5F_ACC_RDONLY_F, file_id, hdferr)
       CALL h5gopen_f(file_id,trim(gname),group_id,hdferr)
       do i=1, num_dset
         CALL h5dopen_f(group_id, dname(i), dset_id, hdferr)
         CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, vars(:,:,:,:,i), dimsf, hdferr)
         CALL h5dclose_f(dset_id, hdferr)
       enddo
       CALL h5gclose_f(group_id, hdferr)
       CALL h5fclose_f(file_id, hdferr)

       deallocate(dname)
       deallocate(dimsf)
!       print *, 'min(rho) = ', minval(vars(:,:,:,:,4)/(287.0*vars(:,:,:,:,5)))
     end subroutine Readtsflow

     subroutine InitHDF5()
      ! Initialize HDF5 library and Fortran interfaces.
       CALL h5open_f(hdferr)
       IsHDF5Initialized = .True.
     end subroutine InitHDF5

     subroutine FinalizeHDF5()
      ! Close FORTRAN interfaces and HDF5 library.
       CALL h5close_f(hdferr)
     end subroutine FinalizeHDF5


!    end module modReadHDF5

!    module modWriteHDF5
!      use modVolume  !!!!!!!
!      use HDF5
!      implicit none

!      integer(HID_T), private :: file_id, group_id, dset_id, dspace_id
!      integer, private :: hdferr

!   contains

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


      ! write 3d single variables
      subroutine Write3dHDF5_svariable(fname,rank,dimsf,dname,vars,iexist)
        character(*), intent(in) :: fname
        character(*), intent(in) :: dname
        integer, intent(in) :: rank
        integer(HSIZE_T), intent(in) :: dimsf(rank)   !!!!!!!!!!!!!!!!
        real, intent(in) :: vars(:,:,:)
        integer, intent(in) :: iexist

        if(iexist.eq.0) then
          CALL h5fcreate_f(trim(fname), H5F_ACC_TRUNC_F, file_id, hdferr)
        else
          CALL h5fopen_f(trim(fname), H5F_ACC_RDWR_F, file_id, hdferr )
        endif

        CALL h5screate_simple_f(rank, dimsf, dspace_id, hdferr)
        CALL h5dcreate_f(file_id, trim(dname), H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)
        CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, vars, dimsf, hdferr)
        CALL h5sclose_f(dspace_id, hdferr)
        CALL h5dclose_f(dset_id, hdferr)

        CALL h5fclose_f(file_id, hdferr)

      end subroutine Write3dHDF5_svariable


!    end module modWriteHDF5

    end module modRWHDF5

