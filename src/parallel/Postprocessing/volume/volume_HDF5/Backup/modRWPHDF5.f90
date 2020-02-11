

    module modPRWHDF5
      use HDF5
      implicit none
      include 'mpif.h'

      integer(HID_T), private :: file_id, group_id, dset_id, dspace_id
      integer(HID_T), private :: memspace, plist_id
      integer(HSIZE_T), dimension(3) :: count, offset, offset2        !!!!!!
      integer, private :: hdferr
      logical :: IsHDF5Initialized = .false.

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



!P(fname,gname,nt,dim1,dim2,dim3,dim1st,dim2st,dim3st,vars,n)

















    end module modPRWHDF5
