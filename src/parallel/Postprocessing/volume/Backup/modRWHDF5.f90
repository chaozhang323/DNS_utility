


    module modReadHDF5
      use HDF5
      implicit none

      integer(HID_T), private :: file_id, group_id, dset_id, dspace_id
      integer, private :: hdferr
      logical :: IsHDF5Initialized = .false.

  contains

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


    end module modReadHDF5

    module modWriteHDF5
      use modVolume  !!!!!!!
      use HDF5
      implicit none

      integer(HID_T), private :: file_id, group_id, dset_id, dspace_id
      integer, private :: hdferr

   contains

      subroutine Write3DHDF5(fname,rank,dimsf,dname,vars,iexit)
        character(*), intent(in) :: fname
        character(*), intent(in) :: dname
        integer, intent(in) :: rank
        integer(HSIZE_T), intent(in) :: dimsf(rank)
        real, intent(in) :: vars(:,:,:)
        integer, intent(in) :: iexit


!        print *, 'vars(1,1,1) = ', real(vars(1,1,1))
!        print *, 'vars(2,1,1) = ', vars(2,1,1)
!        print *, 'dimsf(1) = ', dimsf(1)
!        print *, 'dimsf(2) = ', dimsf(2)
!        print *, 'dimsf(3) = ', dimsf(3)

        if(iexit.eq.0) then
          call h5fcreate_f(trim(fname), H5F_ACC_TRUNC_F, file_id, hdferr)
        else
          call h5fopen_f(trim(fname), H5F_ACC_RDWR_F, file_id, hdferr )
        endif

        call h5screate_simple_f(rank, dimsf, dspace_id, hdferr)
        call h5dcreate_f(file_id, trim(dname), H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, vars, dimsf, hdferr)
        call h5sclose_f(dspace_id, hdferr)
        call h5dclose_f(dset_id, hdferr)

        call h5fclose_f(file_id, hdferr)

      end subroutine Write3DHDF5


    end module modWriteHDF5



