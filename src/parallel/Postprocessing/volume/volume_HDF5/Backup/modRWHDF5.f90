


    module modReadHDF5

      implicit none

  contains

     subroutine readHDF5sol(fname,imax,jmax,kmax,num_dset,vars)
       use HDF5
       character(*), intent(in) :: fname
       integer, intent(in) :: imax, jmax, kmax, num_dset
       real, intent(out) :: vars(kmax,imax,jmax,num_dset)
       integer :: i, j, k, hdferr, rank
       integer(HID_T) :: file_id, group_id, dset_id, dspace_id
       integer(HSIZE_T), dimension(:), allocatable :: dimsf
       character(10), dimension(:), allocatable :: dname

       rank = 3
       allocate(dname(num_dset))
       allocate(dimsf(rank))

       dname(1) = "u"
       dname(2) = "v"
       dname(3) = "w"
       dname(4) = "p"
       dname(5) = "T"
       dimsf(1) = kmax
       dimsf(2) = imax
       dimsf(3) = jmax

       CALL h5open_f(hdferr)

       CALL h5fopen_f(trim(fname), H5F_ACC_RDONLY_F, file_id, hdferr)

       do i=1, num_dset
         CALL h5dopen_f(file_id, dname(i), dset_id, hdferr)
         CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, vars(:,:,:,i), dimsf, hdferr)
         CALL h5dclose_f(dset_id, hdferr)
       enddo
       CALL h5fclose_f(file_id, hdferr)

       CALL h5close_f(hdferr)


     end subroutine readHDF5sol


    end module modReadHDF5

    module modWriteHDF5
      use modVolume  !!!!!!!
      use HDF5
      implicit none



   contains

      subroutine Write3DHDF5(fname,rank,dimsf,dname,vars,iexit)
        character(*), intent(in) :: fname
        character(*), intent(in) :: dname
        integer, intent(in) :: rank
        integer(HSIZE_T), intent(in) :: dimsf(rank)
!        real, intent(in) :: vars(dimsf(1),dimsf(2),dimsf(3))
        real, intent(in) :: vars(:,:,:)
        integer, intent(in) :: iexit
        integer(HID_T) :: file_id, group_id, dset_id, dspace_id
        integer :: hdferr

        call h5open_f(hdferr)

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
        call h5close_f(hdferr)

      end subroutine Write3DHDF5


    end module modWriteHDF5



