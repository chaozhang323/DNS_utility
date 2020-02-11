



     program Cal_kplane
       use HDF5
       implicit none

       character(400) :: filepath
       character(400) :: filepath_output
       character(200) :: fname

       character(10), dimension(:), allocatable :: dname
       character(10), dimension(:), allocatable :: dname_grid
       character(10), dimension(:), allocatable :: gname
       integer(HID_T) :: file_id
       integer(HID_T) :: group_id, dset_id, dspace_id
       integer(HID_T) :: memspace_id
       integer(HSIZE_T), dimension(:), allocatable :: dimsf
       integer(HSIZE_T), dimension(:), allocatable :: dimsf_total

       integer(HSIZE_T), dimension(:), allocatable :: dimsm
       integer(HSIZE_T), dimension(:), allocatable :: count, offset, block, stride


       integer(HSIZE_T), dimension(:), allocatable :: dimsf_grid
       integer(HSIZE_T), dimension(:), allocatable :: dimsf_total_grid

       integer(HSIZE_T), dimension(:), allocatable :: dimsm_grid
       integer(HSIZE_T), dimension(:), allocatable :: count_grid, offset_grid, block_grid, stride_grid


       integer :: rank
       integer :: rank_grid
       integer :: file_be, file_end, file_skip
       character(8) :: fnum_8

       integer :: num_dset, num_group

       integer :: i, j, k, n, num_file, nnn
       integer :: hdferr
       real(8), dimension(:,:,:,:), allocatable :: buffer_kplane, buffer_total
       real(8), dimension(:,:,:,:), allocatable :: buffer_grid

       real(8), dimension(:,:,:), allocatable :: buffer_average_tmp
       real(8), dimension(:,:), allocatable :: buffer_average_1D
       real(8), dimension(:), allocatable :: dudz
       real(8) :: mu
       real(8), dimension(:), allocatable :: Cf
       real(8) :: uinf=1881.88
       real(8) :: rhoinf=0.0172762

       integer :: ntpoint_pfile, nypoint, nxpoint                    !!!!!!!!!!!!!!!


       call Init()

       num_file = (file_end - file_be)/file_skip + 1

       print *, 'Number of Files = ',  num_file
       print *, 'File begin number: ', file_be
       print *, 'File end number: ',   file_end
       print *, 'File skip number: ',  file_skip
       print *, '*********************************************'

       CALL InitHDF5()

       rank = 4
       rank_grid = 3
       num_dset = 2
       num_group = 1

       allocate(dname(num_dset), gname(num_group))
       allocate(dname_grid(num_dset))
       allocate(dimsf(rank), dimsf_total(rank), dimsm(rank))
       allocate(count(rank), offset(rank), block(rank), stride(rank))
       allocate(dimsf_grid(rank_grid), dimsf_total_grid(rank_grid), dimsm_grid(rank_grid))
       allocate(count_grid(rank_grid), offset_grid(rank_grid), block_grid(rank_grid), stride_grid(rank_grid))

       dname(1) = "uk"
       dname_grid(1) = "x"
       dname_grid(2) = "dkdz"
       gname(1) = "/kplane"

       dimsf(1) = ntpoint_pfile
       dimsf(2) = nypoint
       dimsf(3) = nxpoint
       dimsf(4) = 1

       dimsm(1) = dimsf(1)
       dimsm(2) = dimsf(2)
       dimsm(3) = dimsf(3)
       dimsm(4) = dimsf(4)

       block(1) = dimsf(1)
       block(2) = dimsf(2)
       block(3) = dimsf(3)
       block(4) = dimsf(4)

       offset(1) = 0
       offset(2) = 0
       offset(3) = 0
       offset(4) = 0

       stride(1) = 1
       stride(2) = 1
       stride(3) = 1
       stride(4) = 1

       count(1) = 1
       count(2) = 1
       count(3) = 1
       count(4) = 1

       dimsf_grid(1) = nypoint
       dimsf_grid(2) = nxpoint
       dimsf_grid(3) = 1

       dimsm_grid(1) = dimsf_grid(1)
       dimsm_grid(2) = dimsf_grid(2)
       dimsm_grid(3) = dimsf_grid(3)

       block_grid(1) = dimsf_grid(1)
       block_grid(2) = dimsf_grid(2)
       block_grid(3) = dimsf_grid(3)

       offset_grid(1) = 0
       offset_grid(2) = 0
       offset_grid(3) = 0

       stride_grid(1) = 1
       stride_grid(2) = 1
       stride_grid(3) = 1

       count_grid(1) = 1
       count_grid(2) = 1
       count_grid(3) = 1


       allocate(buffer_kplane(ntpoint_pfile, nypoint, nxpoint,1))
       allocate(buffer_total(ntpoint_pfile*num_file, nypoint, nxpoint,1))

       allocate(buffer_average_tmp(ntpoint_pfile*num_file,nxpoint,1))
       allocate(buffer_average_1D(nxpoint,1))

       allocate(buffer_grid(nypoint, nxpoint, 1, num_dset))
       allocate(dudz(nxpoint))
       allocate(Cf(nxpoint))


       nnn = 0
       do n=1, num_file
         write(unit=fnum_8,fmt='(I08.8)') (file_be + (n-1)*file_skip)
         fname = trim(filepath)//fnum_8//'.h5'

         CALL h5fopen_f(trim(fname), H5F_ACC_RDONLY_F, file_id, hdferr)

         CALL h5gopen_f(file_id, trim(gname(1)), group_id, hdferr)
         CALL h5screate_simple_f(rank, dimsm, memspace_id, hdferr)

         CALL h5dopen_f(group_id, trim(dname(1)), dset_id, hdferr)
         CALL h5dget_space_f(dset_id, dspace_id, hdferr)
         CALL h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, &
                                    count, hdferr, stride, block)
         CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, buffer_kplane,dimsf, hdferr, &
                        file_space_id=dspace_id, mem_space_id=memspace_id)
         CALL h5sclose_f(dspace_id, hdferr)
         CALL h5dclose_f(dset_id, hdferr)

         CALL h5gclose_f(group_id, hdferr)
         CALL h5fclose_f(file_id, hdferr)

         buffer_total((ntpoint_pfile*(n-1)+1):ntpoint_pfile*n,1:nypoint, 1:nxpoint,1) = &
         buffer_kplane(1:ntpoint_pfile, 1:nypoint, 1:nxpoint, 1)

       enddo

       print *, 'dimsf(1) = ', size(buffer_total, dim=2)
       print *, 'size = ', size(buffer_average_tmp, dim=1)

 !      print *, 'sum = ', sum(buffer_total, dim=2)

       buffer_average_tmp = sum(buffer_total, dim=2)/size(buffer_total, dim=2)
       buffer_average_1D = sum(buffer_average_tmp, dim=1)/size(buffer_average_tmp, dim=1)
 !      buffer_average_1D = sum(buffer_average_tmp, dim=4)/size(buffer_average_tmp, dim=1)




!       print *, 'buffer_average_1D(1,1) = ', buffer_average_1D(1,1)
!       print *, 'buffer_average_1D(2,1) = ', buffer_average_1D(2,1)
!       print *, 'buffer_average_1D(3,1) = ', buffer_average_1D(3,1)


!       print *, 'buffer_total(1,1,1,1) = ', buffer_total(1,1,1,1)
!       print *, 'buffer_total(1,1,1,1) = ', buffer_total(2,1,1,1)
!       print *, 'buffer_total(1,1,1,1) = ', buffer_total(3,1,1,1)
!       print *, 'buffer_total(1,1,2,1) = ', buffer_total(1,1,2,1)
!       print *, 'buffer_total(1,1,1354,1) = ', buffer_total(1,1,1354,1)


         fname = trim(filepath)//'GridMetrics.h5'

         CALL h5fopen_f(trim(fname), H5F_ACC_RDONLY_F, file_id, hdferr)

         CALL h5gopen_f(file_id, trim(gname(1)), group_id, hdferr)
         CALL h5screate_simple_f(rank_grid, dimsm_grid, memspace_id, hdferr)
         do i=1, num_dset
           CALL h5dopen_f(group_id, trim(dname_grid(i)), dset_id, hdferr)
           CALL h5dget_space_f(dset_id, dspace_id, hdferr)
           CALL h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset_grid, &
                                      count_grid, hdferr, stride_grid, block_grid)
           CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, buffer_grid(:,:,:,i),dimsf_grid, hdferr, &
                          file_space_id=dspace_id, mem_space_id=memspace_id)
           CALL h5sclose_f(dspace_id, hdferr)
           CALL h5dclose_f(dset_id, hdferr)
         enddo
         CALL h5gclose_f(group_id, hdferr)
         CALL h5fclose_f(file_id, hdferr)

!     print *, 'buffer_grid(1,1,1) = ', buffer_grid(1,1,1,2)
!     print *, 'buffer_grid(2,1,1) = ', buffer_grid(2,1,1,2)
!     print *, 'buffer_grid(3,1,1) = ', buffer_grid(3,1,1,2)

!stop

         do i=1, nxpoint
           dudz(i) = buffer_average_1D(i,1)*buffer_grid(1,i,1,2)
         enddo

         mu = (1.458e-6)*(300.**1.5)/dble(300.+110.4)

         do i=1, nxpoint
           Cf(i) = mu*dudz(i)/(0.5*rhoinf*uinf**2)
         enddo

         open(7,file='Cf.dat', status='unknown')
           write(7,*) 'x, Cf'
           do i=1, nxpoint
             write(7,*) buffer_grid(1,i,1,1), Cf(i)

           enddo
         close(7)




!     print *, 'buffer_grid(1,1,1) = ', buffer_grid(1,1,1,2)
!     print *, 'buffer_grid(2,1,1) = ', buffer_grid(2,1,1,2)
!     print *, 'buffer_grid(3,1,1) = ', buffer_grid(3,1,1,2)










       CALL FinalizeHDF5()

  contains


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

         read(*,*)
         read(*,'(a)') filepath
         read(*,*)
         read(*,'(a)') filepath_output
         read(*,*)
         read(*,*) file_be, file_end, file_skip
         read(*,*)
         read(*,*) ntpoint_pfile, nypoint, nxpoint



       end subroutine Init





     end program Cal_kplane
