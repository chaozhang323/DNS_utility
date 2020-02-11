


     program Cal_iplane
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
       real(8), dimension(:,:,:,:), allocatable :: buffer_iplane, buffer_total  !!!!!
       real(8), dimension(:,:,:,:), allocatable :: buffer_igrid

       real(8), dimension(:), allocatable :: pave, p2ave, prms
       real(8) :: uinf
       real(8) :: rhoinf

       integer :: ntpoint_pfile, nypoint, nzpoint
       integer :: numpt


       call Init()

       num_file = (file_end - file_be)/file_skip + 1

       print *, 'Number of Files = ',  num_file
       print *, 'File begin number: ', file_be
       print *, 'File end number: ',   file_end
       print *, 'File skip number: ',  file_skip
       print *, '*********************************************'


       CALL InitHDF5()

       rank = 4
       rank_grid = 4
       num_dset = 3
       num_group = 1

       allocate(dname(num_dset), gname(num_group))
       allocate(dname_grid(num_dset))
       allocate(dimsf(rank), dimsf_total(rank), dimsm(rank))

       allocate(dimsf_grid(rank_grid), dimsf_total_grid(rank_grid), dimsm_grid(rank_grid))


       gname = "/iplane"

       dname(1) = "p"

       dname_grid(1) = "z"


       dimsf(1) = ntpoint_pfile
       dimsf(2) = nypoint
       dimsf(3) = nzpoint
       dimsf(4) = 1

       dimsf_grid(1) = nypoint
       dimsf_grid(2) = nzpoint
       dimsf_grid(3) = 1


       allocate(buffer_iplane(ntpoint_pfile, nypoint, nzpoint, 1))
       allocate(buffer_igrid(nypoint, nzpoint, 1, 1))  !!!!!!!!!!!!!!!!!! "z"
       allocate(pave(nzpoint), p2ave(nzpoint))
       allocate(prms(nzpoint))
!       allocate(dudz(nxpoint))



       pave = 0.
       p2ave = 0.

       do n=1, num_file
         write(unit=fnum_8,fmt='(I08.8)') (file_be + (n-1)*file_skip)
         fname = trim(filepath)//fnum_8//'.h5'
         CALL h5fopen_f(trim(fname), H5F_ACC_RDONLY_F, file_id, hdferr)
         CALL h5gopen_f(file_id, trim(gname(1)), group_id, hdferr)
         CALL h5dopen_f(group_id, trim(dname(1)), dset_id, hdferr)
         CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, buffer_iplane(:,:,:,1),dimsf, hdferr)
         CALL h5dclose_f(dset_id, hdferr)
         CALL h5gclose_f(group_id, hdferr)
         CALL h5fclose_f(file_id, hdferr)


         do k=1, nzpoint
         numpt = 0
           do j=1, nypoint
             do nnn=1, ntpoint_pfile
               pave(k) = pave(k) + buffer_iplane(nnn,j,k,1)
               p2ave(k) = p2ave(k) + buffer_iplane(nnn,j,k,1)**2
               numpt = numpt + 1
             enddo
           enddo
 !          pave(k) = pave(k)/numpt
 !          p2ave(k) = p2ave(k)/numpt
         enddo


       enddo ! end n loop

       pave(1:nzpoint) = pave(1:nzpoint)/(numpt*num_file)
       p2ave(1:nzpoint) = p2ave(1:nzpoint)/(numpt*num_file)

       do i=1, nzpoint
         prms(i) = sqrt(abs(p2ave(i) - pave(i)**2))
       enddo

       print *, 'numpt = ', numpt
!       print *, 'prms(1) = ', prms

         ! read grid info
         fname = trim(filepath)//'GridMetrics.h5'

         CALL h5fopen_f(trim(fname), H5F_ACC_RDONLY_F, file_id, hdferr)
         CALL h5gopen_f(file_id, trim(gname(1)), group_id, hdferr)
         CALL h5dopen_f(group_id, trim(dname_grid(1)), dset_id, hdferr)
         CALL h5dget_space_f(dset_id, dspace_id, hdferr)
         CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, buffer_igrid(:,:,:,1),dimsf_grid, hdferr)
         CALL h5dclose_f(dset_id, hdferr)
         CALL h5gclose_f(group_id, hdferr)
         CALL h5fclose_f(file_id, hdferr)



         open(7,file='prms.dat', status='unknown')
           write(7,*) 'z, prms'
           do i=1, nzpoint
             write(7,*) buffer_igrid(1,i,1,1), prms(i)

           enddo
         close(7)

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
         read(*,*) ntpoint_pfile, nypoint, nzpoint
         read(*,*)
         read(*,*) uinf, rhoinf

 !      real(8) :: uinf=1881.88
 !      real(8) :: rhoinf=0.0172762

       end subroutine Init






     end program Cal_iplane
