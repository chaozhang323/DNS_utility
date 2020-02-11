



     program Cal_jplane
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

       integer :: num_dset, num_dset_grid, num_group

       integer :: i, j, k, n, num_file, nnn
       integer :: kbl
       integer :: hdferr

       real(8), dimension(:,:,:,:,:), allocatable :: buffer_jplane, buffer_total   !!!!!!!
       real(8), dimension(:,:,:,:), allocatable :: buffer_jgrid

       real(8), dimension(:,:), allocatable :: uave2d, rhoave2d  !!!!!!

       real(8), dimension(:), allocatable :: delta, theta, dstar  !!!!!!!!!!!!!

       real(8), dimension(:,:), allocatable :: pave2d, p2, prms      !!!!!!!!!
       integer :: temp(1)
       real(8) :: uinf
       real(8) :: rhoinf, rbar = 287.0
       integer :: numpt

       integer :: ntpoint_pfile, nxpoint, nzpoint

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
       num_dset = 3  !!!!!!!
       num_dset_grid = 3
       num_group = 1

       allocate(dname(num_dset), gname(num_group))
       allocate(dname_grid(num_dset_grid))
       allocate(dimsf(rank), dimsf_total(rank), dimsm(rank))

       allocate(dimsf_grid(rank_grid), dimsf_total_grid(rank_grid), dimsm_grid(rank_grid))

       allocate(delta(nxpoint))
       allocate(theta(nxpoint))
       allocate(dstar(nxpoint))


       allocate(uave2d(nzpoint, nxpoint), rhoave2d(nzpoint, nxpoint))

       allocate(pave2d(nzpoint, nxpoint), p2(nzpoint, nxpoint)) !!!!!!!!!!!!!!!
       allocate(prms(nzpoint, nxpoint))

       gname = "/jplane"

       dname(1) = "u"
       dname(2) = "p"
       dname(3) = "T"

       dname_grid(1) = "x"
       dname_grid(2) = "z"
       dname_grid(3) = "dkdz"

       dimsf(1) = ntpoint_pfile
       dimsf(2) = nxpoint
       dimsf(3) = nzpoint
       dimsf(4) = 1

       dimsf_grid(1) = nxpoint
       dimsf_grid(2) = nzpoint
       dimsf_grid(3) = 1

       allocate(buffer_jplane(ntpoint_pfile, nxpoint, nzpoint, 1, num_dset))
!       allocate(buffer_total(ntpoint_pfile*num_file, nxpoint, nzpoint, 1, num_dset))


       allocate(buffer_jgrid(nxpoint, nzpoint, 1, num_dset_grid))   !!!!!!!!!!!!!!


       uave2d = 0.d0
       rhoave2d = 0.d0
       pave2d = 0.d0
       p2 = 0.d0


       do n=1, num_file
         write(unit=fnum_8,fmt='(I08.8)') (file_be + (n-1)*file_skip)
         fname = trim(filepath)//fnum_8//'.h5'

         CALL h5fopen_f(trim(fname), H5F_ACC_RDONLY_F, file_id, hdferr)
         CALL h5gopen_f(file_id, trim(gname(1)), group_id, hdferr)
         do i=1, 3

           CALL h5dopen_f(group_id, trim(dname(i)), dset_id, hdferr)
           CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, buffer_jplane(:,:,:,1,i),dimsf, hdferr)
           CALL h5dclose_f(dset_id, hdferr)
         enddo
         CALL h5gclose_f(group_id, hdferr)
         CALL h5fclose_f(file_id, hdferr)


         do k=1, nzpoint
           do i=1, nxpoint
             do nnn=1, ntpoint_pfile
               uave2d(k,i) = uave2d(k,i) + buffer_jplane(nnn, i, k, 1, 1)
               rhoave2d(k,i) = rhoave2d(k,i) + buffer_jplane(nnn, i, k, 1, 2)/(rbar*buffer_jplane(nnn, i, k, 1, 3))
               pave2d(k,i) = pave2d(k,i) + buffer_jplane(nnn, i, k, 1, 2)
               p2(k,i) = p2(k,i) + buffer_jplane(nnn, i, k, 1, 2)**2


             enddo
           enddo
         enddo
       enddo ! end n loop

       uave2d = uave2d/(ntpoint_pfile*num_file)
       rhoave2d = rhoave2d/(ntpoint_pfile*num_file)

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       pave2d = pave2d/(ntpoint_pfile*num_file)
       p2 = p2/(ntpoint_pfile*num_file)
       prms = sqrt(abs(pave2d - (p2)**2))


       ! read grid info
       fname = trim(filepath)//'GridMetrics.h5'

       CALL h5fopen_f(trim(fname), H5F_ACC_RDONLY_F, file_id, hdferr)
       CALL h5gopen_f(file_id, trim(gname(1)), group_id, hdferr)
       do i=1, 3
         CALL h5dopen_f(group_id, trim(dname_grid(i)), dset_id, hdferr)
         CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, buffer_jgrid(:,:,1,i),dimsf_grid, hdferr)
         CALL h5dclose_f(dset_id, hdferr)
       enddo
       CALL h5gclose_f(group_id, hdferr)
       CALL h5fclose_f(file_id, hdferr)


       call FinalizeHDF5()

       open(7,file='uave.dat', status='unknown')
       rewind(7)
       write(7,'(4a)') 'variables=x,z,uave'
       write(7,'("Zone T=Ave, I=", I4, ",J=", I4)') nxpoint, nzpoint
       do k=1, nzpoint
          do i=1, nxpoint
             write(7,*) buffer_jgrid(i,k,1,1), buffer_jgrid(i,k,1,2), uave2d(k,i)
           enddo
         enddo
       close(7)

       open(7,file='rhoave2d.dat', status='unknown')
       rewind(7)
       write(7,'(4a)') 'variables=x,z,rhoave'
       write(7,'("Zone T=Ave, I=", I4, ",J=", I4)') nxpoint, nzpoint
       do k=1, nzpoint
          do i=1, nxpoint
             write(7,*) buffer_jgrid(i,k,1,1), buffer_jgrid(i,k,1,2), rhoave2d(k,i)
           enddo
         enddo
       close(7)

       open(7,file='prms_jplane.dat', status='unknown')
       rewind(7)
       write(7,'(4a)') 'variables=x,z,rhoave'
       write(7,'("Zone T=Ave, I=", I4, ",J=", I4)') nxpoint, nzpoint
       do k=1, nzpoint
          do i=1, nxpoint
             write(7,*) buffer_jgrid(i,k,1,1), buffer_jgrid(i,k,1,2), prms(k,i)
           enddo
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
         read(*,*) ntpoint_pfile, nxpoint, nzpoint
!         read(*,*)
!         read(*,*) uinf, rhoinf



       end subroutine Init










     end program
