


     program Validate_Acoustic
       use HDF5
       use modCalInt
       use MTecplotIO
       implicit none

       character(400) :: filepath
       character(400) :: gridpath
       character(400) :: fname

       character(10), dimension(:), allocatable :: dname
       character(10), dimension(:), allocatable :: dname_grid
       character(10) :: gname
       character(10) :: gname_grid

       integer(HID_T) :: file_id
       integer(HID_T) :: group_id, dset_id, dspace_id
       integer(HID_T) :: memspace_id

       integer(HSIZE_T), dimension(:), allocatable :: dimsf
       integer(HSIZE_T), dimension(:), allocatable :: dimsf_total

       integer(HSIZE_T), dimension(:), allocatable :: dimsm
       integer(HSIZE_T), dimension(:), allocatable :: count, offset, block, stride

       integer(HSIZE_T), dimension(:), allocatable :: dimsf_grid
       integer(HSIZE_T), dimension(:), allocatable :: dimsf_total_grid

       integer :: rank, rank_grid
       integer :: num_dset, num_dset_grid  !, num_group

       integer :: temp(1)
       real(8) :: uinf, rhoinf, rbar=287.0

       integer :: num_kinf
       integer, dimension(:), allocatable :: kinf
       real(8), dimension(:,:), allocatable :: delta, theta, dstar, Retau       !!!!!!!!!!
       real(8), dimension(:,:,:), allocatable :: buffer
       real(8), dimension(:,:,:), allocatable :: buffer_grid  !!!!!!!!!!

       integer :: i, j, k, kk, n
       integer :: kbl
       integer :: hdferr
       integer:: nxpoint, nzpoint, nxpoint_grid, nzpoint_grid
       integer :: ibe, iend, kbe, kend, ibe_grid, iend_grid, kbe_grid, kend_grid

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       integer :: TotZone
       type(Pointer3D), dimension(:), allocatable :: vars
       real(8), dimension(:,:,:,:), allocatable, target :: vartmp !!!!!!
       integer, dimension(:), allocatable :: varidx
       logical :: IsFileNew, IsFileClose
       integer :: ilen, jlen, klen, nvarout, numshare
       integer :: nvars
       integer, dimension(:), allocatable :: ShareVar  !!!!!!!!!!!!
       real(8) :: time
       character(600) :: varoutname
       real(8) :: aa

       call InitHDF5()

       call Input()

       nxpoint = iend - ibe + 1
       nzpoint = kend - kbe + 1
       nxpoint_grid = iend_grid - ibe_grid + 1
       nzpoint_grid = kend_grid - kbe_grid + 1

       print *, '*************************************************************'
       print *, 'average file info'
       print *, 'nxpoint =       ', nxpoint,      '  nzpoint =      ', nzpoint
       print *, 'jplane grid file info'
       print *, 'nxpoint_grid =  ', nxpoint_grid, '  nzpoint_grid = ', nzpoint_grid
       print *, '*************************************************************'


       rank = 2
       rank_grid = 2
       num_dset = 4
       num_dset_grid = 2
       gname_grid = "/"

       allocate(dname(num_dset))
       allocate(dname_grid(num_dset_grid))
       allocate(dimsf(rank), dimsf_total(rank), dimsm(rank))
       allocate(count(rank), offset(rank), block(rank), stride(rank))
       allocate(dimsf_grid(rank_grid))

       allocate(buffer(nzpoint_grid, nxpoint_grid, num_dset))
       allocate(buffer_grid(nzpoint_grid, nxpoint_grid,num_dset_grid))
       allocate(delta(nxpoint_grid,num_kinf), theta(nxpoint_grid,num_kinf), dstar(nxpoint_grid,num_kinf), Retau(nxpoint_grid,num_kinf))

       dname(1) = "uave"
       dname(2) = "ru"
       dname(3) = "rhoave"
       dname(4) = "utau"

       dname_grid(1) = "x"
       dname_grid(2) = "z"
!       dname_grid(3) = "dkdz"

       dimsf(1) = nzpoint_grid !!!!!!
       dimsf(2) = nxpoint_grid !!!!!!

       dimsm(1) = dimsf(1)
       dimsm(2) = dimsf(2)

       dimsf_grid(1) = nzpoint_grid
       dimsf_grid(2) = nxpoint_grid
!       dimsf_grid(3) = 1

       block(1) = dimsf(1)
       block(2) = dimsf(2)
       offset(2) = ibe_grid - ibe
       offset(1) = kbe_grid - kbe
       stride(1) = 1
       stride(2) = 1
       count(1) = 1
       count(2) = 1

       ! read file
       fname = trim(filepath)
       print *, 'Reading file ... '
       print *, 'fname = ', trim(fname)

       CALL h5fopen_f(trim(fname), H5F_ACC_RDONLY_F, file_id, hdferr)
       CALL h5gopen_f(file_id, trim(adjustl(gname)), group_id, hdferr)
       CALL h5screate_simple_f(rank, dimsm, memspace_id, hdferr)
       do i=1, num_dset
         CALL h5dopen_f(group_id, trim(dname(i)), dset_id, hdferr)
         CALL h5dget_space_f(dset_id, dspace_id, hdferr)
         CALL h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, &
                                    count, hdferr, stride, block)
         CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, buffer(:,:,i),dimsf, hdferr, &
                        file_space_id=dspace_id, mem_space_id=memspace_id)
         CALL h5sclose_f(dspace_id, hdferr)
         CALL h5dclose_f(dset_id, hdferr)
       enddo
       CALL h5gclose_f(group_id, hdferr)
       CALL h5fclose_f(file_id, hdferr)

       ! read grid
       fname = trim(gridpath)
       print *, 'Reading grid file ...'
       print *, 'fname = ', trim(fname)

       CALL h5fopen_f(trim(fname), H5F_ACC_RDONLY_F, file_id, hdferr)
       CALL h5gopen_f(file_id, trim(gname_grid), group_id, hdferr)
       do i=1, num_dset_grid
         CALL h5dopen_f(group_id, trim(dname_grid(i)), dset_id, hdferr)
         CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, buffer_grid(:,:,i),dimsf_grid, hdferr)
         CALL h5dclose_f(dset_id, hdferr)
       enddo
       CALL h5gclose_f(group_id, hdferr)
       CALL h5fclose_f(file_id, hdferr)

       call FinalizeHDF5()

       do kk=1, num_kinf
         do i=1, nxpoint_grid

           call  Cal_Int(buffer_grid(1:nzpoint_grid,i,2), nzpoint_grid, kinf(kk), buffer(1:nzpoint_grid, i, 1), buffer(1:nzpoint_grid, i, 2),buffer(1:nzpoint_grid, i, 3), buffer(1:nzpoint_grid, i, 4), delta(i,kk), theta(i,kk), dstar(i,kk), Retau(i,kk) )

         enddo
       enddo


       TotZone = num_kinf  !!!!!!!!!
       nvarout = 5       !!!!!
       allocate(varidx(nvarout))
       do i=1, nvarout
         varidx(i) = i
       enddo

       print *, 'nvarout = ', nvarout

       ilen = nxpoint_grid
       jlen = 1
       klen = 1  !!!!!!

       nvars = nvarout
       numshare = 1   !!!!!!!!

        allocate(vartmp(ilen,jlen,klen,nvars))
        allocate(vars(nvarout))
        allocate(ShareVar(nvarout))
        ShareVar = 0
        do n=1, nvarout
          vars(n)%pt => vartmp(1:ilen,1:jlen,1:klen,varidx(n))
          if(n.le.numshare) ShareVar(n) = 1                      !!!!!!!!!
        enddo


        time = 0.

        do kk=1, num_kinf

          do k=1, klen
            do i=1, ilen
              vartmp(i,1,k,1) = buffer_grid(1,i,1)
              vartmp(i,1,k,2) = delta(i,kk)
              vartmp(i,1,k,3) = theta(i,kk)
              vartmp(i,1,k,4) = dstar(i,kk)
              vartmp(i,1,k,5) = Retau(i,kk)
            enddo
          enddo

          varoutname = " x  delta  theta  dstar Retau "

          print *, 'varoutname = ', trim(varoutname)
          time = time + 1
          if(kk.eq.1) then
            IsFileNew=.true.; IsFileClose=.false.
            call WriteTecBin('validate_acoustic.dat',trim(varoutname),ilen,jlen,klen,nvarout,time,vars,IsFileNew,IsFileClose)
          elseif(kk.lt.Totzone) then
            IsFileNew=.false.; IsFileClose=.false.
            call WriteTecBin('validate_acoustic.dat',trim(varoutname),ilen,jlen,klen,nvarout-numshare,time,vars(numshare+1:nvarout),IsFileNew,IsFileClose,ShareVar)
          else
            IsFileNew=.false.; IsFileClose=.true.
            call WriteTecBin('validate_acoustic.dat',trim(varoutname),ilen,jlen,klen,nvarout-numshare,time,vars(numshare+1:nvarout),IsFileNew,IsFileClose,ShareVar)
          endif


        enddo ! end kk loop







  contains

       subroutine InitHDF5()
        ! Initialize HDF5 library and Fortran interfaces.
         CALL h5open_f(hdferr)
       end subroutine InitHDF5

       subroutine FinalizeHDF5()
         ! Close FORTRAN interfaces and HDF5 library.
         CALL h5close_f(hdferr)
       end subroutine FinalizeHDF5


       subroutine Input()
         implicit none
         integer :: i

         read(*,*)
         read(*,*) uinf, rhoinf
         read(*,*)
         read(*,'(a)') filepath
         read(*,*)
         read(*,'(a)') gridpath
         read(*,*)
         read(*,*) ibe, iend, kbe, kend
         read(*,*)
         read(*,*) ibe_grid, iend_grid, kbe_grid, kend_grid
         read(*,*)
         read(*,*) num_kinf
         allocate(kinf(num_kinf))
         read(*,*) (kinf(i), i=1,num_kinf)
         read(*,*)
         read(*,'(a)') gname


       end subroutine Input






     end program Validate_Acoustic
