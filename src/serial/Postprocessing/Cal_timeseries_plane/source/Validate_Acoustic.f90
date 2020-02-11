


     program Validate_Acoustic
       use modRWHDF5
       use modCalInt
       use MTecplotIO
       implicit none

       character(400) :: filepath
       character(400) :: gridpath
       character(400) :: fname

       integer :: temp(1)
       real(8) :: uinf, rhoinf, rbar=287.0
       real(8) :: Tinf, Ma

       integer :: num_kinf
       integer, dimension(:), allocatable :: kinf
!       real(8), dimension(:,:), allocatable :: delta, theta, dstar, Cf
       real(8), dimension(:,:,:), allocatable :: buffer
       real(8), dimension(:,:,:), allocatable :: buffer_grid
       real(8), dimension(:,:,:), allocatable :: buffer_output
       real(8) :: muinf, muw

       integer :: i, j, k, kk, n
       integer :: kbl
       integer :: hdferr
       integer:: nxpoint, nzpoint, nxpoint_grid, nzpoint_grid
       integer :: ibe, iend, kbe, kend, ibe_grid, iend_grid, kbe_grid, kend_grid

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       integer :: TotZone
       type(Pointer3D), dimension(:), allocatable :: vars
       real(8), dimension(:,:,:,:), allocatable, target :: vartmp
       integer, dimension(:), allocatable :: varidx
       logical :: IsFileNew, IsFileClose
       integer :: ilen, jlen, klen, nvarout, numshare
       integer :: nvars
       integer, dimension(:), allocatable :: ShareVar
       real(8) :: time
       character(600) :: varoutname
       real(8) :: aa
       character(4) :: fnum

       type(tp_rdwt_hdf5) :: grd, fsol

       call InitHDF5()
       call Input()

       nxpoint = iend - ibe + 1
       nzpoint = kend - kbe + 1
       nxpoint_grid = iend_grid - ibe_grid + 1
       nzpoint_grid = kend_grid - kbe_grid + 1

       print *, '*************************************************************'
       print *, 'average file info'
       print *, 'nxpoint =       ', nxpoint,      '  nzpoint =      ', nzpoint
       print *, '*************************************************************'

       fsol%rank = 2
       fsol%dnum = 6
       fsol%gname = 'Stat2d'
       allocate(fsol%dname(fsol%dnum))
       allocate(fsol%dimsf(fsol%rank) )
       fsol%dname(1) = 'uave'
       fsol%dname(2) = 'ru'
       fsol%dname(3) = 'rhoave'
       fsol%dname(4) = 'dudk'
       fsol%dname(5) = 'tave'
       fsol%dname(6) = 'dtdk'

       grd%rank = 2
       grd%dnum = 3
       grd%gname = '/'
       allocate(grd%dname(grd%dnum))
       allocate(grd%dimsf(grd%rank))
       grd%dname(1) = 'x'
       grd%dname(2) = 'z'
       grd%dname(3) = 'dkdz'
!       grd%dname(3) = 'dzdz'

       ! nvarout = 7 ! 'x  delta  theta  dstar Cf, qw, Ch'
       nvarout = 12 ! 'x  delta  theta  dstar Cf qw Ch Cf_inf Ch_inf Retheta Redelta2 H'

       allocate(buffer(nzpoint_grid, nxpoint_grid, fsol%dnum))
       allocate(buffer_grid(nzpoint_grid, nxpoint_grid,grd%dnum))
!       allocate(delta(nxpoint_grid,num_kinf), theta(nxpoint_grid,num_kinf), dstar(nxpoint_grid,num_kinf), Cf(nxpoint_grid,num_kinf))
       allocate(buffer_output(nxpoint_grid,num_kinf,nvarout))

       fsol%dimsf = (/nzpoint_grid,nxpoint_grid/)
       fsol%IsHSInitialized = .true.

       grd%dimsf = (/nzpoint_grid,nxpoint_grid/)
       grd%IsHSInitialized = .true.

       ! read file
       fsol%fname = trim(filepath)
       print *, 'Reading file ... '
       print *, 'fname = ', trim(fsol%fname)

       call ReadHDF5_2D(fsol,buffer)

       ! read grid
       grd%fname = trim(gridpath)
       print *, 'Reading grid file ...'
       print *, 'fname = ', trim(grd%fname)

       call ReadHDF5_2D(grd,buffer_grid)
       call FinalizeHDF5()

       do kk=1, num_kinf
         do i=1, nxpoint_grid
!           call  Cal_Int(buffer_grid(1:nzpoint_grid,i,2), nzpoint_grid, kinf(kk), buffer(1:nzpoint_grid, i, 1), &
!                         buffer(1:nzpoint_grid, i, 2),buffer(1:nzpoint_grid, i, 3), buffer(1:nzpoint_grid, i, 4), &
!                         buffer(1:nzpoint_grid, i, 5),buffer_grid(1:nzpoint,i,3), delta(i,kk), theta(i,kk), dstar(i,kk), Cf(i,kk) )
           call Cal_Int(buffer_grid(1:nzpoint_grid,i,2:3), nzpoint_grid, kinf(kk), fsol%dnum, buffer(1:nzpoint_grid,i,1:fsol%dnum),nvarout,buffer_output(i,kk,2:nvarout), Tinf, Ma)

         enddo
       enddo

       TotZone = num_kinf  !!!!!!!!!
       allocate(varidx(nvarout))
       do i=1, nvarout
         varidx(i) = i
       enddo

       print *, 'nvarout = ', nvarout

       ilen = nxpoint_grid
       jlen = 1
       klen = 1

       nvars = nvarout
       numshare = 1

        allocate(vartmp(ilen,jlen,klen,nvars))
        allocate(vars(nvarout))
        allocate(ShareVar(nvarout))
        ShareVar = 0
        do n=1, nvarout
          vars(n)%pt => vartmp(1:ilen,1:jlen,1:klen,varidx(n))
          if(n.le.numshare) ShareVar(n) = 1
        enddo

        time = 0.

        muinf = (1.458E-6)*((tinf)**1.5)/(tinf+110.4)
        muw = (1.458E-6)*((buffer(1,1,5))**1.5)/(buffer(1,1,5)+110.4)

        do kk=1, num_kinf

          do k=1, klen
            do i=1, ilen
              vartmp(i,1,k,1) = buffer_grid(1,i,1)
              vartmp(i,1,k,2) = buffer_output(i,kk,2)  ! delta
              vartmp(i,1,k,3) = buffer_output(i,kk,3)  ! theta
              vartmp(i,1,k,4) = buffer_output(i,kk,4)  ! dstar
              vartmp(i,1,k,5) = buffer_output(i,kk,5)  ! Cf
              vartmp(i,1,k,6) = buffer_output(i,kk,6)  ! qw
              vartmp(i,1,k,7) = buffer_output(i,kk,7)  ! Ch
              vartmp(i,1,k,8) = buffer_output(i,kk,8)/(0.5*rhoinf*uinf**2) ! Cf_inf
              vartmp(i,1,k,9) = buffer_output(i,kk,9)/(rhoinf*uinf)        ! Ch_inf
              vartmp(i,1,k,10) = rhoinf*uinf*buffer_output(i,kk,3)/muinf   ! Retheta
              vartmp(i,1,k,11) = rhoinf*uinf*buffer_output(i,kk,3)/muw     ! Redelta2
              vartmp(i,1,k,12) = buffer_output(i,kk,4)/buffer_output(i,kk,3) ! H
            enddo
          enddo

          ! 'x  delta  theta  dstar Cf qw Ch Cf_inf Ch_inf'
          varoutname = " x  delta  theta  dstar Cf qw Ch Cf_inf Ch_inf Retheta Redelta2 H"

          print *, 'varoutname = ', trim(varoutname)
          time = time + 1
          if(kk.eq.1) then
            IsFileNew=.true.; IsFileClose=.false.
            call WriteTecBin('validate_acoustic.plt',trim(varoutname),ilen,jlen,klen,nvarout,time,vars,IsFileNew,IsFileClose)
          elseif(kk.lt.Totzone) then
            IsFileNew=.false.; IsFileClose=.false.
            call WriteTecBin('validate_acoustic.plt',trim(varoutname),ilen,jlen,klen,nvarout-numshare,time,vars(numshare+1:nvarout),IsFileNew,IsFileClose,ShareVar)
          else
            IsFileNew=.false.; IsFileClose=.true.
            call WriteTecBin('validate_acoustic.plt',trim(varoutname),ilen,jlen,klen,nvarout-numshare,time,vars(numshare+1:nvarout),IsFileNew,IsFileClose,ShareVar)
          endif

          write(unit=fnum,fmt='(I04.4)') kk
          open(7,file='x_theta'//fnum//'.dat',status='unknown')
            write(7,*) 'x, theta'
            do i=1, ilen
              write(7,*) buffer_grid(1,i,1), buffer_output(i,kk,3)
            enddo
          close(7)
          open(7,file='x_Cf'//fnum//'.dat',status='unknown')
            write(7,*) 'x, Cf'
            do i=1, ilen
              write(7,*) buffer_grid(1,i,1), vartmp(i,1,1,8)
            enddo
          close(7)


        enddo ! end kk loop

  contains

       subroutine Input()
         implicit none
         integer :: i

         read(*,*)
         read(*,*) uinf, rhoinf, Tinf
         read(*,*)
         read(*,'(a)') filepath
         read(*,*)
         read(*,'(a)') gridpath
         read(*,*)
         read(*,*) ibe, iend, kbe, kend
         read(*,*)
         read(*,*) num_kinf
         allocate(kinf(num_kinf))
         read(*,*) (kinf(i), i=1,num_kinf)

         Ma = uinf/(sqrt(1.4*287.0*Tinf))

         ibe_grid = ibe; iend_grid = iend
         kbe_grid = kbe; kend_grid = kend

       end subroutine Input






     end program Validate_Acoustic
