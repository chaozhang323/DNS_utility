
!> Using Averaged AveAcousticStat file to calculate the following variables.
!> "Uvd "
!> Andrew Trettel and Johan Larsson, Physics of fluids, 2016, "Mean velocity scaling for compressible wall turbulence with heat transfer"
!>
!> You-Sheng Zhang et, al. JFM, 2013, "A generalized Reynolds analogy for compressible wall-bounded turbulent flows"
!>
!> Rodney D. W. Bowersox, JFM, 2009, "Extension of equilibrium turbulent heat flux models to high-speed shear flows"
!>
!> Jonathan Poggie, AIAA 2015-1983, "Compressible turbulent boundary layer simulations: resolution effects and turbulence modeling"
!>
!> Bin Wu et, al. Journal of Turbulence, "On the invariant mean velocity profile for compressible turbulent boundary layer"

     program Cal_Uvd
       use modRWHDF5
       use modCalInt
       use MTecplotIO
       implicit none

       real(8), parameter :: R=8314.3D0
       integer, parameter :: nvar_dat = 7 ! z, rho, u, v, w, t, p
       character(400) :: gridpath, filepath, datapath
       character(400) :: fname

       integer :: num_kinf
       integer, dimension(:), allocatable :: kinf
       real(8), dimension(:,:,:), allocatable :: buffer
       real(8), dimension(:,:,:), allocatable :: buffer_grid
       real(8), dimension(:,:,:,:), allocatable :: buffer_output
       real(8), dimension(:,:), allocatable :: buffer_dat
       real(8), dimension(:,:,:), allocatable :: buffer_output_dat
       real(8) :: muinf, muw

       integer :: i, j, k, kk, n, nvarout
       integer :: nxpoint, nzpoint, nxpoint_grid, nzpoint_grid
       integer :: ibe, iend, kbe, kend, ibe_grid, iend_grid, kbe_grid, kend_grid

       type(tp_rdwt_hdf5) :: grd, fsol
       integer :: iHDF5, iDat, nskip


       call InitHDF5()
       call Input()

       if(iHDF5.gt.0) then
         call InitReadInfo()
                       !  1  2    3     4      5     6         7          8           9
         nvarout = 39+5  ! 'x, z, utau, ztau, Retau, Retau_star, Uvd, Uvd_Trettel, zplus_Trettel,
                       ! Prt, Prt_ave, Pre, Prt_Bowersox, tau_xz, tau_zz, theta1, tehta2, theta3,
                       ! x2_smooth, delta, Walz, Walz_g, u_ue, T_Te, ZII_plus, Uvd_Wu, SRA_Huang, SRA_Zhang, tau, tau_xz_2, tau_total, ustar, ystar, dkdz'

         allocate(buffer(nzpoint_grid, nxpoint_grid, fsol%dnum))
         allocate(buffer_grid(nzpoint_grid, nxpoint_grid,grd%dnum))
         allocate(buffer_output(nzpoint_grid,nxpoint_grid,num_kinf,nvarout))
         call ReadFile()

         !! calculate Int variables
         do kk=1, num_kinf
           do i=1, nxpoint_grid ! not include the inflow boundary
             call Cal_Int(buffer_grid(1:nzpoint_grid,i,2:3), nzpoint_grid, kinf(kk), fsol%dnum, buffer(1:nzpoint_grid,i,1:fsol%dnum),nvarout,buffer_output(1:nzpoint_grid,i,kk,3:nvarout))
           enddo
         enddo

         do k=kbe, kend
           do i=ibe, iend
             buffer_output(k,i,:,1) = buffer_grid(1,i,1)
             buffer_output(k,i,:,2) = buffer_grid(k,1,2)
           enddo
         enddo

         call WriteFile()
       elseif(iDat.gt.0) then
         nvarout = 14  ! 'z, utau, ztau, Retau, Retau_star, Uvd, Uvd_Trettel, zplus_Trettel, Pre, delta, Walz, Walz_g, u_ue, T_Te'

         allocate(buffer_dat(nzpoint_grid,nvar_dat))
         allocate(buffer_output_dat(nzpoint_grid,num_kinf,nvarout))
         call ReadDat(trim(datapath),nskip,nzpoint_grid,nvar_dat,buffer_dat)

         do kk=1, num_kinf
           call Cal_Int_dat(nzpoint_grid,kinf(kk),nvar_dat,buffer_dat(1:nzpoint_grid,1:nvar_dat),nvarout,buffer_output_dat(1:nzpoint_grid,kk,1:nvarout) )
         enddo ! end kk loop
         call WriteFile_dat(buffer_output_dat)
       endif

       call FinalizeHDF5()

  contains

       subroutine ReadDat(fn,np,nz,nvar,bufferout)
         character(*), intent(in) :: fn
         integer, intent(in) :: np, nz, nvar
         real(8), intent(out) :: bufferout(nz,nvar)
         integer :: k, n

         open(7,file=trim(fn),status='unknown')
           do n=1, np
             read(7,*)
           enddo
           do k=1, nz
             read(7,*) bufferout(k,1:nvar)
           enddo
         close(7)

         ! for test
         print *, 'bufferout(nz,1:nvar) = ', bufferout(nz,1:nvar)

       end subroutine ReadDat

       subroutine WriteFile_dat(bufferin)
         real(8), intent(in) :: bufferin(:,:,:)
         integer :: i, k, kk
         character(400) :: fn

         fn = 'output.dat'
         print *, 'writting file ', trim(fn)
         open(11,file=trim(fn),status='unknown')
           rewind(11)
           write(11,'(4a)') 'variables=z, utau, ztau, Retau, Retau_star, Uvd, Uvd_Trettel, zplus_Trettel, Pre, delta, Walz, Walz_g, u_uinf, T_Tinf'
           do kk=1, num_kinf
             write(11,'("Zone T=kref_",I04.4, ",J=",I4 )') kk, nzpoint_grid
             do k=1, nzpoint_grid
               !do i=1, nxpoint_grid
                 write(11,*) bufferin(k,kk,1:nvarout)
               !enddo
             enddo
           enddo ! end kk loop

         close(11)

       end subroutine WriteFile_dat

       subroutine WriteFile()
         integer :: i, k, kk
         character(400) :: fn

         fn = 'output.dat'
         print *, 'writting file ', trim(fn)
         open(11,file=trim(fn),status='unknown')
           rewind(11)
           write(11,'(4a)') 'variables=x, z, utau, ztau, Retau, Retau_star, Uvd, Uvd_Trettel, zplus_Trettel, Prt, Prt_ave, Pre, Prt_Bowersox, tau_xz, tau_zz, theta1, theta2, theta3, x2_smooth, delta, Walz, Walz_g, u_uinf, T_Tinf, ZII_plus, Uvd_Wu, SRA_Huang, SRA_Zhang, tau,tau_xz_2, tau_total,dkdz,Retaust_Patel,zstar,ustar_patel,l12,l12_star,lm,lm_star, theta, dstar, H, Retheta, Redelta2'
           do kk=1, num_kinf
             write(11,'("Zone T=kref_",I04.4, ",I=",I4, ",J=",I4 )') kk, nxpoint_grid, nzpoint_grid
             do k=1, nzpoint_grid
               do i=1, nxpoint_grid
                 write(11,*) buffer_output(k,i,kk,1:nvarout)
               enddo
             enddo
           enddo ! end kk loop

         close(11)

       end subroutine WriteFile

       subroutine ReadFile()
         implicit none

         ! read data file
         fsol%fname = trim(filepath)
         print *, 'Reading data file ... '
         print *, 'fname = ', trim(fsol%fname)
         call ReadHDF5_2D(fsol,buffer)

         ! read grid
         grd%fname = trim(gridpath)
         print *, 'Reading grid file ...'
         print *, 'fname = ', trim(grd%fname)

         call ReadHDF5_2D(grd,buffer_grid)

       end subroutine ReadFile

       !! select variables to read
       subroutine InitReadInfo()
         implicit none

         fsol%rank = 2
         fsol%dnum = 17
         fsol%gname = 'Stat2d'
         allocate(fsol%dname(fsol%dnum))
         allocate(fsol%dimsf(fsol%rank) )
         fsol%dname(1)  = 'uave'
         fsol%dname(2)  = 'ru'
         fsol%dname(3)  = 'rhoave'
         fsol%dname(4)  = 'dudk'
         fsol%dname(5)  = 'tave'
         fsol%dname(6)  = 'dtdk'
         fsol%dname(7)  = 'divave'

         fsol%dname(8)  = 'wave'
         fsol%dname(9)  = 'ruw'
         fsol%dname(10) = 'rw'
         fsol%dname(11) = 'rwt'
         fsol%dname(12) = 'pave'
         fsol%dname(13) = 'uw'
         fsol%dname(14) = 'w2'
         fsol%dname(15) = 'u2'
         fsol%dname(16) = 't2'
         fsol%dname(17) = 't0ave'


         grd%rank = 2
         grd%dnum = 3
         grd%gname = '/'
         allocate(grd%dname(grd%dnum))
         allocate(grd%dimsf(grd%rank))
         grd%dname(1) = 'x'
         grd%dname(2) = 'z'
         grd%dname(3) = 'dkdz'
         ! grd%dname(3) = 'dzdz'

         fsol%dimsf = (/nzpoint_grid,nxpoint_grid/)
         fsol%IsHSInitialized = .true.

         grd%dimsf = (/nzpoint_grid,nxpoint_grid/)
         grd%IsHSInitialized = .true.

       end subroutine InitReadInfo

       subroutine Input()
         implicit none
         integer :: i

         print *, 'reading parameters ... '
         read(*,*)
         read(*,*) uinf, rhoinf, Tinf, iCalMu, zplus_ref, Pr
         read(*,*)
         read(*,*) Rm
         read(*,*)
         read(*,*)
         read(*,*) iHDF5
         if(iHDF5.gt.0) then
           read(*,*)
           read(*,'(a)') gridpath
           read(*,*)
           read(*,'(a)') filepath
           read(*,*)
           read(*,*) ibe, iend, kbe, kend
           read(*,*)
           read(*,*) nxpoint_grid, nzpoint_grid
           read(*,*)
           read(*,*) num_kinf
           allocate(kinf(num_kinf))
           read(*,*) (kinf(i), i=1,num_kinf)
         else
           do i=1, 9
             read(*,*)
           enddo
         endif
         read(*,*)
         read(*,*)
         read(*,*) iDat, nskip
         if(iDat.gt.0) then
           read(*,*)
           read(*,'(a)') datapath
           read(*,*)
           read(*,*) kbe, kend
           read(*,*)
           read(*,*) num_kinf
           allocate(kinf(num_kinf))
           read(*,*) (kinf(i), i=1,num_kinf)
           ibe = 0; iend = 0
         endif

         !ibe_grid = ibe; iend_grid = iend
         !kbe_grid = kbe; kend_grid = kend

         nxpoint = iend - ibe + 1
         nzpoint = kend - kbe + 1
         !nxpoint_grid = iend_grid - ibe_grid + 1
         !nzpoint_grid = kend_grid - kbe_grid + 1

         print *, '*************************************************************'
         print *, 'reading file info'
         print *, 'nxpoint =       ', nxpoint,      '  nzpoint =      ', nzpoint
         print *, '*************************************************************'

         rbar = R/Rm

         Cp = gamma*rbar/(gamma-1.d0)
         Cv = 2.5*rbar
         print *, 'rbar = ', rbar

       end subroutine Input


     end program Cal_Uvd

       !integer :: temp(1)
       !real(8) :: uinf, rhoinf, rbar=287.0
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !integer :: TotZone
       !type(Pointer3D), dimension(:), allocatable :: vars
       !real(8), dimension(:,:,:,:), allocatable, target :: vartmp
       !integer, dimension(:), allocatable :: varidx
       !logical :: IsFileNew, IsFileClose
       !integer :: ilen, jlen, klen, , numshare
       !integer :: nvars
       !integer, dimension(:), allocatable :: ShareVar
       !real(8) :: time
       !character(600) :: varoutname
       !character(4) :: fnum
