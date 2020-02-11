program CalStatistics
! This code is used for re-calculating boundary layer statistics given the 2D Acoustic HDF5 file:
! eg. AveAcousticStat_timeave00056000-00094000_Stat2d.h5
! Write out new 2D acoustic file containing the new delta
use modRWHDF5
use modTecbin
use MFLOW
implicit none

real(8), parameter :: R = 8314.3d0, gamma=1.4, Pr=0.71, recover = 0.9
! real(8), parameter :: CSutherland = 1.458d-6, TmpSutherland=110.4d0
real(8) ::  rbar, pi, cp, Ma, Mw

! input parameter
integer :: icylin,iwall,iflow,ivis
real(8) :: uinf,tinf,zfreestream
character(400) :: flowfile,grdfile,fileout

integer :: imax,kmax,kw,kinf,ksp,num_append1d,num_append2d
character(100) :: fn_grp1d,varname_append1d,varname_append2d

real(8),dimension(:,:),allocatable :: buffer_rd1d
real(8),dimension(:,:,:),allocatable :: buffer_grd,buffer_rd2d
real(8),dimension(:,:),allocatable :: buffer_stat1d
real(8),dimension(:,:,:),allocatable :: buffer_stat2d

call InitHDF5()
! initialize the input
call Input()
! read the acoustic grid
call InitGrid()
! read the subset of the acoustic data
call ReadAcoustic_Subset()
! calculate the statistics
call DoStatistics()
! append the new calculated statistics to the original HDF5, write new plt file
call WriteStatistics()

! re-write acoustic2D.h5 data file and plt visualization file
!call ReWriteAcoustic()
call FinalizeHDF5()
if(allocated(buffer_grd)) deallocate(buffer_grd)

contains
    ! provided uave,ru,rhoave,tave,utau, calculate delta,dstar,theta,x_tilda,xd_tilda
    ! x: buffer_grd(k,i,1)
    ! z: buffer_grd(k,i,2)
    ! uave:   buffer_rd2d(k,i,1)
    ! ru:     buffer_rd2d(k,i,2)
    ! rhoave: buffer_rd2d(k,i,3)
    ! tave:   buffer_rd2d(k,i,4)
    ! wave:   buffer_rd2d(k,i,5)
    ! uw:     buffer_rd2d(k,i,6)
    ! utau:   buffer_rd1d(i,1)
    ! tauw:   buffer_rd1d(i,2)
    ! buffer_stat1d(i,:): delta,dstar_c,theta_c,delta_i,theta_i,x_tilda,xd_tilda
    ! buffer_uvd(k,i):  uvd
    subroutine DoStatistics()
        real(8) :: taw, twall,muw,muinf, uinf_local,ufs,rufs,zfs, rhow, rhoinf
        real(8) :: delta,dstar_c,theta_c,dstar_i,theta_i,uw_max
        integer :: i,j,k,kbl,kfs, kinf_local
        integer,dimension(1) :: ktmp
        real(8),dimension(:),allocatable :: vdstar_c,vtheta_c,vdstar_i,vtheta_i,x_tilda,xd_tilda,uvand,radius,uw_rs

        !---------------------
        ! determin wall type
        if(icylin .eq. 1) iwall = 2
        kw = 1; kinf = kmax; ksp = 1
        if(iwall .eq. 2) then  ! top wall
            kw = kmax; kinf = 1; ksp = -1
        elseif(iwall .eq. 3) then
            kinf = kmax/2 + mod(kmax,2)
        endif
        kinf_local = kinf - 4*ksp

        num_append1d = 9; num_append2d = 1
        varname_append1d = "delta_c dstar_c theta_c Redelta_c Redelta2_c Redstar_c Retau_c Retheta_c Retau_star"
!        varname_append1d = "delta_c dstar_c theta_c dstar_i theta_i"! x_tilda xd_tilda uw_max"
!        varname_append2d = "uvd_c"

        call InitViscosity(ivis)
!        twall=buffer_rd2d(kw,1,4)      ! wall temperature
!        muw = calMu(twall)             ! viscosity at wall
!        muinf = calMu(tinf)            ! viscosity at freestream
!        taw = tinf*(1.d0+recover*(gamma-1)/2.0*Ma**2)  ! adiabatic wall temperature

        zfs = zfreestream

        allocate(vdstar_c(kmax),vtheta_c(kmax),vdstar_i(kmax),vtheta_i(kmax),uvand(kmax),radius(kmax),uw_rs(kmax))
        allocate(x_tilda(imax),xd_tilda(imax))
        allocate(buffer_stat1d(imax,num_append1d))
        allocate(buffer_stat2d(kmax,imax,num_append2d))

        x_tilda = 0.d0; xd_tilda = 0.d0
        buffer_stat1d = 0.d0
        do i = 1,imax
            ! calculate delta,dstar,theta
            !----------------------------
            if(i.eq.1) print*,'Calculating delta,dstar,theta ...'
            uinf_local = buffer_rd2d(kinf_local,i,1)  ! local u_inf
            ktmp = minloc( abs(buffer_rd2d(kw:kinf:ksp,i,1)/uinf_local-0.99) )    ! 0.99*uinf
            kbl = kw+(ktmp(1)-1)*ksp
            delta = abs( buffer_grd(kbl,i,2)-buffer_grd(kw,i,2) )   ! delta
            ktmp = minloc( abs(abs(buffer_grd(kw:kinf:ksp,i,2)-buffer_grd(kw,i,2))-zfs*delta) )  ! freestream location zfs*delta
            kfs = kw+(ktmp(1)-1)*ksp
            ufs = buffer_rd2d(kfs,i,1)     ! u_{freestream}
            rufs = buffer_rd2d(kfs,i,2)    ! rho*u_{freestream}
!            if(i.eq.imax/2) print*,'i, kfs, uinf, delta, rhou:',i,kfs,uinf_local, delta,rufs
            ! dstar and theta
            radius = 1.0
            if(icylin.eq.1) then
                radius = buffer_grd(:,i,2)/buffer_grd(kw,i,2)
            endif
            do k = 1,kmax
                vdstar_c(k) = radius(k)*(1.-buffer_rd2d(k,i,2)/rufs)
!                vdstar_i(k) = radius(k)*(1.-buffer_rd2d(k,i,1)/buffer_rd2d(kfs,i,1))
                vtheta_c(k) = radius(k)*buffer_rd2d(k,i,2)*(1.-buffer_rd2d(k,i,1)/buffer_rd2d(kfs,i,1))/rufs
!                vtheta_i(k) = radius(k)*buffer_rd2d(k,i,1)*(1.-buffer_rd2d(k,i,1)/buffer_rd2d(kfs,i,1))/buffer_rd2d(kfs,i,1)
            enddo
            if(iwall.eq.2) then     ! upper wall (k=kmax)
                call trapezoid_integration(kmax,vdstar_c,buffer_grd(:,i,2),kfs,kw,dstar_c)
                call trapezoid_integration(kmax,vtheta_c,buffer_grd(:,i,2),kfs,kw,theta_c)
!                call trapezoid_integration(kmax,vdstar_i,buffer_grd(:,i,2),kfs,kw,dstar_i)
!                call trapezoid_integration(kmax,vtheta_i,buffer_grd(:,i,2),kfs,kw,theta_i)
            else                    ! lower wall (k=1)
                call trapezoid_integration(kmax,vdstar_c,buffer_grd(:,i,2),kw,kfs,dstar_c)
                call trapezoid_integration(kmax,vtheta_c,buffer_grd(:,i,2),kw,kfs,theta_c)
!                call trapezoid_integration(kmax,vdstar_i,buffer_grd(:,i,2),kw,kfs,dstar_i)
!                call trapezoid_integration(kmax,vtheta_i,buffer_grd(:,i,2),kw,kfs,theta_i)
            endif
            buffer_stat1d(i,1) = delta
            buffer_stat1d(i,2) = dstar_c
            buffer_stat1d(i,3) = theta_c
!            buffer_stat1d(i,4) = dstar_i
!            buffer_stat1d(i,5) = theta_i

            muw = calMu(buffer_rd2d(kw,i,4))             ! viscosity at wall
            muinf = calMu(buffer_rd2d(kinf,i,4))            ! viscosity at freestream
            rhow = buffer_rd2d(kw,i,3)  ! density at wall
            rhoinf = buffer_rd2d(kinf,i,3)  ! density at freestream
            ! Redelta_c Redelta2_c Redstar Retau Retheta Retau_star
            buffer_stat1d(i,4) = rhoinf*ufs*delta/muinf     ! Redelta
            buffer_stat1d(i,5) = rhoinf*ufs*theta_c/muw     ! Redelta2
            buffer_stat1d(i,6) = rhoinf*ufs*dstar_c/muinf   ! Redstar
            buffer_stat1d(i,7) = rhow*buffer_rd1d(i,1)*delta/muw   ! Retau
            buffer_stat1d(i,8) = rhoinf*ufs*theta_c/muinf   ! Retheta
            buffer_stat1d(i,9) = rhoinf*sqrt(buffer_rd1d(i,2)/rhoinf)*delta/muinf   ! Retau_star

            ! calculate x_tilda,xd_tilda
            !---------------------------
!            if(num_append1d.gt.5) then
!                if(i.eq.1) print*,'Calculating x_tilda,xd_tilda ...'
!                if(i.ge.2) then
!                    x_tilda(i) = (buffer_grd(kw,i,1)-buffer_grd(kw,i-1,1))*buffer_rd1d(i,1)/delta/ufs   ! dx*utau/delta/ufs
!                    xd_tilda(i) = (buffer_grd(kw,i,1)-buffer_grd(kw,i-1,1))/delta                   ! dx/delta
!                    buffer_stat1d(i,6) = sum(x_tilda(1:i))
!                    buffer_stat1d(i,7) = sum(xd_tilda(1:i))
!                endif
!            endif

            ! calculate maximum uw reynolds stress
!            if(num_append1d.gt.7) then
!                forall(k=1:kmax) uw_rs(k) = buffer_rd2d(k,i,1)*buffer_rd2d(k,i,5)-buffer_rd2d(k,i,6)
!                buffer_stat1d(i,8) = -maxval(abs(uw_rs))
!            endif

            ! van Driest transformed velocity
            !--------------------------------
!            if(i.eq.1) print*,'Calculating van Driest transformed velocity ...'
!            uvand = 0.
!            do k = kw+ksp,kinf,ksp
!                uvand(k) = sqrt(buffer_rd2d(kw,i,4)/buffer_rd2d(k,i,4))*(buffer_rd2d(k,i,1)-buffer_rd2d(k-ksp,i,1))
!                buffer_stat2d(k,i,1) = sum(uvand(kw:k:ksp))
!            enddo
!            forall(k=1:kmax) buffer_stat2d(k,i,1) = buffer_stat2d(k,i,1)/buffer_rd1d(i,1)
        enddo

        deallocate(vdstar_c,vtheta_c,vdstar_i,vtheta_i,uvand,radius,x_tilda,xd_tilda)
        deallocate(buffer_rd1d,buffer_rd2d)

    end subroutine DoStatistics

    ! **** read in subset of 2D acoustic data (uave,ru,rhoave,tave,utau) ****
    subroutine ReadAcoustic_Subset()
        type(tp_rdwt_hdf5) :: fsol2d,fsol1d
        integer :: i,num_rd2d,num_rd1d

        num_rd2d = 6
        fsol2d%fname = trim(flowfile)
!        fsol2d%fname = trim(flowfile)//'_Stat2d.h5'
        fsol2d%gname = '/Stat2d'
        call DetectHDF5(fsol2d)
        !     print *, 'fsol2d%gname = ', trim(fsol2d%gname)
        !     print *, 'fsol2d%rank  = ', fsol2d%rank
        !     print *, 'fsol2d%dimsf = ', fsol2d%dimsf
        !     print *, 'fsol2d%dnum  = ', fsol2d%dnum
        !     print *, 'fsol2d%dname = ', (trim(fsol2d%dname(i))//', ',i=1,fsol2d%dnum)
        !     print *, '*********************************************'
        if(fsol2d%dimsf(1) .ne. kmax .or. fsol2d%dimsf(2) .ne. imax) then
            print*,'Dimension inconsistent between grid file and flow file... STOP!'
            print*,'Dimension in grid file:',imax,kmax
            print*,'Dimension in flow file:',fsol2d%dimsf(2),fsol2d%dimsf(1)
            stop
        endif
        fsol2d%rank = 2
        fsol2d%dnum = num_rd2d
        if(allocated(fsol2d%dname)) deallocate(fsol2d%dname)
        if(allocated(fsol2d%dimsf)) deallocate(fsol2d%dimsf)
        allocate(fsol2d%dname(fsol2d%dnum))
        allocate(fsol2d%dimsf(fsol2d%rank) )
        fsol2d%dimsf=(/kmax,imax/)
        fsol2d%dname(1) = 'uave'
        fsol2d%dname(2) = 'ru'
        fsol2d%dname(3) = 'rhoave'
        fsol2d%dname(4) = 'tave'
        fsol2d%dname(5) = 'wave'
        fsol2d%dname(6) = 'uw'
        fsol2d%IsHSInitialized = .true.
        allocate(buffer_rd2d(kmax,imax,fsol2d%dnum))
        print*,'Reading ',('"'//trim(fsol2d%dname(i))//'" ',i=1,fsol2d%dnum)
!        print*,'  from file: ',trim(fsol2d%fname)
        call ReadHDF5_2D(fsol2d,buffer_rd2d)

        fsol1d%fname = trim(fsol2d%fname)
!        fsol1d%gname = '/Int_BotWall.h5'
!        if(iwall.eq.2) '/Int_TopWall.h5'     ! top wall case

        fsol1d%gname = trim(fn_grp1d)
        call DetectHDF5(fsol1d)
!        print *, 'fsol1d%gname = ', trim(fsol1d%gname)
!        print *, 'fsol1d%dnum  = ', fsol1d%dnum
!        print *, 'fsol1d%dname = ', (trim(fsol1d%dname(i))//', ',i=1,fsol1d%dnum)
!        print *, 'fsol1d%rank  = ', fsol1d%rank
!        print *, 'fsol1d%dimsf = ', fsol1d%dimsf
!        print *, '*********************************************'
        fsol1d%dnum = 2
        if(allocated(fsol1d%dname)) deallocate(fsol1d%dname)
        allocate(fsol1d%dname(fsol1d%dnum))
        fsol1d%dname(1) = 'utau'
        fsol1d%dname(2) = 'tauw'
        allocate(buffer_rd1d(imax,fsol1d%dnum))
        print*,'Reading ',('"'//trim(fsol1d%dname(i))//'" ',i=1,fsol1d%dnum)
        call ReadHDF5_1D(fsol1d,buffer_rd1d)

        ! debug
!        open(33,file='check_AcousticRD.dat',status='unknown')
!        write(33,*) 'Variables = x z uave ru rhoave tave utau'
!        write(33,*) 'Zone I = ',imax,', J = ',kmax,', K = 1, F=POINT'
!        do k = 1,kmax
!        do i = 1,imax
!            write(33,*) buffer_grd(k,i,1:2),buffer_rd2d(k,i,1:4),buffer_rd1d(i,1)
!        enddo
!        enddo
!        close(33)
    end subroutine ReadAcoustic_Subset
! ***********************************************

    subroutine WriteStatistics()
        type(tp_rdwt_hdf5) :: fsol1d, fsol2d
        integer :: i,k,n
        real(8),dimension(:,:),allocatable :: mean1d
        real(8),dimension(:,:,:),allocatable :: mean2d
        character(800) :: varname1d,varname2d,varnametec

        ! initialize Acoustic 2D HDF5
        fsol2d%fname = trim(flowfile)
!        fsol2d%fname = trim(flowfile)//'_Stat2d.h5'
        fsol2d%gname = "/Stat2d"
        call DetectHDF5(fsol2d)
        varname2d=''
        do i=1,fsol2d%dnum
            varname2d=trim(adjustl(varname2d))//' '//trim(adjustl(fsol2d%dname(i)))
        enddo
        print*,'var2d: ',fsol2d%dnum,trim(varname2d)

        ! initialize Acoustic 1D HDF5
!        fsol1d%fname = trim(flowfile)//'_Int_BotWall.h5'
!        if(iwall.eq.2) fsol1d%fname = trim(flowfile)//'_Int_TopWall.h5'
!        fsol1d%gname = trim(fn_grp1d)
!        call DetectHDF5(fsol1d)
!        varname1d=''
!        do i=1,fsol1d%dnum
!            varname1d=trim(adjustl(varname1d))//' '//trim(adjustl(fsol1d%dname(i)))
!        enddo
!        print*,'var1d: ',fsol1d%dnum,trim(varname1d)

!        print*,'var2d_append: ',num_append2d,trim(varname_append2d)
        print*,'var1d_append: ',num_append1d,trim(varname_append1d)

        varnametec = 'x z '//trim(varname2d)//' '//' '//' '//trim(varname_append1d)
        call InitTec(1,imax,1,kmax,fsol2d%dnum+num_append1d+2,varnametec,0)
!        varnametec = 'x z '//trim(varname2d)//' '//trim(varname1d)//' '//trim(varname_append2d)//' '//trim(varname_append1d)
!        call InitTec(1,imax,1,kmax,fsol1d%dnum+fsol2d%dnum+2+num_append1d+num_append2d,varnametec,0)
!        varnametec = 'x z '//trim(varname2d)//' '//trim(varname1d)
!        call InitTec(1,imax,1,kmax,fsol1d%dnum+fsol2d%dnum+2,varnametec,0)
        print*,'Output variables: ',trim(varnametec)

        forall(k=1:kmax,i=1:imax)
            vartmp(i,1,k,1)=buffer_grd(k,i,1)
            vartmp(i,1,k,2)=buffer_grd(k,i,2)
        end forall
        deallocate(buffer_grd)

        print*,'Read file: ',trim(fsol2d%fname)
        allocate(mean2d(kmax,imax,fsol2d%dnum))
        call ReadHDF5_2D(fsol2d,mean2d)
        forall(k=1:kmax,i=1:imax)
            vartmp(i,1,k,3:fsol2d%dnum+2)=mean2d(k,i,1:fsol2d%dnum)
        end forall
        deallocate(mean2d)

!        print*,'Read file: ',trim(fsol1d%fname)
!        allocate(mean1d(imax,fsol1d%dnum))
!        call ReadHDF5_1D(fsol1d,mean1d)
!        forall(k=1:kmax,i=1:imax)
!            vartmp(i,1,k,fsol2d%dnum+3:fsol1d%dnum+fsol2d%dnum+2)=mean1d(i,1:fsol1d%dnum)
!        end forall
!        deallocate(mean1d)

!        forall(k=1:kmax,i=1:imax)
!            vartmp(i,1,k,fsol1d%dnum+fsol2d%dnum+3:fsol1d%dnum+fsol2d%dnum+2+num_append2d)=buffer_stat2d(k,i,1:num_append2d)
!        end forall
!        deallocate(buffer_stat2d)

        forall(k=1:kmax,i=1:imax)
            vartmp(i,1,k,fsol2d%dnum+3:fsol2d%dnum+2+num_append1d)=buffer_stat1d(i,1:num_append1d)
        end forall
        deallocate(buffer_stat1d)

        print*,'Writing file: ',trim(fileout)//'-modify.plt'
        call WriteTec(trim(fileout)//'-modify.plt')

    end subroutine WriteStatistics

    ! **** read in grid data 'x','z' ****
    subroutine InitGrid()
        type(tp_rdwt_hdf5) :: grd
        integer :: i

        grd%fname = trim(grdfile)
        grd%gname='/'
        call DetectHDF5(grd)
!     print *, 'grd%gname = ', trim(grd%gname)
!     print *, 'grd%rank  = ', grd%rank
!     print *, 'grd%dimsf = ', grd%dimsf
!     print *, 'grd%dnum  = ', grd%dnum
!     print *, 'grd%dname = ', (trim(grd%dname(i))//', ',i=1,grd%dnum)
!     print *, '*********************************************'
        kmax=grd%dimsf(1); imax=grd%dimsf(2)
        grd%rank=2
        grd%dnum=2
        if(allocated(grd%dname)) deallocate(grd%dname)
        allocate(grd%dname(grd%dnum))
        grd%dname(1)='x'
        grd%dname(2)='z'
        grd%IsHSInitialized = .true.
        allocate(buffer_grd(kmax,imax,grd%dnum))
        call ReadHDF5_2D(grd,buffer_grd)

    end subroutine InitGrid
    ! ***********************************************

    subroutine Input()
        implicit none

        read(*,*)
        read(*,*) icylin, iwall, iflow, ivis
        read(*,*)
        read(*,*) uinf, tinf, zfreestream
        read(*,*)
        read(*,'(a)') grdfile
        read(*,*)
        read(*,'(a)') flowfile
        read(*,*)
        read(*,'(a)') fileout

        pi = 4.d0*atan(1.d0)
        Mw = 28.97
        if(iflow .eq. 2) Mw=28.01
        rbar = R/Mw
        cp = gamma*rbar/(gamma-1.)
        Ma = uinf/(sqrt(gamma*rbar*tinf))

        fn_grp1d = '/Int_BotWall'
        if(iwall.eq.2) then   ! top wall
            fn_grp1d = '/Int_TopWall'
        endif
    end subroutine Input

    subroutine trapezoid_integration(ndim, func_int, y_int, nbeg, nend, res_int)
        ! integration of func_int over y_int between nbeg and nend using trapezoid rule
        ! result is stored in res_int
        implicit none
        integer,intent(in) :: ndim,nbeg,nend
        real(8),dimension(ndim),intent(in) :: y_int,func_int
        real(8),intent(out) :: res_int

        integer :: n,n1

        if(nbeg .lt. 1 .or. nend .gt. nend) then
           print*,'Integration Limit Exceed!!! '
           stop
        endif

        res_int = 0.0
        do n = nbeg,nend-1
          n1 = n+1
          res_int = res_int+(func_int(n1)+func_int(n))*(y_int(n1)-y_int(n))
        enddo
        res_int = 0.5*res_int
    end subroutine trapezoid_integration
end program CalStatistics
