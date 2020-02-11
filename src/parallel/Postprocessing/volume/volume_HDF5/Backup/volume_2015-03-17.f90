program volume
  use MFileIO
  use modVolume
  implicit none
  character(200) :: fname, finfo
  character(200) :: datapath
  character(20) :: surffix
  character(8) :: fnum 
  character(4) :: fnum1
  integer :: imax, jmax, kmax
  integer :: iil, iih, stride
  integer :: ivolume, iplane_xz, iplane_yz, iplane_xy
  type(index_output) :: idx_vol, idx_xyp, idx_yzp, idx_xzp
  integer, dimension(:),  allocatable :: planes_xy, planes_yz, planes_xz
  integer, parameter :: ns=1
  real, dimension(:,:,:), allocatable :: x, y, z
  real, dimension(:,:,:), allocatable :: u, v, w, p, t
  real, dimension(:,:,:,:), allocatable :: vars,rhon
  real, dimension(:,:,:), allocatable :: grad_p, grad_rho
  real, dimension(:,:,:), allocatable :: omx, omy, omz
  real, dimension(:,:,:), allocatable :: swirl
  real, dimension(:,:,:), allocatable :: div
  real, dimension(:,:,:), allocatable :: dudn, dtdn, ds
  integer :: n,i,j,k,m
  real(8) :: uinf, delta, utau, ztau
!  real(8) :: pi, dudy_ana, dudy_num ! for test
  
  call Input()
  do n=iil,iih,stride
    write(unit=fnum,fmt='(I08.8)')n
    fname = trim(datapath)//'flowdata_'//fnum//'.sol'
    call  ReadPlot3DSolPlane(fname,imax,jmax,kmax,ns+5,vars)
    u = vars(:,:,:,1)
    v = vars(:,:,:,2)
    w = vars(:,:,:,3)
    p = vars(:,:,:,4)
    t = vars(:,:,:,5)
    rhon = vars(:,:,:,6:5+ns)

!*********************************************
!   ! For debug purposes
!        pi = 4*atan(1.)
!        forall(i=1:imax,j=1:jmax,k=1:kmax)
!           u(i,j,k) = cos(2.*pi*x(i,j,k))*sin(2.*pi*y(i,j,k))*cos(2.*pi*z(i,j,k))
!        end forall
!        i = 220; j=2; k=2
!        dudy_ana = 2.*pi*cos(2.*pi*x(i,j,k))*cos(2.*pi*y(i,j,k))*cos(2.*pi*z(i,j,k))
!        dudy_num = ddgrid(i,j,k,x,y,z,u,2)
!        print *, 'dudy_ana =', dudy_ana
!        print *, 'dudy_num =', dudy_num
   !*******************************

    if (ivolume.gt.0) then
      fname='volume_'//fnum//'.dat'
      finfo = ' '
      write(finfo,'("DATASETAUXDATA IndexRange = ""ist,iend,isp,jst,jend,jsp,kst,kend,ksp:",I4,",",I4,",",I4,",",I4,",",I4,",",I4,",",I4,",",I4,",",I4,"""")') idx_vol%ist, idx_vol%iend, idx_vol%isp, idx_vol%jst, idx_vol%jend, idx_vol%jsp, idx_vol%kst, idx_vol%kend, idx_vol%ksp
      call CalGradient(idx_vol,x,y,z,sum(rhon(:,:,:,:),dim=4),grad_rho)
      call CalGradient(idx_vol,x,y,z,p,grad_p)
      call CalVorticity(idx_vol,x,y,z,u,v,w,omx,omy,omz)
      call CalDivergence(idx_vol,x,y,z,u,v,w,div)
      call CalSwirl(idx_vol,x,y,z,u,v,w,swirl)
      call OutputVolume(fname,idx_vol)
    endif

    if (iplane_xy.gt.0) then
      write(finfo,'("DATASETAUXDATA IndexRange = ""ist,iend,isp,jst,jend,jsp:",I4,",",I4,",",I4,",",I4,",",I4,",",I4,"""")') idx_xyp%ist, idx_xyp%iend, idx_xyp%isp, idx_xyp%jst, idx_xyp%jend, idx_xyp%jsp
      do m=1,iplane_xy
        write(unit=fnum1,fmt='(I04.4)')planes_xy(m)
        fname='plane_xy_'//fnum//'_k'//fnum1//'.dat'
        idx_xyp%kst=planes_xy(m); idx_xyp%kend=planes_xy(m); idx_xyp%ksp=1
        call CalGradient(idx_xyp,x,y,z,sum(rhon(:,:,:,:),dim=4),grad_rho)
        call CalGradient(idx_xyp,x,y,z,p,grad_p)
        call CalVorticity(idx_xyp,x,y,z,u,v,w,omx,omy,omz)
        call CalDivergence(idx_xyp,x,y,z,u,v,w,div)
        call CalSwirl(idx_xyp,x,y,z,u,v,w,swirl)
        if(planes_xy(m).eq.1) call CalWallStat(idx_xyp,x,y,z,u,t,dudn,dtdn,ds)
        call OutputPlaneXY(planes_xy(m),fname,idx_xyp)
      end do
    end if

    if (iplane_yz.gt.0) then
      write(finfo,'("DATASETAUXDATA IndexRange = ""jst,jend,jsp,kst,kend,ksp:",I4,",",I4,",",I4,",",I4,",",I4,",",I4,"""")') idx_yzp%jst, idx_yzp%jend, idx_yzp%jsp, idx_yzp%kst, idx_yzp%kend, idx_yzp%ksp
      do m=1,iplane_yz
        write(unit=fnum1,fmt='(I04.4)')planes_yz(m)
        fname='plane_yz_'//fnum//'_i'//fnum1//'.dat'
        idx_yzp%ist=planes_yz(m); idx_yzp%iend=planes_yz(m); idx_yzp%isp=1
        call CalGradient(idx_yzp,x,y,z,sum(rhon(:,:,:,:),dim=4),grad_rho)
        call CalGradient(idx_yzp,x,y,z,p,grad_p)
        call CalVorticity(idx_yzp,x,y,z,u,v,w,omx,omy,omz)
        call CalDivergence(idx_yzp,x,y,z,u,v,w,div)
        call CalSwirl(idx_yzp,x,y,z,u,v,w,swirl)
        call OutputPlaneYZ(planes_yz(m),fname,idx_yzp)
      end do
    end if

    if (iplane_xz.gt.0) then
      write(finfo,'("DATASETAUXDATA IndexRange = ""ist,iend,isp,kst,kend,ksp:",I4,",",I4,",",I4,",",I4,",",I4,",",I4,"""")') idx_xzp%ist, idx_xzp%iend, idx_xzp%isp, idx_xzp%kst, idx_xzp%kend, idx_xzp%ksp
      do m=1,iplane_xz
        write(unit=fnum1,fmt='(I04.4)')planes_xz(m)
        fname='plane_xz_'//fnum//'_j'//fnum1//'.dat'
        idx_xzp%jst=planes_xz(m); idx_xzp%jend=planes_xz(m); idx_xzp%jsp=1
        call CalGradient(idx_xzp,x,y,z,sum(rhon(:,:,:,:),dim=4),grad_rho)
        call CalGradient(idx_xzp,x,y,z,p,grad_p)
        call CalVorticity(idx_xzp,x,y,z,u,v,w,omx,omy,omz)
        call CalDivergence(idx_xzp,x,y,z,u,v,w,div)
        call CalSwirl(idx_xzp,x,y,z,u,v,w,swirl)
        call OutputPlaneXZ(planes_xz(m),fname,idx_xzp)
      end do
    end if

    if(allocated(grad_rho)) deallocate(grad_rho)
    if(allocated(grad_p)) deallocate(grad_p)
    if(allocated(omx)) deallocate(omx)
    if(allocated(omy)) deallocate(omy)
    if(allocated(omz)) deallocate(omz)
    if(allocated(div)) deallocate(div)
    if(allocated(swirl)) deallocate(swirl)
    if(allocated(dudn)) deallocate(dudn)
    if(allocated(dtdn)) deallocate(dtdn)
    if(allocated(ds)) deallocate(ds)
  end do  ! end ii loop

  contains
    subroutine Input()
      integer :: i
      read(*,*)
      read(*,'(A)')datapath  !path where data files are stored
      read(*,*)
      read(*,*)uinf,delta,utau,ztau
      read(*,*)
      read(*,*)iil,iih,stride
      read(*,*)
      read(*,*)ivolume, idx_vol%ist, idx_vol%iend, idx_vol%isp, idx_vol%jst, idx_vol%jend, idx_vol%jsp, idx_vol%kst, idx_vol%kend, idx_vol%ksp
!     xy planes
      read(*,*)
      read(*,*)iplane_xy,idx_xyp%ist, idx_xyp%iend, idx_xyp%isp, idx_xyp%jst, idx_xyp%jend, idx_xyp%jsp
      if (iplane_xy.gt.0) then
        allocate(planes_xy(1:iplane_xy))
        read(*,*)(planes_xy(i),i=1,iplane_xy)
      else
        read(*,*)
      end if
!     yz planes
      read(*,*)
      read(*,*)iplane_yz, idx_yzp%jst, idx_yzp%jend, idx_yzp%jsp, idx_yzp%kst, idx_yzp%kend, idx_yzp%ksp
      if (iplane_yz.gt.0) then
        allocate(planes_yz(1:iplane_yz))
        read(*,*)(planes_yz(i),i=1,iplane_yz)
      else
        read(*,*)
      end if
!     xz planes
      read(*,*)
      read(*,*)iplane_xz, idx_xzp%ist, idx_xzp%iend, idx_xzp%isp, idx_xzp%kst, idx_xzp%kend, idx_xzp%ksp
      if (iplane_xz.gt.0) then
        allocate(planes_xz(1:iplane_xz))
        read(*,*)(planes_xz(i),i=1,iplane_xz)
      else
        read(*,*)
      end if
      call ReadPlot3DGridPlaneGen(trim(datapath)//'gridp3d.grd',imax,jmax,kmax,x,y,z,0)
      allocate(vars(imax,jmax,kmax,ns+5))
      allocate(u(imax,jmax,kmax), v(imax,jmax,kmax), w(imax,jmax,kmax))
      allocate(p(imax,jmax,kmax), t(imax,jmax,kmax), rhon(imax,jmax,kmax,ns))
      return
    end subroutine Input

    subroutine OutputVolume(fn,vol)
      character(*), intent(in) :: fn
      type(index_output), intent(in) :: vol
      integer :: i,j,k
      integer :: ilen, jlen, klen
      
      ilen = (vol%iend-vol%ist)/vol%isp+1
      jlen = (vol%jend-vol%jst)/vol%jsp+1
      klen = (vol%kend-vol%kst)/vol%ksp+1
      open(27,file=fn,status='unknown')
      rewind 27
      write(27,'(A)')'variables = "x","y","z","u","v","w","p","t","rho","grad_p","grad_rho","omx", "omy", "omz", "div","swirl"'
      write(27,'(A)') trim(finfo)
!      write(27,*)'zone i=',ilen,' j=',jlen,' k=',klen,' f=block'
      write(27,'("Zone T=Volume, I=", I4, ", J=",I4, ", K=",I4, ", F=BLOCK, AUXDATA uinf=""", E15.9, """, AUXDATA delta0=""", E15.9, """, AUXDATA utau=""", E15.9, """, AUXDATA ztau=""", E15.9,"""")') ilen, jlen, klen, uinf, delta, utau, ztau
      write(27,*)(((x(i,j,k)           ,i=vol%ist,vol%iend,vol%isp),j=vol%jst,vol%jend,vol%jsp),k=vol%kst,vol%kend,vol%ksp)
      write(27,*)(((y(i,j,k)           ,i=vol%ist,vol%iend,vol%isp),j=vol%jst,vol%jend,vol%jsp),k=vol%kst,vol%kend,vol%ksp)
      write(27,*)(((z(i,j,k)           ,i=vol%ist,vol%iend,vol%isp),j=vol%jst,vol%jend,vol%jsp),k=vol%kst,vol%kend,vol%ksp)
      write(27,*)(((u(i,j,k) 	       ,i=vol%ist,vol%iend,vol%isp),j=vol%jst,vol%jend,vol%jsp),k=vol%kst,vol%kend,vol%ksp)
      write(27,*)(((v(i,j,k)           ,i=vol%ist,vol%iend,vol%isp),j=vol%jst,vol%jend,vol%jsp),k=vol%kst,vol%kend,vol%ksp)
      write(27,*)(((w(i,j,k)           ,i=vol%ist,vol%iend,vol%isp),j=vol%jst,vol%jend,vol%jsp),k=vol%kst,vol%kend,vol%ksp)
      write(27,*)(((p(i,j,k)           ,i=vol%ist,vol%iend,vol%isp),j=vol%jst,vol%jend,vol%jsp),k=vol%kst,vol%kend,vol%ksp)
      write(27,*)(((t(i,j,k)           ,i=vol%ist,vol%iend,vol%isp),j=vol%jst,vol%jend,vol%jsp),k=vol%kst,vol%kend,vol%ksp)
      write(27,*)(((sum(rhon(i,j,k,:)) ,i=vol%ist,vol%iend,vol%isp),j=vol%jst,vol%jend,vol%jsp),k=vol%kst,vol%kend,vol%ksp)
      write(27,*)(((grad_p(i,j,k)      ,i=1,ilen),j=1,jlen),k=1,klen)
      write(27,*)(((grad_rho(i,j,k)    ,i=1,ilen),j=1,jlen),k=1,klen)
      write(27,*)(((omx(i,j,k)    ,i=1,ilen),j=1,jlen),k=1,klen)
      write(27,*)(((omy(i,j,k)    ,i=1,ilen),j=1,jlen),k=1,klen)
      write(27,*)(((omz(i,j,k)    ,i=1,ilen),j=1,jlen),k=1,klen)
      write(27,*)(((div(i,j,k)         ,i=1,ilen),j=1,jlen),k=1,klen)
      write(27,*)(((swirl(i,j,k)       ,i=1,ilen),j=1,jlen),k=1,klen)
      close(27)
    end subroutine OutputVolume

    subroutine OutputPlaneXY(kplane,fn,xyp)
      integer, intent(in) :: kplane
      character(*), intent(in) :: fn
      type(index_output), intent(in) :: xyp
      integer :: i,j
      integer :: ilen, jlen
      
      ilen = (xyp%iend-xyp%ist)/xyp%isp+1
      jlen = (xyp%jend-xyp%jst)/xyp%jsp+1
      open(30, file=fn, status='unknown')
      if(kplane.eq.1) then
          write(30,'(A)')'variables = "x","y","z","dudn","dtdn","ds","p","t","rho","grad_p","grad_rho","omz","swirl","div"'
          write(30,'(A)') trim(finfo)
!          write(30,*)'zone i=',ilen,' j=',jlen,' f=block'
          write(30,'("Zone T=XYPLANE, I=", I4, ", J=",I4, ", F=BLOCK, AUXDATA uinf=""", E15.9, """, AUXDATA delta0=""", E15.9, """, AUXDATA utau=""", E15.9, """, AUXDATA ztau=""", E15.9,"""")') ilen, jlen, uinf, delta, utau, ztau
          write(30,*)((x(i,j,kplane)   ,i=xyp%ist,xyp%iend,xyp%isp),j=xyp%jst,xyp%jend,xyp%jsp)
          write(30,*)((y(i,j,kplane)   ,i=xyp%ist,xyp%iend,xyp%isp),j=xyp%jst,xyp%jend,xyp%jsp)
          write(30,*)((z(i,j,kplane)   ,i=xyp%ist,xyp%iend,xyp%isp),j=xyp%jst,xyp%jend,xyp%jsp)
          write(30,*)((dudn(i,j,kplane)   ,i=1,ilen),j=1,jlen)
          write(30,*)((dtdn(i,j,kplane)   ,i=1,ilen),j=1,jlen)
          write(30,*)((ds(i,j,kplane)     ,i=1,ilen),j=1,jlen)
          write(30,*)((p(i,j,kplane)   ,i=xyp%ist,xyp%iend,xyp%isp),j=xyp%jst,xyp%jend,xyp%jsp)
          write(30,*)((t(i,j,kplane)   ,i=xyp%ist,xyp%iend,xyp%isp),j=xyp%jst,xyp%jend,xyp%jsp)
          write(30,*)((sum(rhon(i,j,kplane,:)) ,i=xyp%ist,xyp%iend,xyp%isp),j=xyp%jst,xyp%jend,xyp%jsp)
          write(30,*)((grad_p(i,j,1) ,i=1,ilen),j=1,jlen)
          write(30,*)((grad_rho(i,j,1) ,i=1,ilen),j=1,jlen)
          write(30,*)((omz(i,j,1) ,i=1,ilen),j=1,jlen)
          write(30,*)((swirl(i,j,1) ,i=1,ilen),j=1,jlen)
          write(30,*)((div(i,j,1) ,i=1,ilen),j=1,jlen)
      else
          write(30,'(A)')'variables = "x","y","z","u","v","w","p","t","rho","grad_p","grad_rho","omz","swirl","div"'
          write(30,'(A)') trim(finfo)
!          write(30,*)'zone i=',ilen,' j=',jlen,' f=block'
          write(30,'("Zone T=XYPLANE, I=", I4, ", J=",I4, ", F=BLOCK, AUXDATA uinf=""", E15.9, """, AUXDATA delta0=""", E15.9, """, AUXDATA utau=""", E15.9, """, AUXDATA ztau=""", E15.9,"""")') ilen, jlen, uinf, delta, utau, ztau
          write(30,*)((x(i,j,kplane)   ,i=xyp%ist,xyp%iend,xyp%isp),j=xyp%jst,xyp%jend,xyp%jsp)
          write(30,*)((y(i,j,kplane)   ,i=xyp%ist,xyp%iend,xyp%isp),j=xyp%jst,xyp%jend,xyp%jsp)
          write(30,*)((z(i,j,kplane)   ,i=xyp%ist,xyp%iend,xyp%isp),j=xyp%jst,xyp%jend,xyp%jsp)
          write(30,*)((u(i,j,kplane)   ,i=xyp%ist,xyp%iend,xyp%isp),j=xyp%jst,xyp%jend,xyp%jsp)
          write(30,*)((v(i,j,kplane)   ,i=xyp%ist,xyp%iend,xyp%isp),j=xyp%jst,xyp%jend,xyp%jsp)
          write(30,*)((w(i,j,kplane)   ,i=xyp%ist,xyp%iend,xyp%isp),j=xyp%jst,xyp%jend,xyp%jsp)
          write(30,*)((p(i,j,kplane)   ,i=xyp%ist,xyp%iend,xyp%isp),j=xyp%jst,xyp%jend,xyp%jsp)
          write(30,*)((t(i,j,kplane)   ,i=xyp%ist,xyp%iend,xyp%isp),j=xyp%jst,xyp%jend,xyp%jsp)
          write(30,*)((sum(rhon(i,j,kplane,:)) ,i=xyp%ist,xyp%iend,xyp%isp),j=xyp%jst,xyp%jend,xyp%jsp)
          write(30,*)((grad_p(i,j,1) ,i=1,ilen),j=1,jlen)
          write(30,*)((grad_rho(i,j,1) ,i=1,ilen),j=1,jlen)
          write(30,*)((omz(i,j,1) ,i=1,ilen),j=1,jlen)
          write(30,*)((swirl(i,j,1) ,i=1,ilen),j=1,jlen)
          write(30,*)((div(i,j,1) ,i=1,ilen),j=1,jlen)
      endif
      close(30)

      return
    end subroutine OutputPlaneXY

    subroutine OutputPlaneYZ(iplane,fn,yzp)
      integer, intent(in) :: iplane
      character(*), intent(in) :: fn
      type(index_output), intent(in) :: yzp
      integer :: j,k
      integer :: jlen, klen
      
      jlen = (yzp%jend-yzp%jst)/yzp%jsp+1
      klen = (yzp%kend-yzp%kst)/yzp%ksp+1
      open(30, file=fn, status='unknown')
      write(30,'(A)')'variables = "y","z","u","v","w","p","t","rho","grad_p","grad_rho","omx","swirl","div"'
      write(30,'(A)') trim(finfo)
!      write(30,*)'zone i=',jlen,' j=',klen,' f=block'
      write(30,'("Zone T=YZPLANE, I=", I4, ", J=",I4, ", F=BLOCK, AUXDATA uinf=""", E15.9, """, AUXDATA delta0=""", E15.9, """, AUXDATA utau=""", E15.9, """, AUXDATA ztau=""", E15.9,"""")') jlen, klen, uinf, delta, utau, ztau
      write(30,*)((y(iplane,j,k)   ,j=yzp%jst,yzp%jend,yzp%jsp),k=yzp%kst,yzp%kend,yzp%ksp)
      write(30,*)((z(iplane,j,k)   ,j=yzp%jst,yzp%jend,yzp%jsp),k=yzp%kst,yzp%kend,yzp%ksp)
      write(30,*)((u(iplane,j,k)   ,j=yzp%jst,yzp%jend,yzp%jsp),k=yzp%kst,yzp%kend,yzp%ksp)
      write(30,*)((v(iplane,j,k)   ,j=yzp%jst,yzp%jend,yzp%jsp),k=yzp%kst,yzp%kend,yzp%ksp)
      write(30,*)((w(iplane,j,k)   ,j=yzp%jst,yzp%jend,yzp%jsp),k=yzp%kst,yzp%kend,yzp%ksp)
      write(30,*)((p(iplane,j,k)   ,j=yzp%jst,yzp%jend,yzp%jsp),k=yzp%kst,yzp%kend,yzp%ksp)
      write(30,*)((t(iplane,j,k)   ,j=yzp%jst,yzp%jend,yzp%jsp),k=yzp%kst,yzp%kend,yzp%ksp)
      write(30,*)((sum(rhon(iplane,j,k,:)) ,j=yzp%jst,yzp%jend,yzp%jsp),k=yzp%kst,yzp%kend,yzp%ksp)
      write(30,*)((grad_p(1,j,k) ,j=1,jlen), k=1,klen)
      write(30,*)((grad_rho(1,j,k) ,j=1,jlen), k=1,klen)
      write(30,*)((omx(1,j,k) ,j=1,jlen), k=1,klen)
      write(30,*)((swirl(1,j,k) ,j=1,jlen), k=1,klen)
      write(30,*)((div(1,j,k) ,j=1,jlen), k=1,klen)
      close(30)
      return
    end subroutine OutputPlaneYZ
      
    subroutine OutputPlaneXZ(jplane,fn,xzp)
      integer, intent(in) :: jplane
      character(*) :: fn
      type(index_output), intent(in) :: xzp
      integer :: i,k
      integer :: ilen, klen
      
      ilen = (xzp%iend-xzp%ist)/xzp%isp+1
      klen = (xzp%kend-xzp%kst)/xzp%ksp+1
      open(30, file=fn, status='unknown')
      rewind(30)
      write(30,'(A)')'variables = "x","z","u","v","w","p","t","rho","grad_p","grad_rho","omy","swirl","div"'
      write(30,'(A)') trim(finfo)
!      write(30,*)'zone i=',ilen,' j=',klen,' f=block '
      write(30,'("Zone T=XZPLANE, I=", I4, ", J=",I4, ", F=BLOCK, AUXDATA uinf=""", E15.9, """, AUXDATA delta0=""", E15.9, """, AUXDATA utau=""", E15.9, """, AUXDATA ztau=""", E15.9,"""")') ilen, klen, uinf, delta, utau, ztau
      write(30,*)((x(i,jplane,k)   ,i=xzp%ist,xzp%iend,xzp%isp),k=xzp%kst,xzp%kend,xzp%ksp)
      write(30,*)((z(i,jplane,k)   ,i=xzp%ist,xzp%iend,xzp%isp),k=xzp%kst,xzp%kend,xzp%ksp)
      write(30,*)((u(i,jplane,k)   ,i=xzp%ist,xzp%iend,xzp%isp),k=xzp%kst,xzp%kend,xzp%ksp)
      write(30,*)((v(i,jplane,k)   ,i=xzp%ist,xzp%iend,xzp%isp),k=xzp%kst,xzp%kend,xzp%ksp)
      write(30,*)((w(i,jplane,k)   ,i=xzp%ist,xzp%iend,xzp%isp),k=xzp%kst,xzp%kend,xzp%ksp)
      write(30,*)((p(i,jplane,k)   ,i=xzp%ist,xzp%iend,xzp%isp),k=xzp%kst,xzp%kend,xzp%ksp)
      write(30,*)((t(i,jplane,k)   ,i=xzp%ist,xzp%iend,xzp%isp),k=xzp%kst,xzp%kend,xzp%ksp)
      write(30,*)((sum(rhon(i,jplane,k,:)) ,i=xzp%ist,xzp%iend,xzp%isp),k=xzp%kst,xzp%kend,xzp%ksp)
      write(30,*)((grad_p(i,1,k) ,i=1,ilen), k=1,klen)
      write(30,*)((grad_rho(i,1,k) ,i=1,ilen), k=1,klen)
      write(30,*)((omy(i,1,k) ,i=1,ilen), k=1,klen)
      write(30,*)((swirl(i,1,k) ,i=1,ilen), k=1,klen)
      write(30,*)((div(i,1,k) ,i=1,ilen), k=1,klen)
      close(30)
      return
    end subroutine OutputPlaneXZ

end program volume    
