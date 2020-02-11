program volume
  use MFileIO
  use modVolume
  use modReadHDF5
  use modWriteHDF5
!  use HDF5
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
  real, dimension(:,:,:,:), allocatable :: buffer_grid
!  real, dimension(:,:,:), allocatable :: u, v, w, p, t, rhon
  real, dimension(:,:,:,:), allocatable :: buffer_flow
  real, dimension(:,:,:,:,:), allocatable :: buffer_flow_tmp
!  real, dimension(:,:,:,:), allocatable :: vars
  real, dimension(:,:,:), allocatable :: grad_p, grad_rho
  real, dimension(:,:,:), allocatable :: omx, omy, omz
  real, dimension(:,:,:), allocatable :: swirl
  real, dimension(:,:,:), allocatable :: div
  real, dimension(:,:,:), allocatable :: Q
  real, dimension(:,:,:), allocatable :: dudn, dtdn, ds
  integer :: n,i,j,k,m, nn
  real(8) :: uinf, delta, utau, ztau
  integer :: iformat
!  real(8) :: pi, dudy_ana, dudy_num ! for test
  integer, dimension(:), allocatable :: varindex_vol_output
  character(10), dimension(:), allocatable :: dname_vol_output

  integer :: nvar_vol, nvar_vol_output, nvar_xyp, nvar_xyp_output
  parameter(nvar_vol=17)
  character(10) :: dname_vol(nvar_vol)
  parameter(dname_vol = (/'x','y','z','u','v','w','p','t','rho','grad_p','grad_rho','omx','omy','omz','div','swil','Q'/))
 ! integer(HID_T) :: file_id, group_id, dset_id, dspace_id
  integer(HSIZE_T), dimension(:), allocatable :: dimsf
  integer :: hdferr, rank
  integer :: ilen, jlen, klen

  integer :: ist_tsjplane, iend_tsjplane, isp_tsjplane, kst_tsjplane, kend_tsjplane, ksp_tsjplane, nt


  call Input()
  call InitHDF5()

  do n=iil,iih,stride
    write(unit=fnum,fmt='(I08.8)')n

    fname = trim(datapath)//'timeseries_'//fnum//'.h5'
    call readHDF5sol(fname,imax,jmax,kmax,5,buffer_flow_tmp)

!    call Readtsflow_jplane(fname,nt,iend_tsjplane-ist_tsjplane+1,kend_tsjplane-kst_tsjplane+1,10,buffer_flow_tmp)  !!!!!!!

 !   print *, 'min(T) = ', minval(buffer_flow_tmp(:,:,:,:,5))
!    print *, '', size(buffer_flow_tmp, dim=1)
!    print *, '', size(buffer_flow_tmp, dim=2)
!    print *, '', size(buffer_flow_tmp, dim=3)
!    print *, '', size(buffer_flow_tmp, dim=4)
!    print *, '', size(buffer_flow_tmp, dim=5)

!   print *, 'buffer_flow_tmp(1,10,2,1,5) = ', buffer_flow_tmp(1,10,2,1,5)
!   print *, 'buffer_flow_tmp(1,10,3,1,5) = ', buffer_flow_tmp(1,10,3,1,5)
!   print *, 'buffer_flow_tmp(1,10,2,1,1) = ', buffer_flow_tmp(1,10,2,1,1)
!   print *, 'buffer_flow_tmp(1,10,3,1,1) = ', buffer_flow_tmp(1,10,3,1,1)
!   print *, 'buffer_flow_tmp(1,10,2,1,2) = ', buffer_flow_tmp(1,10,2,1,2)
!   print *, 'buffer_flow_tmp(1,10,3,1,2) = ', buffer_flow_tmp(1,10,3,1,2)
!   print *, 'buffer_flow_tmp(1,10,2,1,3) = ', buffer_flow_tmp(1,10,2,1,3)
!   print *, 'buffer_flow_tmp(1,10,3,1,3) = ', buffer_flow_tmp(1,10,3,1,3)

!  print *, buffer_flow_tmp(1,1,3,1,1)
!  print *, buffer_flow_tmp(1,10,3,1,1)

! stop


    ! # of time steps
    do nn=1, nt
      do i=ist_tsjplane, iend_tsjplane, isp_tsjplane
        do j=1, 1
          do k=kst_tsjplane, kend_tsjplane, ksp_tsjplane
            buffer_flow(i,j,k,1) = buffer_flow_tmp(nn,i,k,1,1)
            buffer_flow(i,j,k,2) = buffer_flow_tmp(nn,i,k,1,2)
            buffer_flow(i,j,k,3) = buffer_flow_tmp(nn,i,k,1,3)
            buffer_flow(i,j,k,4) = buffer_flow_tmp(nn,i,k,1,4)
            buffer_flow(i,j,k,5) = buffer_flow_tmp(nn,i,k,1,5)
            buffer_flow(i,j,k,7) = buffer_flow_tmp(nn,i,k,1,6)  ! dudj
          enddo
        enddo
      enddo

 ! print *, '', buffer_flow(10,1,)



      buffer_flow(ist_tsjplane:iend_tsjplane,1,kst_tsjplane:kend_tsjplane,6) = buffer_flow(ist_tsjplane:iend_tsjplane,1,kst_tsjplane:kend_tsjplane,4)/(287.0*buffer_flow(ist_tsjplane:iend_tsjplane,1,kst_tsjplane:kend_tsjplane,5))

!      print *, 'min(rho) = ', minval(buffer_flow(:,:,:,6))

      if (iplane_xz.gt.0) then

        do m=1,iplane_xz
          write(unit=fnum1,fmt='(I04.4)')planes_xz(m)
          fname='plane_xz_'//fnum//'_j'//fnum1//'.h5'
          idx_xzp%jst=planes_xz(m); idx_xzp%jend=planes_xz(m); idx_xzp%jsp=1

          ilen = (idx_xzp%iend - idx_xzp%ist)/idx_xzp%isp + 1
          jlen = 1
          klen = (idx_xzp%kend - idx_xzp%kst)/idx_xzp%ksp + 1

          rank = 3
          allocate(dimsf(rank))
          dimsf(1) = ilen
          dimsf(2) = jlen
          dimsf(3) = klen

        ! write x
        call Write3DHDF5(fname, rank, dimsf,"x",real(x(idx_xzp%ist:idx_xzp%iend:idx_xzp%isp,idx_xzp%jst:idx_xzp%jend:idx_xzp%jsp, &
                       idx_xzp%kst:idx_xzp%kend:idx_xzp%ksp)),0)
        ! write z
        call Write3DHDF5(fname, rank, dimsf,"z",real(z(idx_xzp%ist:idx_xzp%iend:idx_xzp%isp,idx_xzp%jst:idx_xzp%jend:idx_xzp%jsp, &
                       idx_xzp%kst:idx_xzp%kend:idx_xzp%ksp)),1)
        ! write u
          call Write3DHDF5(fname, rank,dimsf,"u",buffer_flow(idx_xzp%ist:idx_xzp%iend:idx_xzp%isp,1:1:1, &
                       idx_xzp%kst:idx_xzp%kend:idx_xzp%ksp,1),1)
        ! write v
          call Write3DHDF5(fname, rank,dimsf,"v",buffer_flow(idx_xzp%ist:idx_xzp%iend:idx_xzp%isp,1:1:1, &
                       idx_xzp%kst:idx_xzp%kend:idx_xzp%ksp,2),1)
        ! write w
          call Write3DHDF5(fname, rank,dimsf,"w",buffer_flow(idx_xzp%ist:idx_xzp%iend:idx_xzp%isp,1:1:1, &
                       idx_xzp%kst:idx_xzp%kend:idx_xzp%ksp,3),1)
        ! write p
          call Write3DHDF5(fname, rank,dimsf,"p",buffer_flow(idx_xzp%ist:idx_xzp%iend:idx_xzp%isp,1:1:1, &
                       idx_xzp%kst:idx_xzp%kend:idx_xzp%ksp,4),1)
        ! write t
          call Write3DHDF5(fname, rank,dimsf,"T",buffer_flow(idx_xzp%ist:idx_xzp%iend:idx_xzp%isp,1:1:1, &
                       idx_xzp%kst:idx_xzp%kend:idx_xzp%ksp,5),1)
        ! write rho
          call Write3DHDF5(fname, rank,dimsf,"rho",buffer_flow(idx_xzp%ist:idx_xzp%iend:idx_xzp%isp,1:1:1, &
                       idx_xzp%kst:idx_xzp%kend:idx_xzp%ksp,6),1)

      print *, 'buffer_flow(10,1,1,6) = ', buffer_flow(10,1,1,6)

        call CalGradient_timeseries_jplane(idx_xzp,x,y,z,buffer_flow(ist_tsjplane:iend_tsjplane:1,1:1:1,kst_tsjplane:kend_tsjplane:1,6), &
                                           buffer_flow(ist_tsjplane:iend_tsjplane:1,1:1:1,kst_tsjplane:kend_tsjplane:1,7),grad_rho)  !!!!!!!! ddj
        call Write3DHDF5(fname,rank,dimsf,"grad_rho",grad_rho,1)


!        call CalGradient(idx_xzp,x,y,z,p,grad_p)
!        call CalVorticity(idx_xzp,x,y,z,u,v,w,omx,omy,omz)
!        call CalDivergence(idx_xzp,x,y,z,u,v,w,div)
!        call CalSwirl(idx_xzp,x,y,z,u,v,w,swirl)
!        call OutputPlaneXZ(planes_xz(m),fname,idx_xzp)
        end do
      end if

    enddo

  end do  ! end ii loop

  call FinalizeHDF5()

  contains
    subroutine Input()
      integer :: i
      read(*,*)
      read(*,'(A)')datapath  !path where data files are stored
      read(*,*)
      read(*,*)uinf,delta,utau,ztau
      read(*,*)
      read(*,*)iil,iih,stride, iformat
      read(*,*)
      read(*,*)ivolume, idx_vol%ist, idx_vol%iend, idx_vol%isp, idx_vol%jst, idx_vol%jend, idx_vol%jsp, idx_vol%kst, idx_vol%kend, idx_vol%ksp
      read(*,*)
      read(*,*) nvar_vol_output
      allocate(varindex_vol_output(nvar_vol_output), dname_vol_output(nvar_vol_output))
      read(*,*) (varindex_vol_output(i), i=1,nvar_vol_output)
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
      read(*,*)
      read(*,*) ist_tsjplane, iend_tsjplane, isp_tsjplane, kst_tsjplane, kend_tsjplane, ksp_tsjplane, nt

      call ReadPlot3DGridPlaneGen(trim(datapath)//'gridp3d.grd',imax,jmax,kmax,x,y,z,0)

      allocate(buffer_flow(ist_tsjplane:iend_tsjplane,1,kst_tsjplane:kend_tsjplane,10))  !!!!!!!!
      allocate(buffer_flow_tmp(nt,ist_tsjplane:iend_tsjplane,kst_tsjplane:kend_tsjplane,1,10))



!      allocate(buffer_flow(imax,jmax,kmax,6), buffer_flow_tmp(kmax,imax,jmax,5))

!      if(iformat.eq.0) then
!        allocate(vars(imax,jmax,kmax,ns+5))
!      else
!        allocate(vars(kmax,imax,jmax,5))
!      endif
!      allocate(u(imax,jmax,kmax), v(imax,jmax,kmax), w(imax,jmax,kmax))
!      allocate(p(imax,jmax,kmax), t(imax,jmax,kmax), rhon(imax,jmax,kmax))
      return
    end subroutine Input


end program volume    
