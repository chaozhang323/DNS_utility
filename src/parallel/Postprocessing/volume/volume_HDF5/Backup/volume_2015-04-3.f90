program volume
  use MFileIO
  use modVolume
  use modRWHDF5
  use modTecbin
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
  real, dimension(:,:,:,:), allocatable :: buffer_flow, buffer_flow_tmp
!  real, dimension(:,:,:,:), allocatable :: vars
  real, dimension(:,:,:), allocatable :: grad_p, grad_rho
  real, dimension(:,:,:), allocatable :: omx, omy, omz
  real, dimension(:,:,:), allocatable :: swirl
  real, dimension(:,:,:), allocatable :: div
  real, dimension(:,:,:), allocatable :: Q
  real, dimension(:,:,:), allocatable :: dudn, dtdn, ds
  integer :: n, i, j, k, m, nn
  real(8) :: uinf, delta, utau, ztau
  integer :: iformat
  integer, dimension(:), allocatable :: varindex_vol_output
  character(10), dimension(:), allocatable :: dname_vol_output

  integer :: nvar_vol, nvar_vol_output, nvar_xyp, nvar_xyp_output
  parameter(nvar_vol=17)
  character(10) :: dname_vol(nvar_vol)
  parameter(dname_vol = (/'x','y','z','u','v','w','p','t','rho','grad_p','grad_rho','omx','omy','omz','div','swil','Q'/))

  character(200) :: varname_output

  integer, dimension(:), allocatable :: dims
  integer :: rank
  integer :: ilen, jlen, klen

  call Input()
  call InitHDF5()

  do n=iil,iih,stride

!    if(allocated(buffer_flow_tmp)) deallocate(buffer_flow_tmp)
    allocate(buffer_flow_tmp(kmax,imax,jmax,5))

    write(unit=fnum,fmt='(I08.8)') n
    fname = trim(datapath)//'flowdata_'//fnum//'.h5'
    print *, 'Reading file: ', trim(fname)
    call ReadHDF5sol(fname,imax,jmax,kmax,buffer_flow_tmp)

    do i=1, imax
      do j=1, jmax
        do k=1, kmax
          buffer_flow(i,j,k,1) = buffer_flow_tmp(k,i,j,1)
          buffer_flow(i,j,k,2) = buffer_flow_tmp(k,i,j,2)
          buffer_flow(i,j,k,3) = buffer_flow_tmp(k,i,j,3)
          buffer_flow(i,j,k,4) = buffer_flow_tmp(k,i,j,4)
          buffer_flow(i,j,k,5) = buffer_flow_tmp(k,i,j,5)
        enddo
      enddo
    enddo
    deallocate(buffer_flow_tmp)

    buffer_flow(:,:,:,6) = buffer_flow(:,:,:,4)/(287.0*buffer_flow(:,:,:,5))

    if (ivolume.gt.0) then

      ilen = (idx_vol%iend - idx_vol%ist)/idx_vol%isp + 1
      jlen = (idx_vol%jend - idx_vol%jst)/idx_vol%jsp + 1
      klen = (idx_vol%kend - idx_vol%kst)/idx_vol%ksp + 1

      rank = 3
      allocate(dims(rank))
      dims(1) = ilen
      dims(2) = jlen
      dims(3) = klen

      varname_output = ''
      do nn=1, nvar_vol_output
        dname_vol_output(nn) = dname_vol(varindex_vol_output(nn))
        varname_output = (trim(varname_output)//' ', nn=1, nvar_vol_output) !!!!!!!!!!!!!!!!!!!!!!!!!
      enddo
      print *, 'dname_vol_output = ' , (trim(dname_vol_output(nn))//',', nn=1, nvar_vol_output)

      fname='volume_'//fnum//'.h5'

!      call InitTec(1,ilen,jlen,klen,12,"x y z u v w p t rho grad_p grad_rho div",0)


      ! write x
!      call Write3DHDF5(fname, rank, dims,"x",real(x(idx_vol%ist:idx_vol%iend:idx_vol%isp,idx_vol%jst:idx_vol%jend:idx_vol%jsp, &
!                       idx_vol%kst:idx_vol%kend:idx_vol%ksp)),0)
      ! write y
!      call Write3DHDF5(fname, rank, dims,"y",y(idx_vol%ist:idx_vol%iend:idx_vol%isp,idx_vol%jst:idx_vol%jend:idx_vol%jsp, &
!                       idx_vol%kst:idx_vol%kend:idx_vol%ksp),1)
      ! write z
!      call Write3DHDF5(fname, rank, dims,"z",z(idx_vol%ist:idx_vol%iend:idx_vol%isp,idx_vol%jst:idx_vol%jend:idx_vol%jsp, &
!                       idx_vol%kst:idx_vol%kend:idx_vol%ksp),1)

      do nn=1, nvar_vol_output
!        if(varindex_vol_output(nn).eq.4) then
          ! write u
!          call Write3DHDF5(fname, rank,dims,"u",buffer_flow(idx_vol%ist:idx_vol%iend:idx_vol%isp,idx_vol%jst:idx_vol%jend:idx_vol%jsp, &
!                       idx_vol%kst:idx_vol%kend:idx_vol%ksp,1),1)
!        endif
!        if(varindex_vol_output(nn).eq.5) then
          ! write v
!          call Write3DHDF5(fname, rank,dims,"v",buffer_flow(idx_vol%ist:idx_vol%iend:idx_vol%isp,idx_vol%jst:idx_vol%jend:idx_vol%jsp, &
!                       idx_vol%kst:idx_vol%kend:idx_vol%ksp,2),1)
!        endif
!        if(varindex_vol_output(nn).eq.6) then
          ! write w
!          call Write3DHDF5(fname, rank,dims,"w",buffer_flow(idx_vol%ist:idx_vol%iend:idx_vol%isp,idx_vol%jst:idx_vol%jend:idx_vol%jsp, &
!                       idx_vol%kst:idx_vol%kend:idx_vol%ksp,3),1)
!        endif
!        if(varindex_vol_output(nn).eq.7) then
          ! write p
!          call Write3DHDF5(fname, rank,dims,"p",buffer_flow(idx_vol%ist:idx_vol%iend:idx_vol%isp,idx_vol%jst:idx_vol%jend:idx_vol%jsp, &
!                       idx_vol%kst:idx_vol%kend:idx_vol%ksp,4),1)
!        endif
!        if(varindex_vol_output(nn).eq.8) then
          ! write t
!          call Write3DHDF5(fname, rank,dims,"t",buffer_flow(idx_vol%ist:idx_vol%iend:idx_vol%isp,idx_vol%jst:idx_vol%jend:idx_vol%jsp, &
!                       idx_vol%kst:idx_vol%kend:idx_vol%ksp,5),1)
!        endif
!        if(varindex_vol_output(nn).eq.9) then
          ! write rho
!          call Write3DHDF5(fname, rank,dims,"rho",buffer_flow(idx_vol%ist:idx_vol%iend:idx_vol%isp,idx_vol%jst:idx_vol%jend:idx_vol%jsp, &
!                       idx_vol%kst:idx_vol%kend:idx_vol%ksp,6),1)
!        endif
        if(varindex_vol_output(nn).eq.10) then
          call CalGradient(idx_vol,x,y,z,buffer_flow(:,:,:,4),grad_p)
          ! write grad_p
!          call Write3DHDF5(fname, rank,dims,"grad_p",grad_p,1)
        endif

        if(varindex_vol_output(nn).eq.11) then
          call CalGradient(idx_vol,x,y,z,buffer_flow(:,:,:,6),grad_rho)
          ! write grad_rho
!          call Write3DHDF5(fname, rank,dims,"grad_rho",grad_rho,1)
        endif
        if(varindex_vol_output(nn).eq.12) then
          call CalVorticity(idx_vol,x,y,z,buffer_flow(:,:,:,1),buffer_flow(:,:,:,2),buffer_flow(:,:,:,3),omx,omy,omz)
          ! write omx
!          call Write3DHDF5(fname, rank,dims,"omx",omx,1)
          ! write omy
!          call Write3DHDF5(fname, rank,dims,"omy",omy,1)
          ! write omz
!          call Write3DHDF5(fname, rank,dims,"omz",omz,1)
        endif
        if(varindex_vol_output(nn).eq.15) then
          call CalDivergence(idx_vol,x,y,z,buffer_flow(:,:,:,1),buffer_flow(:,:,:,2),buffer_flow(:,:,:,3),div)
          ! write div
!          call Write3DHDF5(fname, rank,dims,"div",div,1)
        endif
        if(varindex_vol_output(nn).eq.16) then
          call CalSwirl(idx_vol,x,y,z,buffer_flow(:,:,:,1),buffer_flow(:,:,:,2),buffer_flow(:,:,:,3),swirl)
          ! write swirl
 !         call Write3DHDF5(fname, rank,dims,"swirl",swirl,1)
        endif
        if(varindex_vol_output(nn).eq.17) then
          call Callambda2(idx_vol,x,y,z,buffer_flow(:,:,:,1),buffer_flow(:,:,:,2),buffer_flow(:,:,:,3),Q)
          ! write Q
 !         call Write3DHDF5(fname, rank,dims,"Q",Q,1)
        endif

      enddo


 !     call WriteTec12("volume.plt",x(idx_vol%ist:idx_vol%iend:idx_vol%isp,idx_vol%jst:idx_vol%jend:idx_vol%jsp, &
 !                      idx_vol%kst:idx_vol%kend:idx_vol%ksp), &
 !                      y(idx_vol%ist:idx_vol%iend:idx_vol%isp,idx_vol%jst:idx_vol%jend:idx_vol%jsp, &
 !                      idx_vol%kst:idx_vol%kend:idx_vol%ksp), &
 !                      z(idx_vol%ist:idx_vol%iend:idx_vol%isp,idx_vol%jst:idx_vol%jend:idx_vol%jsp, &
 !                      idx_vol%kst:idx_vol%kend:idx_vol%ksp), &
 !                      buffer_flow(idx_vol%ist:idx_vol%iend:idx_vol%isp,idx_vol%jst:idx_vol%jend:idx_vol%jsp, &
 !                      idx_vol%kst:idx_vol%kend:idx_vol%ksp,1), &
 !                      buffer_flow(idx_vol%ist:idx_vol%iend:idx_vol%isp,idx_vol%jst:idx_vol%jend:idx_vol%jsp, &
 !                      idx_vol%kst:idx_vol%kend:idx_vol%ksp,2), &
 !                      buffer_flow(idx_vol%ist:idx_vol%iend:idx_vol%isp,idx_vol%jst:idx_vol%jend:idx_vol%jsp, &
 !                      idx_vol%kst:idx_vol%kend:idx_vol%ksp,3), &
 !                      buffer_flow(idx_vol%ist:idx_vol%iend:idx_vol%isp,idx_vol%jst:idx_vol%jend:idx_vol%jsp, &
 !                      idx_vol%kst:idx_vol%kend:idx_vol%ksp,4), &
 !                      buffer_flow(idx_vol%ist:idx_vol%iend:idx_vol%isp,idx_vol%jst:idx_vol%jend:idx_vol%jsp, &
 !                      idx_vol%kst:idx_vol%kend:idx_vol%ksp,5), &
 !                      buffer_flow(idx_vol%ist:idx_vol%iend:idx_vol%isp,idx_vol%jst:idx_vol%jend:idx_vol%jsp, &
 !                      idx_vol%kst:idx_vol%kend:idx_vol%ksp,6), &
 !                      grad_p, grad_rho, div )



    endif




  end do  ! end ii loop

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
      call ReadPlot3DGridPlaneGen(trim(datapath)//'gridp3d.grd',imax,jmax,kmax,x,y,z,0)

      allocate(buffer_flow(imax,jmax,kmax,6))


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
