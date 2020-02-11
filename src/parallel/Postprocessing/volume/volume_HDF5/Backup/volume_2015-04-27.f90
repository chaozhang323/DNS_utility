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
  real, dimension(:,:,:,:), allocatable :: buffer_flow, buffer_flow_tmp
  real, dimension(:,:,:), allocatable :: grad_p, grad_rho
  real, dimension(:,:,:), allocatable :: omx, omy, omz
  real, dimension(:,:,:), allocatable :: swirl
  real, dimension(:,:,:), allocatable :: div
  real, dimension(:,:,:), allocatable :: Q
  real, dimension(:,:,:), allocatable :: dudn, dtdn, ds
  integer :: n, i, j, k, m, nn
  real(8) :: uinf, delta, utau, ztau
  integer, dimension(:), allocatable :: varindex_vol_output, varindex_xyp_output, varindex_yzp_output, varindex_xzp_output
  character(10), dimension(:), allocatable :: dname_vol_output, dname_xyp_output, dname_yzp_output, dname_xzp_output

  integer :: nvar, nvar_vol_output, nvar_xyp_output, nvar_yzp_output, nvar_xzp_output
  parameter(nvar=17)
  character(10) :: dname(nvar)
  parameter(dname = (/'x','y','z','u','v','w','p','t','rho','grad_p','grad_rho','omx','omy','omz','div','swil','Q'/))
  character(200) :: varname_output
  integer, dimension(:), allocatable :: dims
  integer :: rank
  integer :: ilen, jlen, klen
  real(8), dimension(:,:,:,:), allocatable :: buffer_tec

  call Input()
  call InitHDF5()

  do n=iil,iih,stride

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
        dname_vol_output(nn) = dname(varindex_vol_output(nn))
        varname_output = (trim(dname_vol_output(nn))//' ') !!!!!!!!!!!!!!!!!!!!!!!!!
      enddo
      print *, 'dname_output = ' , (trim(dname_vol_output(nn))//',', nn=1, nvar_vol_output)

      fname='volume_'//fnum//'.h5'

!      call InitTec(1,ilen,jlen,klen,12,"x y z u v w p t rho grad_p grad_rho div",0)


      ! write x, y, z
      call Write3DHDF5_svariable(fname, dims,"x",real(x(idx_vol%ist:idx_vol%iend:idx_vol%isp,idx_vol%jst:idx_vol%jend:idx_vol%jsp, &
                       idx_vol%kst:idx_vol%kend:idx_vol%ksp)),0)
      call Write3DHDF5_svariable(fname, dims,"y",y(idx_vol%ist:idx_vol%iend:idx_vol%isp,idx_vol%jst:idx_vol%jend:idx_vol%jsp, &
                       idx_vol%kst:idx_vol%kend:idx_vol%ksp),1)
      call Write3DHDF5_svariable(fname, dims,"z",z(idx_vol%ist:idx_vol%iend:idx_vol%isp,idx_vol%jst:idx_vol%jend:idx_vol%jsp, &
                       idx_vol%kst:idx_vol%kend:idx_vol%ksp),1)

      do nn=1, nvar_vol_output
         ! write u
        if(varindex_vol_output(nn).eq.4) then
          call Write3DHDF5_svariable(fname, dims,"u",buffer_flow(idx_vol%ist:idx_vol%iend:idx_vol%isp,idx_vol%jst:idx_vol%jend:idx_vol%jsp, &
                       idx_vol%kst:idx_vol%kend:idx_vol%ksp,1),1)
        endif
        ! write v
        if(varindex_vol_output(nn).eq.5) then
          call Write3DHDF5_svariable(fname, dims,"v",buffer_flow(idx_vol%ist:idx_vol%iend:idx_vol%isp,idx_vol%jst:idx_vol%jend:idx_vol%jsp, &
                       idx_vol%kst:idx_vol%kend:idx_vol%ksp,2),1)
        endif
        ! write w
        if(varindex_vol_output(nn).eq.6) then
          call Write3DHDF5_svariable(fname, dims,"w",buffer_flow(idx_vol%ist:idx_vol%iend:idx_vol%isp,idx_vol%jst:idx_vol%jend:idx_vol%jsp, &
                       idx_vol%kst:idx_vol%kend:idx_vol%ksp,3),1)
        endif
        ! write p
        if(varindex_vol_output(nn).eq.7) then
          call Write3DHDF5_svariable(fname, dims,"p",buffer_flow(idx_vol%ist:idx_vol%iend:idx_vol%isp,idx_vol%jst:idx_vol%jend:idx_vol%jsp, &
                       idx_vol%kst:idx_vol%kend:idx_vol%ksp,4),1)
        endif
        ! write t
        if(varindex_vol_output(nn).eq.8) then

          call Write3DHDF5_svariable(fname, dims,"t",buffer_flow(idx_vol%ist:idx_vol%iend:idx_vol%isp,idx_vol%jst:idx_vol%jend:idx_vol%jsp, &
                       idx_vol%kst:idx_vol%kend:idx_vol%ksp,5),1)
        endif
        ! write rho
        if(varindex_vol_output(nn).eq.9) then
          call Write3DHDF5_svariable(fname, dims,"rho",buffer_flow(idx_vol%ist:idx_vol%iend:idx_vol%isp,idx_vol%jst:idx_vol%jend:idx_vol%jsp, &
                       idx_vol%kst:idx_vol%kend:idx_vol%ksp,6),1)
        endif
        ! write grad_p
        if(varindex_vol_output(nn).eq.10) then
          call CalGradient(idx_vol,x,y,z,buffer_flow(:,:,:,4),grad_p)
          call Write3DHDF5_svariable(fname, dims,"grad_p",grad_p,1)
        endif
        ! write grad_rho
        if(varindex_vol_output(nn).eq.11) then
          call CalGradient(idx_vol,x,y,z,buffer_flow(:,:,:,6),grad_rho)
          call Write3DHDF5_svariable(fname, dims,"grad_rho",grad_rho,1)
        endif
        ! write omx, omy, omz
        if(varindex_vol_output(nn).eq.12) then
          call CalVorticity(idx_vol,x,y,z,buffer_flow(:,:,:,1),buffer_flow(:,:,:,2),buffer_flow(:,:,:,3),omx,omy,omz)
          call Write3DHDF5_svariable(fname, dims,"omx",omx,1)
          call Write3DHDF5_svariable(fname, dims,"omy",omy,1)
          call Write3DHDF5_svariable(fname, dims,"omz",omz,1)
        endif
        ! write div
        if(varindex_vol_output(nn).eq.15) then
          call CalDivergence(idx_vol,x,y,z,buffer_flow(:,:,:,1),buffer_flow(:,:,:,2),buffer_flow(:,:,:,3),div)
          call Write3DHDF5_svariable(fname, dims,"div",div,1)
        endif
        ! write swirl
        if(varindex_vol_output(nn).eq.16) then
          call CalSwirl(idx_vol,x,y,z,buffer_flow(:,:,:,1),buffer_flow(:,:,:,2),buffer_flow(:,:,:,3),swirl)
          call Write3DHDF5_svariable(fname, dims,"swirl",swirl,1)
        endif
        ! write Q
        if(varindex_vol_output(nn).eq.17) then
          call Callambda2(idx_vol,x,y,z,buffer_flow(:,:,:,1),buffer_flow(:,:,:,2),buffer_flow(:,:,:,3),Q)
          call Write3DHDF5_svariable(fname, dims,"Q",Q,1)
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


    if(iplane_xy.gt.0) then

      ilen = (idx_xyp%iend - idx_xyp%ist)/idx_xyp%isp + 1
      jlen = (idx_xyp%jend - idx_xyp%jst)/idx_xyp%jsp + 1
      klen = 1

      rank = 3
      if(allocated(dims)) deallocate(dims)
      allocate(dims(rank))
      dims(1) = ilen
      dims(2) = jlen
      dims(3) = klen

      varname_output = ''
      do nn=1, nvar_xyp_output
        dname_xyp_output(nn) = dname(varindex_xyp_output(nn))
        varname_output = (trim(dname_xyp_output(nn))//' ')
      enddo
      print *, 'dname_output = ' , (trim(dname_xyp_output(nn))//',', nn=1, nvar_xyp_output)

      call InitTec(1,ilen,jlen,klen,nvar_xyp_output,varname_output,0)
      if(allocated(buffer_tec)) deallocate(buffer_tec)
      allocate(buffer_tec(ilen,jlen,klen,nvar_xyp_output))

      do m=1, iplane_xy

        write(unit=fnum1,fmt='(I04.4)') planes_xy(m)
        fname = 'plane_xy_'//fnum//'_k'//fnum1//'.h5'

        idx_xyp%kst = planes_xy(m); idx_xyp%kend = planes_xy(m); idx_xyp%ksp = 1
        ! write x, y
        call Write3DHDF5_svariable(fname, dims,"x",x(idx_xyp%ist:idx_xyp%iend:idx_xyp%isp,idx_xyp%jst:idx_xyp%jend:idx_xyp%jsp, &
                                   idx_xyp%kst:idx_xyp%kend:idx_xyp%ksp),0)
        call Write3DHDF5_svariable(fname, dims,"y",y(idx_xyp%ist:idx_xyp%iend:idx_xyp%isp,idx_xyp%jst:idx_xyp%jend:idx_xyp%jsp, &
                                   idx_xyp%kst:idx_xyp%kend:idx_xyp%ksp),1)
        buffer_tec(1:ilen,1:jlen,1:klen,1) =  x(idx_xyp%ist:idx_xyp%iend:idx_xyp%isp,idx_xyp%jst:idx_xyp%jend:idx_xyp%jsp, &
                                   idx_xyp%kst:idx_xyp%kend:idx_xyp%ksp)
        buffer_tec(1:ilen,1:jlen,1:klen,2) =  y(idx_xyp%ist:idx_xyp%iend:idx_xyp%isp,idx_xyp%jst:idx_xyp%jend:idx_xyp%jsp, &
                                   idx_xyp%kst:idx_xyp%kend:idx_xyp%ksp)
        buffer_tec(1:ilen,1:jlen,1:klen,3) =  z(idx_xyp%ist:idx_xyp%iend:idx_xyp%isp,idx_xyp%jst:idx_xyp%jend:idx_xyp%jsp, &
                                   idx_xyp%kst:idx_xyp%kend:idx_xyp%ksp)
        do nn=1, nvar_xyp_output
          ! write u
          if(varindex_xyp_output(nn).eq.4) then
            call Write3DHDF5_svariable(fname, dims,"u",buffer_flow(idx_xyp%ist:idx_xyp%iend:idx_xyp%isp,idx_xyp%jst:idx_xyp%jend:idx_xyp%jsp, &
                                       idx_xyp%kst:idx_xyp%kend:idx_xyp%ksp,1),1)
            buffer_tec(1:ilen,1:jlen,1:klen,nn) = buffer_flow(idx_xyp%ist:idx_xyp%iend:idx_xyp%isp,idx_xyp%jst:idx_xyp%jend:idx_xyp%jsp, &
                                       idx_xyp%kst:idx_xyp%kend:idx_xyp%ksp,1)
          endif
          ! write v
          if(varindex_xyp_output(nn).eq.5) then
            call Write3DHDF5_svariable(fname, dims,"v",buffer_flow(idx_xyp%ist:idx_xyp%iend:idx_xyp%isp,idx_xyp%jst:idx_xyp%jend:idx_xyp%jsp, &
                                       idx_xyp%kst:idx_xyp%kend:idx_xyp%ksp,2),1)
            buffer_tec(1:ilen,1:jlen,1:klen,nn) =  buffer_flow(idx_xyp%ist:idx_xyp%iend:idx_xyp%isp,idx_xyp%jst:idx_xyp%jend:idx_xyp%jsp, &
                                       idx_xyp%kst:idx_xyp%kend:idx_xyp%ksp,2)
          endif
          ! write w
          if(varindex_xyp_output(nn).eq.6) then
            call Write3DHDF5_svariable(fname, dims,"w",buffer_flow(idx_xyp%ist:idx_xyp%iend:idx_xyp%isp,idx_xyp%jst:idx_xyp%jend:idx_xyp%jsp, &
                                       idx_xyp%kst:idx_xyp%kend:idx_xyp%ksp,3),1)
            buffer_tec(1:ilen,1:jlen,1:klen,nn) = buffer_flow(idx_xyp%ist:idx_xyp%iend:idx_xyp%isp,idx_xyp%jst:idx_xyp%jend:idx_xyp%jsp, &
                                       idx_xyp%kst:idx_xyp%kend:idx_xyp%ksp,3)
          endif
          ! write p
          if(varindex_xyp_output(nn).eq.7) then
            call Write3DHDF5_svariable(fname, dims,"p",buffer_flow(idx_xyp%ist:idx_xyp%iend:idx_xyp%isp,idx_xyp%jst:idx_xyp%jend:idx_xyp%jsp, &
                                       idx_xyp%kst:idx_xyp%kend:idx_xyp%ksp,4),1)
            buffer_tec(1:ilen,1:jlen,1:klen,nn) = buffer_flow(idx_xyp%ist:idx_xyp%iend:idx_xyp%isp,idx_xyp%jst:idx_xyp%jend:idx_xyp%jsp, &
                                       idx_xyp%kst:idx_xyp%kend:idx_xyp%ksp,4)
          endif
          ! write T
          if(varindex_xyp_output(nn).eq.8) then
            call Write3DHDF5_svariable(fname, dims,"t",buffer_flow(idx_xyp%ist:idx_xyp%iend:idx_xyp%isp,idx_xyp%jst:idx_xyp%jend:idx_xyp%jsp, &
                                       idx_xyp%kst:idx_xyp%kend:idx_xyp%ksp,5),1)
            buffer_tec(1:ilen,1:jlen,1:klen,nn) = buffer_flow(idx_xyp%ist:idx_xyp%iend:idx_xyp%isp,idx_xyp%jst:idx_xyp%jend:idx_xyp%jsp, &
                                       idx_xyp%kst:idx_xyp%kend:idx_xyp%ksp,5)
          endif
          ! write rho
          if(varindex_xyp_output(nn).eq.9) then
            call Write3DHDF5_svariable(fname, dims,"rho",buffer_flow(idx_xyp%ist:idx_xyp%iend:idx_xyp%isp,idx_xyp%jst:idx_xyp%jend:idx_xyp%jsp, &
                                       idx_xyp%kst:idx_xyp%kend:idx_xyp%ksp,6),1)
            buffer_tec(1:ilen,1:jlen,1:klen,nn) = buffer_flow(idx_xyp%ist:idx_xyp%iend:idx_xyp%isp,idx_xyp%jst:idx_xyp%jend:idx_xyp%jsp, &
                                       idx_xyp%kst:idx_xyp%kend:idx_xyp%ksp,6)
          endif

          ! write grad_p
          if(varindex_xyp_output(nn).eq.10) then
            call CalGradient(idx_xyp,x,y,z,buffer_flow(:,:,:,4),grad_p)
            call Write3DHDF5_svariable(fname, dims,"grad_p",grad_p,1)
            buffer_tec(1:ilen,1:jlen,1:klen,nn) = grad_p(1:ilen,1:jlen,1:klen)
          endif
          ! write grad_rho
          if(varindex_xyp_output(nn).eq.11) then
            call CalGradient(idx_xyp,x,y,z,buffer_flow(:,:,:,6),grad_rho)
            call Write3DHDF5_svariable(fname, dims,"grad_rho",grad_rho,1)
            buffer_tec(1:ilen,1:jlen,1:klen,nn-1) = grad_rho(1:ilen,1:jlen,1:klen)
          endif
          ! write omx, omy, omz
          if(varindex_xyp_output(nn).eq.12) then
            call CalVorticity(idx_xyp,x,y,z,buffer_flow(:,:,:,1),buffer_flow(:,:,:,2),buffer_flow(:,:,:,3),omx,omy,omz)
            call Write3DHDF5_svariable(fname, dims,"omx",omx,1)
            call Write3DHDF5_svariable(fname, dims,"omx",omy,1)
            call Write3DHDF5_svariable(fname, dims,"omx",omz,1)
            buffer_tec(1:ilen,1:jlen,1:klen,nn) = omx(1:ilen,1:jlen,1:klen)
            buffer_tec(1:ilen,1:jlen,1:klen,nn+1) = omy(1:ilen,1:jlen,1:klen)
            buffer_tec(1:ilen,1:jlen,1:klen,nn+2) = omz(1:ilen,1:jlen,1:klen)
          endif
          ! write div
          if(varindex_xyp_output(nn).eq.15) then
            call CalDivergence(idx_xyp,x,y,z,buffer_flow(:,:,:,1),buffer_flow(:,:,:,2),buffer_flow(:,:,:,3),div)
            call Write3DHDF5_svariable(fname, dims,"div",div,1)
            buffer_tec(1:ilen,1:jlen,1:klen,nn) = div(1:ilen,1:jlen,1:klen)
          endif
          ! write swirl
          if(varindex_xyp_output(nn).eq.16) then
            call CalSwirl(idx_xyp,x,y,z,buffer_flow(:,:,:,1),buffer_flow(:,:,:,2),buffer_flow(:,:,:,3),swirl)
            call Write3DHDF5_svariable(fname, dims,"swirl",swirl,1)
            buffer_tec(1:ilen,1:jlen,1:klen,nn) = swirl(1:ilen,1:jlen,1:klen)
          endif
          ! write Q
          if(varindex_xyp_output(nn).eq.17) then
            call Callambda2(idx_xyp,x,y,z,buffer_flow(:,:,:,1),buffer_flow(:,:,:,2),buffer_flow(:,:,:,3),Q)
            call Write3DHDF5_svariable(fname, dims,"Q",Q,1)
            buffer_tec(1:ilen,1:jlen,1:klen,nn) = Q(1:ilen,1:jlen,1:klen)
          endif

        enddo ! end nn loop
        fname = 'plane_xy_'//fnum//'_k'//fnum1//'.plt'
!        call WriteTec3D(trim(fname),buffer_tec)
      enddo ! end m loop

    endif

    if(iplane_yz.gt.0) then

      ilen = 1
      jlen = (idx_yzp%jend - idx_xyp%jst)/idx_xyp%jsp + 1
      klen = (idx_yzp%kend - idx_yzp%kst)/idx_yzp%ksp + 1

      rank = 3
      if(allocated(dims)) deallocate(dims)
      allocate(dims(rank))
      dims(1) = ilen
      dims(2) = jlen
      dims(3) = klen

      varname_output = ''
      do nn=1, nvar_yzp_output
        dname_yzp_output(nn) = dname(varindex_yzp_output(nn))
        varname_output = (trim(dname_yzp_output(nn))//' ')
      enddo
      print *, 'dname_output = ' , (trim(dname_yzp_output(nn))//',', nn=1, nvar_yzp_output)

     do m=1,iplane_yz
        write(unit=fnum1,fmt='(I04.4)')planes_yz(m)
        fname='plane_yz_'//fnum//'_i'//fnum1//'.h5'
        idx_yzp%ist=planes_yz(m); idx_yzp%iend=planes_yz(m); idx_yzp%isp=1
        ! write y, z
        call Write3DHDF5_svariable(fname, dims,"y",y(idx_yzp%ist:idx_yzp%iend:idx_yzp%isp,idx_yzp%jst:idx_yzp%jend:idx_yzp%jsp, &
                                   idx_yzp%kst:idx_yzp%kend:idx_yzp%ksp),0)
        call Write3DHDF5_svariable(fname, dims,"z",z(idx_yzp%ist:idx_yzp%iend:idx_yzp%isp,idx_yzp%jst:idx_yzp%jend:idx_yzp%jsp, &
                                   idx_yzp%kst:idx_yzp%kend:idx_yzp%ksp),1)
        do nn=1, nvar_yzp_output
          ! write u
          if(varindex_yzp_output(nn).eq.4) then
            call Write3DHDF5_svariable(fname, dims,"u",buffer_flow(idx_yzp%ist:idx_yzp%iend:idx_yzp%isp,idx_yzp%jst:idx_yzp%jend:idx_yzp%jsp, &
                                       idx_yzp%kst:idx_yzp%kend:idx_yzp%ksp,1),1)
          endif
          ! write v
          if(varindex_yzp_output(nn).eq.5) then
            call Write3DHDF5_svariable(fname, dims,"v",buffer_flow(idx_yzp%ist:idx_yzp%iend:idx_yzp%isp,idx_yzp%jst:idx_yzp%jend:idx_yzp%jsp, &
                                       idx_yzp%kst:idx_yzp%kend:idx_yzp%ksp,2),1)
          endif
          ! write w
          if(varindex_yzp_output(nn).eq.6) then
            call Write3DHDF5_svariable(fname, dims,"w",buffer_flow(idx_yzp%ist:idx_yzp%iend:idx_yzp%isp,idx_yzp%jst:idx_yzp%jend:idx_yzp%jsp, &
                                       idx_yzp%kst:idx_yzp%kend:idx_yzp%ksp,3),1)
          endif
          ! write p
          if(varindex_yzp_output(nn).eq.7) then
            call Write3DHDF5_svariable(fname, dims,"p",buffer_flow(idx_yzp%ist:idx_yzp%iend:idx_yzp%isp,idx_yzp%jst:idx_yzp%jend:idx_yzp%jsp, &
                                       idx_yzp%kst:idx_yzp%kend:idx_yzp%ksp,4),1)
          endif
          ! write T
          if(varindex_yzp_output(nn).eq.8) then
            call Write3DHDF5_svariable(fname, dims,"t",buffer_flow(idx_yzp%ist:idx_yzp%iend:idx_yzp%isp,idx_yzp%jst:idx_yzp%jend:idx_yzp%jsp, &
                                       idx_yzp%kst:idx_yzp%kend:idx_yzp%ksp,5),1)
          endif
          ! write rho
          if(varindex_yzp_output(nn).eq.9) then
            call Write3DHDF5_svariable(fname, dims,"rho",buffer_flow(idx_yzp%ist:idx_yzp%iend:idx_yzp%isp,idx_yzp%jst:idx_yzp%jend:idx_yzp%jsp, &
                                       idx_yzp%kst:idx_yzp%kend:idx_yzp%ksp,6),1)
          endif

          ! write grad_p
          if(varindex_yzp_output(nn).eq.10) then
            call CalGradient(idx_yzp,x,y,z,buffer_flow(:,:,:,4),grad_p)
            call Write3DHDF5_svariable(fname, dims,"grad_p",grad_p,1)
          endif
          ! write grad_rho
          if(varindex_yzp_output(nn).eq.11) then
            call CalGradient(idx_yzp,x,y,z,buffer_flow(:,:,:,6),grad_rho)
            call Write3DHDF5_svariable(fname, dims,"grad_rho",grad_rho,1)
          endif
          ! write omx, omy, omz
          if(varindex_yzp_output(nn).eq.12) then
            call CalVorticity(idx_yzp,x,y,z,buffer_flow(:,:,:,1),buffer_flow(:,:,:,2),buffer_flow(:,:,:,3),omx,omy,omz)
            call Write3DHDF5_svariable(fname, dims,"omx",omx,1)
            call Write3DHDF5_svariable(fname, dims,"omx",omy,1)
            call Write3DHDF5_svariable(fname, dims,"omx",omz,1)
          endif
          ! write div
          if(varindex_yzp_output(nn).eq.15) then
            call CalDivergence(idx_yzp,x,y,z,buffer_flow(:,:,:,1),buffer_flow(:,:,:,2),buffer_flow(:,:,:,3),div)
            call Write3DHDF5_svariable(fname, dims,"div",div,1)
          endif
          ! write swirl
          if(varindex_yzp_output(nn).eq.16) then
            call CalSwirl(idx_yzp,x,y,z,buffer_flow(:,:,:,1),buffer_flow(:,:,:,2),buffer_flow(:,:,:,3),swirl)
            call Write3DHDF5_svariable(fname, dims,"swirl",swirl,1)
          endif
          ! write Q
          if(varindex_yzp_output(nn).eq.17) then
            call Callambda2(idx_yzp,x,y,z,buffer_flow(:,:,:,1),buffer_flow(:,:,:,2),buffer_flow(:,:,:,3),Q)
            call Write3DHDF5_svariable(fname, dims,"Q",Q,1)
          endif

        enddo ! end nn loop
     enddo ! end m loop
    endif


    if(iplane_xz.gt.0) then

      ilen = (idx_xzp%iend - idx_xzp%ist)/idx_xzp%isp + 1
      jlen = 1
      klen = (idx_xzp%kend - idx_xzp%kst)/idx_xzp%ksp + 1

      rank = 3
      if(allocated(dims)) deallocate(dims)
      allocate(dims(rank))
      dims(1) = ilen
      dims(2) = jlen
      dims(3) = klen

      varname_output = ''
      do nn=1, nvar_xzp_output
        dname_xzp_output(nn) = dname(varindex_xzp_output(nn))
        varname_output = (trim(dname_xzp_output(nn))//' ')
      enddo
      print *, 'dname_output = ' , (trim(dname_xzp_output(nn))//',', nn=1, nvar_xzp_output)

      do m=1,iplane_xz
        write(unit=fnum1,fmt='(I04.4)') planes_xz(m)
        fname='plane_xz_'//fnum//'_j'//fnum1//'.h5'
        idx_xzp%jst=planes_xz(m); idx_xzp%jend=planes_xz(m); idx_xzp%jsp=1

        ! write x, z
        call Write3DHDF5_svariable(fname, dims,"x",x(idx_xzp%ist:idx_xzp%iend:idx_xzp%isp,idx_xzp%jst:idx_xzp%jend:idx_xzp%jsp, &
                                   idx_xzp%kst:idx_xzp%kend:idx_xzp%ksp),0)
        call Write3DHDF5_svariable(fname, dims,"y",z(idx_xzp%ist:idx_xzp%iend:idx_xzp%isp,idx_xzp%jst:idx_xzp%jend:idx_xzp%jsp, &
                                   idx_xzp%kst:idx_xzp%kend:idx_xzp%ksp),1)
        do nn=1, nvar_xzp_output
          ! write u
          if(varindex_xzp_output(nn).eq.4) then
            call Write3DHDF5_svariable(fname, dims,"u",buffer_flow(idx_xzp%ist:idx_xzp%iend:idx_xzp%isp,idx_xzp%jst:idx_xzp%jend:idx_xzp%jsp, &
                                       idx_xzp%kst:idx_xzp%kend:idx_xzp%ksp,1),1)
          endif
          ! write v
          if(varindex_xyp_output(nn).eq.5) then
            call Write3DHDF5_svariable(fname, dims,"v",buffer_flow(idx_xzp%ist:idx_xzp%iend:idx_xzp%isp,idx_xzp%jst:idx_xzp%jend:idx_xzp%jsp, &
                                       idx_xzp%kst:idx_xzp%kend:idx_xzp%ksp,2),1)
          endif
          ! write w
          if(varindex_xyp_output(nn).eq.6) then
            call Write3DHDF5_svariable(fname, dims,"w",buffer_flow(idx_xzp%ist:idx_xzp%iend:idx_xzp%isp,idx_xzp%jst:idx_xzp%jend:idx_xzp%jsp, &
                                       idx_xzp%kst:idx_xzp%kend:idx_xzp%ksp,3),1)
          endif
          ! write p
          if(varindex_xyp_output(nn).eq.7) then
            call Write3DHDF5_svariable(fname, dims,"p",buffer_flow(idx_xzp%ist:idx_xzp%iend:idx_xzp%isp,idx_xzp%jst:idx_xzp%jend:idx_xzp%jsp, &
                                       idx_xzp%kst:idx_xzp%kend:idx_xzp%ksp,4),1)
          endif
          ! write T
          if(varindex_xyp_output(nn).eq.8) then
            call Write3DHDF5_svariable(fname, dims,"t",buffer_flow(idx_xzp%ist:idx_xzp%iend:idx_xzp%isp,idx_xzp%jst:idx_xzp%jend:idx_xzp%jsp, &
                                       idx_xzp%kst:idx_xzp%kend:idx_xzp%ksp,5),1)
          endif
          ! write rho
          if(varindex_xyp_output(nn).eq.9) then
            call Write3DHDF5_svariable(fname, dims,"rho",buffer_flow(idx_xzp%ist:idx_xzp%iend:idx_xzp%isp,idx_xzp%jst:idx_xzp%jend:idx_xzp%jsp, &
                                       idx_xzp%kst:idx_xzp%kend:idx_xzp%ksp,6),1)
          endif

          ! write grad_p
          if(varindex_xzp_output(nn).eq.10) then
            call CalGradient(idx_xzp,x,y,z,buffer_flow(:,:,:,4),grad_p)
            call Write3DHDF5_svariable(fname, dims,"grad_p",grad_p,1)
          endif
          ! write grad_rho
          if(varindex_xzp_output(nn).eq.11) then
            call CalGradient(idx_xzp,x,y,z,buffer_flow(:,:,:,6),grad_rho)
            call Write3DHDF5_svariable(fname, dims,"grad_rho",grad_rho,1)
          endif
          ! write omx, omy, omz
          if(varindex_xzp_output(nn).eq.12) then
            call CalVorticity(idx_xzp,x,y,z,buffer_flow(:,:,:,1),buffer_flow(:,:,:,2),buffer_flow(:,:,:,3),omx,omy,omz)
            call Write3DHDF5_svariable(fname, dims,"omx",omx,1)
            call Write3DHDF5_svariable(fname, dims,"omx",omy,1)
            call Write3DHDF5_svariable(fname, dims,"omx",omz,1)
          endif
          ! write div
          if(varindex_xzp_output(nn).eq.15) then
            call CalDivergence(idx_xzp,x,y,z,buffer_flow(:,:,:,1),buffer_flow(:,:,:,2),buffer_flow(:,:,:,3),div)
            call Write3DHDF5_svariable(fname, dims,"div",div,1)
          endif
          ! write swirl
          if(varindex_xzp_output(nn).eq.16) then
            call CalSwirl(idx_xzp,x,y,z,buffer_flow(:,:,:,1),buffer_flow(:,:,:,2),buffer_flow(:,:,:,3),swirl)
            call Write3DHDF5_svariable(fname, dims,"swirl",swirl,1)
          endif
          ! write Q
          if(varindex_xzp_output(nn).eq.17) then
            call Callambda2(idx_xzp,x,y,z,buffer_flow(:,:,:,1),buffer_flow(:,:,:,2),buffer_flow(:,:,:,3),Q)
            call Write3DHDF5_svariable(fname, dims,"Q",Q,1)
          endif
         enddo ! end nn loop
      enddo ! end m loop
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
      read(*,*)iil,iih,stride
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
      read(*,*)
      read(*,*) nvar_xyp_output
      allocate(varindex_xyp_output(nvar_xyp_output), dname_xyp_output(nvar_xyp_output))
      read(*,*) (varindex_xyp_output(i), i=1,nvar_xyp_output)

!     yz planes
      read(*,*)
      read(*,*)iplane_yz, idx_yzp%jst, idx_yzp%jend, idx_yzp%jsp, idx_yzp%kst, idx_yzp%kend, idx_yzp%ksp
      if (iplane_yz.gt.0) then
        allocate(planes_yz(1:iplane_yz))
        read(*,*)(planes_yz(i),i=1,iplane_yz)
      else
        read(*,*)
      end if
      read(*,*)
      read(*,*) nvar_yzp_output
      allocate(varindex_yzp_output(nvar_yzp_output), dname_yzp_output(nvar_yzp_output))
      read(*,*) (varindex_yzp_output(i), i=1,nvar_yzp_output)
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
      read(*,*) nvar_xzp_output
      allocate(varindex_xzp_output(nvar_xzp_output), dname_xzp_output(nvar_xzp_output))
      read(*,*) (varindex_xzp_output(i), i=1,nvar_xzp_output)

      call ReadPlot3DGridPlaneGen(trim(datapath)//'gridp3d.grd',imax,jmax,kmax,x,y,z,0)

      allocate(buffer_flow(imax,jmax,kmax,6))


      return
    end subroutine Input


end program volume    
