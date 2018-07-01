program volume
  use MFileIO
  use modVolume
  use modRWHDF5
  use modTecbin
  implicit none

  real(8), parameter :: R=8314.3D0
  character(400) :: fname
  character(400) :: datapath
  character(8) :: fnum 
  character(4) :: fnum1, fnum4_1, fnum4_2
  integer :: imax, jmax, kmax
  integer :: file_be, file_end, file_skip, num_file, icylin, iSpanave
  integer :: ivolume, iplane_xz, iplane_yz, iplane_xy
  type(index_output) :: idx_vol, idx_xyp, idx_yzp, idx_xzp
  integer, dimension(:),  allocatable :: planes_xy, planes_yz, planes_xz
  real, dimension(:,:,:), allocatable :: x, y, z
  real, dimension(:,:,:,:), allocatable :: buffer_grid
  real, dimension(:,:,:,:), allocatable :: buffer_flow, buffer_flow_tmp, buffer_flow_tmp2
  real, dimension(:,:,:), allocatable :: grad_p, grad_rho
  real, dimension(:,:,:), allocatable :: omx, omy, omz
  real, dimension(:,:,:), allocatable :: swirl
  real, dimension(:,:,:), allocatable :: div
  real, dimension(:,:,:), allocatable :: Q, src_Philip, p0
  real, dimension(:,:,:), allocatable :: dudn, dtdn, ds
  integer :: n, i, j, k, m, nn
  real(8) :: uinf, delta, utau, ztau
  real(8) :: rbar, Rm
  integer, dimension(:), allocatable :: varindex_vol_output, varindex_xyp_output, varindex_yzp_output, varindex_xzp_output
  character(10), dimension(:), allocatable :: dname_vol_output, dname_xyp_output, dname_yzp_output, dname_xzp_output

  integer :: nvar, nvar_vol_output, nvar_xyp_output, nvar_yzp_output, nvar_xzp_output
  parameter(nvar=22)
  character(10) :: dname(nvar)
  parameter(dname = (/'x','y','z','u','v','w','p','t','rho','grad_p','grad_rho','omx','omy','omz','div','swil','Q', 'src_Philip', 'p0', 'ru', 'rv', 'rw'/))
  character(400) :: varname_output, varname_dat
  integer, dimension(:), allocatable :: dims
  integer :: rank
  integer :: ilen, jlen, klen
  real(8), dimension(:,:,:,:), allocatable :: buffer_tec

  type(tp_rdwt_hdf5) :: fsol, grd
  integer :: iHDF5, iAppend, iBaseflowSubtract
  integer :: tmp, cal_2D
  real(8), dimension(:,:,:), allocatable :: var_dat, var2_dat



  call InitHDF5()
  call Input()

  call InitFlowHDF5(fsol)
  fsol%dimsf = (/kmax,imax,jmax/)

  num_file = (file_end-file_be)/file_skip + 1
  write(unit=fnum4_1,fmt='(I04.4)') file_be/1000
  write(unit=fnum4_2,fmt='(I04.4)') file_end/1000

  if(num_file.eq.1.and.iAppend.eq.1) then
    print *, 'num_file cannot be 1 if iAppend=1 ... Stop! '
    print *, 'num_file = ', num_file, 'iAppend = ', iAppend
    stop
  endif

  if(iBaseflowSubtract.eq.1) then
     if(iHDF5.eq.0) then
        allocate(buffer_flow_tmp(imax,jmax,kmax,6))
        fsol%fname = trim(datapath)//'baseflow.sol'
        print *, 'Reading baseflow file: ', trim(fsol%fname)        
        call ReadPlot3DSolPlane(trim(fsol%fname),imax,jmax,kmax,6,buffer_flow_tmp)
      else
        allocate(buffer_flow_tmp(kmax,imax,jmax,5))
        fsol%fname = trim(datapath)//'baseflow.h5'
        print *, 'Reading file: ', trim(fsol%fname)
        call ReadHDF5_3D(fsol,buffer_flow_tmp)
      endif  
  endif
  


  do n=1, num_file
     write(unit=fnum,fmt='(I08.8)') file_be + (n-1)*file_skip
     if(iHDF5.eq.0) then
       fsol%fname = trim(datapath)//'flowdata_'//fnum//'.sol'
       print *, 'Reading file: ', trim(fsol%fname)
       call ReadPlot3DSolPlane(trim(fsol%fname),imax,jmax,kmax,6,buffer_flow)
       if(iBaseflowSubtract.eq.1) then
!           allocate(buffer_flow_tmp(imax,jmax,kmax,6))
!           fsol%fname = trim(datapath)//'baseflow.sol'
!            print *, 'Reading baseflow file: ', trim(fsol%fname)        
!            call ReadPlot3DSolPlane(trim(fsol%fname),imax,jmax,kmax,6,buffer_flow_tmp)
!            print *, 'Subtracting baseflow' 
           buffer_flow = buffer_flow - buffer_flow_tmp
!           deallocate(buffer_flow_tmp)            
       endif   
     else
       allocate(buffer_flow_tmp2(kmax,imax,jmax,5))
       fsol%fname = trim(datapath)//'flowdata_'//fnum//'.h5'
       print *, 'Reading file: ', trim(fsol%fname)
       call ReadHDF5_3D(fsol,buffer_flow_tmp2)

       do i=1, imax
         do j=1, jmax
           do k=1, kmax
             buffer_flow(i,j,k,1) = buffer_flow_tmp2(k,i,j,1)
             buffer_flow(i,j,k,2) = buffer_flow_tmp2(k,i,j,2)
             buffer_flow(i,j,k,3) = buffer_flow_tmp2(k,i,j,3)
             buffer_flow(i,j,k,4) = buffer_flow_tmp2(k,i,j,4)
             buffer_flow(i,j,k,5) = buffer_flow_tmp2(k,i,j,5)
           enddo
         enddo
       enddo
       buffer_flow(:,:,:,6) = buffer_flow(:,:,:,4)/(rbar*buffer_flow(:,:,:,5))
       deallocate(buffer_flow_tmp2)
       
       if(iBaseflowSubtract.eq.1) then
!          fsol%fname = trim(datapath)//'baseflow.h5'
!          print *, 'Reading baseflow file: ', trim(fsol%fname)        
!          call ReadHDF5_3D(fsol,buffer_flow_tmp)
          print *, 'Subtracting baseflow' 
          do i=1, imax
             do j=1, jmax
                do k=1, kmax
                   buffer_flow(i,j,k,1) = buffer_flow(i,j,k,1) - buffer_flow_tmp(k,i,j,1)
                   buffer_flow(i,j,k,2) = buffer_flow(i,j,k,2) - buffer_flow_tmp(k,i,j,2)
                   buffer_flow(i,j,k,3) = buffer_flow(i,j,k,3) - buffer_flow_tmp(k,i,j,3)
                   buffer_flow(i,j,k,4) = buffer_flow(i,j,k,4) - buffer_flow_tmp(k,i,j,4)
                   buffer_flow(i,j,k,5) = buffer_flow(i,j,k,5) - buffer_flow_tmp(k,i,j,5)
                   buffer_flow(i,j,k,6) = buffer_flow(i,j,k,6) &
                                        - buffer_flow_tmp(k,i,j,4)/(rbar*buffer_flow_tmp(k,i,j,5))
                enddo
             enddo
           enddo      
       endif                             
!       deallocate(buffer_flow_tmp)
            
     endif
    ! if(iBaseflowSubtract.eq.1) deallocate(buffer_flow_tmp)
 
     if(ivolume.gt.0) then

      ilen = (idx_vol%iend - idx_vol%ist)/idx_vol%isp + 1
      jlen = (idx_vol%jend - idx_vol%jst)/idx_vol%jsp + 1
      klen = (idx_vol%kend - idx_vol%kst)/idx_vol%ksp + 1

      varname_output = ''
      do nn=1, nvar_vol_output
        dname_vol_output(nn) = dname(varindex_vol_output(nn))
        varname_output = trim(adjustl(varname_output))//' '//trim(adjustl(dname_vol_output(nn)))
      enddo
      print *, 'Calculate volume'
      print *, 'varname_output = ', trim(varname_output)

      if(iSpanave.eq.1.and.n.eq.1) allocate(var_dat(ilen,klen,nvar_vol_output),var2_dat(ilen,klen,nvar_vol_output))

      if(iAppend.eq.0) then
        if(icylin .eq. 1 .and. idx_vol%jend .eq. jmax) then
          call InitTec(1,ilen,jlen+1,klen,nvar_vol_output,trim(varname_output),0)
        else
          call InitTec(1,ilen,jlen,klen,nvar_vol_output,trim(varname_output),0)
        endif
      elseif(iAppend.eq.1.and.n.eq.1) then
        do i=1, 3
          if(varindex_vol_output(i).le.3)then
            tmp = i
          endif
        enddo
        print *, ' number of share variable = ', tmp
        if(icylin .eq. 1 .and. idx_vol%jend .eq. jmax) then
          call InitTec(num_file,ilen,jlen+1,klen,nvar_vol_output,trim(varname_output),tmp)
        else
          call InitTec(num_file,ilen,jlen,klen,nvar_vol_output,trim(varname_output),tmp)
        endif
      endif

      fname='volume_'//fnum//'.plt'

      do nn=1, nvar_vol_output
        ! write x, y, z
        if(varindex_vol_output(nn).eq.1) then
          vartmp(1:ilen,1:jlen,1:klen,nn) = x(idx_vol%ist:idx_vol%iend:idx_vol%isp,idx_vol%jst:idx_vol%jend:idx_vol%jsp, &
                                              idx_vol%kst:idx_vol%kend:idx_vol%ksp)
        endif
        if(varindex_vol_output(nn).eq.2) then
          if(icylin.eq.0) then
            vartmp(1:ilen,1:jlen,1:klen,nn) = y(idx_vol%ist:idx_vol%iend:idx_vol%isp,idx_vol%jst:idx_vol%jend:idx_vol%jsp, &
                                                idx_vol%kst:idx_vol%kend:idx_vol%ksp)
          else
            vartmp(1:ilen,1:jlen,1:klen,nn) = z(idx_vol%ist:idx_vol%iend:idx_vol%isp, idx_vol%jst:idx_vol%jend:idx_vol%jsp, &
                                                idx_vol%kst:idx_vol%kend:idx_vol%ksp)*cos(y(idx_vol%ist:idx_vol%iend:idx_vol%isp, &
                                                idx_vol%jst:idx_vol%jend:idx_vol%jsp, idx_vol%kst:idx_vol%kend:idx_vol%ksp))
          endif
        endif
        if(varindex_vol_output(nn).eq.3) then
          if(icylin.eq.0) then
            vartmp(1:ilen,1:jlen,1:klen,nn) = z(idx_vol%ist:idx_vol%iend:idx_vol%isp,idx_vol%jst:idx_vol%jend:idx_vol%jsp, &
                                                idx_vol%kst:idx_vol%kend:idx_vol%ksp)
          else
            vartmp(1:ilen,1:jlen,1:klen,nn) = z(idx_vol%ist:idx_vol%iend:idx_vol%isp,idx_vol%jst:idx_vol%jend:idx_vol%jsp, &
                                                idx_vol%kst:idx_vol%kend:idx_vol%ksp)*sin(y(idx_vol%ist:idx_vol%iend:idx_vol%isp, &
                                                idx_vol%jst:idx_vol%jend:idx_vol%jsp, idx_vol%kst:idx_vol%kend:idx_vol%ksp))
          endif
        endif

         ! write u
        if(varindex_vol_output(nn).eq.4) then
          vartmp(1:ilen,1:jlen,1:klen,nn) = buffer_flow(idx_vol%ist:idx_vol%iend:idx_vol%isp,idx_vol%jst:idx_vol%jend:idx_vol%jsp, &
                                                        idx_vol%kst:idx_vol%kend:idx_vol%ksp,1)
        endif
        ! write v
        if(varindex_vol_output(nn).eq.5) then
          vartmp(1:ilen,1:jlen,1:klen,nn) = buffer_flow(idx_vol%ist:idx_vol%iend:idx_vol%isp,idx_vol%jst:idx_vol%jend:idx_vol%jsp, &
                                                        idx_vol%kst:idx_vol%kend:idx_vol%ksp,2)
        endif
        ! write w
        if(varindex_vol_output(nn).eq.6) then
          vartmp(1:ilen,1:jlen,1:klen,nn) = buffer_flow(idx_vol%ist:idx_vol%iend:idx_vol%isp,idx_vol%jst:idx_vol%jend:idx_vol%jsp, &
                                                        idx_vol%kst:idx_vol%kend:idx_vol%ksp,3)
        endif
        ! write p
        if(varindex_vol_output(nn).eq.7) then
          vartmp(1:ilen,1:jlen,1:klen,nn) = buffer_flow(idx_vol%ist:idx_vol%iend:idx_vol%isp,idx_vol%jst:idx_vol%jend:idx_vol%jsp, &
                                                        idx_vol%kst:idx_vol%kend:idx_vol%ksp,4)
        endif
        ! write t
        if(varindex_vol_output(nn).eq.8) then
          vartmp(1:ilen,1:jlen,1:klen,nn) = buffer_flow(idx_vol%ist:idx_vol%iend:idx_vol%isp,idx_vol%jst:idx_vol%jend:idx_vol%jsp, &
                                                        idx_vol%kst:idx_vol%kend:idx_vol%ksp,5)
        endif
        ! write rho
        if(varindex_vol_output(nn).eq.9) then
          vartmp(1:ilen,1:jlen,1:klen,nn) = buffer_flow(idx_vol%ist:idx_vol%iend:idx_vol%isp,idx_vol%jst:idx_vol%jend:idx_vol%jsp, &
                                                        idx_vol%kst:idx_vol%kend:idx_vol%ksp,6)
        endif
        ! write grad_p
        if(varindex_vol_output(nn).eq.10) then
          if(icylin.eq.0) then
            call CalGradient(idx_vol,x,y,z,buffer_flow(:,:,:,4),grad_p)
          else
            call CalGradient_cylin(idx_vol,x,y,z,buffer_flow(:,:,:,4),grad_p)
            grad_p(1:ilen,1:jlen,1) = grad_p(1:ilen,1:jlen,2)
          endif
          vartmp(1:ilen,1:jlen,1:klen,nn) = grad_p(1:ilen,1:jlen,1:klen)
        endif
        ! write grad_rho
        if(varindex_vol_output(nn).eq.11) then
          if(icylin.eq.0) then
            call CalGradient(idx_vol,x,y,z,buffer_flow(:,:,:,6),grad_rho)
          else
            call CalGradient_cylin(idx_vol,x,y,z,buffer_flow(:,:,:,6),grad_rho)
            grad_rho(1:ilen,1:jlen,1) = grad_rho(1:ilen,1:jlen,2)
          endif
          vartmp(1:ilen,1:jlen,1:klen,nn) = grad_rho(1:ilen,1:jlen,1:klen)
        endif
        ! write omx, omy, omz
        if(varindex_vol_output(nn).eq.12) then
          if(icylin.eq.0) then
            call CalVorticity(idx_vol,x,y,z,buffer_flow(:,:,:,1),buffer_flow(:,:,:,2),buffer_flow(:,:,:,3),omx,omy,omz)
          else
            call CalVorticity_cylin(idx_vol,x,y,z,buffer_flow(:,:,:,1),buffer_flow(:,:,:,2),buffer_flow(:,:,:,3),omx,omy,omz)
            omx(1:ilen,1:jlen,1) = omx(1:ilen,1:jlen,2)
            omy(1:ilen,1:jlen,1) = omy(1:ilen,1:jlen,2)
            omz(1:ilen,1:jlen,1) = omz(1:ilen,1:jlen,2)
          endif
          vartmp(1:ilen,1:jlen,1:klen,nn)   = omx(1:ilen,1:jlen,1:klen)
          vartmp(1:ilen,1:jlen,1:klen,nn+1) = omy(1:ilen,1:jlen,1:klen)
          vartmp(1:ilen,1:jlen,1:klen,nn+2) = omz(1:ilen,1:jlen,1:klen)
        endif
        ! write div
        if(varindex_vol_output(nn).eq.15) then
          if(icylin.eq.0) then
            call CalDivergence(idx_vol,x,y,z,buffer_flow(:,:,:,1),buffer_flow(:,:,:,2),buffer_flow(:,:,:,3),div)
          else
            call CalDivergence_cylin(idx_vol,x,y,z,buffer_flow(:,:,:,1),buffer_flow(:,:,:,2),buffer_flow(:,:,:,3),div)
            div(1:ilen,1:jlen,1) = div(1:ilen,1:jlen,2)
          endif
          vartmp(1:ilen,1:jlen,1:klen,nn) = div(1:ilen,1:jlen,1:klen)
        endif
        ! write swirl
        if(varindex_vol_output(nn).eq.16) then
          call CalSwirl(idx_vol,x,y,z,buffer_flow(:,:,:,1),buffer_flow(:,:,:,2),buffer_flow(:,:,:,3),swirl)
          vartmp(1:ilen,1:jlen,1:klen,nn) = swirl(1:ilen,1:jlen,1:klen)
        endif
        ! write Q
        if(varindex_vol_output(nn).eq.17) then
          call Callambda2(idx_vol,x,y,z,buffer_flow(:,:,:,1),buffer_flow(:,:,:,2),buffer_flow(:,:,:,3),Q)
          vartmp(1:ilen,1:jlen,1:klen,nn) = Q(1:ilen,1:jlen,1:klen)
        endif

        ! write src_Philip
        if(varindex_vol_output(nn).eq.18) then
          call Calsrc_Philip(idx_vol,x,y,z,buffer_flow(:,:,:,1),buffer_flow(:,:,:,2),buffer_flow(:,:,:,3),src_Philip)
          vartmp(1:ilen,1:jlen,1:klen,nn) = src_Philip(1:ilen,1:jlen,1:klen)
        endif

        ! write p0
        if(varindex_vol_output(nn).eq.19) then
          call Calp0(idx_vol,buffer_flow(:,:,:,1),buffer_flow(:,:,:,2),buffer_flow(:,:,:,3),buffer_flow(:,:,:,4),buffer_flow(:,:,:,5),p0)
          vartmp(1:ilen,1:jlen,1:klen,nn) = p0(1:ilen,1:jlen,1:klen)
        endif

        ! write ru
        !if(varindex_vol_output(nn).eq.20) then
        !  vartmp(1:ilen,1:jlen,1:klen,nn) = buffer_flow(idx_vol%ist:idx_vol%iend:idx_vol%isp,idx_vol%jst:idx_vol%jend:idx_vol%jsp, &
        !                                                idx_vol%kst:idx_vol%kend:idx_vol%ksp,1)
        !endif

        ! write rv
        !if(varindex_vol_output(nn).eq.4) then
        !  vartmp(1:ilen,1:jlen,1:klen,nn) = buffer_flow(idx_vol%ist:idx_vol%iend:idx_vol%isp,idx_vol%jst:idx_vol%jend:idx_vol%jsp, &
        !                                                idx_vol%kst:idx_vol%kend:idx_vol%ksp,1)
        !endif

        ! write rw
        !if(varindex_vol_output(nn).eq.4) then
        !  vartmp(1:ilen,1:jlen,1:klen,nn) = buffer_flow(idx_vol%ist:idx_vol%iend:idx_vol%isp,idx_vol%jst:idx_vol%jend:idx_vol%jsp, &
        !                                                idx_vol%kst:idx_vol%kend:idx_vol%ksp,1)
        !endif

      enddo

      if(icylin .eq. 1 .and. idx_vol%jend .eq. jmax) vartmp(:,jlen+1,:,:) = vartmp(:,1,:,:)

      if(iAppend.eq.0) then
        call WriteTec(fname)
      else
        fname = 'volume_'//fnum4_1//'-'//fnum4_2//'.plt'
        print *, 'Writing file: ', trim(fname)
        call WriteTec(fname,n)
      endif

      if(iSpanave.eq.1) then
        fname='volume_'//fnum//'.dat'
        varname_dat = 'variables='//trim(varname_output)
        print *, 'varname_dat = ', trim(varname_dat)
        var_dat  = sum(vartmp,dim=2)/dble(jlen)
        var2_dat = sum(vartmp**2,dim=2)/dble(jlen)
        open(7,file=trim(fname),status='unknown')
          write(7,'(2a)') 'variables=x,z,omx,omy,omz,div,omx2,omy2,omz2,div2'
          write(7,'("Zone T=AveSpanStat, I=", I4, ",J=",I4)') klen, ilen
          do i=1, ilen
            do k=1, klen
              write(7,*) var_dat(i,k,:), var2_dat(i,k,3:)
            enddo
          enddo
        close(7)

      endif



    endif ! end ivolume.gt.0


    if(iplane_xy.gt.0) then

      ilen = (idx_xyp%iend - idx_xyp%ist)/idx_xyp%isp + 1
      jlen = (idx_xyp%jend - idx_xyp%jst)/idx_xyp%jsp + 1
      klen = 1

      varname_output = ''
      do nn=1, nvar_xyp_output
        dname_xyp_output(nn) = dname(varindex_xyp_output(nn))
        varname_output = trim(adjustl(varname_output))//' '//trim(adjustl(dname_xyp_output(nn)))
      enddo

      print *, 'Calculate x-y plane'
      print *, 'varname_output = ', trim(varname_output)

      if(iAppend.eq.0) then
        call InitTec(1,ilen,jlen,klen,nvar_xyp_output,varname_output,0)
      elseif(iAppend.eq.1.and.n.eq.1) then
        do i=1, 3
          if(varindex_xyp_output(i).le.3)then
            tmp = i
          endif
        enddo
        print *, ' number of share variable = ', tmp
        call InitTec(num_file,ilen,jlen,klen,nvar_xyp_output,trim(varname_output),tmp)
      endif

      do m=1, iplane_xy

        write(unit=fnum1,fmt='(I04.4)') planes_xy(m)
        if(iAppend.eq.0) then
          fname = 'plane_xy_'//fnum//'_k'//fnum1//'.plt'
        endif
        idx_xyp%kst = planes_xy(m); idx_xyp%kend = planes_xy(m); idx_xyp%ksp = 1

        do nn=1, nvar_xyp_output
           ! write x, y
           if(varindex_xyp_output(nn).eq.1) then
             vartmp(1:ilen,1:jlen,1:klen,nn) = x(idx_xyp%ist:idx_xyp%iend:idx_xyp%isp,idx_xyp%jst:idx_xyp%jend:idx_xyp%jsp, &
                                                idx_xyp%kst:idx_xyp%kend:idx_xyp%ksp)
           endif
           if(varindex_xyp_output(nn).eq.2) then
             vartmp(1:ilen,1:jlen,1:klen,nn) = y(idx_xyp%ist:idx_xyp%iend:idx_xyp%isp,idx_xyp%jst:idx_xyp%jend:idx_xyp%jsp, &
                                                idx_xyp%kst:idx_xyp%kend:idx_xyp%ksp)
           endif
           if(varindex_xyp_output(nn).eq.3) then
             vartmp(1:ilen,1:jlen,1:klen,nn) = z(idx_xyp%ist:idx_xyp%iend:idx_xyp%isp,idx_xyp%jst:idx_xyp%jend:idx_xyp%jsp, &
                                                idx_xyp%kst:idx_xyp%kend:idx_xyp%ksp)
           endif
          ! write u
          if(varindex_xyp_output(nn).eq.4) then
            vartmp(1:ilen,1:jlen,1:klen,nn) = buffer_flow(idx_xyp%ist:idx_xyp%iend:idx_xyp%isp,idx_xyp%jst:idx_xyp%jend:idx_xyp%jsp, &
                                                          idx_xyp%kst:idx_xyp%kend:idx_xyp%ksp,1)
          endif
          ! write v
          if(varindex_xyp_output(nn).eq.5) then
            vartmp(1:ilen,1:jlen,1:klen,nn) = buffer_flow(idx_xyp%ist:idx_xyp%iend:idx_xyp%isp,idx_xyp%jst:idx_xyp%jend:idx_xyp%jsp, &
                                                          idx_xyp%kst:idx_xyp%kend:idx_xyp%ksp,2)
          endif
          ! write w
          if(varindex_xyp_output(nn).eq.6) then
            vartmp(1:ilen,1:jlen,1:klen,nn) = buffer_flow(idx_xyp%ist:idx_xyp%iend:idx_xyp%isp,idx_xyp%jst:idx_xyp%jend:idx_xyp%jsp, &
                                                          idx_xyp%kst:idx_xyp%kend:idx_xyp%ksp,3)
          endif
          ! write p
          if(varindex_xyp_output(nn).eq.7) then
            vartmp(1:ilen,1:jlen,1:klen,nn) = buffer_flow(idx_xyp%ist:idx_xyp%iend:idx_xyp%isp,idx_xyp%jst:idx_xyp%jend:idx_xyp%jsp, &
                                                          idx_xyp%kst:idx_xyp%kend:idx_xyp%ksp,4)

          endif
          ! write T
          if(varindex_xyp_output(nn).eq.8) then
            vartmp(1:ilen,1:jlen,1:klen,nn) = buffer_flow(idx_xyp%ist:idx_xyp%iend:idx_xyp%isp,idx_xyp%jst:idx_xyp%jend:idx_xyp%jsp, &
                                                          idx_xyp%kst:idx_xyp%kend:idx_xyp%ksp,5)
          endif
          ! write rho
          if(varindex_xyp_output(nn).eq.9) then
            vartmp(1:ilen,1:jlen,1:klen,nn) = buffer_flow(idx_xyp%ist:idx_xyp%iend:idx_xyp%isp,idx_xyp%jst:idx_xyp%jend:idx_xyp%jsp, &
                                                          idx_xyp%kst:idx_xyp%kend:idx_xyp%ksp,6)
          endif

          ! write grad_p
          if(varindex_xyp_output(nn).eq.10) then
            call CalGradient(idx_xyp,x,y,z,buffer_flow(:,:,:,4),grad_p)
            vartmp(1:ilen,1:jlen,1:klen,nn) = grad_p(1:ilen,1:jlen,1:klen)
          endif
          ! write grad_rho
          if(varindex_xyp_output(nn).eq.11) then
            call CalGradient(idx_xyp,x,y,z,buffer_flow(:,:,:,6),grad_rho)
            vartmp(1:ilen,1:jlen,1:klen,nn) = grad_rho(1:ilen,1:jlen,1:klen)
          endif
          ! write omx, omy, omz
          if(varindex_xyp_output(nn).eq.12) then
            call CalVorticity(idx_xyp,x,y,z,buffer_flow(:,:,:,1),buffer_flow(:,:,:,2),buffer_flow(:,:,:,3),omx,omy,omz)
            vartmp(1:ilen,1:jlen,1:klen,nn) = omx(1:ilen,1:jlen,1:klen)
            vartmp(1:ilen,1:jlen,1:klen,nn+1) = omy(1:ilen,1:jlen,1:klen)
            vartmp(1:ilen,1:jlen,1:klen,nn+2) = omz(1:ilen,1:jlen,1:klen)
          endif
          ! write div
          if(varindex_xyp_output(nn).eq.15) then
            call CalDivergence(idx_xyp,x,y,z,buffer_flow(:,:,:,1),buffer_flow(:,:,:,2),buffer_flow(:,:,:,3),div)
            vartmp(1:ilen,1:jlen,1:klen,nn) = div(1:ilen,1:jlen,1:klen)
          endif
          ! write swirl
          if(varindex_xyp_output(nn).eq.16) then
            call CalSwirl(idx_xyp,x,y,z,buffer_flow(:,:,:,1),buffer_flow(:,:,:,2),buffer_flow(:,:,:,3),swirl)
            vartmp(1:ilen,1:jlen,1:klen,nn) = swirl(1:ilen,1:jlen,1:klen)
          endif
          ! write Q
          if(varindex_xyp_output(nn).eq.17) then
            call Callambda2(idx_xyp,x,y,z,buffer_flow(:,:,:,1),buffer_flow(:,:,:,2),buffer_flow(:,:,:,3),Q)
            vartmp(1:ilen,1:jlen,1:klen,nn) = Q(1:ilen,1:jlen,1:klen)
          endif
        enddo ! end nn loop
        if(iAppend.eq.0) then
          call WriteTec(fname)
        else
          fname = 'plane_xy_'//fnum4_1//'-'//fnum4_2//'_k'//fnum1//'.plt'
          print *, 'Writing file: ', trim(fname)
          call WriteTec(fname,n)
        endif
      enddo ! end m loop

    endif ! end iplane_xy.gt.0

    if(iplane_yz.gt.0) then

      ilen = 1
      jlen = (idx_yzp%jend - idx_yzp%jst)/idx_yzp%jsp + 1
      klen = (idx_yzp%kend - idx_yzp%kst)/idx_yzp%ksp + 1

      varname_output = ''
      do nn=1, nvar_yzp_output
        dname_yzp_output(nn) = dname(varindex_yzp_output(nn))
        varname_output = trim(adjustl(varname_output))//' '//trim(adjustl(dname_yzp_output(nn)))
      enddo
      print *, 'Calculate y-z plane'
      print *, 'varname_output = ', trim(varname_output)

      if(iAppend.eq.0) then
        if(icylin .eq. 1 .and. idx_yzp%jend .eq. jmax) then
          call InitTec(1,ilen,jlen+1,klen,nvar_yzp_output,trim(varname_output),0)
        else
          call InitTec(1,ilen,jlen,klen,nvar_yzp_output,trim(varname_output),0)
        endif
      elseif(iAppend.eq.1.and.n.eq.1) then
        do i=1, 3
          if(varindex_yzp_output(i).le.3)then
            tmp = i
          endif
        enddo
        print *, ' number of share variable = ', tmp
        if(icylin .eq. 1 .and. idx_yzp%jend .eq. jmax) then
          call InitTec(num_file,ilen,jlen+1,klen,nvar_yzp_output,trim(varname_output),tmp)
        else
          call InitTec(num_file,ilen,jlen,klen,nvar_yzp_output,trim(varname_output),tmp)
        endif
      endif

     do m=1,iplane_yz
        write(unit=fnum1,fmt='(I04.4)')planes_yz(m)
        if(iAppend.eq.0) then
          fname='plane_yz_'//fnum//'_i'//fnum1//'.plt'
        endif
        idx_yzp%ist=planes_yz(m); idx_yzp%iend=planes_yz(m); idx_yzp%isp=1

        do nn=1, nvar_yzp_output
           ! write y, z
           if(varindex_yzp_output(nn).eq.1) then
             vartmp(1:ilen,1:jlen,1:klen,nn) = x(idx_yzp%ist:idx_yzp%iend:idx_yzp%isp,idx_yzp%jst:idx_yzp%jend:idx_yzp%jsp, &
                                                 idx_yzp%kst:idx_yzp%kend:idx_yzp%ksp)
           endif
           if(varindex_yzp_output(nn).eq.2) then
             if(icylin.eq.0) then
               vartmp(1:ilen,1:jlen,1:klen,nn) = y(idx_yzp%ist:idx_yzp%iend:idx_yzp%isp,idx_yzp%jst:idx_yzp%jend:idx_yzp%jsp, &
                                                   idx_yzp%kst:idx_yzp%kend:idx_yzp%ksp)
             else
               vartmp(1:ilen,1:jlen,1:klen,nn) = z(idx_yzp%ist:idx_yzp%iend:idx_yzp%isp,idx_yzp%jst:idx_yzp%jend:idx_yzp%jsp, &
                                                   idx_yzp%kst:idx_yzp%kend:idx_yzp%ksp)*cos(y(idx_yzp%ist:idx_yzp%iend:idx_yzp%isp, &
                                                   idx_yzp%jst:idx_yzp%jend:idx_yzp%jsp, idx_yzp%kst:idx_yzp%kend:idx_yzp%ksp))
             endif
           endif
           if(varindex_yzp_output(nn).eq.3) then
             if(icylin.eq.0) then
               vartmp(1:ilen,1:jlen,1:klen,nn) = z(idx_yzp%ist:idx_yzp%iend:idx_yzp%isp,idx_yzp%jst:idx_yzp%jend:idx_yzp%jsp, &
                                                   idx_yzp%kst:idx_yzp%kend:idx_yzp%ksp)
             else
               vartmp(1:ilen,1:jlen,1:klen,nn) = z(idx_yzp%ist:idx_yzp%iend:idx_yzp%isp,idx_yzp%jst:idx_yzp%jend:idx_yzp%jsp, &
                                                   idx_yzp%kst:idx_yzp%kend:idx_yzp%ksp)*sin(y(idx_yzp%ist:idx_yzp%iend:idx_yzp%isp, &
                                                   idx_yzp%jst:idx_yzp%jend:idx_yzp%jsp, idx_yzp%kst:idx_yzp%kend:idx_yzp%ksp))
             endif
           endif
          ! write u
          if(varindex_yzp_output(nn).eq.4) then
             vartmp(1:ilen,1:jlen,1:klen,nn) = buffer_flow(idx_yzp%ist:idx_yzp%iend:idx_yzp%isp,idx_yzp%jst:idx_yzp%jend:idx_yzp%jsp, &
                                                           idx_yzp%kst:idx_yzp%kend:idx_yzp%ksp,1)
          endif
          ! write v
          if(varindex_yzp_output(nn).eq.5) then
             if(icylin.eq.0) then
                vartmp(1:ilen,1:jlen,1:klen,nn) = buffer_flow(idx_yzp%ist:idx_yzp%iend:idx_yzp%isp,idx_yzp%jst:idx_yzp%jend:idx_yzp%jsp, &
                                                           idx_yzp%kst:idx_yzp%kend:idx_yzp%ksp,2)
             else
                ! convert to u_y: u_y = w*cos(theta)-v*sin(theta)
                vartmp(1:ilen,1:jlen,1:klen,nn) = buffer_flow(idx_yzp%ist:idx_yzp%iend:idx_yzp%isp,idx_yzp%jst:idx_yzp%jend:idx_yzp%jsp, &
                                                           idx_yzp%kst:idx_yzp%kend:idx_yzp%ksp,3)*cos(y(idx_yzp%ist:idx_yzp%iend:idx_yzp%isp, &
                                                           idx_yzp%jst:idx_yzp%jend:idx_yzp%jsp, idx_yzp%kst:idx_yzp%kend:idx_yzp%ksp)) - &
                                                  buffer_flow(idx_yzp%ist:idx_yzp%iend:idx_yzp%isp,idx_yzp%jst:idx_yzp%jend:idx_yzp%jsp, &
                                                           idx_yzp%kst:idx_yzp%kend:idx_yzp%ksp,2)*sin(y(idx_yzp%ist:idx_yzp%iend:idx_yzp%isp, &
                                                           idx_yzp%jst:idx_yzp%jend:idx_yzp%jsp, idx_yzp%kst:idx_yzp%kend:idx_yzp%ksp))
             endif
          endif
          ! write w
          if(varindex_yzp_output(nn).eq.6) then
             if(icylin.eq.0) then
             vartmp(1:ilen,1:jlen,1:klen,nn) = buffer_flow(idx_yzp%ist:idx_yzp%iend:idx_yzp%isp,idx_yzp%jst:idx_yzp%jend:idx_yzp%jsp, &
                                                           idx_yzp%kst:idx_yzp%kend:idx_yzp%ksp,3)
             else
                ! convert to u_z: u_z = w*sin(theta)+v*cos(theta)
                vartmp(1:ilen,1:jlen,1:klen,nn) = buffer_flow(idx_yzp%ist:idx_yzp%iend:idx_yzp%isp,idx_yzp%jst:idx_yzp%jend:idx_yzp%jsp, &
                                                           idx_yzp%kst:idx_yzp%kend:idx_yzp%ksp,3)*sin(y(idx_yzp%ist:idx_yzp%iend:idx_yzp%isp, &
                                                           idx_yzp%jst:idx_yzp%jend:idx_yzp%jsp, idx_yzp%kst:idx_yzp%kend:idx_yzp%ksp)) + &
                                                  buffer_flow(idx_yzp%ist:idx_yzp%iend:idx_yzp%isp,idx_yzp%jst:idx_yzp%jend:idx_yzp%jsp, &
                                                           idx_yzp%kst:idx_yzp%kend:idx_yzp%ksp,2)*cos(y(idx_yzp%ist:idx_yzp%iend:idx_yzp%isp, &
                                                           idx_yzp%jst:idx_yzp%jend:idx_yzp%jsp, idx_yzp%kst:idx_yzp%kend:idx_yzp%ksp))
             endif
          endif
          ! write p
          if(varindex_yzp_output(nn).eq.7) then
             vartmp(1:ilen,1:jlen,1:klen,nn) = buffer_flow(idx_yzp%ist:idx_yzp%iend:idx_yzp%isp,idx_yzp%jst:idx_yzp%jend:idx_yzp%jsp, &
                                                           idx_yzp%kst:idx_yzp%kend:idx_yzp%ksp,4)
          endif
          ! write T
          if(varindex_yzp_output(nn).eq.8) then
             vartmp(1:ilen,1:jlen,1:klen,nn) = buffer_flow(idx_yzp%ist:idx_yzp%iend:idx_yzp%isp,idx_yzp%jst:idx_yzp%jend:idx_yzp%jsp, &
                                                           idx_yzp%kst:idx_yzp%kend:idx_yzp%ksp,5)
          endif
          ! write rho
          if(varindex_yzp_output(nn).eq.9) then
             vartmp(1:ilen,1:jlen,1:klen,nn) = buffer_flow(idx_yzp%ist:idx_yzp%iend:idx_yzp%isp,idx_yzp%jst:idx_yzp%jend:idx_yzp%jsp, &
                                                           idx_yzp%kst:idx_yzp%kend:idx_yzp%ksp,6)
          endif

          ! write grad_p
          if(varindex_yzp_output(nn).eq.10) then
            if(icylin.eq.0) then
              call CalGradient(idx_yzp,x,y,z,buffer_flow(:,:,:,4),grad_p)
            else
              call CalGradient_cylin(idx_yzp,x,y,z,buffer_flow(:,:,:,4),grad_p)
              grad_p(1:ilen,1:jlen,1) = grad_p(1:ilen,1:jlen,2)
            endif
            vartmp(1:ilen,1:jlen,1:klen,nn) = grad_p(1:ilen,1:jlen,1:klen)
          endif
          ! write grad_rho
          if(varindex_yzp_output(nn).eq.11) then
            if(icylin.eq.0) then
              call CalGradient(idx_yzp,x,y,z,buffer_flow(:,:,:,6),grad_rho)
            else
              call CalGradient_cylin(idx_yzp,x,y,z,buffer_flow(:,:,:,6),grad_rho)
              grad_rho(1:ilen,1:jlen,1) = grad_rho(1:ilen,1:jlen,2)
            endif
            vartmp(1:ilen,1:jlen,1:klen,nn) = grad_rho(1:ilen,1:jlen,1:klen)
          endif
          ! write omx, omy, omz
          if(varindex_yzp_output(nn).eq.12) then
            if(icylin.eq.0) then
              call CalVorticity(idx_yzp,x,y,z,buffer_flow(:,:,:,1),buffer_flow(:,:,:,2),buffer_flow(:,:,:,3),omx,omy,omz)
            else
              call CalVorticity_cylin(idx_yzp,x,y,z,buffer_flow(:,:,:,1),buffer_flow(:,:,:,2),buffer_flow(:,:,:,3),omx,omy,omz)
              omx(1:ilen,1:jlen,1) = omx(1:ilen,1:jlen,2)
              omy(1:ilen,1:jlen,1) = omy(1:ilen,1:jlen,2)
              omz(1:ilen,1:jlen,1) = omz(1:ilen,1:jlen,2)
            endif
            vartmp(1:ilen,1:jlen,1:klen,nn)   = omx(1:ilen,1:jlen,1:klen)
            vartmp(1:ilen,1:jlen,1:klen,nn+1) = omy(1:ilen,1:jlen,1:klen)
            vartmp(1:ilen,1:jlen,1:klen,nn+2) = omz(1:ilen,1:jlen,1:klen)
          endif
          ! write div
          if(varindex_yzp_output(nn).eq.15) then
            if(icylin.eq.0) then
              call CalDivergence(idx_yzp,x,y,z,buffer_flow(:,:,:,1),buffer_flow(:,:,:,2),buffer_flow(:,:,:,3),div)
            else
              call CalDivergence_cylin(idx_yzp,x,y,z,buffer_flow(:,:,:,1),buffer_flow(:,:,:,2),buffer_flow(:,:,:,3),div)
              div(1:ilen,1:jlen,1) = div(1:ilen,1:jlen,2)
            endif
            vartmp(1:ilen,1:jlen,1:klen,nn) = div(1:ilen,1:jlen,1:klen)
          endif
          ! write swirl
          if(varindex_yzp_output(nn).eq.16) then
            call CalSwirl(idx_yzp,x,y,z,buffer_flow(:,:,:,1),buffer_flow(:,:,:,2),buffer_flow(:,:,:,3),swirl)
            vartmp(1:ilen,1:jlen,1:klen,nn) = swirl(1:ilen,1:jlen,1:klen)
          endif
          ! write Q
          if(varindex_yzp_output(nn).eq.17) then
            call Callambda2(idx_yzp,x,y,z,buffer_flow(:,:,:,1),buffer_flow(:,:,:,2),buffer_flow(:,:,:,3),Q)
            vartmp(1:ilen,1:jlen,1:klen,nn) = Q(1:ilen,1:jlen,1:klen)
          endif
          ! write src_Philip
          if(varindex_yzp_output(nn).eq.18) then
            call Calsrc_Philip(idx_yzp,x,y,z,buffer_flow(:,:,:,1),buffer_flow(:,:,:,2),buffer_flow(:,:,:,3),src_Philip)
            vartmp(1:ilen,1:jlen,1:klen,nn) = src_Philip(1:ilen,1:jlen,1:klen)
          endif

        enddo ! end nn loop

        if(icylin .eq. 1 .and. idx_yzp%jend .eq. jmax) vartmp(:,jlen+1,:,:) = vartmp(:,1,:,:)

        if(iAppend.eq.0) then
          call WriteTec(fname)
        else
           fname = 'plane_yz_'//fnum4_1//'-'//fnum4_2//'_i'//fnum1//'.plt'
           print *, 'Writing file: ', trim(fname)
           call WriteTec(fname,n)
        endif
      enddo ! end m loop
    endif ! end iplane_yz.gt.0


    if(iplane_xz.gt.0) then

      ilen = (idx_xzp%iend - idx_xzp%ist)/idx_xzp%isp + 1
      jlen = 1
      klen = (idx_xzp%kend - idx_xzp%kst)/idx_xzp%ksp + 1

      varname_output = ''
      do nn=1, nvar_xzp_output
        dname_xzp_output(nn) = dname(varindex_xzp_output(nn))
        varname_output = trim(adjustl(varname_output))//' '//trim(adjustl(dname_xzp_output(nn)))
      enddo
      print *, 'Calculate x-z plane'
      print *, 'varname_output = ', trim(varname_output)

      if(iAppend.eq.0) then
        call InitTec(1,ilen,jlen,klen,nvar_xzp_output,trim(varname_output),0)
      elseif(iAppend.eq.1.and.n.eq.1) then
        do i=1, 3
          if(varindex_xzp_output(i).le.3)then
            tmp = i
          endif
        enddo
        print *, ' number of share variable = ', tmp
        call InitTec(num_file,ilen,jlen,klen,nvar_xzp_output,trim(varname_output),tmp)
      endif

      do m=1,iplane_xz
        write(unit=fnum1,fmt='(I04.4)') planes_xz(m)
        if(iAppend.eq.0) then
          fname='plane_xz_'//fnum//'_j'//fnum1//'.plt'
        endif
        idx_xzp%jst=planes_xz(m); idx_xzp%jend=planes_xz(m); idx_xzp%jsp=1

        do nn=1, nvar_xzp_output
          if(varindex_xzp_output(nn).eq.1) then
            vartmp(1:ilen,1:jlen,1:klen,nn) = x(idx_xzp%ist:idx_xzp%iend:idx_xzp%isp,idx_xzp%jst:idx_xzp%jend:idx_xzp%jsp, &
                                                idx_xzp%kst:idx_xzp%kend:idx_xzp%ksp)
          endif
          if(varindex_xzp_output(nn).eq.2) then
            vartmp(1:ilen,1:jlen,1:klen,nn) = y(idx_xzp%ist:idx_xzp%iend:idx_xzp%isp,idx_xzp%jst:idx_xzp%jend:idx_xzp%jsp, &
                                                  idx_xzp%kst:idx_xzp%kend:idx_xzp%ksp)
          endif
          if(varindex_xzp_output(nn).eq.3) then
            vartmp(1:ilen,1:jlen,1:klen,nn) = z(idx_xzp%ist:idx_xzp%iend:idx_xzp%isp,idx_xzp%jst:idx_xzp%jend:idx_xzp%jsp, &
                                                idx_xzp%kst:idx_xzp%kend:idx_xzp%ksp)
          endif

          ! write u
          if(varindex_xzp_output(nn).eq.4) then
            vartmp(1:ilen,1:jlen,1:klen,nn) = buffer_flow(idx_xzp%ist:idx_xzp%iend:idx_xzp%isp,idx_xzp%jst:idx_xzp%jend:idx_xzp%jsp, &
                                                          idx_xzp%kst:idx_xzp%kend:idx_xzp%ksp,1)
          endif
          ! write v
          if(varindex_xzp_output(nn).eq.5) then
            vartmp(1:ilen,1:jlen,1:klen,nn) = buffer_flow(idx_xzp%ist:idx_xzp%iend:idx_xzp%isp,idx_xzp%jst:idx_xzp%jend:idx_xzp%jsp, &
                                                          idx_xzp%kst:idx_xzp%kend:idx_xzp%ksp,2)
          endif
          ! write w
          if(varindex_xzp_output(nn).eq.6) then
            vartmp(1:ilen,1:jlen,1:klen,nn) = buffer_flow(idx_xzp%ist:idx_xzp%iend:idx_xzp%isp,idx_xzp%jst:idx_xzp%jend:idx_xzp%jsp, &
                                                          idx_xzp%kst:idx_xzp%kend:idx_xzp%ksp,3)
          endif
          ! write p
          if(varindex_xzp_output(nn).eq.7) then
            vartmp(1:ilen,1:jlen,1:klen,nn) = buffer_flow(idx_xzp%ist:idx_xzp%iend:idx_xzp%isp,idx_xzp%jst:idx_xzp%jend:idx_xzp%jsp, &
                                                          idx_xzp%kst:idx_xzp%kend:idx_xzp%ksp,4)
          endif
          ! write T
          if(varindex_xzp_output(nn).eq.8) then
            vartmp(1:ilen,1:jlen,1:klen,nn) = buffer_flow(idx_xzp%ist:idx_xzp%iend:idx_xzp%isp,idx_xzp%jst:idx_xzp%jend:idx_xzp%jsp, &
                                                          idx_xzp%kst:idx_xzp%kend:idx_xzp%ksp,5)
          endif
          ! write rho
          if(varindex_xzp_output(nn).eq.9) then
            vartmp(1:ilen,1:jlen,1:klen,nn) = buffer_flow(idx_xzp%ist:idx_xzp%iend:idx_xzp%isp,idx_xzp%jst:idx_xzp%jend:idx_xzp%jsp, &
                                                          idx_xzp%kst:idx_xzp%kend:idx_xzp%ksp,6)
          endif

          ! write grad_p
          if(varindex_xzp_output(nn).eq.10) then
            if(icylin.eq.0) then
              call CalGradient(idx_xzp,x,y,z,buffer_flow(:,:,:,4),grad_p)
            else
              call CalGradient_cylin(idx_xzp,x,y,z,buffer_flow(:,:,:,4),grad_p)
              grad_p(1:ilen,1:jlen,1) = grad_p(1:ilen,1:jlen,2)
            endif
            vartmp(1:ilen,1:jlen,1:klen,nn) = grad_p(1:ilen,1:jlen,1:klen)
            deallocate(grad_p)
          endif
          ! write grad_rho
          if(varindex_xzp_output(nn).eq.11) then
            if(icylin.eq.0) then
              call CalGradient(idx_xzp,x,y,z,buffer_flow(:,:,:,6),grad_rho)
            else
              call CalGradient_cylin(idx_xzp,x,y,z,buffer_flow(:,:,:,6),grad_rho)
              grad_rho(1:ilen,1:jlen,1) = grad_rho(1:ilen,1:jlen,2)
            endif
            vartmp(1:ilen,1:jlen,1:klen,nn) = grad_rho(1:ilen,1:jlen,1:klen)
            deallocate(grad_rho)
          endif
          ! write omx, omy, omz
          if(varindex_xzp_output(nn).eq.12) then
            if(icylin.eq.0) then
              call CalVorticity(idx_xzp,x,y,z,buffer_flow(:,:,:,1),buffer_flow(:,:,:,2),buffer_flow(:,:,:,3),omx,omy,omz)
            else
              call CalVorticity_cylin(idx_xzp,x,y,z,buffer_flow(:,:,:,1),buffer_flow(:,:,:,2),buffer_flow(:,:,:,3),omx,omy,omz)
              omx(1:ilen,1:jlen,1) = omx(1:ilen,1:jlen,2)
              omy(1:ilen,1:jlen,1) = omy(1:ilen,1:jlen,2)
              omz(1:ilen,1:jlen,1) = omz(1:ilen,1:jlen,2)
            endif
            vartmp(1:ilen,1:jlen,1:klen,nn)   = omx(1:ilen,1:jlen,1:klen)
            vartmp(1:ilen,1:jlen,1:klen,nn+1) = omy(1:ilen,1:jlen,1:klen)
            vartmp(1:ilen,1:jlen,1:klen,nn+2) = omz(1:ilen,1:jlen,1:klen)
          endif
          ! write div
          if(varindex_xzp_output(nn).eq.15) then
            if(icylin.eq.0) then
              call CalDivergence(idx_xzp,x,y,z,buffer_flow(:,:,:,1),buffer_flow(:,:,:,2),buffer_flow(:,:,:,3),div)
            else
              call CalDivergence_cylin(idx_xzp,x,y,z,buffer_flow(:,:,:,1),buffer_flow(:,:,:,2),buffer_flow(:,:,:,3),div)
              div(1:ilen,1:jlen,1) = div(1:ilen,1:jlen,2)
            endif
            vartmp(1:ilen,1:jlen,1:klen,nn) = div(1:ilen,1:jlen,1:klen)
          endif
          ! write swirl
          if(varindex_xzp_output(nn).eq.16) then
            call CalSwirl(idx_xzp,x,y,z,buffer_flow(:,:,:,1),buffer_flow(:,:,:,2),buffer_flow(:,:,:,3),swirl)
             vartmp(1:ilen,1:jlen,1:klen,nn) = swirl(1:ilen,1:jlen,1:klen)
          endif
          ! write Q
          if(varindex_xzp_output(nn).eq.17) then
            call Callambda2(idx_xzp,x,y,z,buffer_flow(:,:,:,1),buffer_flow(:,:,:,2),buffer_flow(:,:,:,3),Q)
             vartmp(1:ilen,1:jlen,1:klen,nn) = Q(1:ilen,1:jlen,1:klen)
          endif

          ! write src_Philip
          if(varindex_xzp_output(nn).eq.18) then
            call Calsrc_Philip(idx_xzp,x,y,z,buffer_flow(:,:,:,1),buffer_flow(:,:,:,2),buffer_flow(:,:,:,3),src_Philip)
             vartmp(1:ilen,1:jlen,1:klen,nn) = src_Philip(1:ilen,1:jlen,1:klen)
          endif


          ! write rho*u
          if(varindex_xzp_output(nn).eq.20) then
            vartmp(1:ilen,1:jlen,1:klen,nn) = buffer_flow(idx_xzp%ist:idx_xzp%iend:idx_xzp%isp,idx_xzp%jst:idx_xzp%jend:idx_xzp%jsp, &
                                                          idx_xzp%kst:idx_xzp%kend:idx_xzp%ksp,6)*buffer_flow(idx_xzp%ist:idx_xzp%iend:idx_xzp%isp,idx_xzp%jst:idx_xzp%jend:idx_xzp%jsp, &
                                                          idx_xzp%kst:idx_xzp%kend:idx_xzp%ksp,1)
          endif

          ! write rho*v
          if(varindex_xzp_output(nn).eq.21) then
            vartmp(1:ilen,1:jlen,1:klen,nn) = buffer_flow(idx_xzp%ist:idx_xzp%iend:idx_xzp%isp,idx_xzp%jst:idx_xzp%jend:idx_xzp%jsp, &
                                                          idx_xzp%kst:idx_xzp%kend:idx_xzp%ksp,6)*buffer_flow(idx_xzp%ist:idx_xzp%iend:idx_xzp%isp,idx_xzp%jst:idx_xzp%jend:idx_xzp%jsp, &
                                                          idx_xzp%kst:idx_xzp%kend:idx_xzp%ksp,2)
          endif

          ! write rho*w
          if(varindex_xzp_output(nn).eq.22) then
            vartmp(1:ilen,1:jlen,1:klen,nn) = buffer_flow(idx_xzp%ist:idx_xzp%iend:idx_xzp%isp,idx_xzp%jst:idx_xzp%jend:idx_xzp%jsp, &
                                                          idx_xzp%kst:idx_xzp%kend:idx_xzp%ksp,6)*buffer_flow(idx_xzp%ist:idx_xzp%iend:idx_xzp%isp,idx_xzp%jst:idx_xzp%jend:idx_xzp%jsp, &
                                                          idx_xzp%kst:idx_xzp%kend:idx_xzp%ksp,3)
          endif

         enddo ! end nn loop
         if(iAppend.eq.0) then
           call WriteTec(fname)
         else
           fname = 'plane_xz_'//fnum4_1//'-'//fnum4_2//'_j'//fnum1//'.plt'
           print *, 'Writing file: ', trim(fname)
           call WriteTec(fname,n)
         endif
      enddo ! end m loop
    endif ! end iplane_xz.gt.0

  end do  ! end n loop

  contains
    subroutine Input()
      integer :: i
      read(*,*)
      read(*,'(A)')datapath  !path where data files are stored
      read(*,*)
      read(*,*) Rm
      read(*,*)
      read(*,*) file_be, file_end, file_skip, iHDF5, iAppend, icylin, iBaseflowSubtract, iSpanave
      read(*,*)
      read(*,*)
      read(*,*)ivolume, idx_vol%ist, idx_vol%iend, idx_vol%isp, idx_vol%jst, idx_vol%jend, idx_vol%jsp, idx_vol%kst, idx_vol%kend, idx_vol%ksp
      read(*,*)
      read(*,*) nvar_vol_output
      allocate(varindex_vol_output(nvar_vol_output), dname_vol_output(nvar_vol_output))
      read(*,*) (varindex_vol_output(i), i=1,nvar_vol_output)
!     xy planes
      read(*,*)
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

      if(iAppend.gt.0) then
        if((ivolume+iplane_xy+iplane_yz+iplane_xz).gt.1) then
          print *, 'Only one of ivolume, iplane_xy, iplane_yz and iplane_xz can be 1. '
          stop
        endif
      endif

      if(iHDF5.eq.0) then
        print *, 'Reading file: ', trim(datapath)//'gridp3d.grd'
        call ReadPlot3DGridPlaneGen(trim(datapath)//'gridp3d.grd',imax,jmax,kmax,x,y,z,0)
      else
        grd%fname = trim(datapath)//'grid.h5'
        grd%gname = '/'
        print *, 'Reading file: ', trim(grd%fname)
        call DetectHDF5(grd)
        allocate(buffer_grid(grd%dimsf(1), grd%dimsf(2),grd%dimsf(3),3))
        grd%dname(1) = 'x'
        grd%dname(2) = 'y'
        grd%dname(3) = 'z'
        call ReadHDF5_3D(grd,buffer_grid)
        imax = grd%dimsf(2)
        jmax = grd%dimsf(3)
        kmax = grd%dimsf(1)
        allocate(x(imax,jmax,kmax),y(imax,jmax,kmax),z(imax,jmax,kmax))
        !if(icylin .ne. 1) then
          do k=1, kmax
            do j=1, jmax
              do i=1, imax
                x(i,j,k) = buffer_grid(k,i,j,1)
                y(i,j,k) = buffer_grid(k,i,j,2)
                z(i,j,k) = buffer_grid(k,i,j,3)
              enddo
            enddo
          enddo
        !else
        !  do k=1, kmax
        !    do j=1, jmax
        !      do i=1, imax
        !        x(i,j,k) = buffer_grid(k,i,j,1)
        !        y(i,j,k) = buffer_grid(k,i,j,3)*cos(buffer_grid(k,i,j,2))
        !        z(i,j,k) = buffer_grid(k,i,j,3)*sin(buffer_grid(k,i,j,2))
        !      enddo
        !    enddo
        !  enddo
        !endif

        deallocate(buffer_grid)
      endif
      allocate(buffer_flow(imax,jmax,kmax,6))
      print *, 'imax = ', imax, 'jmax = ', jmax, 'kmax = ', kmax

      rbar = R/Rm
      return
    end subroutine Input


end program volume    
