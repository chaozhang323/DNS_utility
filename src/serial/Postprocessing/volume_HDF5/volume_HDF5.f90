program volume
  use MFileIO
  use modVolume
  use modRWHDF5
  use modTecbin
  implicit none

  real(8), parameter :: R=8314.3D0
  real(8) :: rbar, Rm
  character(400) :: fname, datapath, datapath_RM, gridpath_RM
  character(8) :: fnum 
  character(4) :: fnum1, fnum4_1, fnum4_2, fibe, fiend, fjbe, fjend, fkbe, fkend
  integer :: n, i, j, k, m, nn, k2
  integer :: imax, jmax, kmax, file_be, file_end, file_skip, num_file, icylin, isym
  type(index_output) :: idx_vol, idx_vol2
  logical :: IsCylPeriod=.false., IsAxiSym=.false.

  real, dimension(:,:,:), allocatable :: x, y, z, x2, y2, z2
  real, dimension(:,:,:,:), allocatable :: buffer_grid, buffer_grid2, buffer_RM, buffer_RM_rd
  real, dimension(:,:,:,:), allocatable :: buffer_flow, buffer_flow_tmp, buffer_flow_tmp2, buffer_flow2
  real, dimension(:,:,:), allocatable :: grad_p, grad_rho, omx, omy, omz
  real, dimension(:,:,:), allocatable :: swirl, div, Q, src_Philip, p0, lap, shockloc
  real, dimension(:,:,:), allocatable :: dudn, dtdn, ds

  integer, dimension(:), allocatable :: varindex_vol_output
  character(10), dimension(:), allocatable :: dname_vol_output
  integer :: nvar, nvar_vol_output
  parameter(nvar=21)
  character(10) :: dname(nvar)
  parameter(dname = (/'x','y','z','u','v','w','p','t','rho','grad_p','grad_rho','omx','omy','omz','div','swil','Q', 'src_Philip', 'p0', 'Lap_rho', 'shockloc'/))
  character(400) :: varname_output, varname_dat
  integer, dimension(:), allocatable :: dims
  integer :: rank, ilen, jlen, klen, ilen_wt, jlen_wt, klen_wt
  type(tp_rdwt_hdf5) :: fsol, grd, fsol2, grd2, fRM
  integer :: iAppend, iBaseflowSubtract
  integer :: tmp, cal_2D
  integer :: iFlowData, iRescaleMean, istencilsize, kdim_RM
  real(8) :: uinf, rhoinf

  call InitHDF5()
  call Input()

  num_file = (file_end-file_be)/file_skip + 1
  write(unit=fnum4_1,fmt='(I04.4)') file_be/1000
  write(unit=fnum4_2,fmt='(I04.4)') file_end/1000

if(iFlowData.gt.0) then

  if(num_file.eq.1.and.iAppend.eq.1) then
    print *, 'num_file cannot be 1 if iAppend=1 ... Stop! '
    print *, 'num_file = ', num_file, 'iAppend = ', iAppend
    stop
  endif

  if(iBaseflowSubtract.eq.1) then
    allocate(buffer_flow_tmp(klen,ilen,jlen,5))
    fsol%fname = trim(datapath)//'baseflow.h5'
    print *, 'Reading file: ', trim(fsol%fname)
    call ReadHDF5_3D_test(fsol,buffer_flow_tmp)
  endif
  
  do n=1, num_file
     write(unit=fnum,fmt='(I08.8)') file_be + (n-1)*file_skip

       allocate(buffer_flow_tmp2(klen,ilen,jlen,5))
       fsol%fname = trim(datapath)//'flowdata_'//fnum//'.h5'
       print *, 'Reading file: ', trim(fsol%fname)
       call ReadHDF5_3D_test(fsol,buffer_flow_tmp2)

       do i=1, ilen
         do j=1, jlen
           do k=1, klen
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
          print *, 'Subtracting baseflow' 
          do i=1, ilen
             do j=1, jlen
                do k=1, klen
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

       if(IsAxiSym) then
           allocate(buffer_flow_tmp2(klen,ilen,jlen,5))
           fsol2%fname = trim(datapath)//'flowdata_'//fnum//'.h5'
!           print *, 'Reading file: ', trim(fsol2%fname)
           call ReadHDF5_3D_test(fsol2,buffer_flow_tmp2)

           do i=1, ilen
             do j=1, jlen
               do k=1, klen
                 buffer_flow2(i,j,k,1) = buffer_flow_tmp2(k,i,j,1)
                 buffer_flow2(i,j,k,2) = buffer_flow_tmp2(k,i,j,2)
                 buffer_flow2(i,j,k,3) = buffer_flow_tmp2(k,i,j,3)
                 buffer_flow2(i,j,k,4) = buffer_flow_tmp2(k,i,j,4)
                 buffer_flow2(i,j,k,5) = buffer_flow_tmp2(k,i,j,5)
               enddo
             enddo
           enddo
           buffer_flow2(:,:,:,6) = buffer_flow2(:,:,:,4)/(rbar*buffer_flow2(:,:,:,5))
           deallocate(buffer_flow_tmp2)

           if(iBaseflowSubtract.eq.1) then
              print *, 'Subtracting baseflow'
              do i=1, ilen
                 do j=1, jlen
                    do k=1, klen
                       buffer_flow2(i,j,k,1) = buffer_flow2(i,j,k,1) - buffer_flow_tmp(k,i,j,1)
                       buffer_flow2(i,j,k,2) = buffer_flow2(i,j,k,2) - buffer_flow_tmp(k,i,j,2)
                       buffer_flow2(i,j,k,3) = buffer_flow2(i,j,k,3) - buffer_flow_tmp(k,i,j,3)
                       buffer_flow2(i,j,k,4) = buffer_flow2(i,j,k,4) - buffer_flow_tmp(k,i,j,4)
                       buffer_flow2(i,j,k,5) = buffer_flow2(i,j,k,5) - buffer_flow_tmp(k,i,j,5)
                       buffer_flow2(i,j,k,6) = buffer_flow2(i,j,k,6) &
                                            - buffer_flow_tmp(k,i,j,4)/(rbar*buffer_flow_tmp(k,i,j,5))
                    enddo
                 enddo
               enddo
           endif
       endif    ! endif IsAxiSym

      varname_output = ''
      do nn=1, nvar_vol_output
        dname_vol_output(nn) = dname(varindex_vol_output(nn))
        varname_output = trim(adjustl(varname_output))//' '//trim(adjustl(dname_vol_output(nn)))
      enddo
      print *, 'Calculate volume'
      print *, 'varname_output = ', trim(varname_output)

      if(iAppend.eq.0) then
        if(IsCylPeriod) then
          call InitTec(1,ilen_wt,jlen_wt+1,klen_wt,nvar_vol_output,trim(varname_output),0)
        else
          if(IsAxiSym) then ! axisymmetric plane
            call InitTec(1,ilen_wt,jlen_wt,2*klen_wt,nvar_vol_output,trim(varname_output),0)
          else
            call InitTec(1,ilen_wt,jlen_wt,klen_wt,nvar_vol_output,trim(varname_output),0)
          endif
        endif
      elseif(iAppend.eq.1.and.n.eq.1) then
        do i=1, 3
          if(varindex_vol_output(i).le.3)then
            tmp = i
          endif
        enddo
        print *, ' number of share variable = ', tmp
        if(IsCylPeriod) then
          call InitTec(num_file,ilen_wt,jlen_wt+1,klen_wt,nvar_vol_output,trim(varname_output),tmp)
        else
          if(IsAxiSym) then ! axisymmetric plane
            call InitTec(num_file,ilen_wt,jlen_wt,2*klen_wt,nvar_vol_output,trim(varname_output),tmp)
          else
            call InitTec(num_file,ilen_wt,jlen_wt,klen_wt,nvar_vol_output,trim(varname_output),tmp)
          endif
        endif
      endif

      !fname='volume_'//fnum//'.plt'
      fname = 'volume_i'//fibe//'-'//fiend//'_j'//fjbe//'-'//fjend//'_k'//fkbe//'-'//fkend//'_'//fnum//'.plt'

      do nn=1, nvar_vol_output
        ! write x, y, z
        if(varindex_vol_output(nn).eq.1) then
          if(IsAxiSym) then
            vartmp(1:ilen_wt,1:jlen_wt,klen_wt+1:2*klen_wt,nn) = x(idx_vol%ist:idx_vol%iend:idx_vol%isp,idx_vol%jst:idx_vol%jend:idx_vol%jsp, &
                                              idx_vol%kst:idx_vol%kend:idx_vol%ksp)
            do k = 1,klen_wt
                k2 = idx_vol2%kend-k+1
                vartmp(1:ilen_wt,1:jlen_wt,k,nn) = x2(idx_vol2%ist:idx_vol2%iend:idx_vol2%isp,idx_vol2%jst:idx_vol2%jend:idx_vol2%jsp,k2)
            enddo
          else
            vartmp(1:ilen_wt,1:jlen_wt,1:klen_wt,nn) = x(idx_vol%ist:idx_vol%iend:idx_vol%isp,idx_vol%jst:idx_vol%jend:idx_vol%jsp, &
                                              idx_vol%kst:idx_vol%kend:idx_vol%ksp)
          endif
        endif
        if(varindex_vol_output(nn).eq.2) then
!          if(icylin.eq.0) then
            vartmp(1:ilen_wt,1:jlen_wt,1:klen_wt,nn) = y(idx_vol%ist:idx_vol%iend:idx_vol%isp,idx_vol%jst:idx_vol%jend:idx_vol%jsp, &
                                                idx_vol%kst:idx_vol%kend:idx_vol%ksp)
!          else
!            vartmp(1:ilen_wt,1:jlen_wt,1:klen_wt,nn) = z(idx_vol%ist:idx_vol%iend:idx_vol%isp, idx_vol%jst:idx_vol%jend:idx_vol%jsp, &
!                                                idx_vol%kst:idx_vol%kend:idx_vol%ksp)*cos(y(idx_vol%ist:idx_vol%iend:idx_vol%isp, &
!                                                idx_vol%jst:idx_vol%jend:idx_vol%jsp, idx_vol%kst:idx_vol%kend:idx_vol%ksp))
!          endif
        endif
        if(varindex_vol_output(nn).eq.3) then
!          if(icylin.eq.0) then
          if(IsAxiSym) then
            vartmp(1:ilen_wt,1:jlen_wt,klen_wt+1:2*klen_wt,nn) = z(idx_vol%ist:idx_vol%iend:idx_vol%isp,idx_vol%jst:idx_vol%jend:idx_vol%jsp, &
                                              idx_vol%kst:idx_vol%kend:idx_vol%ksp)
            do k = 1,klen_wt
                k2 = idx_vol2%kend-k+1
                vartmp(1:ilen_wt,1:jlen_wt,k,nn) = -z2(idx_vol2%ist:idx_vol2%iend:idx_vol2%isp,idx_vol2%jst:idx_vol2%jend:idx_vol2%jsp,k2)
            enddo
          else
            vartmp(1:ilen_wt,1:jlen_wt,1:klen_wt,nn) = z(idx_vol%ist:idx_vol%iend:idx_vol%isp,idx_vol%jst:idx_vol%jend:idx_vol%jsp, &
                                                idx_vol%kst:idx_vol%kend:idx_vol%ksp)
!          else
!            vartmp(1:ilen_wt,1:jlen_wt,1:klen_wt,nn) = z(idx_vol%ist:idx_vol%iend:idx_vol%isp,idx_vol%jst:idx_vol%jend:idx_vol%jsp, &
!                                                idx_vol%kst:idx_vol%kend:idx_vol%ksp)*sin(y(idx_vol%ist:idx_vol%iend:idx_vol%isp, &
!                                                idx_vol%jst:idx_vol%jend:idx_vol%jsp, idx_vol%kst:idx_vol%kend:idx_vol%ksp))
          endif
        endif

         ! write u
        if(varindex_vol_output(nn).eq.4) then
          if(IsAxiSym) then
            vartmp(1:ilen_wt,1:jlen_wt,klen_wt+1:2*klen_wt,nn) = buffer_flow(idx_vol%ist:idx_vol%iend:idx_vol%isp,&
                                              idx_vol%jst:idx_vol%jend:idx_vol%jsp, idx_vol%kst:idx_vol%kend:idx_vol%ksp,1)
            do k = 1,klen_wt
                k2 = idx_vol2%kend-k+1
                vartmp(1:ilen_wt,1:jlen_wt,k,nn) = buffer_flow2(idx_vol2%ist:idx_vol2%iend:idx_vol2%isp,&
                                                                idx_vol2%jst:idx_vol2%jend:idx_vol2%jsp,k2,1)
            enddo
          else
            vartmp(1:ilen_wt,1:jlen_wt,1:klen_wt,nn) = buffer_flow(idx_vol%ist:idx_vol%iend:idx_vol%isp,idx_vol%jst:idx_vol%jend:idx_vol%jsp, &
                                                        idx_vol%kst:idx_vol%kend:idx_vol%ksp,1)
          endif     ! endif IsAxiSym
        endif
        ! write v
        if(varindex_vol_output(nn).eq.5) then
          if(IsAxiSym) then
            vartmp(1:ilen_wt,1:jlen_wt,klen_wt+1:2*klen_wt,nn) = buffer_flow(idx_vol%ist:idx_vol%iend:idx_vol%isp,&
                                              idx_vol%jst:idx_vol%jend:idx_vol%jsp, idx_vol%kst:idx_vol%kend:idx_vol%ksp,2)
            do k = 1,klen_wt
                k2 = idx_vol2%kend-k+1
                vartmp(1:ilen_wt,1:jlen_wt,k,nn) = buffer_flow2(idx_vol2%ist:idx_vol2%iend:idx_vol2%isp,&
                                                                idx_vol2%jst:idx_vol2%jend:idx_vol2%jsp,k2,2)
            enddo
          else
            vartmp(1:ilen_wt,1:jlen_wt,1:klen_wt,nn) = buffer_flow(idx_vol%ist:idx_vol%iend:idx_vol%isp,idx_vol%jst:idx_vol%jend:idx_vol%jsp, &
                                                        idx_vol%kst:idx_vol%kend:idx_vol%ksp,2)
          endif     ! endif IsAxiSym
        endif
        ! write w
        if(varindex_vol_output(nn).eq.6) then
          if(IsAxiSym) then
            vartmp(1:ilen_wt,1:jlen_wt,klen_wt+1:2*klen_wt,nn) = buffer_flow(idx_vol%ist:idx_vol%iend:idx_vol%isp,&
                                              idx_vol%jst:idx_vol%jend:idx_vol%jsp, idx_vol%kst:idx_vol%kend:idx_vol%ksp,3)
            do k = 1,klen_wt
                k2 = idx_vol2%kend-k+1
                vartmp(1:ilen_wt,1:jlen_wt,k,nn) = -buffer_flow2(idx_vol2%ist:idx_vol2%iend:idx_vol2%isp,&
                                                                idx_vol2%jst:idx_vol2%jend:idx_vol2%jsp,k2,3)
            enddo
          else
            vartmp(1:ilen_wt,1:jlen_wt,1:klen_wt,nn) = buffer_flow(idx_vol%ist:idx_vol%iend:idx_vol%isp,idx_vol%jst:idx_vol%jend:idx_vol%jsp, &
                                                        idx_vol%kst:idx_vol%kend:idx_vol%ksp,3)
          endif     ! endif IsAxiSym
        endif
        ! write p
        if(varindex_vol_output(nn).eq.7) then
          if(IsAxiSym) then
            vartmp(1:ilen_wt,1:jlen_wt,klen_wt+1:2*klen_wt,nn) = buffer_flow(idx_vol%ist:idx_vol%iend:idx_vol%isp,&
                                              idx_vol%jst:idx_vol%jend:idx_vol%jsp, idx_vol%kst:idx_vol%kend:idx_vol%ksp,4)
            do k = 1,klen_wt
                k2 = idx_vol2%kend-k+1
                vartmp(1:ilen_wt,1:jlen_wt,k,nn) = buffer_flow2(idx_vol2%ist:idx_vol2%iend:idx_vol2%isp,&
                                                                idx_vol2%jst:idx_vol2%jend:idx_vol2%jsp,k2,4)
            enddo
          else
            vartmp(1:ilen_wt,1:jlen_wt,1:klen_wt,nn) = buffer_flow(idx_vol%ist:idx_vol%iend:idx_vol%isp,idx_vol%jst:idx_vol%jend:idx_vol%jsp, &
                                                        idx_vol%kst:idx_vol%kend:idx_vol%ksp,4)
          endif     ! endif IsAxiSym
        endif
        ! write t
        if(varindex_vol_output(nn).eq.8) then
          if(IsAxiSym) then
            vartmp(1:ilen_wt,1:jlen_wt,klen_wt+1:2*klen_wt,nn) = buffer_flow(idx_vol%ist:idx_vol%iend:idx_vol%isp,&
                                              idx_vol%jst:idx_vol%jend:idx_vol%jsp, idx_vol%kst:idx_vol%kend:idx_vol%ksp,5)
            do k = 1,klen_wt
                k2 = idx_vol2%kend-k+1
                vartmp(1:ilen_wt,1:jlen_wt,k,nn) = buffer_flow2(idx_vol2%ist:idx_vol2%iend:idx_vol2%isp,&
                                                                idx_vol2%jst:idx_vol2%jend:idx_vol2%jsp,k2,5)
            enddo
          else
            vartmp(1:ilen_wt,1:jlen_wt,1:klen_wt,nn) = buffer_flow(idx_vol%ist:idx_vol%iend:idx_vol%isp,idx_vol%jst:idx_vol%jend:idx_vol%jsp, &
                                                        idx_vol%kst:idx_vol%kend:idx_vol%ksp,5)
          endif     ! endif IsAxiSym
        endif
        ! write rho
        if(varindex_vol_output(nn).eq.9) then
          if(IsAxiSym) then
            vartmp(1:ilen_wt,1:jlen_wt,klen_wt+1:2*klen_wt,nn) = buffer_flow(idx_vol%ist:idx_vol%iend:idx_vol%isp,&
                                              idx_vol%jst:idx_vol%jend:idx_vol%jsp, idx_vol%kst:idx_vol%kend:idx_vol%ksp,6)
            do k = 1,klen_wt
                k2 = idx_vol2%kend-k+1
                vartmp(1:ilen_wt,1:jlen_wt,k,nn) = buffer_flow2(idx_vol2%ist:idx_vol2%iend:idx_vol2%isp,&
                                                                idx_vol2%jst:idx_vol2%jend:idx_vol2%jsp,k2,6)
            enddo
          else
            vartmp(1:ilen_wt,1:jlen_wt,1:klen_wt,nn) = buffer_flow(idx_vol%ist:idx_vol%iend:idx_vol%isp,idx_vol%jst:idx_vol%jend:idx_vol%jsp, &
                                                        idx_vol%kst:idx_vol%kend:idx_vol%ksp,6)
          endif     ! endif IsAxiSym
        endif
        ! write grad_p
        if(varindex_vol_output(nn).eq.10) then
          if(IsAxiSym) then
              call CalGradient_cylin(idx_vol,x,y,z,buffer_flow(:,:,:,4),grad_p)
              grad_p(1:ilen_wt,1:jlen_wt,1) = grad_p(1:ilen_wt,1:jlen_wt,2)
              vartmp(1:ilen_wt,1:jlen_wt,klen_wt+1:2*klen_wt,nn) = grad_p(1:ilen_wt,1:jlen_wt,1:klen_wt)
              call CalGradient_cylin(idx_vol2,x2,y2,z2,buffer_flow2(:,:,:,4),grad_p)
              grad_p(1:ilen_wt,1:jlen_wt,1) = grad_p(1:ilen_wt,1:jlen_wt,2)
              do k = 1,klen_wt
                k2 = idx_vol2%kend-k+1
                vartmp(1:ilen_wt,1:jlen_wt,k,nn) = grad_p(1:ilen_wt,1:jlen_wt,k2)
              enddo
              deallocate(grad_p)
          else
              if(icylin.eq.0) then
                call CalGradient(idx_vol,x,y,z,buffer_flow(:,:,:,4),grad_p)
              else
                call CalGradient_cylin(idx_vol,x,y,z,buffer_flow(:,:,:,4),grad_p)
                grad_p(1:ilen_wt,1:jlen_wt,1) = grad_p(1:ilen_wt,1:jlen_wt,2)
              endif
              vartmp(1:ilen_wt,1:jlen_wt,1:klen_wt,nn) = grad_p(1:ilen_wt,1:jlen_wt,1:klen_wt)
              deallocate(grad_p)
          endif     ! endif IsAxiSym
        endif
        ! write grad_rho
        if(varindex_vol_output(nn).eq.11) then
          if(IsAxiSym) then
              call CalGradient_cylin(idx_vol,x,y,z,buffer_flow(:,:,:,6),grad_rho)
              grad_rho(1:ilen_wt,1:jlen_wt,1) = grad_rho(1:ilen_wt,1:jlen_wt,2)
              vartmp(1:ilen_wt,1:jlen_wt,klen_wt+1:2*klen_wt,nn) = grad_rho(1:ilen_wt,1:jlen_wt,1:klen_wt)
              call CalGradient_cylin(idx_vol2,x2,y2,z2,buffer_flow2(:,:,:,6),grad_rho)
              grad_rho(1:ilen_wt,1:jlen_wt,1) = grad_rho(1:ilen_wt,1:jlen_wt,2)
              do k = 1,klen_wt
                k2 = idx_vol2%kend-k+1
                vartmp(1:ilen_wt,1:jlen_wt,k,nn) = grad_rho(1:ilen_wt,1:jlen_wt,k2)
              enddo
              deallocate(grad_rho)
          else
              if(icylin.eq.0) then
                call CalGradient(idx_vol,x,y,z,buffer_flow(:,:,:,6),grad_rho)
              else
                call CalGradient_cylin(idx_vol,x,y,z,buffer_flow(:,:,:,6),grad_rho)
                grad_rho(1:ilen_wt,1:jlen_wt,1) = grad_rho(1:ilen_wt,1:jlen_wt,2)
              endif
              vartmp(1:ilen_wt,1:jlen_wt,1:klen_wt,nn) = grad_rho(1:ilen_wt,1:jlen_wt,1:klen_wt)
              deallocate(grad_rho)
          endif     ! endif IsAxiSym
        endif
        ! write omx, omy, omz
        if(varindex_vol_output(nn).eq.12.or.varindex_vol_output(nn).eq.13.or.varindex_vol_output(nn).eq.14) then
          if(IsAxiSym) then
              call CalVorticity_cylin(idx_vol,x,y,z,buffer_flow(:,:,:,1),buffer_flow(:,:,:,2),buffer_flow(:,:,:,3),omx,omy,omz)
              omx(1:ilen_wt,1:jlen_wt,1) = omx(1:ilen_wt,1:jlen_wt,2)
              omy(1:ilen_wt,1:jlen_wt,1) = omy(1:ilen_wt,1:jlen_wt,2)
              omz(1:ilen_wt,1:jlen_wt,1) = omz(1:ilen_wt,1:jlen_wt,2)
              if(varindex_vol_output(nn).eq.12) vartmp(1:ilen_wt,1:jlen_wt,klen_wt+1:2*klen_wt,nn) = omx(1:ilen_wt,1:jlen_wt,1:klen_wt)
              if(varindex_vol_output(nn).eq.13) vartmp(1:ilen_wt,1:jlen_wt,klen_wt+1:2*klen_wt,nn) = omy(1:ilen_wt,1:jlen_wt,1:klen_wt)
              if(varindex_vol_output(nn).eq.14) vartmp(1:ilen_wt,1:jlen_wt,klen_wt+1:2*klen_wt,nn) = omz(1:ilen_wt,1:jlen_wt,1:klen_wt)
              call CalVorticity_cylin(idx_vol2,x2,y2,z2,buffer_flow2(:,:,:,1),buffer_flow2(:,:,:,2),buffer_flow2(:,:,:,3),omx,omy,omz)
              omx(1:ilen_wt,1:jlen_wt,1) = omx(1:ilen_wt,1:jlen_wt,2)
              omy(1:ilen_wt,1:jlen_wt,1) = omy(1:ilen_wt,1:jlen_wt,2)
              omz(1:ilen_wt,1:jlen_wt,1) = omz(1:ilen_wt,1:jlen_wt,2)
              if(varindex_vol_output(nn).eq.12) then
                  do k = 1,klen_wt
                    k2 = idx_vol2%kend-k+1
                    vartmp(1:ilen_wt,1:jlen_wt,k,nn) = omx(1:ilen_wt,1:jlen_wt,k2)
                  enddo
              endif
              if(varindex_vol_output(nn).eq.13) then
                  do k = 1,klen_wt
                    k2 = idx_vol2%kend-k+1
                    vartmp(1:ilen_wt,1:jlen_wt,k,nn) = omy(1:ilen_wt,1:jlen_wt,k2)
                  enddo
              endif
              if(varindex_vol_output(nn).eq.14) then
                  do k = 1,klen_wt
                    k2 = idx_vol2%kend-k+1
                    vartmp(1:ilen_wt,1:jlen_wt,k,nn) = -omz(1:ilen_wt,1:jlen_wt,k2)
                  enddo
              endif
          else
              if(icylin.eq.0) then
                call CalVorticity(idx_vol,x,y,z,buffer_flow(:,:,:,1),buffer_flow(:,:,:,2),buffer_flow(:,:,:,3),omx,omy,omz)
              else
                call CalVorticity_cylin(idx_vol,x,y,z,buffer_flow(:,:,:,1),buffer_flow(:,:,:,2),buffer_flow(:,:,:,3),omx,omy,omz)
                omx(1:ilen_wt,1:jlen_wt,1) = omx(1:ilen_wt,1:jlen_wt,2)
                omy(1:ilen_wt,1:jlen_wt,1) = omy(1:ilen_wt,1:jlen_wt,2)
                omz(1:ilen_wt,1:jlen_wt,1) = omz(1:ilen_wt,1:jlen_wt,2)
              endif
              if(varindex_vol_output(nn).eq.12) vartmp(1:ilen_wt,1:jlen_wt,1:klen_wt,nn)   = omx(1:ilen_wt,1:jlen_wt,1:klen_wt)
              if(varindex_vol_output(nn).eq.13) vartmp(1:ilen_wt,1:jlen_wt,1:klen_wt,nn) = omy(1:ilen_wt,1:jlen_wt,1:klen_wt)
              if(varindex_vol_output(nn).eq.14) vartmp(1:ilen_wt,1:jlen_wt,1:klen_wt,nn) = omz(1:ilen_wt,1:jlen_wt,1:klen_wt)
          endif     ! endif IsAxiSym
          deallocate(omx,omy,omz)
        endif
        ! write div
        if(varindex_vol_output(nn).eq.15) then
          if(IsAxiSym) then
              call CalDivergence_cylin(idx_vol,x,y,z,buffer_flow(:,:,:,1),buffer_flow(:,:,:,2),buffer_flow(:,:,:,3),div)
              div(1:ilen_wt,1:jlen_wt,1) = div(1:ilen_wt,1:jlen_wt,2)
              vartmp(1:ilen_wt,1:jlen_wt,klen_wt+1:2*klen_wt,nn) = div(1:ilen_wt,1:jlen_wt,1:klen_wt)
              call CalDivergence_cylin(idx_vol2,x2,y2,z2,buffer_flow2(:,:,:,1),buffer_flow2(:,:,:,2),buffer_flow2(:,:,:,3),div)
              div(1:ilen_wt,1:jlen_wt,1) = div(1:ilen_wt,1:jlen_wt,2)
              do k = 1,klen_wt
                k2 = idx_vol2%kend-k+1
                vartmp(1:ilen_wt,1:jlen_wt,k,nn) = div(1:ilen_wt,1:jlen_wt,k2)
              enddo
          else
              if(icylin.eq.0) then
                call CalDivergence(idx_vol,x,y,z,buffer_flow(:,:,:,1),buffer_flow(:,:,:,2),buffer_flow(:,:,:,3),div)
              else
                call CalDivergence_cylin(idx_vol,x,y,z,buffer_flow(:,:,:,1),buffer_flow(:,:,:,2),buffer_flow(:,:,:,3),div)
                div(1:ilen_wt,1:jlen_wt,1) = div(1:ilen_wt,1:jlen_wt,2)
              endif
              vartmp(1:ilen_wt,1:jlen_wt,1:klen_wt,nn) = div(1:ilen_wt,1:jlen_wt,1:klen_wt)
          endif     ! endif IsAxiSym
          deallocate(div)
        endif
        ! write swirl
        if(varindex_vol_output(nn).eq.16) then
          call CalSwirl(idx_vol,x,y,z,buffer_flow(:,:,:,1),buffer_flow(:,:,:,2),buffer_flow(:,:,:,3),swirl)
          vartmp(1:ilen_wt,1:jlen_wt,1:klen_wt,nn) = swirl(1:ilen_wt,1:jlen_wt,1:klen_wt)
        endif
        ! write Q
        if(varindex_vol_output(nn).eq.17) then
          call Callambda2(idx_vol,x,y,z,buffer_flow(:,:,:,1),buffer_flow(:,:,:,2),buffer_flow(:,:,:,3),Q)
          vartmp(1:ilen_wt,1:jlen_wt,1:klen_wt,nn) = Q(1:ilen_wt,1:jlen_wt,1:klen_wt)
        endif

        ! write src_Philip
        if(varindex_vol_output(nn).eq.18) then
          call Calsrc_Philip(idx_vol,x,y,z,buffer_flow(:,:,:,1),buffer_flow(:,:,:,2),buffer_flow(:,:,:,3),src_Philip)
          vartmp(1:ilen_wt,1:jlen_wt,1:klen_wt,nn) = src_Philip(1:ilen_wt,1:jlen_wt,1:klen_wt)
        endif

        ! write p0
        if(varindex_vol_output(nn).eq.19) then
          call Calp0(idx_vol,buffer_flow(:,:,:,1),buffer_flow(:,:,:,2),buffer_flow(:,:,:,3),buffer_flow(:,:,:,4),buffer_flow(:,:,:,5),p0)
          vartmp(1:ilen_wt,1:jlen_wt,1:klen_wt,nn) = p0(1:ilen_wt,1:jlen_wt,1:klen_wt)
        endif

        ! write laplacian of density
        if(varindex_vol_output(nn).eq.20) then
          call CalLaplacianGCC(idx_vol,x,y,z,buffer_flow(:,:,:,6),buffer_flow(:,:,:,1),buffer_flow(:,:,:,2),buffer_flow(:,:,:,3),lap,shockloc)
          vartmp(1:ilen_wt,1:jlen_wt,1:klen_wt,nn)   = lap(1:ilen_wt,1:jlen_wt,1:klen_wt)
          vartmp(1:ilen_wt,1:jlen_wt,1:klen_wt,nn+1) = shockloc(1:ilen_wt,1:jlen_wt,1:klen_wt)

          if(n.eq.1) then
            open(22,file='ShockLoc_volume.dat',status='unknown')
            write(22,'(a)') 'variables=z,np_shock,np_shock_percentage'
            close(22)
          endif

          open(22,file='ShockLoc_volume.dat',position='append',status='old')
            write(22,'("Zone T=ndv_", I08.8, ", I=", I4)') file_be + (n-1)*file_skip, klen_wt
            do k=1, klen_wt
              write(22,*) z(1,1,k), sum(shockloc(:,:,k)), sum(shockloc(:,:,k))/dble(ilen_wt*jlen_wt)
            enddo
          close(16)

        endif

      enddo ! end nn loop

      if(IsCylPeriod) vartmp(:,jlen_wt+1,:,:) = vartmp(:,1,:,:)

      if(iAppend.eq.0) then
        call WriteTec(fname)
      else
        !fname = 'volume_'//fnum4_1//'-'//fnum4_2//'.plt'
        fname = 'volume_i'//fibe//'-'//fiend//'_j'//fjbe//'-'//fjend//'_k'//fkbe//'-'//fkend//'-'//fnum4_1//'-'//fnum4_2//'.plt'
        if(n.eq.num_file) print *, 'Writing file: ', trim(fname)
        call WriteTec(fname,n)
      endif

  enddo  ! end n loop

  endif ! iFlowData.gt.0

  if(iRescaleMean.gt.0) then
    varname_output = 'z uave wave rhoave tave pave'
    call InitTec(num_file,istencilsize,1,kdim_RM,6,trim(varname_output),0)
    vartmp(:,1,:,1) = buffer_RM(:,1,:,1)
    do n=1, num_file
      write(unit=fnum,fmt='(I08.8)') file_be + (n-1)*file_skip
      fRM%fname = trim(datapath_RM)//'rescalemean_'//fnum//'.h5'
      print *, 'reading file: ', trim(fRM%fname)
      call ReadHDF5_3D(fRM,buffer_RM_rd(:,:,:,1:4))
      do k=1, kdim_RM
        vartmp(:,1,k,2:3) = buffer_RM_rd(k,:,1,1:2)*uinf
        vartmp(:,1,k,4)   = buffer_RM_rd(k,:,1,3)*rhoinf
        vartmp(:,1,k,5)   = buffer_RM_rd(k,:,1,4)*uinf**2
        vartmp(:,1,k,6)   = vartmp(:,1,k,4)*rbar*vartmp(:,1,k,5)
      enddo
      fname = 'rescalemean_'//fnum4_1//'-'//fnum4_2//'.plt'
      if(n.eq.num_file) print *, 'Writing file: ', trim(fname)
      call WriteTec(fname,n)
    enddo ! end n loop
  endif ! iRescaleMean.gt.0


  contains
    subroutine Input()
      integer :: i
      integer :: ist, iend, isp, jst, jend, jsp, kst, kend, ksp, jst2

      read(*,*)
      read(*,'(A)')datapath  !path where data files are stored
      read(*,*)
      read(*,*) Rm
      read(*,*)
      read(*,*) file_be, file_end, file_skip, iAppend, icylin, iBaseflowSubtract, isym
      read(*,*)
      read(*,*) iFlowData, iRescaleMean
      read(*,*)
      read(*,*)
      read(*,*) ist, iend, isp, jst, jend, jsp, kst, kend, ksp
      read(*,*)
      read(*,*) nvar_vol_output
      allocate(varindex_vol_output(nvar_vol_output), dname_vol_output(nvar_vol_output))
      read(*,*) (varindex_vol_output(i), i=1,nvar_vol_output)
      read(*,*)
      read(*,*)
      read(*,'(A)') datapath_RM
      read(*,*)
      read(*,'(A)') gridpath_RM
      read(*,*)
      read(*,*) istencilsize, kdim_RM
      read(*,*)
      read(*,*) uinf, rhoinf

      if(iFlowData.gt.0) then
        if(iend.lt.ist.or.jend.lt.jst.or.kend.lt.kst) then
          print *, 'selected region out of dimension'
          stop
        endif

        write(unit=fibe, fmt='(I04.4)') ist
        write(unit=fiend,fmt='(I04.4)') iend
        write(unit=fjbe, fmt='(I04.4)') jst
        write(unit=fjend,fmt='(I04.4)') jend
        write(unit=fkbe, fmt='(I04.4)') kst
        write(unit=fkend,fmt='(I04.4)') kend

        grd%fname = trim(datapath)//'grid.h5'
        grd%gname = '/'
        call DetectHDF5(grd)
        kmax = grd%dimsf(1); imax = grd%dimsf(2); jmax = grd%dimsf(3)
        print *, 'File dimension: '
        print *, 'imax = ', imax, 'jmax = ', jmax, 'kmax = ', kmax

        if(jend.eq.jmax .and. icylin.eq.1) IsCylPeriod = .true.

        ilen_wt = (iend-ist)/isp + 1;jlen_wt = (jend-jst)/jsp + 1;klen_wt = (kend-kst)/ksp + 1

        if(jst.eq.jend .and. jsp.eq.1 .and. icylin.eq.1 .and. isym.eq.1) IsAxiSym = .true.

        if(ilen_wt.eq.1) then
          ilen = 5
        else
          ilen = ilen_wt
        endif
        if(jlen_wt.eq.1) then
          jlen = 5
        else
          jlen = jlen_wt
        endif
        if(klen_wt.eq.1) then
          klen = 5
        else
          klen = klen_wt
        endif

        iend = ist + (ilen-1)*isp; jend = jst + (jlen-1)*jsp; kend = kst + (klen-1)*ksp

        if(iend.gt.grd%dimsf(2).or.jend.gt.grd%dimsf(3).or.kend.gt.grd%dimsf(1)) then
          print *, 'selected region out of dimension ... '
          stop
        endif
        print *, '##########################################################'
        print *, 'selected region: '
        print *, 'ist = ', ist, 'iend = ', iend, 'isp = ', isp
        print *, 'jst = ', jst, 'jend = ', jend, 'jsp = ', jsp
        print *, 'kst = ', kst, 'kend = ', kend, 'ksp = ', ksp
        print *, '##########################################################'
        allocate(buffer_grid(klen,ilen,jlen,3))
        grd%dname(1) = 'x'; grd%dname(2) = 'y'; grd%dname(3) = 'z'
        grd%block  = 1
        grd%dimsm  = (/klen,ilen,jlen/)
        grd%offset = (/kst-1,ist-1,jst-1/)
        grd%stride = (/ksp,isp,jsp/)
        grd%count  = (/klen,ilen,jlen/)
        print *, 'Reading file: ', trim(grd%fname)
        call ReadHDF5_3D_test(grd,buffer_grid)
        allocate(x(ilen,jlen,klen),y(ilen,jlen,klen),z(ilen,jlen,klen))

        do k=1, klen
        do j=1, jlen
          do i=1, ilen
            x(i,j,k) = buffer_grid(k,i,j,1)
            y(i,j,k) = buffer_grid(k,i,j,2)
            z(i,j,k) = buffer_grid(k,i,j,3)
          enddo
        enddo
        enddo
        deallocate(buffer_grid)

        call InitFlowHDF5(fsol)
        fsol%dimsf  = (/kmax,imax,jmax/)
        fsol%block  = 1
        fsol%dimsm  = grd%dimsm
        fsol%offset = grd%offset
        fsol%stride = grd%stride
        fsol%count  = grd%count
        allocate(buffer_flow(ilen,jlen,klen,6))

        idx_vol%ist = 1; idx_vol%iend = ilen_wt; idx_vol%isp = 1
        idx_vol%jst = 1; idx_vol%jend = jlen_wt; idx_vol%jsp = 1
        idx_vol%kst = 1; idx_vol%kend = klen_wt; idx_vol%ksp = 1

        if(IsAxiSym) then
          jst2 = mod(jst+jmax/2,jmax)
          print*,'Axisymmetric 2D cut is turned on!'
          print*,'jst: ',jst,' , jst_sym:',jst2
          grd2 = grd

          grd2%block  = 1
          grd2%dimsm  = (/klen,ilen,jlen/)
          grd2%offset = (/kst-1,ist-1,jst2-1/)
          grd2%stride = (/ksp,isp,jsp/)
          grd2%count  = (/klen,ilen,jlen/)
          allocate(buffer_grid2(klen,ilen,jlen,3))
          call ReadHDF5_3D_test(grd2,buffer_grid2)
          allocate(x2(ilen,jlen,klen),y2(ilen,jlen,klen),z2(ilen,jlen,klen))

          do k=1, klen
            do j=1, jlen
              do i=1, ilen
                x2(i,j,k) = buffer_grid2(k,i,j,1)
                y2(i,j,k) = buffer_grid2(k,i,j,2)
                z2(i,j,k) = buffer_grid2(k,i,j,3)
              enddo
            enddo
          enddo
          deallocate(buffer_grid2)

          call InitFlowHDF5(fsol2)
          fsol2%dimsf  = (/kmax,imax,jmax/)
          fsol2%block  = 1
          fsol2%dimsm  = grd2%dimsm
          fsol2%offset = grd2%offset
          fsol2%stride = grd2%stride
          fsol2%count  = grd2%count
          allocate(buffer_flow2(ilen,jlen,klen,6))

          idx_vol2%ist = 1; idx_vol2%iend = ilen_wt; idx_vol2%isp = 1
          idx_vol2%jst = 1; idx_vol2%jend = jlen_wt; idx_vol2%jsp = 1
          idx_vol2%kst = 1; idx_vol2%kend = klen_wt; idx_vol2%ksp = 1
        endif ! endif IsAxiSym
      endif ! end iFlowData.gt.0

      if(iRescaleMean.gt.0) then
        call InitFlowHDF5(fRM)
        fRM%dnum = 4
        fRM%dname(1) = 'uave_recycle'
        fRM%dname(2) = 'wave_recycle'
        fRM%dname(3) = 'rhoave_recycle'
        fRM%dname(4) = 'tave_recycle'
        fRM%dimsf  = (/kdim_RM,istencilsize,1/)
        fRM%block  = fRM%dimsf
        fRM%dimsm  = fRM%dimsf
        fRM%offset = 0
        fRM%stride = 1
        fRM%count  = 1
        allocate(buffer_RM(istencilsize,1,kdim_RM,6))
        allocate(buffer_RM_rd(kdim_RM,istencilsize,1,4)) ! u, w, rho, T

        grd%fname = trim(gridpath_RM)//'grid.h5'
        grd%gname = '/'
        call DetectHDF5(grd)
        grd%dnum = 1
        grd%dname(1) = 'z'
        grd%block  = 1
        grd%dimsm  = (/kdim_RM,1,1/)
        grd%offset = (/0,0,0/)
        grd%stride = (/1,1,1/)
        grd%count  = (/kdim_RM,1,1/)
        print *, 'Reading file: ', trim(grd%fname)
        allocate(buffer_grid(kdim_RM,1,1,1))
        call ReadHDF5_3D_test(grd,buffer_grid)
        do k=1, kdim_RM
          buffer_RM(:,1,k,1)   = buffer_grid(k,1,1,1)
        enddo
        deallocate(buffer_grid)
      endif

        rbar = R/Rm
        print *, 'rbar = ', rbar
      return
    end subroutine Input


end program volume    
