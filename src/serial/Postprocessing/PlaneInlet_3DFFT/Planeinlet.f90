
   ! convert timesereis volume data files to plane inlet files

   Program planeinlet
     use modRWHDF5
     use modTecbin
     use modHelmholtzDecomp
     implicit none

     real(8), parameter :: R=8314.3D0, iwrite3D = 0
     integer :: istencilsize = 4
     integer :: i, j, k, kk, n, jj, nn
     integer :: imax, jmax, kmax, imax_n, jmax_n, kmax_n
     integer :: file_be, file_end, file_skip, num_file
     real(8), dimension(:,:,:,:,:), allocatable :: varstmp, varsnew
     real(8), dimension(:,:,:,:), allocatable :: vars_rd, vars_wt, buffer_grd, buffer_pgrd
     real(8), dimension(:,:,:,:), allocatable :: buffer_grd_tmp

     real(8), dimension(:), allocatable :: ttime, ttime_total
     character(4) :: fnum4, fnum_iloc
     character(8) :: fnum8

     integer :: iplanes, nbuffer, nperfile
     character(400) :: fname, filepath, filepath_output, gridfile, planegrid

     real(8), dimension(5) :: varave  ! u, v, w, p, T
     integer :: nmax, iwindow, ntrans, nperiod, isubave, ifilter_np, Pfile_skip
     real(8) :: pi, zlen, ave1, ave2, ave3, ave4
     complex, dimension(:), allocatable :: vars1d
     real(8) :: tmean, tmeanenergy, tmean2, tmeanenergy2

     type(tp_rdwt_hdf5) :: fsol, f1D, grd, pgrd, psol, fsol_wt, grd_wt
     type(tp_DNSIndex)  :: fvol, pvol, tsplane

     character(200) :: varname_output, varname_plane
     real(8) :: xbeg, xend, ybeg, yend, zbeg, zend
     real(8), dimension(:,:), allocatable :: vars_ave
     integer :: count_file

     real(8) :: dx, dy, dz, Lx, Ly, Lz, rbar, Rm, didx, djdy, dkdz, xloc_output
     real(8), dimension(:,:,:), allocatable :: dudi, dudj, dudk, dvdi, dvdj, dvdk, dwdi, dwdj, dwdk
     real(8), dimension(:,:,:,:), allocatable :: buffer_var, buffer_sole, buffer_dila, buffer_rd, buffer_rd_tmp
     integer :: num_iplane, iloc_output
     real(8), dimension(:), allocatable :: xgrid, ygrid, zgrid, xgrido, ygrido, zgrido

     varname_output = 'x y z u v w omx omy omz div' ! p T'
     pi = 4.d0*atan(1.)

     call InitHDF5()
     call InitGridHDF5(grd)
     call InitGridHDF5(grd_wt)
     call InitFlowHDF5(fsol)
     call InitFlowHDF5(fsol_wt)
     call InitFlow_PlaneInlet(psol)
     call Input()

     call ReadGridfiles()

     ! init 1D time data
     call InitFlow_1D(f1D)
     f1D%dimsf = (/nperfile/)
     allocate(ttime(nperfile))

     ! initialize 3D FFT
     iwindow = 0
     call InitHDecomp3d(imax, jmax_n, kmax_n, iwindow, dx, dy, dz)

     write(30,*) '###############################################'
     write(30,*) '# of data points per segments in space domain kx = ', imax
     write(30,*) '# of segments in space domain = 1 '
     write(30,*) 'Interval of sampling (m) =', dx
     write(30,*) 'Length per segments (m) =', dble(imax)*dx
     write(30,*) '###############################################'
     write(30,*) '# of data points per segments in space domain ky = ', jmax_n
     write(30,*) '# of segments in space domain = 1 '
     write(30,*) 'Interval of sampling (m) =', dy
     write(30,*) 'Length per segments (m) =', dble(jmax_n)*dy
     write(30,*) '###############################################'
     write(30,*) '# of data points per segments in space domain kz = ', kmax_n
     write(30,*) '# of segments in space domain = 1 '
     write(30,*) 'Interval of sampling (m) =', dz
     write(30,*) 'Length per segments (m) =', dble(kmax_n)*dz

     fsol%dimsf = (/kmax,imax,jmax/)
     fsol%offset(1) = pvol%kbe -1; fsol%offset(2) = pvol%ibe -1; fsol%offset(3) = pvol%jbe -1;

     allocate(buffer_rd(kmax,imax,jmax,1),buffer_rd_tmp(kmax_n,imax,jmax_n,1))
     allocate(buffer_var(imax,jmax_n,kmax_n,5),buffer_sole(imax,jmax_n,kmax_n,3),buffer_dila(imax,jmax_n,kmax_n,3))

     ! information for writing new 3D dataset
     fsol_wt%dimsf = (/kmax_n,imax,jmax_n/)
     fsol_wt%dimsm = fsol_wt%dimsf
     fsol_wt%block = fsol_wt%dimsf
     fsol_wt%offset = 0

     open(7,file='planeindex.dat',status='unknown')
       fsol%dimsf  = grd%dimsf
       count_file = 0
       fsol%dimsm = fsol%dimsf
       fsol%block = fsol%dimsf

     ! begin the main code
     do n=1+Pfile_skip*nperfile, num_file, nperfile

       do jj=1, nperfile
         write(unit=fnum8,fmt='(I08.8)') (file_be + (jj+n-1-1)*file_skip)
         fsol%fname = trim(filepath)//'timeseriesVol_'//fnum8//'.h5'
         print *, 'Reading file: ', trim(fsol%fname)
         fsol%gname = 'vol'

         ! information for new 3D dataset
         fsol_wt%fname = 'timeseriesVol_'//fnum8//'.h5'
         fsol_wt%gname = 'vol'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! writing Debug file
         if(n.eq.nperfile.and.jj.eq.1) then
           varname_plane = 'x, y, z, u, v, w, p, T, pave, tave, p2, t2'
           call InitTec(1,imax,jmax,kmax, 12,trim(varname_plane),0)
           do i=1, imax
             vartmp(i,:,:,1) = xgrido(i)
           enddo
           do j=1, jmax
             vartmp(:,j,:,2) = ygrido(j)
           enddo
           do k=1, kmax
             vartmp(:,:,k,3) = zgrido(k)
           enddo
         endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         do nn=1, 5
           call ReadHDF5_3D_1V(fsol,buffer_rd(:,:,:,1),nn)
           call InterpSol(xgrido,ygrido,zgrido,buffer_rd(:,:,:,:), &
                          xgrido,ygrid,zgrid,buffer_rd_tmp(:,:,:,:))
           do k=1, kmax_n
             buffer_var(1:imax,1:jmax_n,k,nn) = buffer_rd_tmp(k,1:imax,1:jmax_n,1)
           enddo
           if(n.eq.nperfile.and.jj.eq.1) then !! debug
             do k=1, kmax
               vartmp(1:imax,1:jmax,k,3+nn) = buffer_rd(k,1:imax,1:jmax,1)
             enddo
           endif
         enddo ! end nn loop
         fsol%gname = '/'
         call ReadHDF5_scalar(fsol,ttime(jj))
         print *, 'reading file info: time = ', ttime(jj)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!! debug file !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! writing Debug file
         if(n.eq.nperfile.and.jj.eq.1) then
           fname = 'original_data_beforeInterp.plt'

           do i=1, imax
             do k=1, kmax
               vartmp(i,:,k,9)  = sum(vartmp(i,:,k,7))/dble(jmax)
               vartmp(i,:,k,10) = sum(vartmp(i,:,k,8))/dble(jmax)
               vartmp(i,:,k,11) = sum(vartmp(i,:,k,7)**2)/dble(jmax)
               vartmp(i,:,k,12) = sum(vartmp(i,:,k,8)**2)/dble(jmax)
             enddo
           enddo
           print *, 'Writing file: ', trim(fname)
           call WriteTec(fname)

           allocate(dudi(imax,jmax_n,kmax_n),dudj(imax,jmax_n,kmax_n),dudk(imax,jmax_n,kmax_n))
           allocate(dvdi(imax,jmax_n,kmax_n),dvdj(imax,jmax_n,kmax_n),dvdk(imax,jmax_n,kmax_n))
           allocate(dwdi(imax,jmax_n,kmax_n),dwdj(imax,jmax_n,kmax_n),dwdk(imax,jmax_n,kmax_n))
           didx = 1./dx
           djdy = 1./dy
           dkdz = 1./dz
           call Calddijk(imax,jmax_n,kmax_n,buffer_var(:,:,:,1:3),dudi,dudj,dudk,dvdi,dvdj,dvdk,dwdi,dwdj,dwdk)
           call InitTec(1,imax,jmax_n,kmax_n, 12,trim(varname_plane),0)
           do i=1, imax
             vartmp(i,:,:,1) = xgrido(i)
           enddo
           do j=1, jmax_n
             vartmp(:,j,:,2) = ygrid(j)
           enddo
           do k=1, kmax_n
             vartmp(:,:,k,3) = zgrid(k)
           enddo
           do k=1, kmax_n
             vartmp(1:imax,1:jmax_n,k,4) = buffer_var(1:imax,1:jmax_n,k,1)
             vartmp(1:imax,1:jmax_n,k,5) = buffer_var(1:imax,1:jmax_n,k,2)
             vartmp(1:imax,1:jmax_n,k,6) = buffer_var(1:imax,1:jmax_n,k,3)
             vartmp(1:imax,1:jmax_n,k,7) = buffer_var(1:imax,1:jmax_n,k,4)
             vartmp(1:imax,1:jmax_n,k,8) = buffer_var(1:imax,1:jmax_n,k,5)
!             vartmp(1:imax,1:jmax_n,k,7) = dwdj(1:imax,1:jmax_n,k)*djdy - dvdk(1:imax,1:jmax_n,k)*dkdz
!             vartmp(1:imax,1:jmax_n,k,8) = dudk(1:imax,1:jmax_n,k)*dkdz - dwdi(1:imax,1:jmax_n,k)*didx
!             vartmp(1:imax,1:jmax_n,k,9) = dvdi(1:imax,1:jmax_n,k)*didx - dudj(1:imax,1:jmax_n,k)*djdy
!             vartmp(1:imax,1:jmax_n,k,10)= dudi(1:imax,1:jmax_n,k)*didx + dvdj(1:imax,1:jmax_n,k)*djdy + dwdk(1:imax,1:jmax_n,k)*dkdz
           enddo
           do i=1, imax
             do k=1, kmax_n
               vartmp(i,:,k,9)  = sum(vartmp(i,:,k,7))/dble(jmax_n)
               vartmp(i,:,k,10) = sum(vartmp(i,:,k,8))/dble(jmax_n)
               vartmp(i,:,k,11) = sum(vartmp(i,:,k,7)**2)/dble(jmax_n)
               vartmp(i,:,k,12) = sum(vartmp(i,:,k,8)**2)/dble(jmax_n)
             enddo
           enddo

           fname = 'original_data_afterInterp.plt'
           print *, 'Writing file: ', trim(fname)
           call WriteTec(fname)
         endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         ! FFT filter
         do nn=1, 5
           call FFTFilter3d_2(buffer_var(:,:,:,nn),buffer_var(:,:,:,nn),1./Lx,1./(dble(ifilter_np)*dx),1./Ly,1./(dble(ifilter_np)*dy),1./Lz,1./(dble(ifilter_np)*dz))
          ! call FFTFilter3d_3(buffer_var(:,:,:,nn),buffer_var(:,:,:,nn))
         enddo

         ! Decomp
         call DoHDecomp3d(buffer_var(:,:,:,1:3),buffer_sole,buffer_dila)

         do nn=1, 5
           if(nn.le.3) then
             do k=1, kmax_n
               buffer_rd_tmp(k,1:imax,1:jmax_n,1) = buffer_dila(1:imax,1:jmax_n,k,nn)
             enddo
           else
             do k=1, kmax_n
               buffer_rd_tmp(k,1:imax,1:jmax_n,1) = buffer_var(1:imax,1:jmax_n,k,nn)
             enddo
           endif
           call InterpSol(xgrido, ygrid, zgrid, buffer_rd_tmp(:,:,:,1:1), &
                 xgrid(3:imax_n), ygrid, zgrid, varsnew(:,3:pgrd%dimsf(2),:,nn:nn,jj) )
           ! writing new 3D dataset
           if(iwrite3D.eq.1) then
             if(nn.eq.1) then
               print *, 'writing file: ', trim(fsol_wt%fname)
               call WriteHDF5_3D_1v(fsol_wt,buffer_rd_tmp(:,:,:,1)+varave(1),nn)
             else
               call WriteHDF5_3D_1v(fsol_wt,buffer_rd_tmp(:,:,:,1),nn)
             endif
             if(nn.eq.5) then
               fsol_wt%gname = '/'
               call WriteHDF5_scalar(fsol_wt,ttime(jj))
             endif
           endif
         enddo ! end nn loop

       enddo ! end jj loop

       varsnew(:,:,:,4,:) = varsnew(:,:,:,4,:)/(rbar*varsnew(:,:,:,5,:))

       ! rescale Uave
       varsnew(:,3:pgrd%dimsf(2),:,1,:) = varsnew(:,3:pgrd%dimsf(2),:,1,:) + varave(1)

       ! write planeinlet files
       write(unit=fnum4,fmt='(I04.4)') n/nperfile
       psol%fname = trim(filepath_output)//'planedata_'//fnum4//'.h5'
       psol%dimsf = (/int(pgrd%dimsf(1)),int(pgrd%dimsf(3)),nperfile/)
       print *, 'Writing file: ', trim(psol%fname)
       do i=1, num_iplane
         write(unit=fnum_iloc,fmt='(I04.4)') iloc_output + i + 2
         psol%gname = 'i'//trim(fnum_iloc)
         call WriteHDF5_planeinlet(psol,varsnew(:,i+2,:,:,:),i)
       enddo ! end i loop
       f1D%fname = trim(psol%fname)
       call WriteHDF5_1D_1V(f1D,ttime)
       write(7,'(I8,2E23.15)') n/nperfile, ttime(1), ttime(nperfile)

       ! write 3D debug *.plt files
       if(n.eq.nperfile) then
         varname_plane = 'x, y, z, u, v, w, rho, T, rhoave, tave, rho2, t2'
         call InitTec(1,int(pgrd%dimsf(2)-2),int(pgrd%dimsf(3)),int(pgrd%dimsf(1)), 12,trim(varname_plane),0)
         do k=1, pgrd%dimsf(1)
           vartmp(1:pgrd%dimsf(2)-2,1:pgrd%dimsf(3),k,1) = buffer_pgrd(k,3:pgrd%dimsf(2),1:pgrd%dimsf(3),1)
           vartmp(1:pgrd%dimsf(2)-2,1:pgrd%dimsf(3),k,2) = buffer_pgrd(k,3:pgrd%dimsf(2),1:pgrd%dimsf(3),2)
           vartmp(1:pgrd%dimsf(2)-2,1:pgrd%dimsf(3),k,3) = buffer_pgrd(k,3:pgrd%dimsf(2),1:pgrd%dimsf(3),3)
           vartmp(1:pgrd%dimsf(2)-2,1:pgrd%dimsf(3),k,4) = varsnew(k,3:pgrd%dimsf(2),1:pgrd%dimsf(3),1,1)
           vartmp(1:pgrd%dimsf(2)-2,1:pgrd%dimsf(3),k,5) = varsnew(k,3:pgrd%dimsf(2),1:pgrd%dimsf(3),2,1)
           vartmp(1:pgrd%dimsf(2)-2,1:pgrd%dimsf(3),k,6) = varsnew(k,3:pgrd%dimsf(2),1:pgrd%dimsf(3),3,1)
           vartmp(1:pgrd%dimsf(2)-2,1:pgrd%dimsf(3),k,7) = varsnew(k,3:pgrd%dimsf(2),1:pgrd%dimsf(3),4,1)
           vartmp(1:pgrd%dimsf(2)-2,1:pgrd%dimsf(3),k,8) = varsnew(k,3:pgrd%dimsf(2),1:pgrd%dimsf(3),5,1)
         enddo

         do i=1, pgrd%dimsf(2)-2
           do k=1, pgrd%dimsf(1)
             vartmp(i,:,k,9)  = sum(vartmp(i,:,k,7))/dble(pgrd%dimsf(3))
             vartmp(i,:,k,10) = sum(vartmp(i,:,k,8))/dble(pgrd%dimsf(3))
             vartmp(i,:,k,11) = sum(vartmp(i,:,k,7)**2)/dble(pgrd%dimsf(3))
             vartmp(i,:,k,12) = sum(vartmp(i,:,k,8)**2)/dble(pgrd%dimsf(3))
           enddo
         enddo

         fname = 'planeinlet.plt'
         print *, 'Writing file: ', trim(fname)
         call WriteTec(fname)

         varname_plane = 'x, y, z, u, v, w, p, T, pave, tave, p2, t2'
         call Calddijk(imax,jmax_n,kmax_n,buffer_var(:,:,:,1:3),dudi,dudj,dudk,dvdi,dvdj,dvdk,dwdi,dwdj,dwdk)
         call InitTec(1,imax,jmax_n,kmax_n, 12,trim(varname_plane),0)
         do i=1, imax
           vartmp(i,:,:,1) = xgrido(i)
         enddo
         do j=1, jmax_n
           vartmp(:,j,:,2) = ygrid(j)
         enddo
         do k=1, kmax_n
           vartmp(:,:,k,3) = zgrid(k)
         enddo
         do k=1, kmax_n
           vartmp(1:imax,1:jmax_n,k,4) = buffer_var(1:imax,1:jmax_n,k,1)
           vartmp(1:imax,1:jmax_n,k,5) = buffer_var(1:imax,1:jmax_n,k,2)
           vartmp(1:imax,1:jmax_n,k,6) = buffer_var(1:imax,1:jmax_n,k,3)
           vartmp(1:imax,1:jmax_n,k,7) = buffer_var(1:imax,1:jmax_n,k,4)
           vartmp(1:imax,1:jmax_n,k,8) = buffer_var(1:imax,1:jmax_n,k,5)
!           vartmp(1:imax,1:jmax_n,k,9) = dvdi(1:imax,1:jmax_n,k)*didx - dudj(1:imax,1:jmax_n,k)*djdy
!           vartmp(1:imax,1:jmax_n,k,10)= dudi(1:imax,1:jmax_n,k)*didx + dvdj(1:imax,1:jmax_n,k)*djdy + dwdk(1:imax,1:jmax_n,k)*dkdz
         enddo
         do i=1, imax
           do k=1, kmax_n
             vartmp(i,:,k,9)  = sum(vartmp(i,:,k,7))/dble(jmax_n)
             vartmp(i,:,k,10) = sum(vartmp(i,:,k,8))/dble(jmax_n)
             vartmp(i,:,k,11) = sum(vartmp(i,:,k,7)**2)/dble(jmax_n)
             vartmp(i,:,k,12) = sum(vartmp(i,:,k,8)**2)/dble(jmax_n)
           enddo
         enddo

         fname = 'original_data_FFT.plt'
         print *, 'Writing file: ', trim(fname)
         call WriteTec(fname)

         call Calddijk(imax,jmax_n,kmax_n,buffer_dila,dudi,dudj,dudk,dvdi,dvdj,dvdk,dwdi,dwdj,dwdk)
         do i=1, imax
           vartmp(i,:,:,1) = xgrido(i)
         enddo
         do j=1, jmax_n
           vartmp(:,j,:,2) = ygrid(j)
         enddo
         do k=1, kmax_n
           vartmp(:,:,k,3) = zgrid(k)
         enddo
         do k=1, kmax_n
           vartmp(1:imax,1:jmax_n,k,4) = buffer_dila(1:imax,1:jmax_n,k,1)
           vartmp(1:imax,1:jmax_n,k,5) = buffer_dila(1:imax,1:jmax_n,k,2)
           vartmp(1:imax,1:jmax_n,k,6) = buffer_dila(1:imax,1:jmax_n,k,3)
           vartmp(1:imax,1:jmax_n,k,7) = dwdj(1:imax,1:jmax_n,k)*djdy - dvdk(1:imax,1:jmax_n,k)*dkdz
           vartmp(1:imax,1:jmax_n,k,8) = dudk(1:imax,1:jmax_n,k)*dkdz - dwdi(1:imax,1:jmax_n,k)*didx
           vartmp(1:imax,1:jmax_n,k,9) = dvdi(1:imax,1:jmax_n,k)*didx - dudj(1:imax,1:jmax_n,k)*djdy
           vartmp(1:imax,1:jmax_n,k,10)= dudi(1:imax,1:jmax_n,k)*didx + dvdj(1:imax,1:jmax_n,k)*djdy + dwdk(1:imax,1:jmax_n,k)*dkdz
         enddo
         fname = 'new_data_dila.plt'
         print *, 'Writing file: ', trim(fname)
         call WriteTec(fname)

         call Calddijk(imax,jmax_n,kmax_n,buffer_sole,dudi,dudj,dudk,dvdi,dvdj,dvdk,dwdi,dwdj,dwdk)
         do i=1, imax
           vartmp(i,:,:,1) = xgrido(i)
         enddo
         do j=1, jmax_n
           vartmp(:,j,:,2) = ygrid(j)
         enddo
         do k=1, kmax_n
           vartmp(:,:,k,3) = zgrid(k)
         enddo
         do k=1, kmax_n
           vartmp(1:imax,1:jmax_n,k,4) = buffer_sole(1:imax,1:jmax_n,k,1)
           vartmp(1:imax,1:jmax_n,k,5) = buffer_sole(1:imax,1:jmax_n,k,2)
           vartmp(1:imax,1:jmax_n,k,6) = buffer_sole(1:imax,1:jmax_n,k,3)
           vartmp(1:imax,1:jmax_n,k,7) = dwdj(1:imax,1:jmax_n,k)*djdy - dvdk(1:imax,1:jmax_n,k)*dkdz
           vartmp(1:imax,1:jmax_n,k,8) = dudk(1:imax,1:jmax_n,k)*dkdz - dwdi(1:imax,1:jmax_n,k)*didx
           vartmp(1:imax,1:jmax_n,k,9) = dvdi(1:imax,1:jmax_n,k)*didx - dudj(1:imax,1:jmax_n,k)*djdy
           vartmp(1:imax,1:jmax_n,k,10)= dudi(1:imax,1:jmax_n,k)*didx + dvdj(1:imax,1:jmax_n,k)*djdy + dwdk(1:imax,1:jmax_n,k)*dkdz
         enddo
         fname = 'new_data_sole.plt'
         print *, 'Writing file: ', trim(fname)
         call WriteTec(fname)
       endif !

     enddo ! end n loop

     close(7)
 
     call Writeinput()
     call FinalizeHDF5()

contains

     subroutine Calddijk(idim,jdim,kdim,varin,dudi,dudj,dudk,dvdi,dvdj,dvdk,dwdi,dwdj,dwdk)
       implicit none
       integer, intent(in) :: idim, jdim, kdim
       real(8), intent(in) :: varin(idim,jdim,kdim,3)
       real(8), intent(out) :: dudi(idim,jdim,kdim), dudj(idim,jdim,kdim), dudk(idim,jdim,kdim), dvdi(idim,jdim,kdim), dvdj(idim,jdim,kdim), dvdk(idim,jdim,kdim), dwdi(idim,jdim,kdim), dwdj(idim,jdim,kdim), dwdk(idim,jdim,kdim)
       integer :: i, j, k

         do j=1, jdim
           do k=1, kdim
             do i=1, idim
               if(i.eq.1) then
                 dudi(i,j,k) = (varin(idim-1,j,k,1)-8.*(varin(idim,j,k,1)-varin(2,j,k,1))-varin(3,j,k,1))/12.d0
                 dvdi(i,j,k) = (varin(idim-1,j,k,2)-8.*(varin(idim,j,k,2)-varin(2,j,k,2))-varin(3,j,k,2))/12.d0
                 dwdi(i,j,k) = (varin(idim-1,j,k,3)-8.*(varin(idim,j,k,3)-varin(2,j,k,3))-varin(3,j,k,3))/12.d0
               elseif(i.eq.2) then
                 dudi(i,j,k) = (varin(idim,j,k,1)-8.*(varin(1,j,k,1)-varin(3,j,k,1))-varin(4,j,k,1))/12.d0
                 dvdi(i,j,k) = (varin(idim,j,k,2)-8.*(varin(1,j,k,2)-varin(3,j,k,2))-varin(4,j,k,2))/12.d0
                 dwdi(i,j,k) = (varin(idim,j,k,3)-8.*(varin(1,j,k,3)-varin(3,j,k,3))-varin(4,j,k,3))/12.d0
               elseif(i.eq.idim-1) then
                 dudi(i,j,k) = (varin(idim-3,j,k,1)-8.*(varin(idim-2,j,k,1)-varin(idim,j,k,1))-varin(1,j,k,1))/12.d0
                 dvdi(i,j,k) = (varin(idim-3,j,k,2)-8.*(varin(idim-2,j,k,2)-varin(idim,j,k,2))-varin(1,j,k,2))/12.d0
                 dwdi(i,j,k) = (varin(idim-3,j,k,3)-8.*(varin(idim-2,j,k,3)-varin(idim,j,k,3))-varin(1,j,k,3))/12.d0
               elseif(i.eq.idim) then
                 dudi(i,j,k) = (varin(idim-2,j,k,1)-8.*(varin(idim-1,j,k,1)-varin(1,j,k,1))-varin(2,j,k,1))/12.d0
                 dvdi(i,j,k) = (varin(idim-2,j,k,2)-8.*(varin(idim-1,j,k,2)-varin(1,j,k,2))-varin(2,j,k,2))/12.d0
                 dwdi(i,j,k) = (varin(idim-2,j,k,3)-8.*(varin(idim-1,j,k,3)-varin(1,j,k,3))-varin(2,j,k,3))/12.d0
               else
                 dudi(i,j,k) = (varin(i-2,j,k,1)-8.*(varin(i-1,j,k,1)-varin(i+1,j,k,1))-varin(i+2,j,k,1))/12.d0
                 dvdi(i,j,k) = (varin(i-2,j,k,2)-8.*(varin(i-1,j,k,2)-varin(i+1,j,k,2))-varin(i+2,j,k,2))/12.d0
                 dwdi(i,j,k) = (varin(i-2,j,k,3)-8.*(varin(i-1,j,k,3)-varin(i+1,j,k,3))-varin(i+2,j,k,3))/12.d0
               endif
             enddo
           enddo
         enddo

         do i=1, idim
           do k=1, kdim
             do j=1, jdim
               if(j.eq.1) then
                 dudj(i,j,k) = (varin(i,jdim-1,k,1)-8.*(varin(i,jdim,k,1)-varin(i,2,k,1))-varin(i,3,k,1))/12.d0
                 dvdj(i,j,k) = (varin(i,jdim-1,k,2)-8.*(varin(i,jdim,k,2)-varin(i,2,k,2))-varin(i,3,k,2))/12.d0
                 dwdj(i,j,k) = (varin(i,jdim-1,k,3)-8.*(varin(i,jdim,k,3)-varin(i,2,k,3))-varin(i,3,k,3))/12.d0
               elseif(j.eq.2) then
                 dudj(i,j,k) = (varin(i,jdim,k,1)-8.*(varin(i,1,k,1)-varin(i,3,k,1))-varin(i,4,k,1))/12.d0
                 dvdj(i,j,k) = (varin(i,jdim,k,2)-8.*(varin(i,1,k,2)-varin(i,3,k,2))-varin(i,4,k,2))/12.d0
                 dwdj(i,j,k) = (varin(i,jdim,k,3)-8.*(varin(i,1,k,3)-varin(i,3,k,3))-varin(i,4,k,3))/12.d0
               elseif(j.eq.jdim-1) then
                 dudj(i,j,k) = (varin(i,jdim-3,k,1)-8.*(varin(i,jdim-2,k,1)-varin(i,jdim,k,1))-varin(i,1,k,1))/12.d0
                 dvdj(i,j,k) = (varin(i,jdim-3,k,2)-8.*(varin(i,jdim-2,k,2)-varin(i,jdim,k,2))-varin(i,1,k,2))/12.d0
                 dwdj(i,j,k) = (varin(i,jdim-3,k,3)-8.*(varin(i,jdim-2,k,3)-varin(i,jdim,k,3))-varin(i,1,k,3))/12.d0
               elseif(j.eq.jdim) then
                 dudj(i,j,k) = (varin(i,jdim-2,k,1)-8.*(varin(i,jdim-1,k,1)-varin(i,1,k,1))-varin(i,2,k,1))/12.d0
                 dvdj(i,j,k) = (varin(i,jdim-2,k,2)-8.*(varin(i,jdim-1,k,2)-varin(i,1,k,2))-varin(i,2,k,2))/12.d0
                 dwdj(i,j,k) = (varin(i,jdim-2,k,3)-8.*(varin(i,jdim-1,k,3)-varin(i,1,k,3))-varin(i,2,k,3))/12.d0
               else
                 dudj(i,j,k) = (varin(i,j-2,k,1)-8.*(varin(i,j-1,k,1)-varin(i,j+1,k,1))-varin(i,j+2,k,1))/12.d0
                 dvdj(i,j,k) = (varin(i,j-2,k,2)-8.*(varin(i,j-1,k,2)-varin(i,j+1,k,2))-varin(i,j+2,k,2))/12.d0
                 dwdj(i,j,k) = (varin(i,j-2,k,3)-8.*(varin(i,j-1,k,3)-varin(i,j+1,k,3))-varin(i,j+2,k,3))/12.d0
               endif
             enddo
           enddo
         enddo

         do i=1, idim
           do j=1, jdim
             do k=1, kdim
               if(k.eq.1) then
                 dudk(i,j,k) = (varin(i,j,kdim-1,1)-8.*(varin(i,j,kdim,1)-varin(i,j,2,1))-varin(i,j,3,1))/12.d0
                 dvdk(i,j,k) = (varin(i,j,kdim-1,2)-8.*(varin(i,j,kdim,2)-varin(i,j,2,2))-varin(i,j,3,2))/12.d0
                 dwdk(i,j,k) = (varin(i,j,kdim-1,3)-8.*(varin(i,j,kdim,3)-varin(i,j,2,3))-varin(i,j,3,3))/12.d0
               elseif(k.eq.2) then
                 dudk(i,j,k) = (varin(i,j,kdim,1)-8.*(varin(i,j,1,1)-varin(i,j,3,1))-varin(i,j,4,1))/12.d0
                 dvdk(i,j,k) = (varin(i,j,kdim,2)-8.*(varin(i,j,1,2)-varin(i,j,3,2))-varin(i,j,4,2))/12.d0
                 dwdk(i,j,k) = (varin(i,j,kdim,3)-8.*(varin(i,j,1,3)-varin(i,j,3,3))-varin(i,j,4,3))/12.d0
               elseif(k.eq.kdim-1) then
                 dudk(i,j,k) = (varin(i,j,kdim-3,1)-8.*(varin(i,j,kdim-2,1)-varin(i,j,kdim,1))-varin(i,j,1,1))/12.d0
                 dvdk(i,j,k) = (varin(i,j,kdim-3,2)-8.*(varin(i,j,kdim-2,2)-varin(i,j,kdim,2))-varin(i,j,1,2))/12.d0
                 dwdk(i,j,k) = (varin(i,j,kdim-3,3)-8.*(varin(i,j,kdim-2,3)-varin(i,j,kdim,3))-varin(i,j,1,3))/12.d0
               elseif(k.eq.kdim) then
                 dudk(i,j,k) = (varin(i,j,kdim-2,1)-8.*(varin(i,j,kdim-1,1)-varin(i,j,1,1))-varin(i,j,2,1))/12.d0
                 dvdk(i,j,k) = (varin(i,j,kdim-2,2)-8.*(varin(i,j,kdim-1,2)-varin(i,j,1,2))-varin(i,j,2,2))/12.d0
                 dwdk(i,j,k) = (varin(i,j,kdim-2,3)-8.*(varin(i,j,kdim-1,3)-varin(i,j,1,3))-varin(i,j,2,3))/12.d0
               else
                 dudk(i,j,k) = (varin(i,j,k-2,1)-8.*(varin(i,j,k-1,1)-varin(i,j,k+1,1))-varin(i,j,k+2,1))/12.d0
                 dvdk(i,j,k) = (varin(i,j,k-2,2)-8.*(varin(i,j,k-1,2)-varin(i,j,k+1,2))-varin(i,j,k+2,2))/12.d0
                 dwdk(i,j,k) = (varin(i,j,k-2,3)-8.*(varin(i,j,k-1,3)-varin(i,j,k+1,3))-varin(i,j,k+2,3))/12.d0
               endif
             enddo
           enddo
         enddo

     end subroutine Calddijk

     subroutine Input()
       implicit none
       integer :: i, iloc_tmp1, iloc_tmp2

       read(*,*)
       read(*,*) varave(1:5), Rm
       read(*,*)
       read(*,'(a)') filepath
       read(*,*)
       read(*,'(a)') filepath_output
       read(*,*)
       read(*,*) file_be, file_end, file_skip
       read(*,*)
       read(*,*)
       read(*,*) xbeg, xend, ybeg, yend, zbeg, zend
       read(*,*)
       read(*,*)
       read(*,*) nperfile
       read(*,*)
       read(*,'(a)') planegrid
       read(*,*)
       read(*,*) ifilter_np, Pfile_skip
       read(*,*)
       read(*,*) num_iplane, xloc_output
       !allocate(iplane(num_iplane))
       !read(*,*) (iplane(i),i=1,num_iplane)

       rbar = R/Rm
       print *, 'rbar = ', rbar
       print *, 'ifilter_np = ', ifilter_np

       if(xloc_output.lt.xbeg.or.xloc_output.gt.xend) then
         print *, 'the output location xloc_output is out of the selected region'
         stop
       endif

       open(30,file='debug_info.dat',status='unknown')

       ! read timesereis volume
         num_file = (file_end - file_be)/file_skip + 1
         grd%fname = trim(filepath)//'timeseriesVol_grid.h5'
         call DetectHDF5(grd)
         print *, '*************************'
         write(*,'(a,I8,I8,I8)') 'Timeseries data dimension, imaxo, jmaxo, kmaxo: ', grd%dimsf(2), grd%dimsf(3), grd%dimsf(1)

         allocate(buffer_grd(grd%dimsf(1),grd%dimsf(2),grd%dimsf(3),3))
         call ReadHDF5_3D(grd,buffer_grd)

         call FindIndex(buffer_grd(1,:,1,1),pvol%ibe,pvol%iend,xbeg,xend)
         print *, 'i index for reading: ibe = ', pvol%ibe, 'iend = ', pvol%iend
         call FindIndex(buffer_grd(1,1,:,2),pvol%jbe,pvol%jend,ybeg,yend)
         print *, 'j index for reading: jbe = ', pvol%jbe, 'jend = ', pvol%jend
         call FindIndex(buffer_grd(:,1,1,3),pvol%kbe,pvol%kend,zbeg,zend)
         print *, 'k index for reading: kbe = ', pvol%kbe, 'kend = ', pvol%kend
         print *, '*************************'
         call FindIndex(buffer_grd(1,:,1,1),iloc_tmp1,iloc_tmp2,xloc_output,xloc_output)
         print *, 'the selected index for output the plane inlet data files. i = ', iloc_tmp1
         if(iloc_tmp1.lt.pvol%ibe.or.iloc_tmp1.gt.pvol%iend) then
           print *, 'the selected i-index for output the plane data files is out of bound'
           stop
         endif
         xloc_output = buffer_grd(1,iloc_tmp1,1,1)
         iloc_output = iloc_tmp1
         print *, '*************************'

         imax = pvol%iend - pvol%ibe + 1
         jmax = pvol%jend - pvol%jbe + 1
         kmax = pvol%kend - pvol%kbe + 1

         print *, 'number of samples per file = ', nperfile
         write(*,'(a,I8,I8,I8)') 'Selected region dimensions for reading, imax, jmax, kmax: ', imax, jmax, kmax
         print *, '*************************'

         write(30,*) '*************************'
         write(30,'(a,I8,I8,I8)') 'Timeseries data dimension, imaxo, jmaxo, kmaxo: ', grd%dimsf(2), grd%dimsf(3), grd%dimsf(1)
         write(30,*) 'i index for reading: ibe = ', pvol%ibe, 'iend = ', pvol%iend
         write(30,*) 'j index for reading: jbe = ', pvol%jbe, 'jend = ', pvol%jend
         write(30,*) 'k index for reading: kbe = ', pvol%kbe, 'kend = ', pvol%kend
         write(30,*) '*************************'
         write(30,*) 'number of samples per file = ', nperfile
         write(30,'(a,I8,I8,I8)') 'Selected region dimensions for reading, imax, jmax, kmax: ', imax, jmax, kmax
         write(30,*) '*************************'

     end subroutine Input

     subroutine Writeinput()
       implicit none
       integer :: tmp=0

       fname = trim(filepath_output)//'planeinlet.inp'
       print *, 'writing file: ', trim(fname)
       print *, '   assuming nbuffer=nperfile/10'
       !> assume nbuffer = nperfile/10
       open(7,file=trim(fname),status='unknown')
         write(7,*) 'nfiles    nperfile    nbuffer   iperiod'
         write(7,'(I8,I8,I8,I8)') num_file/nperfile, nperfile, nperfile/10, tmp       
         write(7,*) 'nplanes'
         write(7,*) '4'
         write(7,'(I8,I8,I8,I8)') pvol%ibe, pvol%ibe+1, pvol%ibe+2, pvol%ibe+3
         write(7,*) 'planedatapath'
         write(7,'(a)') trim(filepath_output)
       close(7)

     end subroutine Writeinput

     subroutine ReadGridfiles()
       implicit none
       integer :: i, j, k
       real(8) :: dx1, dx2, dy1, dy2, dz1, dz2

       ! read grid.h5 file
       grd%dimsf = (/kmax,imax,jmax/)
       grd%fname = trim(filepath)//'timeseriesVol_grid.h5'
       if(allocated(buffer_grd)) deallocate(buffer_grd)
       allocate(buffer_grd(kmax,imax,jmax,3))
       print *, 'Reading file: ', trim(grd%fname)
       grd%offset(1) = pvol%kbe - 1
       grd%offset(2) = pvol%ibe - 1
       grd%offset(3) = pvol%jbe - 1
       write(*,'(a,I8,I8,I8)') 'The reading dimension for grid.h5 (kmax,imax,jmax): ', grd%dimsf
       call ReadHDF5_3D(grd, buffer_grd)

       ! resize the x, y, z location. make sure y and z begin with 0. "x=0" is the location for matching the inlet grid file
       buffer_grd(:,:,:,1) = buffer_grd(:,:,:,1) - xloc_output
       buffer_grd(:,:,:,2) = buffer_grd(:,:,:,2) - buffer_grd(1,1,1,2)
       buffer_grd(:,:,:,3) = buffer_grd(:,:,:,3) - buffer_grd(1,1,1,3)

       dx1 = buffer_grd(1,2,1,1) - buffer_grd(1,1,1,1)
       dy1 = buffer_grd(1,1,2,2) - buffer_grd(1,1,1,2)
       dz1 = buffer_grd(2,1,1,3) - buffer_grd(1,1,1,3)

       Lx = buffer_grd(1,imax,1,1) - buffer_grd(1,1,1,1) + dx1
       Ly = buffer_grd(1,1,jmax,2) - buffer_grd(1,1,1,2) + dy1
       Lz = buffer_grd(kmax,1,1,3) - buffer_grd(1,1,1,3) + dz1
       print *, 'test ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
       print *, 'Lx = ', Lx, 'Ly = ', Ly, 'Lz = ', Lz
       print *, 'test ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

       ! read PlaneInletgrid.h5
       call InitGridHDF5(pgrd)
       pgrd%fname = trim(planegrid)
       print *, 'Reading file: ', trim(pgrd%fname)
       call DetectHDF5(pgrd)
       write(*,'(a,I8,I8,I8)') 'The dimension for PlaneInletgrid.h5: ', pgrd%dimsf
       allocate(buffer_pgrd(pgrd%dimsf(1),pgrd%dimsf(2),pgrd%dimsf(3),pgrd%dnum))
       call ReadHDF5_3D(pgrd,buffer_pgrd)

       ! write PlaneInletgrid.h5 to the output directory
       pgrd%fname = trim(filepath_output)//'PlaneInletgrid.h5'
       print *, 'Writing file: ', trim(pgrd%fname)
       call WriteHDF5_3D(pgrd,buffer_pgrd)

       allocate(varsnew(pgrd%dimsf(1), 3:pgrd%dimsf(2), pgrd%dimsf(3), 5, nperfile))

       ! for the plane inlet grid file. make sure x=0, y=0 and z=0 at the first(1,1,1) grid point.
       buffer_pgrd(:,:,:,1) = buffer_pgrd(:,:,:,1) - buffer_pgrd(1,1,1,1)
       buffer_pgrd(:,:,:,2) = buffer_pgrd(:,:,:,2) - buffer_pgrd(1,1,1,2)
       buffer_pgrd(:,:,:,3) = buffer_pgrd(:,:,:,3) - buffer_pgrd(1,1,1,3)

       allocate(xgrido(imax),ygrido(jmax),zgrido(kmax))
       xgrido(1:imax) = buffer_grd(1,1:imax,1,1)
       ygrido(1:jmax) = buffer_grd(1,1,1:jmax,2)
       zgrido(1:kmax) = buffer_grd(1:kmax,1,1,3)

       ! reset  imax, jmax, kmax
       imax_n = pgrd%dimsf(2)
       jmax_n = pgrd%dimsf(3)
       kmax_n = pgrd%dimsf(1)
       allocate(xgrid(pgrd%dimsf(2)),ygrid(jmax_n),zgrid(kmax_n))
       xgrid(1:imax_n) = buffer_pgrd(1,1:imax_n,1,1)
       ygrid(1:jmax_n) = buffer_pgrd(1,1,1:jmax_n,2)
       zgrid(1:kmax_n) = buffer_pgrd(1:kmax_n,1,1,3)

       dx = xgrido(2) - xgrido(1)
       dy = ygrid(2)  - ygrid(1)
       dz = zgrid(2)  - zgrid(1)
       Lx = xgrido(imax)  - xgrido(1) + dx
       Ly = ygrid(jmax_n) - ygrid(1)  + dy
       Lz = zgrid(kmax_n) - zgrid(1)  + dz
       print *, 'dx = ', dx, '(m)'
       print *, 'dy = ', dy, '(m)'
       print *, 'dz = ', dz, '(m)'
       print *, 'Lx = ', Lx, '(m)'
       print *, 'Ly = ', Ly, '(m)'
       print *, 'Lz = ', Lz, '(m)'

       write(30,*) 'dx = ', dx, '(m)'
       write(30,*) 'dy = ', dy, '(m)'
       write(30,*) 'dz = ', dz, '(m)'
       write(30,*) 'Lx = ', Lx, '(m)'
       write(30,*) 'Ly = ', Ly, '(m)'
       write(30,*) 'Lz = ', Lz, '(m)'

       ! write new 3D dataset
       if(iwrite3D.eq.1) then
         allocate(buffer_grd_tmp(kmax_n,imax,jmax_n,3))
         do i=1, imax
           buffer_grd_tmp(:,i,:,1) = xgrido(i)
         enddo
         do j=1, jmax_n
           buffer_grd_tmp(:,:,j,2) = ygrid(j)
         enddo
         do k=1, kmax_n
           buffer_grd_tmp(k,:,:,3) = zgrid(k)
         enddo
         grd_wt%fname = 'timeseries_grid.h5'
         grd_wt%gname = '/'
         grd_wt%dimsf = (/kmax_n,imax,jmax_n/)
         grd_wt%dimsm = grd_wt%dimsf
         grd_wt%block = grd_wt%dimsf
         grd_wt%offset = 0
         print *, 'writing file: ', trim(grd_wt%fname)
         call WriteHDF5_3D(grd_wt,buffer_grd_tmp)
         deallocate(buffer_grd_tmp)
       endif
     end subroutine Readgridfiles

     subroutine FindIndex(varinput,index1,index2,val1,val2)
       implicit none
       real(8), dimension(:), intent(in) :: varinput
       real(8), intent(in) :: val1, val2
       integer, intent(out) :: index1, index2
       integer :: i, dimin

       index1 = 0
       index2 = 0
       dimin = size(varinput)
       do i=1, dimin
         if(varinput(i).ge.val1) then
           index1 = i
           exit
         elseif(i.eq.dimin) then
           print *, 'cannot find index'
           stop
         endif
       enddo
       do i=1, dimin
         if(varinput(i).ge.val2) then
           index2 = i
           exit
         elseif(i.eq.dimin) then
           print *, 'cannot find index'
           stop
         endif
       enddo
     end subroutine FindIndex

     subroutine InterpSol(xold,yold,zold,profilevars_old, &
                          xnew,ynew,znew,profilevars_new)
       implicit none
       real(8),dimension(:), intent(in) :: xold, yold, zold
       real(8),dimension(:), intent(in) :: xnew, ynew, znew
       real(8),dimension(:,:,:,:), intent(in) :: profilevars_old
       real(8),dimension(:,:,:,:), intent(out) :: profilevars_new
       integer :: imax_old, jmax_old, kmax_old
       integer :: imax_new, jmax_new, kmax_new
       integer :: nvars
       integer :: i, j, k, n
       real(8),dimension(:,:,:), allocatable :: znew1, znew2
       real(8),dimension(:,:,:,:), allocatable :: profilevars_new1
       real(8),dimension(:,:,:,:), allocatable :: profilevars_new2

       imax_old = size(xold)
       jmax_old = size(yold)
       kmax_old = size(zold)
       imax_new = size(xnew)
       jmax_new = size(ynew)
       kmax_new = size(znew)
       nvars = size(profilevars_old,dim=4)

       allocate(profilevars_new1(kmax_old,imax_new,jmax_old,nvars))

       if(imax_new.ne.imax_old) then
         do j = 1, jmax_old
           do k = 1, kmax_old
             do n = 1, nvars
               call SplineInterp(xold(:),profilevars_old(k,:,j,n), imax_old, &
                                 xnew(:),profilevars_new1(k,:,j,n),imax_new,0)
             enddo
           enddo
         enddo
       else ! imax_new = imax_old
         profilevars_new1(:,:,:,:) = profilevars_old(:,:,:,:)
       endif

       allocate( profilevars_new2(kmax_old,imax_new,jmax_new,nvars) )

       if(jmax_new.eq.1) then
         profilevars_new2(:,:,1,:) = profilevars_new1(:,:,1,:)
       else
         if(jmax_old.ne.jmax_new) then
           do i = 1, imax_new
             do k = 1, kmax_old
               do n = 1, nvars
                 call SplineInterp(yold(:),profilevars_new1(k,i,:,n), jmax_old, &
                                   ynew(:),profilevars_new2(k,i,:,n), jmax_new,2)
                 !call SplineInterp_periodic(yold(:),profilevars_new1(k,i,:,n), jmax_old, &
                 !                           ynew(:),profilevars_new2(k,i,:,n), jmax_new,2)
               enddo
             enddo
           enddo
         else
           profilevars_new2 = profilevars_new1
         endif
       endif
       deallocate( profilevars_new1)

       if(kmax_old.ne.kmax_new) then
         do n = 1, nvars
           do j = 1, jmax_new
             do i = 1, imax_new
               call SplineInterp(zold(:), profilevars_new2(:,i,j,n),kmax_old, &
                                 znew(:), profilevars_new(:,i,j,n) ,kmax_new,2)
               !call SplineInterp_periodic(zold(:), profilevars_new2(:,i,j,n),kmax_old, &
               !                  znew(:), profilevars_new(:,i,j,n) ,kmax_new,2)
             enddo
           enddo
         enddo
       else
         profilevars_new = profilevars_new2
       endif

       deallocate( profilevars_new2 )

     end subroutine InterpSol


   end program planeinlet
