
   ! convert timesereis volume data files to plane inlet files

   Program planeinlet
     use modRWHDF5
     use modTecbin
     use modHelmholtzDecomp
     implicit none

     real(8), parameter :: R=8314.3D0
     integer :: istencilsize = 4
     integer :: i, j, k, kk, n, jj, nn
     integer :: imax, jmax, kmax
     integer :: file_be, file_end, file_skip, num_file
     real(8), dimension(:,:,:,:,:), allocatable :: varstmp, varsnew
     real(8), dimension(:,:,:,:), allocatable :: vars_rd, vars_wt, buffer_grd, buffer_pgrd
     real(8), dimension(:,:,:,:), allocatable :: buffer_grd_tmp

     real(8), dimension(:), allocatable :: ttime, ttime_total
     integer, dimension(:,:), allocatable :: buffer_attr
     character(4) :: fnum4, fnum_iloc
     character(8) :: fnum8

     integer :: iplanes, nbuffer, nperfile
     character(400) :: fname, filepath, filepath_output, gridfile, planegrid

     real(8), dimension(5) :: varave  ! u, v, w, p, T
     integer :: nmax, iwindow, ntrans, nperiod, isubave
     real(8) :: pi, zlen, ave1, ave2, ave3, ave4
     complex, dimension(:), allocatable :: vars1d
     real(8) :: tmean, tmeanenergy, tmean2, tmeanenergy2

     type(tp_rdwt_hdf5) :: fsol, f1D, grd, pgrd, psol
     type(tp_DNSIndex)  :: fvol, pvol, tsplane

     character(200) :: varname_output, varname_plane
     real(8) :: xbeg, xend, ybeg, yend, zbeg, zend
     real(8), dimension(:,:), allocatable :: vars_ave
     integer :: count_file

     type(fprop), dimension(:), allocatable :: fileprop
     integer :: iinterp
     real(8) :: dx, dy, dz, Lx, Ly, Lz, rbar, Rm, didx, djdy, dkdz, xloc_output
     real(8), dimension(:,:,:), allocatable :: dudi, dudj, dudk, dvdi, dvdj, dvdk, dwdi, dwdj, dwdk
     real(8), dimension(:,:,:,:), allocatable :: buffer_var, buffer_sole, buffer_dila, buffer_rd, buffer_rd_tmp
     integer :: num_iplane, iloc_output

     varname_output = 'x y z u v w omx omy omz div' ! p T'
     pi = 4.d0*atan(1.)

     call InitHDF5()
     call InitGridHDF5(grd)
     call InitFlowHDF5(fsol)
     call InitFlow_PlaneInlet(psol)
     call Input()

     call ReadGridfiles()

     ! init 1D time data
     call InitFlow_1D(f1D)
     f1D%dimsf = (/nperfile/)
     allocate(ttime(nperfile))

     ! initialize 3D FFT
     call InitHDecomp3d(imax, jmax, kmax, iwindow, dx, dy, dz)

     write(30,*) '###############################################'
     write(30,*) '# of data points per segments in space domain kx = ', imax
     write(30,*) '# of segments in space domain = 1 '
     write(30,*) 'Interval of sampling (m) =', dx
     write(30,*) 'Length per segments (m) =', dble(imax)*dx
     write(30,*) '###############################################'
     write(30,*) '# of data points per segments in space domain ky = ', jmax
     write(30,*) '# of segments in space domain = 1 '
     write(30,*) 'Interval of sampling (m) =', dy
     write(30,*) 'Length per segments (m) =', dble(jmax)*dy
     write(30,*) '###############################################'
     write(30,*) '# of data points per segments in space domain kz = ', kmax
     write(30,*) '# of segments in space domain = 1 '
     write(30,*) 'Interval of sampling (m) =', dz
     write(30,*) 'Length per segments (m) =', dble(kmax)*dz

     fsol%dimsf = (/kmax,imax,jmax/)
     fsol%offset(1) = pvol%kbe -1; fsol%offset(2) = pvol%ibe -1; fsol%offset(3) = pvol%jbe -1;

     allocate(buffer_rd(kmax,imax,jmax,1),buffer_rd_tmp(kmax,imax,jmax,1))
     allocate(buffer_var(imax,jmax,kmax,5),buffer_sole(imax,jmax,kmax,3),buffer_dila(imax,jmax,kmax,3))

     open(7,file='planeindex.dat',status='unknown')
       fsol%dimsf  = grd%dimsf
       count_file = 0

     ! begin the main code
     do n=1, num_file, nperfile

       do jj=1, nperfile
         write(unit=fnum8,fmt='(I08.8)') (file_be + (jj+n-1-1)*file_skip)
         fsol%fname = trim(filepath)//'timeseriesVol_'//fnum8//'.h5'
         print *, 'Reading file: ', trim(fsol%fname)
         fsol%gname = 'vol'
         do nn=1, 5
           call ReadHDF5_3D_1V(fsol,buffer_rd(:,:,:,1),nn)
           do k=1, kmax
             buffer_var(1:imax,1:jmax,k,nn) = buffer_rd(k,1:imax,1:jmax,1)
           enddo
         enddo ! end nn loop
         fsol%gname = '/'
         call ReadHDF5_scalar(fsol,ttime(jj))
         print *, 'reading file info: time = ', ttime(jj)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!! debug file !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! writing Debug file
         if(n.eq.nperfile) then
           allocate(dudi(imax,jmax,kmax),dudj(imax,jmax,kmax),dudk(imax,jmax,kmax))
           allocate(dvdi(imax,jmax,kmax),dvdj(imax,jmax,kmax),dvdk(imax,jmax,kmax))
           allocate(dwdi(imax,jmax,kmax),dwdj(imax,jmax,kmax),dwdk(imax,jmax,kmax))
           didx = 1./dx
           djdy = 1./dy
           dkdz = 1./dz
           call Calddijk(imax,jmax,kmax,buffer_var(:,:,:,1:3),dudi,dudj,dudk,dvdi,dvdj,dvdk,dwdi,dwdj,dwdk)
           call InitTec(1,imax,jmax,kmax, 10,trim(varname_output),0)
           do k=1, kmax
             if(iinterp.eq.0) then
               vartmp(1:imax,1:jmax,k,1) = buffer_grd(k,1:imax,1:jmax,1)
               vartmp(1:imax,1:jmax,k,2) = buffer_grd(k,1:imax,1:jmax,2)
               vartmp(1:imax,1:jmax,k,3) = buffer_grd(k,1:imax,1:jmax,3)
             else
               vartmp(1:imax,1:jmax,k,1) = buffer_grd_tmp(k,1:imax,1:jmax,1)
               vartmp(1:imax,1:jmax,k,2) = buffer_grd_tmp(k,1:imax,1:jmax,2)
               vartmp(1:imax,1:jmax,k,3) = buffer_grd_tmp(k,1:imax,1:jmax,3)
             endif
             vartmp(1:imax,1:jmax,k,4) = buffer_var(1:imax,1:jmax,k,1)
             vartmp(1:imax,1:jmax,k,5) = buffer_var(1:imax,1:jmax,k,2)
             vartmp(1:imax,1:jmax,k,6) = buffer_var(1:imax,1:jmax,k,3)
             vartmp(1:imax,1:jmax,k,7) = dwdj(1:imax,1:jmax,k)*djdy - dvdk(1:imax,1:jmax,k)*dkdz
             vartmp(1:imax,1:jmax,k,8) = dudk(1:imax,1:jmax,k)*dkdz - dwdi(1:imax,1:jmax,k)*didx
             vartmp(1:imax,1:jmax,k,9) = dvdi(1:imax,1:jmax,k)*didx - dudj(1:imax,1:jmax,k)*djdy
             vartmp(1:imax,1:jmax,k,10)= dudi(1:imax,1:jmax,k)*didx + dvdj(1:imax,1:jmax,k)*djdy + dwdk(1:imax,1:jmax,k)*dkdz
             !vartmp(1:imax,1:jmax,k,11)= dwdk(1:imax,1:jmax,k)
           enddo
           fname = 'original_data.plt'
           print *, 'Writing file: ', trim(fname)
           call WriteTec(fname)
         endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         ! FFT filter
         do nn=1, 5
           call FFTFilter3d_2(buffer_var(:,:,:,nn),buffer_var(:,:,:,nn),1./Lx,1./(12.*dx),1./Ly,1./(12.*dy),1./Lz,1./(12.*dz))
         enddo

         ! Decomp
         call DoHDecomp3d(buffer_var(:,:,:,1:3),buffer_sole,buffer_dila)

         do nn=1, 5
           if(nn.le.3) then
             do k=1, kmax
               buffer_rd(k,1:imax,1:jmax,1) = buffer_dila(1:imax,1:jmax,k,nn)
             enddo
           else
             do k=1, kmax
               buffer_rd(k,1:imax,1:jmax,1) = buffer_var(1:imax,1:jmax,k,nn)
             enddo
           endif
           call InterpSol(buffer_grd(:,:,:,1), buffer_grd(:,:,:,2), buffer_grd(:,:,:,3), buffer_rd(:,:,:,1:1), &
                          buffer_pgrd(:,3:pgrd%dimsf(2),:,1), buffer_pgrd(:,3:pgrd%dimsf(2),:,2), buffer_pgrd(:,3:pgrd%dimsf(2),:,3), varsnew(:,3:pgrd%dimsf(2),:,nn:nn,jj) )
         enddo ! end nn loop

       enddo ! end jj loop

       varsnew(:,:,:,4,:) = varsnew(:,:,:,4,:)/(rbar*varsnew(:,:,:,5,:))

       ! rescale Uave
  !     varsnew(:,3:pgrd%dimsf(2),:,1,:) = varsnew(:,3:pgrd%dimsf(2),:,1,:) + varave(1)

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
         varname_plane = 'x, y, z, u, v, w, rho, T'
         call InitTec(1,int(pgrd%dimsf(2)-2),int(pgrd%dimsf(3)),int(pgrd%dimsf(1)), 8,trim(varname_plane),0)
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
         fname = 'planeinlet.plt'
         print *, 'Writing file: ', trim(fname)
         call WriteTec(fname)

         call Calddijk(imax,jmax,kmax,buffer_var(:,:,:,1:3),dudi,dudj,dudk,dvdi,dvdj,dvdk,dwdi,dwdj,dwdk)
         call InitTec(1,imax,jmax,kmax, 10,trim(varname_output),0)
         do k=1, kmax
           vartmp(1:imax,1:jmax,k,1) = buffer_grd(k,1:imax,1:jmax,1)
           vartmp(1:imax,1:jmax,k,2) = buffer_grd(k,1:imax,1:jmax,2)
           vartmp(1:imax,1:jmax,k,3) = buffer_grd(k,1:imax,1:jmax,3)
           vartmp(1:imax,1:jmax,k,4) = buffer_var(1:imax,1:jmax,k,1)
           vartmp(1:imax,1:jmax,k,5) = buffer_var(1:imax,1:jmax,k,2)
           vartmp(1:imax,1:jmax,k,6) = buffer_var(1:imax,1:jmax,k,3)
           vartmp(1:imax,1:jmax,k,7) = dwdj(1:imax,1:jmax,k)*djdy - dvdk(1:imax,1:jmax,k)*dkdz
           vartmp(1:imax,1:jmax,k,8) = dudk(1:imax,1:jmax,k)*dkdz - dwdi(1:imax,1:jmax,k)*didx
           vartmp(1:imax,1:jmax,k,9) = dvdi(1:imax,1:jmax,k)*didx - dudj(1:imax,1:jmax,k)*djdy
           vartmp(1:imax,1:jmax,k,10)= dudi(1:imax,1:jmax,k)*didx + dvdj(1:imax,1:jmax,k)*djdy + dwdk(1:imax,1:jmax,k)*dkdz
         enddo
         fname = 'original_data_FFT.plt'
         print *, 'Writing file: ', trim(fname)
         call WriteTec(fname)

         call Calddijk(imax,jmax,kmax,buffer_dila,dudi,dudj,dudk,dvdi,dvdj,dvdk,dwdi,dwdj,dwdk)
         do k=1, kmax
           vartmp(1:imax,1:jmax,k,1) = buffer_grd(k,1:imax,1:jmax,1)
           vartmp(1:imax,1:jmax,k,2) = buffer_grd(k,1:imax,1:jmax,2)
           vartmp(1:imax,1:jmax,k,3) = buffer_grd(k,1:imax,1:jmax,3)
           vartmp(1:imax,1:jmax,k,4) = buffer_dila(1:imax,1:jmax,k,1)
           vartmp(1:imax,1:jmax,k,5) = buffer_dila(1:imax,1:jmax,k,2)
           vartmp(1:imax,1:jmax,k,6) = buffer_dila(1:imax,1:jmax,k,3)
           vartmp(1:imax,1:jmax,k,7) = dwdj(1:imax,1:jmax,k)*djdy - dvdk(1:imax,1:jmax,k)*dkdz
           vartmp(1:imax,1:jmax,k,8) = dudk(1:imax,1:jmax,k)*dkdz - dwdi(1:imax,1:jmax,k)*didx
           vartmp(1:imax,1:jmax,k,9) = dvdi(1:imax,1:jmax,k)*didx - dudj(1:imax,1:jmax,k)*djdy
           vartmp(1:imax,1:jmax,k,10)= dudi(1:imax,1:jmax,k)*didx + dvdj(1:imax,1:jmax,k)*djdy + dwdk(1:imax,1:jmax,k)*dkdz
         enddo
         fname = 'new_data_dila.plt'
         print *, 'Writing file: ', trim(fname)
         call WriteTec(fname)

         call Calddijk(imax,jmax,kmax,buffer_sole,dudi,dudj,dudk,dvdi,dvdj,dvdk,dwdi,dwdj,dwdk)
         do k=1, kmax
           vartmp(1:imax,1:jmax,k,1) = buffer_grd(k,1:imax,1:jmax,1)
           vartmp(1:imax,1:jmax,k,2) = buffer_grd(k,1:imax,1:jmax,2)
           vartmp(1:imax,1:jmax,k,3) = buffer_grd(k,1:imax,1:jmax,3)
           vartmp(1:imax,1:jmax,k,4) = buffer_sole(1:imax,1:jmax,k,1)
           vartmp(1:imax,1:jmax,k,5) = buffer_sole(1:imax,1:jmax,k,2)
           vartmp(1:imax,1:jmax,k,6) = buffer_sole(1:imax,1:jmax,k,3)
           vartmp(1:imax,1:jmax,k,7) = dwdj(1:imax,1:jmax,k)*djdy - dvdk(1:imax,1:jmax,k)*dkdz
           vartmp(1:imax,1:jmax,k,8) = dudk(1:imax,1:jmax,k)*dkdz - dwdi(1:imax,1:jmax,k)*didx
           vartmp(1:imax,1:jmax,k,9) = dvdi(1:imax,1:jmax,k)*didx - dudj(1:imax,1:jmax,k)*djdy
           vartmp(1:imax,1:jmax,k,10)= dudi(1:imax,1:jmax,k)*didx + dvdj(1:imax,1:jmax,k)*djdy + dwdk(1:imax,1:jmax,k)*dkdz
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
       read(*,*) iwindow
       read(*,*)
       read(*,*) num_iplane, xloc_output
       !allocate(iplane(num_iplane))
       !read(*,*) (iplane(i),i=1,num_iplane)

       rbar = R/Rm
       print *, 'rbar = ', rbar

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
       allocate(buffer_grd_tmp(kmax,imax,jmax,3))
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

       dx2 = buffer_grd(1,imax,1,1) - buffer_grd(1,imax-1,1,1)
       dy2 = buffer_grd(1,1,jmax,2) - buffer_grd(1,1,jmax-1,2)
       dz2 = buffer_grd(kmax,1,1,3) - buffer_grd(kmax-1,1,1,3)

       Lx = buffer_grd(1,imax,1,1) - buffer_grd(1,1,1,1) + dx1
       Ly = buffer_grd(1,1,jmax,2) - buffer_grd(1,1,1,2) + dy1
       Lz = buffer_grd(kmax,1,1,3) - buffer_grd(1,1,1,3) + dz1

       iinterp = 0
       if(abs(dx2-dx1).gt.1.e-10) then
         print *, 'Interp in x direction using imax'
         do i=1, imax
           buffer_grd_tmp(:,i,:,1) = dble(i-1)*Lx/dble(imax)
         enddo
         iinterp = 1
       endif
       if(abs(dy2-dy1).gt.1.e-10) then
         print *, 'Interp in y direction using jmax'
         do j=1, jmax
           buffer_grd_tmp(:,:,j,2) = dble(j-1)*Ly/dble(jmax)
         enddo
         iinterp = 1
       endif
       if(abs(dz2-dz1).gt.1.e-10) then
         print *, 'Interp in z direction using kmax'
         do k=1, kmax
           buffer_grd_tmp(k,:,:,3) = dble(k-1)*Lz/dble(kmax)
         enddo
         iinterp = 1
       endif
       if(iinterp.eq.1) then
         dx = buffer_grd_tmp(1,2,1,1) - buffer_grd_tmp(1,1,1,1)
         dy = buffer_grd_tmp(1,1,2,2) - buffer_grd_tmp(1,1,1,2)
         dz = buffer_grd_tmp(2,1,1,3) - buffer_grd_tmp(1,1,1,3)
         Lx = buffer_grd_tmp(1,imax,1,1) - buffer_grd_tmp(1,1,1,1)
         Ly = buffer_grd_tmp(1,1,jmax,2) - buffer_grd_tmp(1,1,1,2)
         Lz = buffer_grd_tmp(kmax,1,1,3) - buffer_grd_tmp(1,1,1,3)
       else
         dx = dx1
         dy = dy1
         dz = dz1
       endif
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

       ! for the plane inlet grid file. make sure x=0, y=0 and z=0 at the first grid point.
       buffer_pgrd(:,:,:,1) = buffer_pgrd(:,:,:,1) - buffer_pgrd(1,1,1,1)
       buffer_pgrd(:,:,:,2) = buffer_pgrd(:,:,:,2) - buffer_pgrd(1,1,1,2)
       buffer_pgrd(:,:,:,3) = buffer_pgrd(:,:,:,3) - buffer_pgrd(1,1,1,3)

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
       real(8),dimension(:,:,:), intent(in) :: xold, yold, zold
       real(8),dimension(:,:,:), intent(in) :: xnew, ynew, znew
       real(8),dimension(:,:,:,:), intent(in) :: profilevars_old
       real(8),dimension(:,:,:,:), intent(out) :: profilevars_new
       integer :: imax_old, jmax_old, kmax_old
       integer :: imax_new, jmax_new, kmax_new
       integer :: nvars
       integer :: i, j, k, n
       real(8),dimension(:,:,:), allocatable :: znew1, znew2
       real(8),dimension(:,:,:,:), allocatable :: profilevars_new1
       real(8),dimension(:,:,:,:), allocatable :: profilevars_new2

       imax_old = size(xold,dim=2)
       jmax_old = size(xold,dim=3)
       kmax_old = size(xold,dim=1)
       imax_new = size(xnew,dim=2)
       jmax_new = size(xnew,dim=3)
       kmax_new = size(xnew,dim=1)
       nvars = size(profilevars_old,dim=4)

       allocate(profilevars_new1(kmax_old,imax_new,jmax_old,nvars))

       if(imax_new.ne.imax_old) then
         do j = 1, jmax_old
           do k = 1, kmax_old
             do n = 1, nvars
               call SplineInterp(xold(1,:,1),profilevars_old(k,:,j,n), imax_old, &
                                 xnew(1,:,1),profilevars_new1(k,:,j,n),imax_new,0)
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
                 call SplineInterp(yold(1,1,:),profilevars_new1(k,i,:,n), jmax_old, &
                                   ynew(1,1,:),profilevars_new2(k,i,:,n), jmax_new,2)
                 !call SplineInterp_periodic(yold(1,1,:),profilevars_new1(k,i,:,n), jmax_old, &
                 !                  ynew(1,1,:),profilevars_new2(k,i,:,n), jmax_new,2)
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
               call SplineInterp(zold(:,i,j), profilevars_new2(:,i,j,n),kmax_old, &
                                 znew(:,i,j), profilevars_new(:,i,j,n) ,kmax_new,2)
               !call SplineInterp_periodic(zold(:,i,j), profilevars_new2(:,i,j,n),kmax_old, &
               !                  znew(:,i,j), profilevars_new(:,i,j,n) ,kmax_new,2)
             enddo
           enddo
         enddo
       else
         profilevars_new = profilevars_new2
       endif

       deallocate( profilevars_new2 )

     end subroutine InterpSol


   end program planeinlet
