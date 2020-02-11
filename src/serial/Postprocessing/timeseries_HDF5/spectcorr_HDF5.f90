! program to compute time series spectrum

!  spectrum is computed if ispect=1,
!  nperiod is the section length which should be less than the toltal length of the signal
!  ioverlap=1 means compute spectrum using overlapped sections, ie the actually section will be 2*nperiod
!
! correlation is always computed
! icorrspace & jcorrspace are the spacing for the time correlation,
! icorrave & jcorrave are the average distance, ie skip how many points before doing averaging in space
! ntcorrmax is the largest time interval for the correlation
! iwindow: 0 flat, 1 hann, 2 modified hann
program spectcorr
  use MTSKplane
  use MTSJplane
  use MTSIplane
  use MTSVOL
  use modSpect
  use modSpecty
  implicit none

  integer :: ical_Taylor = 0, ilocal_ave, icalPDF = 0
  integer :: icalPSD, icalcorrxt, icalcorrxy, icalphasespeed, icalWNSD_ky, icalWNSD_kx, icalky_t, icalkx_t, ical3D, ical_skewness
  integer :: icalWNSD_kz, icalkx_ky, icalcorrxz, icalcorrzt
  integer :: nskip, ixwindowl, ixwindowr, jywindow, kzwindowd, kzwindowu
  integer :: nperiod, ioverlap, ntrans, ntrans_kx, nsection, nsection_kx
  integer :: iwindow, ntpoint_total

  !> timeseries planes information
  integer :: ireadkp, ireadip, ireadjp, ireadTSVol
  integer :: num_kplane, num_iplane, num_jplane
  integer, dimension(:), allocatable :: kplane, iplane, jplane
  integer :: ibe_k, iend_k, jbe_k, jend_k, ibe_j, iend_j, kbe_j, kend_j, jbe_i, jend_i, kbe_i, kend_i
  integer :: ibe_Vol, iend_Vol, jbe_Vol, jend_Vol, kbe_Vol, kend_Vol
  integer :: kplane_be, kplane_end, kplane_skip, iplane_be, iplane_end, iplane_skip, jplane_be, jplane_end, jplane_skip
  type(tp_DNSIndex) :: DNSIndex_k, DNSIndex_i, DNSIndex_j, DNSIndex_Vol

  !> temporary variables
  real(8), dimension(:), allocatable :: tsdata1d, tsdata1dx, tsdata1d_tmp, pdf
  real(8), dimension(:,:), allocatable :: tsdata2d, tsdata2d_tmp
  real(8), dimension(:,:,:), allocatable :: tsdata3d
  complex, dimension(:), allocatable :: spect1dtmp, spect1dave, spect1d1, spect1d2
  complex, dimension(:), allocatable :: spect1dxtmp, spect1dxave
  complex, dimension(:,:), allocatable :: spect1dtmp2, spect1dave2
  complex, dimension(:,:), allocatable :: spect1dtmp3, spect1dave3
  complex, dimension(:,:), allocatable :: spect2dtmp, spect2dave
  complex, dimension(:,:,:), allocatable :: spect3dtmp, spect3dave
  complex, dimension(:,:), allocatable :: spectspanave, spectxspanave
  real(8) :: t1, t2, t3, t4, rms_tmp, skewness_tmp, flatness_tmp, moment(4)
!  real(8) :: t1, t2, t3, t4, rms_tmp, skewness_tmp, flatness_tmp, rms, skewness, flatness
!  integer :: subsample, ntpoint_tmp
  !real(8), dimension(:), allocatable :: powerspect

  !> timeseries files information
  integer :: npath
  type(fprop), dimension(:), allocatable :: fileprop
  character(400) :: fileinput, fileoutput, filename
  character(4) :: fnum1, fnum2, fnum3, fnum4, fnum5
  character(8) :: fnum_8
  integer :: nvar_output
  integer, dimension(:), allocatable :: varindex_output

  integer :: ii, jj, i, j, k, n, nn, nnn, m, numpt1, numpt2, kk
  real(8) :: dt_sample, dx_kp, dy_kp, dy_ip, dz_ip, dx_jp, dz_jp, dx_Vol, dy_Vol, dz_Vol
  real(8) :: tmean, tmeanenergy, tmpenergy, ff, fperiod
  integer, dimension(:), allocatable :: numpt
  real(8), dimension(:), allocatable :: phasespeed, R2, phase
  real(8), dimension(:,:,:,:,:), allocatable :: buffer_kplane, buffer_iplane, buffer_jplane, buffer_Vol, buffer_p0

  !> the following variables are for calculating skewness and flatness
  real(8) :: pave, p2, p3, p4, prms, prms_tmp, Sp, Fp, Sp_tmp, Fp_tmp, iradius, iradius_tmp, std, dx_tmp
  real(8) :: Sp_max, Sp_min, Fp_max, Fp_min
  real(8), dimension(:), allocatable :: fl, fu
  integer :: num_freq, ntpoint_tmp1, ntpoint_tmp2, ixw, jyw

  ! debug file
  character(80) :: DebugFileName

  DebugFileName = 'timeseries_hdf5.debug'
  call InitHDF5()
  call Input()

  print *, 'ntpoint_total =', ntpoint_total, 'dt_sample =', dt_sample

  !-----------------------------------------
  ! write DebugFileName
  open(66,file=trim(DebugFileName),position='append')
  write(66,*) 'ntpoint_total =', ntpoint_total, 'dt_sample =', dt_sample
  !-----------------------------------------
  if(icalPSD.eq.1.or.icalcorrxt.eq.1.or.icalPhaseSpeed.eq.1.or.icalcorrzt.eq.1) then
      ! call InitSpect1d(ntpoint_total, nperiod, iwindow, ioverlap, dt_sample, ntrans, nsection)
      call InitSpect1d(ntpoint_total, nsection, iwindow, ioverlap, dt_sample, ntrans, nperiod)
      write(66,*) '# of data points per segments =', ntrans
      write(66,*) '# of segments =', nsection
      write(66,*) 'Sampling frequency (Herz) =', 1./dt_sample
      write(66,*) 'Interval of sampling (sec) =', dt_sample
      write(66,*) 'Length per segments (sec) =', (ntrans-1)*dt_sample
      if(ioverlap.eq.0) write(66,*) 'Segments with no overlap'
      if(ioverlap.eq.1) write(66,*) 'Segments with half overlap'
      if(iwindow.eq.0) write(66,*) 'FFT - window type: flat top window'
      if(iwindow.eq.1) write(66,*) 'FFT - window type: hanning window'
      if(iwindow.eq.2) write(66,*) 'FFT - window type: modified hanning window'
      if(iwindow.eq.3) write(66,*) 'FFT - window type: hamming window'
      allocate(tsdata1d(1:ntpoint_total))
      tsdata1d = 0.d0
  endif



  !!#################################################################################
  if(ireadTSVol.eq.1) then  ! Using 3D timeseries volume dataset
    call InitTSVol(ntpoint_total,npath,DNSIndex_Vol, ibe_Vol, iend_Vol, jbe_Vol, jend_Vol, kbe_Vol, kend_Vol, &
                   ixwindowl, ixwindowr, kzwindowd, kzwindowu, nvar_output, varindex_output, fileprop(1))
    allocate(buffer_Vol(ntpoint_total,kmints_Vol:kmaxts_Vol,imints_Vol:imaxts_Vol,nypoint_Vol,nvar_output))

    write(unit=fnum1,fmt='(I04.4)') ibe_Vol
    write(unit=fnum2,fmt='(I04.4)') iend_Vol
    write(unit=fnum3,fmt='(I04.4)') kbe_Vol
    write(unit=fnum4,fmt='(I04.4)') kend_Vol
    write(unit=fnum5,fmt='(I04.4)') nsection

    call ReadTSVol(buffer_Vol)

    !print *, 'buffer_Vol(1,:,1,1,1) = ', buffer_Vol(1,:,1,1,1)
    if(icalPSD.eq.1) then
      call InitFFT1D(ntrans)
      allocate(spect1dtmp(ntrans), spect1dave(ntrans))
      do m=1, nvar_output
        spect1dave = 0.d0;  numpt1 = 0
        tmean = 0.d0; tmeanenergy = 0.d0
        fileoutput = '_Vol_i'//fnum1//'-'//fnum2//'_k'//fnum3//'-'//fnum4//'_'//trim(varout_Vol(m))//'_nsection'//fnum5
        print *, 'fileoutput =', trim(fileoutput)
        do k=1, nzpoint_Vol
          do i=1, nxpoint_Vol
            do j=1, nypoint_Vol
              tsdata1d(1:ntpoint_total) = buffer_Vol(1:ntpoint_total,k,i,j,m)
              spect1dtmp = autospect1d(tsdata1d(1:ntpoint_total))
              spect1dave = spect1dave + spect1dtmp
              tmean = tmean + sum(tsdata1d(1:ntpoint_total))
              tmeanenergy = tmeanenergy + sum(tsdata1d(1:ntpoint_total)**2)
              numpt1 = numpt1 + 1
            enddo
          enddo
        enddo
        tmean = tmean/dble(ntpoint_total*numpt1)
        tmeanenergy = tmeanenergy/dble(ntpoint_total*numpt1)
        tmeanenergy = tmeanenergy - tmean**2
        spect1dave = spect1dave/dble(numpt1)
        ! print *, 'spect1dave =', real(sum(spect1dave))
        print *, 'tmean =', tmean, 'tmeanenergy =', tmeanenergy
        call CalPowerspect(spect1dave,tmeanenergy,fileoutput)
        call CalCorrt(spect1dave,fileoutput)
      enddo ! end m loop
      deallocate(spect1dtmp, spect1dave,tsdata1d)
    endif ! end icalPSD.eq.1

    if(icalcorrxy) then
      call InitFFT1D(nypoint_Vol)
      call InitCorry(nypoint_Vol, dy_Vol)

      allocate(spect1dtmp3(-ixwindowl:ixwindowr,nypoint_Vol), spect1dave3(-ixwindowl:ixwindowr,nypoint_Vol))
      if(.not.allocated(numpt)) allocate(numpt(-ixwindowl:ixwindowr))
      do m=1, nvar_output
        fileoutput = '_Vol_i'//fnum1//'-'//fnum2//'_k'//fnum3//'-'//fnum4//'_'//trim(varout_Vol(m))
        print *, 'fileoutput =', trim(fileoutput)
        numpt = 0
        spect1dtmp3 = (0.,0.); spect1dave3 = (0.,0.)
        do n=1, ntpoint_total
          do k=1, nzpoint_Vol
            do i=1, nxpoint_Vol
              do ii = -min(ixwindowl,i-imints_Vol), min(ixwindowr,imaxts_Vol-i)
                spect1dtmp3(ii,:) = crossspect1dy(buffer_Vol(n,k,i+ii,1:nypoint_Vol,m),buffer_Vol(n,k,i,1:nypoint_Vol,m))
                spect1dave3(ii,:) = spect1dave3(ii,:) + spect1dtmp3(ii,:)
                numpt(ii) = numpt(ii) + 1
              enddo
            enddo
          enddo
        enddo
        forall(ii=-ixwindowl:ixwindowr,numpt(ii).gt.0)
          spect1dave3(ii,:) = spect1dave3(ii,:)/dble(numpt(ii))
        end forall
        call CalCorrxy(ixwindowl,ixwindowr,dx_Vol,spect1dave3,fileoutput)
      enddo ! end m loop
      deallocate(spect1dtmp3, spect1dave3)
    endif ! end icalcorrxy

    if(icalcorrxt.eq.1) then
         call InitFFT1D(ntrans)
         allocate(spect1dtmp2(ntrans,-ixwindowl:ixwindowr), spect1dave2(ntrans,-ixwindowl:ixwindowr))
         if(.not.allocated(numpt)) allocate(numpt(-ixwindowl:ixwindowr))
         do m = 1, nvar_output
            fileoutput = '_Vol_i'//fnum1//'_i'//fnum2//'_k'//fnum3//'-'//fnum4//'_'//trim(varout_Vol(m))//'_nsection'//fnum5
            numpt = 0
            spect1dtmp2 = (0.,0.); spect1dave2 = (0.,0.)
            do k=1, nzpoint_Vol
              do i = 1, nxpoint_Vol
                do j = 1, nypoint_Vol
                  do ii = -min(ixwindowl,i-imints_Vol), min(ixwindowr,imaxts_Vol-i)
                    spect1dtmp2(:,ii) = crossspect1d(buffer_Vol(1:ntpoint_total,k,i+ii,j,m),buffer_Vol(1:ntpoint_total,k,i,j,m))
                    spect1dave2(:,ii) = spect1dave2(:,ii) + spect1dtmp2(:,ii)
                    numpt(ii) = numpt(ii) + 1
                  enddo
                enddo
              enddo
            enddo
            forall(ii=-ixwindowl:ixwindowr,numpt(ii).gt.0)
               spect1dave2(:,ii) = spect1dave2(:,ii)/dble(numpt(ii))
            end forall
            call CalCorrxt(ixwindowl,ixwindowr,dx_Vol,spect1dave2,fileoutput)
         enddo
         deallocate(spect1dtmp2, spect1dave2)


    endif ! icalcorrxt.eq.1

    if(icalcorrzt) then
      call InitFFT1D(ntrans)
      allocate(spect1dtmp2(ntrans,-kzwindowd:kzwindowu), spect1dave2(ntrans,-kzwindowd:kzwindowu))
      if(.not.allocated(numpt)) allocate(numpt(-kzwindowd:kzwindowu))
      do m=1, nvar_output
        fileoutput = '_Vol_i'//fnum1//'-'//fnum2//'_k'//fnum3//'-'//fnum4//'_'//trim(varout_Vol(m))//'_nsection'//fnum5
        print *, 'fileoutput =', trim(fileoutput)
        numpt = 0
        spect1dtmp2 = (0.,0.); spect1dave2 = (0.,0.)
        do k=1, nzpoint_Vol
          do i=1, nxpoint_Vol
            do j=1, nypoint_Vol
              do kk = -min(kzwindowd,k-kmints_Vol), min(kzwindowu,kmaxts_Vol-k)
                spect1dtmp2(:,kk) = crossspect1d(buffer_Vol(1:ntpoint_total,k+kk,i,j,m),buffer_Vol(1:ntpoint_total,k,i,j,m))
                spect1dave2(:,kk) = spect1dave2(:,kk) + spect1dtmp2(:,kk)
                numpt(kk) = numpt(kk) + 1
              enddo
            enddo
          enddo
        enddo
        forall(kk=-kzwindowd:kzwindowu,numpt(kk).gt.0)
          spect1dave2(:,kk) = spect1dave2(:,kk)/dble(numpt(kk))
         end forall
         call CalCorrzt(kzwindowd,kzwindowu,dz_Vol,spect1dave2,fileoutput)
      enddo ! end m loop
      deallocate(spect1dtmp2, spect1dave2)
    endif ! end icalcorrzt

  endif ! end ireadTSVol.eq.1


  !!#################################################################################
  if(ireadip.eq.1.and.num_iplane.gt.0) then ! Using constant i-plane time series

    if(iFileType.eq.0) then
       call InitTSiplane(ntpoint_total, npath, DNSIndex_i, jbe_i, jend_i, kbe_i, kend_i, &
                              nvar_output, varindex_output, fileprop)
    else
      call InitVOLTSIplane(ntpoint_total, npath, DNSIndex_i, jbe_i, jend_i, kbe_i, kend_i, &
                        ixwindowl, ixwindowr, nvar_output, varindex_output, fileprop)
    endif

    allocate( buffer_iplane(ntpoint_total, nypoint_i, nzpoint_i, 1, nvar_output) )

    do i=1, num_iplane
      write(unit=fnum1,fmt='(I04.4)') iplane(i)
      write(unit=fnum2,fmt='(I04.4)') kbe_i
      write(unit=fnum3,fmt='(I04.4)') kend_i
      write(unit=fnum4,fmt='(I04.4)') nsection
      if(iFileType.eq.0) then
        call ReadTSiplane(iplane(i),buffer_iplane)
      else
        do jj=1, DNSIndex_i%niplane
          if(iplane(i).eq.DNSIndex_i%iplane(jj)) then
            iplane_be = jj
          endif
        enddo
        call ReadVOLTSIplane_files(iplane_be,iplane(i),buffer_iplane)
      endif
      print *, 'streamwise location i =', iplane(i)
      print *, 'wall normal average range k =', kbe_i, kend_i

!print *, 'buffer = ', buffer_iplane(:,1,1,1,1)

      if(icalPSD.eq.1) then
        call InitFFT1D(ntrans)
        allocate(spect1dtmp(ntrans),spect1dave(ntrans))
        allocate(spect3dtmp(ntrans,nypoint_i,nzpoint_i), spect3dave(ntrans,nypoint_i,nzpoint_i))
        allocate(phase(ntrans))
        do m = 1, nvar_output
              spect1dave = 0.d0;  numpt1 = 0
              tmean = 0.d0; tmeanenergy = 0.d0
              fileoutput = '_i'//fnum1//'_k'//fnum2//'-'//fnum3//'_'//trim(varout_i(m))//'_nsection'//fnum4
              print *, 'fileoutput =', trim(fileoutput)
              do k = 1, nzpoint_i
              do j = 1, nypoint_i
                 tsdata1d(1:ntpoint_total) = buffer_iplane(1:ntpoint_total,j,k,1,m)
                 spect1dtmp = autospect1d(tsdata1d(1:ntpoint_total))
                 spect1dave = spect1dave + spect1dtmp
                 tmean = tmean + sum(tsdata1d(1:ntpoint_total))
                 tmeanenergy = tmeanenergy + sum(tsdata1d(1:ntpoint_total)**2)
                 numpt1 = numpt1 + 1
              enddo
              enddo
              tmean = tmean/dble(ntpoint_total*numpt1)
              tmeanenergy = tmeanenergy/dble(ntpoint_total*numpt1)
              tmeanenergy = tmeanenergy - tmean**2
              spect1dave = spect1dave/dble(numpt1)
              ! print *, 'spect1dave =', real(sum(spect1dave))
              print *, 'tmean =', tmean, 'tmeanenergy =', tmeanenergy
              call CalPowerspect(spect1dave,tmeanenergy,fileoutput)
              !call CalCorrt(spect1dave,fileoutput)

              ! calculate the phase
              do k=1, nzpoint_i
                do j=1, nypoint_i
                  tsdata1d(1:ntpoint_total) = buffer_iplane(1:ntpoint_total,j,k,1,m)
                  spect3dtmp(:,j,k) = sumspect1d(tsdata1d(1:ntpoint_total))
                enddo
              enddo
              spect3dave(1,:,:) = spect3dtmp(1,:,:)
              spect3dave(2:ntrans/2,:,:) = spect3dtmp(2:ntrans/2,:,:)*sqrt(2.d0)
              print *, 'fileoutput =', 'phase'//trim(fileoutput)
              open(7,file='phase'//trim(fileoutput)//'.dat',status='unknown')
                rewind(7)
                write(7,'(a)') 'variables=f,spectrum_r,spectrum_i,phase'
                write(7,'("Zone T=PSD, I=", I8.8, ",J=", I8.8, ",K=", I8.8, ", AUXDATA rms2=""", E15.8,"""")') ntrans/2, nypoint_i, nzpoint_i, tmeanenergy
                write(7,'("AUXDATA TotalTime=""", E15.8,"""")') dble(ntrans-1)*dt_sample
                write(7,'(a,I8.8,a,I3.3)') '# Data points per seg = ', ntrans, '   Num of seg = ', nsection
                write(7,'(a,I1.1,a,I1.1,a,E15.8)') '# Window type = ', iwindow, '   ioverlap = ', ioverlap, '   dt_sample = ', dt_sample
                write(7,'(a,E15.8)') '# meanenergy (mean square amplitude) = ', tmeanenergy

                do k=1, nzpoint_i
                  do j=1, nypoint_i
                    do ii=1, ntrans/2
                      phase(ii) = atan2(imag(spect3dtmp(ii,j,k)),real(spect3dtmp(ii,j,k)))
                      if(phase(ii).lt.0) phase(ii) = phase(ii) + 8.d0*atan(1.d0)
                    write(7,'(4E20.12)') real(ii-1)/(dble(ntrans-1)*dt_sample), real(spect3dave(ii,j,k)), imag(spect3dave(ii,j,k)),phase(ii)
                    enddo
                  enddo
                enddo
              close(7)

        enddo ! end m loop
        deallocate(spect1dtmp,spect1dave)
        deallocate(spect3dtmp,spect3dave)
        deallocate(phase)
      endif ! end icalPSD.eq.1

      ! Calculating spanwise wave-number spectrum
      if(icalWNSD_ky.eq.1) then

print *, 'dy_ip = ', dy_ip

        call InitSpect1dy(nypoint_i, nsection, iwindow, ioverlap, dy_ip, ntrans, nperiod)
        allocate(tsdata1d(1:nypoint_i))
        tsdata1d = 0.d0

        allocate(spect1dtmp(ntrans), spect1dave(ntrans))
        do m=1, nvar_output
          spect1dave = 0.d0;  numpt1 = 0
          tmean = 0.d0; tmeanenergy = 0.d0
          fileoutput = '_ky_i'//fnum1//'_k'//fnum2//'-'//fnum3//'_'//trim(varout_i(m))
          print *, 'fileoutput =', trim(fileoutput)
          do k=1, nzpoint_i
            do n=1, ntpoint_total
              tsdata1d(1:nypoint_i) = buffer_iplane(n,1:nypoint_i,k,1,m)
              spect1dtmp = autospect1dty(tsdata1d(1:nypoint_i))
              spect1dave = spect1dave + spect1dtmp
              tmean = tmean + sum(tsdata1d(1:nypoint_i))
              tmeanenergy = tmeanenergy + sum(tsdata1d(1:nypoint_i)**2)
              numpt1 = numpt1 + 1
            enddo ! end n loop
          enddo ! end k loop
          tmean = tmean/dble(nypoint_i*numpt1)
          tmeanenergy = tmeanenergy/dble(nypoint_i*numpt1)
          tmeanenergy = tmeanenergy - tmean**2
          spect1dave = spect1dave/dble(numpt1)
          print *, 'tmean =', tmean, 'tmeanenergy =', tmeanenergy
          call CalPowerspecty(spect1dave,tmeanenergy,fileoutput)
        enddo ! end m loop
        deallocate(spect1dtmp, spect1dave,tsdata1d)

      endif ! end calculating spanwise wave-number spectrum

      ! Calculating wall-normal wave-number spectrum
      if(icalWNSD_kz.eq.1) then

        call InitSpect1dy(nzpoint_i, nsection, iwindow, ioverlap, dz_ip, ntrans, nperiod)
        allocate(tsdata1d(1:nzpoint_i))
        tsdata1d = 0.d0

        allocate(spect1dtmp(ntrans), spect1dave(ntrans))
        do m=1, nvar_output
          spect1dave = 0.d0;  numpt1 = 0
          tmean = 0.d0; tmeanenergy = 0.d0
          fileoutput = '_kz_k'//fnum1//'_k'//fnum2//'-'//fnum3//'_'//trim(varout_i(m))
          print *, 'fileoutput =', trim(fileoutput)
          do j=1, nypoint_i
            do n=1, ntpoint_total
              tsdata1d(1:nzpoint_i) = buffer_iplane(n,j,1:nzpoint_i,1,m)
              spect1dtmp = autospect1dty(tsdata1d(1:nzpoint_i))
              spect1dave = spect1dave + spect1dtmp
              tmean = tmean + sum(tsdata1d(1:nzpoint_i))
              tmeanenergy = tmeanenergy + sum(tsdata1d(1:nzpoint_i)**2)
              numpt1 = numpt1 + 1
            enddo ! end n loop
          enddo ! end j loop
          tmean = tmean/dble(nzpoint_i*numpt1)
          tmeanenergy = tmeanenergy/dble(nzpoint_i*numpt1)
          tmeanenergy = tmeanenergy - tmean**2
          spect1dave = spect1dave/dble(numpt1)
          print *, 'tmean =', tmean, 'tmeanenergy =', tmeanenergy
          call CalPowerspecty(spect1dave,tmeanenergy,fileoutput)
        enddo ! end m loop
        deallocate(spect1dtmp, spect1dave,tsdata1d)

      endif ! end calculating wall-normal wave-number spectrum


     enddo ! end i loop


  endif ! end ireadip



  !!#################################################################################
  if(ireadjp.eq.1.and.num_jplane.gt.0) then ! Using constant jplane time series
    if(iFileType.eq.0) then
      print *, 'iFileType should be set to 1 '
    else
      call InitVOLTSJplane(ntpoint_total, npath, DNSIndex_j, ibe_j, iend_j, kbe_j, kend_j, &
                        ixwindowl, ixwindowr, kzwindowd, kzwindowu, nvar_output, varindex_output, fileprop)
    endif

    allocate( buffer_jplane(ntpoint_total, imints_j:imaxts_j, kmints_j:kmaxts_j, 1, nvar_output) )
    if(nvar_output.eq.5) then
      allocate( buffer_p0(ntpoint_total, imints_j:imaxts_j, kmints_j:kmaxts_j, 1, 1) )
    endif
    do j=1, num_jplane
      write(unit=fnum1,fmt='(I04.4)') jplane(j)
      write(unit=fnum2,fmt='(I04.4)') ibe_j
      write(unit=fnum3,fmt='(I04.4)') iend_j
      write(unit=fnum4,fmt='(I04.4)') nsection

      do jj=1, DNSIndex_j%njplane
        if(jplane(j).eq.DNSIndex_j%jplane(jj)) then
          jplane_be = jj
        endif
      enddo
      call ReadVOLTSJplane_files(jplane_be,jplane(j),buffer_jplane)
      print *, 'spanwise location j =', jplane(j)
      print *, 'streamwise average range i =', ibe_j, iend_j

!print *, '(1): ', buffer_jplane(:,1,1,1,1)
!print *, '(-128):', buffer_jplane(:,1,-128,1,1)


      ! Calculating PSD
      if(icalPSD.eq.1) then
        call InitFFT1D(ntrans)
        allocate(spect1dtmp(ntrans), spect1dave(ntrans))
        do m = 1, nvar_output
          spect1dave = 0.d0;  numpt1 = 0
          tmean = 0.d0; tmeanenergy = 0.d0
          fileoutput = '_j'//fnum1//'_i'//fnum2//'-'//fnum3//'_'//trim(varout_j(m))//'_nsection'//fnum4
          print *, 'fileoutput =', trim(fileoutput)
          do k = 1, nzpoint_j
            do i = 1, nxpoint_j
              tsdata1d(1:ntpoint_total) = buffer_jplane(1:ntpoint_total,i,k,1,m)
              spect1dtmp = autospect1d(tsdata1d(1:ntpoint_total))
              spect1dave = spect1dave + spect1dtmp
              tmean = tmean + sum(tsdata1d(1:ntpoint_total))
              tmeanenergy = tmeanenergy + sum(tsdata1d(1:ntpoint_total)**2)
              numpt1 = numpt1 + 1
            enddo
          enddo
          tmean = tmean/dble(ntpoint_total*numpt1)
          tmeanenergy = tmeanenergy/dble(ntpoint_total*numpt1)
          tmeanenergy = tmeanenergy - tmean**2
          spect1dave = spect1dave/dble(numpt1)
          ! print *, 'spect1dave =', real(sum(spect1dave))
          print *, 'tmean =', tmean, 'tmeanenergy =', tmeanenergy
          call CalPowerspect(spect1dave,tmeanenergy,fileoutput)
          call CalCorrt(spect1dave,fileoutput)
        enddo
        deallocate(spect1dtmp, spect1dave)
      endif

      ! Calculate x-t correlation
      if(icalcorrxt.eq.1) then
         call InitFFT1D(ntrans)
         allocate(spect1dtmp2(ntrans,-ixwindowl:ixwindowr), spect1dave2(ntrans,-ixwindowl:ixwindowr))
         if(.not.allocated(numpt)) allocate(numpt(-ixwindowl:ixwindowr))
         do m = 1, nvar_output
            fileoutput = '_j'//fnum1//'_i'//fnum2//'-'//fnum3//'_'//trim(varout_j(m))//'_nsection'//fnum4
            numpt = 0
            spect1dtmp2 = (0.,0.); spect1dave2 = (0.,0.)
            do k = 1, nzpoint_j
            do i = 1, nxpoint_j
               do ii = -min(ixwindowl,i-imints_j), min(ixwindowr,imaxts_j-i)
                  spect1dtmp2(:,ii) = crossspect1d(buffer_jplane(1:ntpoint_total,i+ii,k,1,m),buffer_jplane(1:ntpoint_total,i,k,1,m))
                  spect1dave2(:,ii) = spect1dave2(:,ii) + spect1dtmp2(:,ii)
                  numpt(ii) = numpt(ii) + 1
               enddo
            enddo
            enddo
            forall(ii=-ixwindowl:ixwindowr,numpt(ii).gt.0)
               spect1dave2(:,ii) = spect1dave2(:,ii)/dble(numpt(ii))
            end forall
            call CalCorrxt(ixwindowl,ixwindowr,dx_jp,spect1dave2,fileoutput)
         enddo
         deallocate(spect1dtmp2, spect1dave2)
      endif

      ! Calculate z-t correlation
      if(icalcorrzt.eq.1) then
         call InitFFT1D(ntrans)
         allocate(spect1dtmp2(ntrans,-kzwindowd:kzwindowu), spect1dave2(ntrans,-kzwindowd:kzwindowu))
         if(.not.allocated(numpt)) allocate(numpt(-kzwindowd:kzwindowu))
         do m = 1, nvar_output
            fileoutput = '_j'//fnum1//'_i'//fnum2//'-'//fnum3//'_'//trim(varout_j(m))//'_nsection'//fnum4
            numpt = 0
            spect1dtmp2 = (0.,0.); spect1dave2 = (0.,0.)
            do k = 1, nzpoint_j
            do i = 1, nxpoint_j
               do kk = -min(kzwindowd,k-kmints_j), min(kzwindowu,kmaxts_j-k)
                  spect1dtmp2(:,kk) = crossspect1d(buffer_jplane(1:ntpoint_total,i,k+kk,1,m),buffer_jplane(1:ntpoint_total,i,k,1,m))
                  spect1dave2(:,kk) = spect1dave2(:,kk) + spect1dtmp2(:,kk)
                  numpt(kk) = numpt(kk) + 1
               enddo
            enddo
            enddo
            forall(kk=-kzwindowd:kzwindowu,numpt(kk).gt.0)
               spect1dave2(:,kk) = spect1dave2(:,kk)/dble(numpt(kk))
            end forall
            call CalCorrzt(kzwindowd,kzwindowu,dz_jp,spect1dave2,fileoutput)
         enddo
         deallocate(spect1dtmp2, spect1dave2)
      endif

      ! Calculate x-z correlation
      if(icalcorrxz.eq.1) then
         call InitFFT1D(nzpoint_j)
         call InitCorry(nzpoint_j, dz_jp)
         allocate(spect1dtmp3(-ixwindowl:ixwindowr,nzpoint_j), spect1dave3(-ixwindowl:ixwindowr,nzpoint_j))
         if(.not.allocated(numpt)) allocate(numpt(-ixwindowl:ixwindowr))
         do m = 1, nvar_output
            fileoutput = '_j'//fnum1//'_i'//fnum2//'-'//fnum3//'_'//trim(varout_j(m))
            numpt = 0
            spect1dtmp3 = (0.,0.); spect1dave3 = (0.,0.)
            do i = 1, nxpoint_j
            do jj = 1, ntpoint_total
               do ii = -min(ixwindowl,i-imints_j), min(ixwindowr,imaxts_j-i)
                  spect1dtmp3(ii,:) = crossspect1dy(buffer_jplane(jj,i+ii,1:nzpoint_j,1,m),buffer_jplane(jj,i,1:nzpoint_j,1,m))
                  spect1dave3(ii,:) = spect1dave3(ii,:) + spect1dtmp3(ii,:)
                  numpt(ii) = numpt(ii) + 1
               enddo
            enddo
            enddo
            forall(ii=-ixwindowl:ixwindowr,numpt(ii).gt.0)
               spect1dave3(ii,:) = spect1dave3(ii,:)/dble(numpt(ii))
            end forall
            call CalCorrxz(ixwindowl,ixwindowr,dx_jp,spect1dave3,fileoutput)
         enddo
         deallocate(spect1dtmp3, spect1dave3)
      endif ! end icalcorrxz

      if(nvar_output.eq.5) then
        buffer_jplane(:,:,:,1,2) = 0.d0
        buffer_p0(:,:,:,1,1) = buffer_jplane(:,:,:,1,4)*( (buffer_jplane(:,:,:,1,5) + 0.5*(1.4-1.0)/1.4/287.0*(buffer_jplane(:,:,:,1,1)**2 + &
                                buffer_jplane(:,:,:,1,2)**2 + buffer_jplane(:,:,:,1,3)**2   ) )/buffer_jplane(:,:,:,1,5)  )**(1.4/(1.4-1.0))
        deallocate(buffer_jplane)
      endif


      ! Calculating streamwise wave-number spectrum for p0
    if(nvar_output.eq.5) then
      if(icalWNSD_kx.eq.1) then
        call InitSpect1dy(nxpoint_j,nsection,iwindow,ioverlap,dx_jp,ntrans,nperiod)
        allocate(tsdata1d(1:nxpoint_j))
        tsdata1d = 0.d0
        allocate(spect1dtmp(ntrans), spect1dave(ntrans))

        spect1dave = 0.d0;  numpt1 = 0
        tmean = 0.d0; tmeanenergy = 0.d0
        fileoutput = '_kx_j'//fnum1//'_i'//fnum2//'-'//fnum3//'_'//trim('p0')
        print *, 'fileoutput =', trim(fileoutput)
        do k=1, nzpoint_j
          do n=1, ntpoint_total
            tsdata1d(1:nxpoint_j) = buffer_p0(n,1:nxpoint_j,k,1,1)
            spect1dtmp = autospect1dty(tsdata1d(1:nxpoint_j))
            spect1dave = spect1dave + spect1dtmp
            tmean = tmean + sum(tsdata1d(1:nxpoint_j))
            tmeanenergy = tmeanenergy + sum(tsdata1d(1:nxpoint_j)**2)
            numpt1 = numpt1 + 1
          enddo
        enddo
        tmean = tmean/dble(nxpoint_j*numpt1)
        tmeanenergy = tmeanenergy/dble(nxpoint_j*numpt1)
        tmeanenergy = tmeanenergy - tmean**2
        spect1dave = spect1dave/dble(numpt1)
        print *, 'tmean =', tmean, 'tmeanenergy =', tmeanenergy
        call CalPowerspecty(spect1dave,tmeanenergy,fileoutput)

        deallocate(spect1dtmp, spect1dave)
      endif ! end icalWNSD_kx.eq.1
    endif ! end nvar_output.eq.5

!      ! Calculating streamwise wave-number spectrum
!      if(icalWNSD_kx.eq.1) then
!        call InitSpect1dy(nxpoint_j,nsection,iwindow,ioverlap,dx_jp,ntrans,nperiod)
!        allocate(tsdata1d(1:nxpoint_j))
!        tsdata1d = 0.d0
!        allocate(spect1dtmp(ntrans), spect1dave(ntrans))
!        do m=1, nvar_output
!          spect1dave = 0.d0;  numpt1 = 0
!          tmean = 0.d0; tmeanenergy = 0.d0
!          fileoutput = '_kx_j'//fnum1//'_i'//fnum2//'-'//fnum3//'_'//trim(varout_j(m))
!          print *, 'fileoutput =', trim(fileoutput)
!          do k=1, nzpoint_j
!            do n=1, ntpoint_total
!              tsdata1d(1:nxpoint_j) = buffer_jplane(n,1:nxpoint_j,k,1,m)
!              spect1dtmp = autospect1dty(tsdata1d(1:nxpoint_j))
!              spect1dave = spect1dave + spect1dtmp
!              tmean = tmean + sum(tsdata1d(1:nxpoint_j))
!              tmeanenergy = tmeanenergy + sum(tsdata1d(1:nxpoint_j)**2)
!              numpt1 = numpt1 + 1
!            enddo
!          enddo
!          tmean = tmean/dble(nxpoint_j*numpt1)
!          tmeanenergy = tmeanenergy/dble(nxpoint_j*numpt1)
!          tmeanenergy = tmeanenergy - tmean**2
!          spect1dave = spect1dave/dble(numpt1)
!          print *, 'tmean =', tmean, 'tmeanenergy =', tmeanenergy
!          call CalPowerspecty(spect1dave,tmeanenergy,fileoutput)
!        enddo ! end m loop
!        deallocate(spect1dtmp, spect1dave)
!      endif ! end icalWNSD_kx.eq.1

    enddo ! end j loop
  endif ! end ireadjp



  !!#################################################################################
  if(ireadkp.eq.1.and.num_kplane.gt.0) then  ! Using constant-k plane time series

    if(iFileType.eq.0) then
      call InitTSkplane(ntpoint_total, npath, DNSIndex_k, ibe_k, iend_k, jbe_k, jend_k, &
                        ixwindowl, ixwindowr, nvar_output, varindex_output, fileprop)
    elseif(iFileType.ge.1) then
      call InitVOLTSkplane(ntpoint_total, npath, DNSIndex_k, ibe_k, iend_k, jbe_k, jend_k, &
                  ixwindowl, ixwindowr, nvar_output, varindex_output, fileprop)
    else
      print *, "iFileType should be 0, 1 or 2"
      stop
    endif

!     if(icalcorrxy.eq.1) then
!         call InitFFT1D(nypoint_k)
!         call InitCorry(nypoint_k, dy_kp)
!     endif
!    if(icalWNSD.eq.1) then
!      call InitSpect1dy(nypoint_k, nsection, iwindow, ioverlap, dy_kp, ntrans, nperiod)
!      allocate(tsdata1d(1:nypoint_k))
!      tsdata1d = 0.d0
!    endif
!    if(icalky_t.eq.1) then
!      call InitSpect2d_ty(ntpoint_total, nypoint_k, nsection, iwindow, ioverlap, dt_sample, dy_kp, ntrans, nperiod)
!      allocate(tsdata2d(1:ntpoint_total, 1:nypoint_k))
!      tsdata2d = 0.d0
!    endif

     allocate( buffer_kplane(ntpoint_total,-jywindow+1:nypoint_k+jywindow, imints_k:imaxts_k, 1, nvar_output) )

    ! Reading k plane data
     do k = 1, num_kplane
       write(unit=fnum1,fmt='(I04.4)') kplane(k)
       write(unit=fnum2,fmt='(I04.4)') ibe_k
       write(unit=fnum3,fmt='(I04.4)') iend_k
       write(unit=fnum4,fmt='(I04.4)') nsection
       if(iFileType.eq.0) then
         call ReadTSkplane(kplane(k),buffer_kplane)
       else
         do kk=1, DNSIndex_k%nkplane
            if(kplane(k).eq.DNSIndex_k%kplane(kk)) then
              kplane_be = kk
            endif
         enddo ! end kk loop
         call ReadVOLTSkplane(kplane_be,kplane(k),buffer_kplane(1:ntpoint_total,1:nypoint_k,imints_k:imaxts_k,1:1,1:nvar_output))
       endif
       print *, 'wall-normal location k =', kplane(k)
       print *, 'streamwise average range i =', ibe_k, iend_k

       !
       if(ical3D.eq.1) then

         call InitSpect3d(ntpoint_total, nypoint_k, nxpoint_k, nsection, nsection_kx, iwindow, ioverlap, dt_sample, dy_kp, dx_kp, ntrans, ntrans_kx)
         allocate(tsdata3d(1:ntpoint_total,1:nypoint_k,1:nxpoint_k))
         tsdata3d = 0.d0

         allocate(spect3dtmp(ntrans,nypoint_k,ntrans_kx),spect3dave(ntrans,nypoint_k,ntrans_kx))
         do m=1, nvar_output
           spect3dave = 0.d0; numpt1 = 0
           tmean = 0.d0; tmeanenergy = 0.d0
           fileoutput = '_k'//fnum1//'_i'//fnum2//'-'//fnum3//'_'//trim(varout_k(m))//'_nsection'//fnum4
           print *, 'fileoutput =', trim(fileoutput)

          ! do i=1, nxpoint_k
           tsdata3d(1:ntpoint_total,1:nypoint_k,1:nxpoint_k) = buffer_kplane(1:ntpoint_total,1:nypoint_k,1:nxpoint_k,1,m)
           spect3dtmp(1:ntrans,1:nypoint_k,1:ntrans_kx) = autospect3d(tsdata3d(1:ntpoint_total,1:nypoint_k,1:nxpoint_k))
          !  tmean = sum(tsdata3d(1:ntpoint_total,1:nypoint_k,i))/dble(ntpoint_total*nypoint_k)
          !  enddo ! end i loop

           call CalPowerspect3D(spect3dtmp,fileoutput)

         enddo ! end m loop

       endif


       ! Calculating ky-t spectrum
       if(icalky_t.eq.1) then
         call InitSpect2d_ty(ntpoint_total, nypoint_k, nsection, iwindow, ioverlap, dt_sample, dy_kp, ntrans, nperiod)
         allocate(tsdata2d(1:ntpoint_total, 1:nypoint_k))
         tsdata2d = 0.d0

         allocate(spect2dtmp(ntrans,nypoint_k), spect2dave(ntrans,nypoint_k))
         do m = 1, nvar_output
           spect2dave = 0.d0; numpt1 = 0
           tmean = 0.d0; tmeanenergy = 0.d0
           fileoutput = 'ky-t_powerspect_k'//fnum1//'_i'//fnum2//'-'//fnum3//'_'//trim(varout_k(m))//'_nsection'//fnum4
           print *, 'fileoutput =', trim(fileoutput)

           do i=1, nxpoint_k
             tsdata2d(1:ntpoint_total,1:nypoint_k) = buffer_kplane(1:ntpoint_total,1:nypoint_k,i,1,m)
             spect2dtmp(:,:) = autospect2d(tsdata2d(1:ntpoint_total,1:nypoint_k))
             spect2dave = spect2dave + spect2dtmp
             tmean = tmean + sum(tsdata2d(1:ntpoint_total,1:nypoint_k))
             tmeanenergy = tmeanenergy + sum(tsdata2d(1:ntpoint_total,1:nypoint_k)**2)
             numpt1 = numpt1 + 1
           enddo ! end i loop
           tmean = tmean/dble(ntpoint_total*nypoint_k*numpt1)
           tmeanenergy = tmeanenergy/dble(ntpoint_total*nypoint_k*numpt1)
           tmeanenergy = tmeanenergy - tmean**2
           spect2dave = spect2dave/dble(numpt1)
           print *, 'tmean =', tmean, 'tmeanenergy =', tmeanenergy
           call CalPowerspect_ky_t(spect2dave,tmeanenergy,fileoutput)
         enddo ! end m loop
         deallocate(spect2dtmp,spect2dave,tsdata2d)
       endif


       ! Calculating kx-t spectrum
       if(icalkx_t.eq.1) then
         call InitSpect2d_ty(ntpoint_total, nxpoint_k, nsection, iwindow, ioverlap, dt_sample, dx_kp, ntrans, nperiod)
         allocate(tsdata2d(1:ntpoint_total, 1:nxpoint_k),tsdata2d_tmp(1:ntpoint_total, 1:nxpoint_k))
         tsdata2d = 0.d0

         allocate(spect2dtmp(ntrans,nxpoint_k), spect2dave(ntrans,nxpoint_k))
         do m = 1, nvar_output
           spect2dave = 0.d0; numpt1 = 0
           tmean = 0.d0; tmeanenergy = 0.d0
           fileoutput = 'kx-t_powerspect_k'//fnum1//'_i'//fnum2//'-'//fnum3//'_'//trim(varout_k(m))//'_nsection'//fnum4
           print *, 'fileoutput =', trim(fileoutput)
           do n=1, ntpoint_total
             do i=1, nxpoint_k
               tsdata2d_tmp(n,i) = sum(buffer_kplane(n,:,i,1,m))/dble(nypoint_k)
             enddo
           enddo

           do i=1, nypoint_k
             tsdata2d(1:ntpoint_total,1:nxpoint_k) = buffer_kplane(1:ntpoint_total,i,1:nxpoint_k,1,m) - tsdata2d_tmp(1:ntpoint_total,1:nxpoint_k)
             spect2dtmp(:,:) = autospect2d(tsdata2d(1:ntpoint_total,1:nxpoint_k))
             spect2dave = spect2dave + spect2dtmp
             tmean = tmean + sum(tsdata2d(1:ntpoint_total,1:nxpoint_k))
             tmeanenergy = tmeanenergy + sum(tsdata2d(1:ntpoint_total,1:nxpoint_k)**2)
             numpt1 = numpt1 + 1
           enddo ! end i loop
           tmean = tmean/dble(ntpoint_total*nxpoint_k*numpt1)
           tmeanenergy = tmeanenergy/dble(ntpoint_total*nxpoint_k*numpt1)
           tmeanenergy = tmeanenergy - tmean**2
           spect2dave = spect2dave/dble(numpt1) 
           print *, 'tmean =', tmean, 'tmeanenergy =', tmeanenergy
           call CalPowerspect_kx_t(spect2dave,tmeanenergy,fileoutput)
         enddo ! end m loop
         deallocate(tsdata2d,tsdata2d_tmp)
         deallocate(spect2dtmp,spect2dave)
       endif

       ! Calculating kx-ky spectrum
       if(icalkx_ky.eq.1) then
         call InitSpect2d_ky_kx(nypoint_k, nxpoint_k, iwindow, dy_kp, dx_kp)
         allocate(tsdata2d(1:nypoint_k,1:nxpoint_k),tsdata2d_tmp(1:ntpoint_total,1:nxpoint_k))
         tsdata2d = 0.d0

         allocate(spect2dtmp(nypoint_k,nxpoint_k), spect2dave(nypoint_k,nxpoint_k))
         do m = 1, nvar_output
           spect2dave = 0.d0; numpt1 = 0
           tmean = 0.d0; tmeanenergy = 0.d0
           fileoutput = 'kx-ky_powerspect_k'//fnum1//'_i'//fnum2//'-'//fnum3//'_'//trim(varout_k(m))
           print *, 'fileoutput =', trim(fileoutput)
           do n=1, ntpoint_total
             do i=1, nxpoint_k
               tsdata2d_tmp(n,i) = sum(buffer_kplane(n,:,i,1,m))/dble(nypoint_k)
             enddo
           enddo

           do n=1, ntpoint_total
             do j=1, nypoint_k
               tsdata2d(j,1:nxpoint_k) = buffer_kplane(n,j,1:nxpoint_k,1,m) - tsdata2d_tmp(n,1:nxpoint_k)
             enddo
             spect2dtmp(:,:) = autospect2d_ky_kx(tsdata2d(1:nypoint_k,1:nxpoint_k))
             spect2dave = spect2dave + spect2dtmp
             tmean = tmean + sum(tsdata2d(1:nypoint_k,1:nxpoint_k))
             tmeanenergy = tmeanenergy + sum(tsdata2d(1:nypoint_k,1:nxpoint_k)**2)
             numpt1 = numpt1 + 1
           enddo ! end i loop
           tmean = tmean/dble(nypoint_k*nxpoint_k*numpt1)
           tmeanenergy = tmeanenergy/dble(nypoint_k*nxpoint_k*numpt1)
           tmeanenergy = tmeanenergy - tmean**2
           spect2dave = spect2dave/dble(numpt1)
           print *, 'tmean =', tmean, 'tmeanenergy =', tmeanenergy
           call CalPowerspect_kx_ky(spect2dave,tmeanenergy,fileoutput)
         enddo ! end m loop
         deallocate(spect2dtmp,spect2dave,tsdata2d,tsdata2d_tmp)
       endif

       ! Calculating PSD
       if(icalPSD.eq.1) then
          call InitFFT1D(ntrans)
          allocate(spect1dtmp(ntrans), spect1dave(ntrans))
          do m = 1, nvar_output
              spect1dave = 0.d0;  numpt1 = 0
              tmean = 0.d0; tmeanenergy = 0.d0
!              subsample = 3; ntpoint_tmp = (ntpoint_total-1)/subsample + 1
!              rms = 0.d0; skewness = 0.d0; flatness = 0.d0
              moment = 0.d0
              fileoutput = '_k'//fnum1//'_i'//fnum2//'-'//fnum3//'_'//trim(varout_k(m))//'_nsection'//fnum4
              print *, 'fileoutput =', trim(fileoutput)
              do i = 1, nxpoint_k
              do j = 1, nypoint_k
                 tsdata1d(1:ntpoint_total) = buffer_kplane(1:ntpoint_total,j,i,1,m)
                 spect1dtmp = autospect1d(tsdata1d(1:ntpoint_total))
                 spect1dave = spect1dave + spect1dtmp
                 tmean = tmean + sum(tsdata1d(1:ntpoint_total))
                 tmeanenergy = tmeanenergy + sum(tsdata1d(1:ntpoint_total)**2)
!                 t1 = sum(tsdata1d(1:ntpoint_total:subsample)   )/ntpoint_tmp
!                 t2 = sum(tsdata1d(1:ntpoint_total:subsample)**2)/ntpoint_tmp
!                 t3 = sum(tsdata1d(1:ntpoint_total:subsample)**3)/ntpoint_tmp
!                 t4 = sum(tsdata1d(1:ntpoint_total:subsample)**4)/ntpoint_tmp
                 t1 = sum(tsdata1d(1:ntpoint_total)   )/ntpoint_total
                 t2 = sum(tsdata1d(1:ntpoint_total)**2)/ntpoint_total
                 t3 = sum(tsdata1d(1:ntpoint_total)**3)/ntpoint_total
                 t4 = sum(tsdata1d(1:ntpoint_total)**4)/ntpoint_total
                 rms_tmp = sqrt(abs(t2 - t1**2))  
                 skewness_tmp = (t3-t1**3-3.0*t1*rms_tmp**2)/(rms_tmp**3+1.e-30)
                 flatness_tmp = (t4-t1**4-6.0*t1**2*rms_tmp**2-4.0*t1*skewness_tmp*rms_tmp**3)/(rms_tmp**4+1.e-30)
                 moment(2) = moment(2) + rms_tmp
                 moment(3) = moment(3) + skewness_tmp
                 moment(4) = moment(4) + flatness_tmp
                 numpt1 = numpt1 + 1
              enddo
              enddo
              tmean = tmean/dble(ntpoint_total*numpt1)
              tmeanenergy = tmeanenergy/dble(ntpoint_total*numpt1)
              tmeanenergy = tmeanenergy - tmean**2
              spect1dave = spect1dave/dble(numpt1)
              moment(1) = tmean
              moment(2:4) = moment(2:4)/dble(numpt1)
!              rms = rms/dble(numpt1)
!              skewness = skewness/dble(numpt1)
!              flatness = flatness/dble(numpt1)
              
              ! print *, 'spect1dave =', real(sum(spect1dave))
              print *, 'tmean =', tmean, 'tmeanenergy =', tmeanenergy
              print *, 'rms =', moment(2), 'skewness =', moment(3), 'flatness =', moment(4)
              call CalPowerspect(spect1dave,tmeanenergy,fileoutput,moment)
              call CalCorrt(spect1dave,fileoutput)

!              ! calculate the phase
!              spect1dave = 0.d0;  numpt1 = 0
!              fileoutput = '_k'//fnum1//'_i'//fnum2//'-'//fnum3//'_'//trim(varout_k(m))//'_nsection'//fnum4
!
!              do i = 1, nxpoint_k
!              do j = 1, nypoint_k
!                 tsdata1d(1:ntpoint_total) = buffer_kplane(1:ntpoint_total,j,i,1,m)
!                 spect1dtmp = sumspect1d(tsdata1d(1:ntpoint_total))
!                 spect1dave = spect1dave + spect1dtmp
!                 numpt1 = numpt1 + 1
!              enddo
!              enddo
!              spect1dave = spect1dave/dble(numpt1)
!              print *, 'fileoutput =', 'phase'//trim(fileoutput)
!              open(7,file='phase'//trim(fileoutput)//'.dat',status='unknown')
!                write(7,'(a)') 'variables=f,phase'
!                do i=1, ntrans/2
!                  write(7,'(2E20.12)') real(i-1)/(dble(ntrans-1)*dt_sample), atan2(imag(spect1dave(i)),real(spect1dave(i)))
!                enddo
!              close(7)
!
!             ! for test
!             open(7,file='pave.dat',status='unknown')
!               write(7,'(a)') 'variables=t,p'
!               do i=1, ntrans
!                 write(7,'(2E20.12)') real(i-1)*dt_sample, tsdata1d(i)
!               enddo
!             close(7)
!
!              call FFT1DB(spect1dave,spect1dtmp)
!               open(7,file='p'//trim(fileoutput)//'.dat',status='unknown')
!                write(7,'(a)') 'variables=t,p'
!                do i=1, ntrans
!                  write(7,'(2E20.12)') real(i-1)*dt_sample, real(spect1dtmp(i))
!                enddo
!              close(7)
          enddo
          deallocate(spect1dtmp, spect1dave)
        endif

    ! Calculating streamwise wave-number spectrum
    if(icalWNSD_kx.eq.1) then

        call InitSpect1dy(nxpoint_k, nsection, iwindow, ioverlap, dx_kp, ntrans, nperiod)
        allocate(tsdata1d(1:nxpoint_k),tsdata1d_tmp(1:nxpoint_k))
        tsdata1d = 0.d0

        allocate(spect1dtmp(ntrans), spect1dave(ntrans))
        do m = 1, nvar_output
          spect1dave = 0.d0;  numpt1 = 0
          tmean = 0.d0; tmeanenergy = 0.d0
          fileoutput = '_kx_k'//fnum1//'_i'//fnum2//'-'//fnum3//'_'//trim(varout_k(m))
          print *, 'fileoutput =', trim(fileoutput)
          do j = 1, nypoint_k
            do n = 1, ntpoint_total
              !tsdata1d(1:nxpoint_k) = buffer_kplane(n,j,1:nxpoint_k,1,m)
              tsdata1d_tmp(1:nxpoint_k) = sum(buffer_kplane(n,1:nypoint_k,1:nxpoint_k,1,m),dim=1)/dble(nypoint_k)
              tsdata1d(1:nxpoint_k) = buffer_kplane(n,j,1:nxpoint_k,1,m) - tsdata1d_tmp(1:nxpoint_k)
              spect1dtmp = autospect1dty(tsdata1d(1:nxpoint_k))
              spect1dave = spect1dave + spect1dtmp
              tmean = tmean + sum(tsdata1d(1:nxpoint_k))
              tmeanenergy = tmeanenergy + sum(tsdata1d(1:nxpoint_k)**2)
              numpt1 = numpt1 + 1
            enddo
          enddo
          tmean = tmean/dble(nxpoint_k*numpt1)
          tmeanenergy = tmeanenergy/dble(nxpoint_k*numpt1)
          tmeanenergy = tmeanenergy - tmean**2
          spect1dave = spect1dave/dble(numpt1)
          print *, 'tmean =', tmean, 'tmeanenergy =', tmeanenergy
          call CalPowerspectx(spect1dave,tmeanenergy,fileoutput)
        enddo ! end m loop
        deallocate(spect1dtmp, spect1dave, tsdata1d,tsdata1d_tmp)

    endif ! end calculating streamwise wave-number spectrum

    ! Calculating spanwise wave-number spectrum
    if(icalWNSD_ky.eq.1) then

        call InitSpect1dy(nypoint_k, nsection, iwindow, ioverlap, dy_kp, ntrans, nperiod)
        allocate(tsdata1d(1:nypoint_k))
        tsdata1d = 0.d0

        allocate(spect1dtmp(ntrans), spect1dave(ntrans))
        do m = 1, nvar_output
          spect1dave = 0.d0;  numpt1 = 0
          tmean = 0.d0; tmeanenergy = 0.d0
          fileoutput = '_ky_k'//fnum1//'_i'//fnum2//'-'//fnum3//'_'//trim(varout_k(m))
          print *, 'fileoutput =', trim(fileoutput)
          do i = 1, nxpoint_k
            do n = 1, ntpoint_total
              tsdata1d(1:nypoint_k) = buffer_kplane(n,1:nypoint_k,i,1,m)
              spect1dtmp = autospect1dty(tsdata1d(1:nypoint_k))
              spect1dave = spect1dave + spect1dtmp
              tmean = tmean + sum(tsdata1d(1:nypoint_k))
              tmeanenergy = tmeanenergy + sum(tsdata1d(1:nypoint_k)**2)
              numpt1 = numpt1 + 1
            enddo
          enddo
          tmean = tmean/dble(nypoint_k*numpt1)
          tmeanenergy = tmeanenergy/dble(nypoint_k*numpt1)
          tmeanenergy = tmeanenergy - tmean**2
          spect1dave = spect1dave/dble(numpt1)
          print *, 'tmean =', tmean, 'tmeanenergy =', tmeanenergy
          call CalPowerspecty(spect1dave,tmeanenergy,fileoutput)
        enddo ! end m loop
        deallocate(spect1dtmp, spect1dave,tsdata1d)


    endif ! end calculating spanwise wave-number spectrum


    ! Reading k plane data
!      do k = 1, num_kplane
!        write(unit=fnum1,fmt='(I04.4)') kplane(k)
!        write(unit=fnum2,fmt='(I04.4)') ibe_k
!        write(unit=fnum3,fmt='(I04.4)') iend_k
!        if(ireadVOL.eq.0) then
!          call ReadTSkplane(kplane(k),buffer_kplane)
!        else
!          do kk=1, DNSIndex_k%nkplane
!            if(kplane(k).eq.DNSIndex_k%kplane(kk)) then
!              kplane_be = kk
!            endif
!          enddo ! end kk loop
!          call ReadVOLTSkplane(kplane_be,kplane(k),buffer_kplane)   !!!!!!!!????????
!        endif
!        print *, 'wall-normal location k =', kplane(k)
!        print *, 'streamwise average range i =', ibe_k, iend_k

!      enddo ! end k loop


       ! Calculate x-t correlation
       if(icalcorrxt.eq.1) then
         call InitFFT1D(ntrans)
         allocate(spect1dtmp2(ntrans,-ixwindowl:ixwindowr), spect1dave2(ntrans,-ixwindowl:ixwindowr))
         if(.not.allocated(numpt)) allocate(numpt(-ixwindowl:ixwindowr))
         do m = 1, nvar_output
            fileoutput = '_k'//fnum1//'_i'//fnum2//'-'//fnum3//'_'//trim(varout_k(m))//'_nsection'//fnum4
            numpt = 0
            spect1dtmp2 = (0.,0.); spect1dave2 = (0.,0.)
            do i = 1, nxpoint_k
            do j = 1, nypoint_k
               do ii = -min(ixwindowl,i-imints_k), min(ixwindowr,imaxts_k-i)
                  spect1dtmp2(:,ii) = crossspect1d(buffer_kplane(1:ntpoint_total,j,i+ii,1,m),buffer_kplane(1:ntpoint_total,j,i,1,m))
                  spect1dave2(:,ii) = spect1dave2(:,ii) + spect1dtmp2(:,ii)
                  numpt(ii) = numpt(ii) + 1
               enddo
            enddo
            enddo
            forall(ii=-ixwindowl:ixwindowr,numpt(ii).gt.0)
               spect1dave2(:,ii) = spect1dave2(:,ii)/dble(numpt(ii))
            end forall
            call CalCorrxt(ixwindowl,ixwindowr,dx_kp,spect1dave2,fileoutput)
         enddo
         deallocate(spect1dtmp2, spect1dave2)
      endif

      ! Calculate x-y correlation
      if(icalcorrxy.eq.1) then

         call InitFFT1D(nypoint_k)
         call InitCorry(nypoint_k, dy_kp)

         allocate(spect1dtmp3(-ixwindowl:ixwindowr,nypoint_k), spect1dave3(-ixwindowl:ixwindowr,nypoint_k))
         if(.not.allocated(numpt)) allocate(numpt(-ixwindowl:ixwindowr))
         do m = 1, nvar_output
            fileoutput = '_k'//fnum1//'_i'//fnum2//'-'//fnum3//'_'//trim(varout_k(m))
            numpt = 0
            spect1dtmp3 = (0.,0.); spect1dave3 = (0.,0.)
            do i = 1, nxpoint_k
            do j = 1, ntpoint_total
               do ii = -min(ixwindowl,i-imints_k), min(ixwindowr,imaxts_k-i)
                  spect1dtmp3(ii,:) = crossspect1dy(buffer_kplane(j,1:nypoint_k,i+ii,1,m),buffer_kplane(j,1:nypoint_k,i,1,m))
                  spect1dave3(ii,:) = spect1dave3(ii,:) + spect1dtmp3(ii,:)
                  numpt(ii) = numpt(ii) + 1
               enddo
            enddo
            enddo
            forall(ii=-ixwindowl:ixwindowr,numpt(ii).gt.0)
               spect1dave3(ii,:) = spect1dave3(ii,:)/dble(numpt(ii))
            end forall
            call CalCorrxy(ixwindowl,ixwindowr,dx_kp,spect1dave3,fileoutput)
         enddo
         deallocate(spect1dtmp3, spect1dave3)
      endif

      ! Calculating Phase speed
      ! Stegen & Van Atta, "A technique for phase speed measurements in turbulent flow", JFM, vol 42, part 4, 689-699
      ! Ref: Juan C. Alamo & Javier Jimenez, JFM, 640, pp 5-26, 2009
      if(icalphasespeed.eq.1) then
         allocate(spectspanave(ntrans,nxpoint_k), spectxspanave(ntrans,nxpoint_k))
         allocate(spect1dtmp(ntrans), spect1dave(ntrans))
         allocate(spect1dxtmp(ntrans), spect1dxave(ntrans))
         allocate(phasespeed(ntrans/2),phase(ntrans/2))

         allocate(spect1d1(ntrans),spect1d2(ntrans))
         allocate(R2(ntrans))

         fperiod = dble(ntrans-1)*dt_sample  ! period per segments
         do m = 1, nvar_output
             spect1dave = 0.d0;  spect1dxave = 0.d0; numpt1 = 0; spectspanave = 0.
             fileoutput = 'k'//fnum1//'_i'//fnum2//'-'//fnum3//'_'//trim(varout_k(m))//'_nsection'//fnum4

             do i = 1, nxpoint_k
                numpt1 = 0
                do j = 1, nypoint_k
                    tsdata1d(1:ntpoint_total) = buffer_kplane(1:ntpoint_total,j,i,1,m)
                    spect1dtmp = sum(spect1d(tsdata1d(1:ntpoint_total)),dim=2)
                    spectspanave(:,i) = spectspanave(:,i) + spect1dtmp
                    numpt1 = numpt1 + 1
                enddo
                numpt1 = numpt1*nsection
                spectspanave(:,i) = spectspanave(:,i)/numpt1
             enddo

             spect1dave = sum(spectspanave(:,3:nxpoint_k-2),dim=2)/dble(nxpoint_k-4)
             spectxspanave = 0.d0
             do i = 3, nxpoint_k-2
                spectxspanave(:,i) = spectspanave(:,i-2)-8.d0*(spectspanave(:,i-1)-spectspanave(:,i+1))-spectspanave(:,i+2)
             enddo
             spectxspanave = spectxspanave/(12.d0*dx_kp)
             spect1dxave = sum(spectxspanave(:,3:nxpoint_k-2),dim=2)/dble(nxpoint_k-4)
             !call CalPhaseSpeed_Jimenez(spect1dave,spect1dxave,fileoutput)

             open(21, file='phasespeed_VanAtta_'//trim(fileoutput)//'.dat')
             write(21,'(a)')'variables=f,dx,phasespeed,phase,coh'
             write(21,'("Zone T=phasespeed, I=", I8, ",J=",I8)') ntrans/2, nxpoint_k/2
             do nn=1, nxpoint_k/2
               spect1dave = 0.
               spect1d1 = 0.; spect1d2 = 0.
               do i=1, nxpoint_k-nn ! nxpoint_k-1
                 spect1dave(:) = spect1dave(:) + spectspanave(:,i)*conjg(spectspanave(:,i+nn))
                 spect1d1(:) = spect1d1(:) + spectspanave(:,i)*conjg(spectspanave(:,i))
                 spect1d2(:) = spect1d2(:) + spectspanave(:,i+nn)*conjg(spectspanave(:,i+nn))
               enddo
               spect1dave = spect1dave/dble(nxpoint_k-nn)
               spect1d1 = spect1d1/dble(nxpoint_k-nn)
               spect1d2 = spect1d2/dble(nxpoint_k-nn)
               R2 = real(spect1dave*conjg(spect1dave)/(spect1d1*spect1d2))

               call CalPhaseSpeed_VanAtta(spect1dave,dx_kp*dble(nn),phasespeed,phase)
               do n=1, ntrans/2
                 write(21,'(3E20.12)') dble(n-1)/fperiod, dx_kp*dble(nn), phasespeed(n), phase(n), R2(n)
               enddo
             enddo ! end nn loop
             close(21)

         enddo  ! end variables loop
         deallocate(spectspanave, spectxspanave)
         deallocate(spect1dtmp, spect1dave)
         deallocate(spect1dxtmp, spect1dxave)
         deallocate(phasespeed,phase)
         deallocate(spect1d1,spect1d2)
         deallocate(R2)
      endif ! icalphasespeed.eq.1

      if(ical_skewness.eq.1) then
          if(allocated(tsdata1d)) deallocate(tsdata1d)
          allocate(tsdata1d(ntpoint_total),tsdata2d(ntpoint_total,-jywindow+1:nypoint_k+jywindow),tsdata3d(ntpoint_total,nypoint_k,nxpoint_k) )
          do m = 1, nvar_output
            if(jywindow.gt.0) then
              ! update j-direction halo
              buffer_kplane(1:ntpoint_total,-jywindow+1:0,imints_k:imaxts_k,1,m) = buffer_kplane(1:ntpoint_total,nypoint_k-jywindow+1:nypoint_k,imints_k:imaxts_k,1,m)
              buffer_kplane(1:ntpoint_total,nypoint_k+1:nypoint_k+jywindow,imints_k:imaxts_k,1,m) = buffer_kplane(1:ntpoint_total,1:jywindow,imints_k:imaxts_k,1,m)
            endif

            dx_tmp = max(dx_kp,dy_kp)
            do nn=1, ixwindowr+1
              iradius_tmp = dble(nn-1)*dx_tmp
              print *, 'radius (m) = ', iradius_tmp
              tsdata3d = 0.; tsdata1d = 0.
              jyw = iradius_tmp/dy_kp; ixw = iradius_tmp/dx_kp
              print *, 'ixw = ', ixw, 'jyw = ', jyw
              if(ilocal_ave.eq.1.and.iradius_tmp.gt.0) then ! circle average
                do i=1, nxpoint_k
                  do j=1, nypoint_k
                    numpt1 = 0; tsdata1d = 0
                    do ii=-ixwindowl, ixwindowr
                      do jj=-jywindow, jywindow
                        tmean = sqrt((dble(ii)*dx_kp)**2 + (dble(jj)*dy_kp)**2)
                        if(tmean.le.iradius_tmp) then
                          tsdata1d(1:ntpoint_total) = tsdata1d(1:ntpoint_total) + buffer_kplane(1:ntpoint_total,j+jj,i+ii,1,m)
                          numpt1 = numpt1 + 1
                        endif
                      enddo
                    enddo
                    tsdata3d(1:ntpoint_total,j,i) = tsdata1d(1:ntpoint_total)/dble(numpt1)
                  enddo
                enddo
                print *, 'numpt1 = ', numpt1
              elseif(ilocal_ave.eq.2.and.iradius_tmp.gt.0) then ! square average

                do i=1, nxpoint_k
                  tsdata2d(1:ntpoint_total,-jyw+1:nypoint_k+jyw) = &
                          sum(buffer_kplane(1:ntpoint_total,-jyw+1:nypoint_k+jyw,i-ixw:i+ixw,1,m),dim=3)/dble(2*ixw+1)
                  do j=1, nypoint_k
                    tsdata1d(1:ntpoint_total) = sum(tsdata2d(1:ntpoint_total,j-jyw:j+jyw),dim=2)/dble(2*jyw+1)
                    tsdata3d(1:ntpoint_total,j,i) = tsdata1d(1:ntpoint_total)
                  enddo ! end j loop
                enddo ! end i loop
              else
                tsdata3d(1:ntpoint_total,1:nypoint_k,1:nxpoint_k) = buffer_kplane(1:ntpoint_total,1:nypoint_k,1:nxpoint_k,1,m)
              endif

              !print *, 'p = ', buffer_kplane(1,-jywindow+1:1,1,1,1)
              !print *, '', buffer_kplane(1,nypoint_k:nypoint_k+jywindow,1,1,1)
              ! for debug
!              write(unit=fnum4,fmt='(I04.4)') nn
!              filename = 'p'//fnum4//'.dat'
!              open(12,file=trim(filename),status='unknown')
!                write(12,'(a)') 'variables=p'
!                write(12,*) 'zone T = k001, i =', ntpoint_total, ' j =', nypoint_k, ' k =', nxpoint_k
!                do i=1, nxpoint_k
!                  do j=1, nypoint_k
!                    do n=1, ntpoint_total
!                      write(12,*) tsdata3d(n,j,i)
!                    enddo
!                  enddo
!                enddo

              write(unit=fnum4,fmt='(I04.4)') nn
              filename = 'Skewness_k'//fnum1//'_i'//fnum2//'-'//fnum3//'_'//trim(varout_k(m))//'_'//fnum4//'.dat'
              open(7,file=trim(filename),status='unknown')
              write(7,*) 'variables = radius, freq_l, freq_u, pave, prms, Sp, Fp, Sp_min, Sp_max, Fp_min, Fp_max'
              if(ilocal_ave.eq.1) write(7,'(a)') '# using local circle average '
              if(ilocal_ave.eq.2) write(7,'(a)') '# using local square average '
              do n=1, num_freq
                pave=0.d0; p2=0.d0; p3=0.d0; p4=0.d0;numpt1=0
                prms = 0.d0; Sp = 0.d0; Fp = 0.d0
                Sp_max = 0.; Sp_min = 10.; Fp_max = 0.; Fp_min = 10.
                do i = 1, nxpoint_k, 2*ixw+1
                  do j = 1, nypoint_k, 2*jyw+1
                    call FFTFilter(ntpoint_total,tsdata3d(1:ntpoint_total,j,i),dt_sample,fl(n),fu(n),tsdata1d(1:ntpoint_total))
                    pave = sum(tsdata1d(1:ntpoint_total))/dble(ntpoint_total)
                    p2  = sum(tsdata1d(1:ntpoint_total)**2)/dble(ntpoint_total)
                    p3  = sum(tsdata1d(1:ntpoint_total)**3)/dble(ntpoint_total)
                    p4  = sum(tsdata1d(1:ntpoint_total)**4)/dble(ntpoint_total)

                    prms_tmp = sqrt(abs(p2-pave**2))
                    Sp_tmp = (p3 - pave**3 - 3.d0*pave*prms_tmp**2 ) /(prms_tmp**3)
                    Fp_tmp = (p4 - pave**4 - 6.d0*pave**2*prms_tmp**2 - 4.d0*pave*Sp_tmp*prms_tmp**3 )/(prms_tmp**4)

                    if(Sp_tmp.lt.Sp_min) Sp_min = Sp_tmp
                    if(Sp_tmp.gt.Sp_max) Sp_max = Sp_tmp
                    if(Fp_tmp.lt.Fp_min) Fp_min = Fp_tmp
                    if(Fp_tmp.gt.Fp_max) Fp_max = Fp_tmp

                    prms = prms + prms_tmp
                    Sp = Sp + Sp_tmp
                    Fp = Fp + Fp_tmp
                    numpt1 = numpt1 + 1
                  enddo
                enddo
                print *, 'numpt = ', numpt1
                prms = prms/dble(numpt1)
                Sp = Sp/dble(numpt1)
                Fp = Fp/dble(numpt1)
                print *, 'pave = ', pave, 'prms = ', prms, 'Sp = ', Sp, 'Fp = ', Fp
                write(7,'(11E20.12)') iradius_tmp, fl(n), fu(n), pave, prms, Sp, Fp, Sp_min, Sp_max, Fp_min, Fp_max
              enddo ! end n loop
              close(7)
            enddo  ! end nn loop
          enddo ! end m loop
      endif

      if(ical_Taylor.eq.1) then
          if(nxpoint_k.ne.1) then
            print *, 'nxpoint_k should be equal to 1 ... '
            stop
          endif
          allocate(tsdata2d(nypoint_k,nxpoint_k))
          do m = 1, nvar_output
            ! original data are plotted for checking
            filename = 'Taylor_k'//fnum1//'_i'//fnum2//'_'//trim(varout_k(m))//'.dat'
            print *, 'writing file: ', trim(filename)
            open(7,file=trim(filename),status='unknown')
            write(7,*) 'variables = t,y,p'
            write(7,'("Zone T=Taylor, I=", I8, ",J=",I4)') ntpoint_total, nypoint_k
            do j=1, nypoint_k
              do nn=1, ntpoint_total
                tmean = sum(buffer_kplane(:,j,1,1,m))/dble(ntpoint_total)
                write(7,*) dble(nn-1)*dt_sample, dble(j-1)*dy_kp, buffer_kplane(nn,j,1,1,m) - tmean
              enddo
            enddo
            close(7)

            do n=1, num_freq
              write(unit=fnum4,fmt='(I04.4)') int(fl(n))/1000
              write(unit=fnum5,fmt='(I04.4)') int(fu(n))/1000
              filename = 'Taylor_k'//fnum1//'_i'//fnum2//'_'//fnum4//'-'//fnum5//'kHz_'//trim(varout_k(m))//'.dat'
              print *, 'writing file: ', trim(filename)
              open(7,file=trim(filename),status='unknown')
              write(7,*) 'variables = t,y,p'
              write(7,'("Zone T=Taylor, I=", I8, ",J=",I4)') ntpoint_total, nypoint_k
              tsdata2d = sum(buffer_kplane(1:ntpoint_total,1:nypoint_k,1:nxpoint_k,1,m),dim=1)/dble(ntpoint_total)

              do i = 1, nxpoint_k
                do j = 1, nypoint_k
                  call FFTFilter(ntpoint_total,buffer_kplane(1:ntpoint_total,j,i,1,m),dt_sample,fl(n),fu(n),buffer_kplane(1:ntpoint_total,j,i,1,m))
                enddo
              enddo
              do j=1, nypoint_k
                do nn=1, ntpoint_total
                  write(7,*) dble(nn-1)*dt_sample, dble(j-1)*dy_kp, buffer_kplane(nn,j,1,1,m)
                enddo
              enddo
              close(7)
            enddo ! end n loop

          enddo ! end m loop
      endif

     ! calculate PDF
      if(icalPDF.eq.1) then
          if(allocated(tsdata1d)) deallocate(tsdata1d)
          if(allocated(tsdata2d)) deallocate(tsdata2d)
          if(allocated(tsdata3d)) deallocate(tsdata3d)

          allocate(tsdata1d(ntpoint_total),tsdata2d(ntpoint_total,-jywindow+1:nypoint_k+jywindow),tsdata3d(ntpoint_total,nypoint_k,nxpoint_k) )

          do m = 1, nvar_output
            if(jywindow.gt.0) then
              ! update j-direction halo
              buffer_kplane(1:ntpoint_total,-jywindow+1:0,imints_k:imaxts_k,1,m) = buffer_kplane(1:ntpoint_total,nypoint_k-jywindow+1:nypoint_k,imints_k:imaxts_k,1,m)
              buffer_kplane(1:ntpoint_total,nypoint_k+1:nypoint_k+jywindow,imints_k:imaxts_k,1,m) = buffer_kplane(1:ntpoint_total,1:jywindow,imints_k:imaxts_k,1,m)
            endif

            dx_tmp = max(dx_kp,dy_kp)
            !do nn=1, ixwindowr+1
            do nn=2, ixwindowr+1
              iradius_tmp = dble(nn-1)*dx_tmp
              print *, 'radius (m) = ', iradius_tmp
              tsdata3d = 0.; tsdata1d = 0.
              jyw = iradius_tmp/dy_kp; ixw = iradius_tmp/dx_kp
              print *, 'ixw = ', ixw, 'jyw = ', jyw
              ntpoint_tmp1 = ntpoint_total*( (nxpoint_k-1)/(2*ixw+1)+1)*( (nypoint_k-1)/(2*jyw+1)+1)
              print *, 'ntpoint_tmp1 = ', ntpoint_tmp1
              if(allocated(tsdata1d_tmp)) deallocate(tsdata1d_tmp,pdf)
              allocate(tsdata1d_tmp(ntpoint_tmp1 ),pdf(ntpoint_tmp1))

              if(ilocal_ave.eq.1.and.iradius_tmp.gt.0) then ! circle average
                do i=1, nxpoint_k
                  do j=1, nypoint_k
                    numpt1 = 0; tsdata1d = 0
                    do ii=-ixwindowl, ixwindowr
                      do jj=-jywindow, jywindow
                        tmean = sqrt((dble(ii)*dx_kp)**2 + (dble(jj)*dy_kp)**2)
                        if(tmean.le.iradius_tmp) then
                          tsdata1d(1:ntpoint_total) = tsdata1d(1:ntpoint_total) + buffer_kplane(1:ntpoint_total,j+jj,i+ii,1,m)
                          numpt1 = numpt1 + 1
                        endif
                      enddo
                    enddo
                    tsdata3d(1:ntpoint_total,j,i) = tsdata1d(1:ntpoint_total)/dble(numpt1)
                  enddo
                enddo
                print *, 'numpt1 = ', numpt1
              elseif(ilocal_ave.eq.2.and.iradius_tmp.gt.0) then ! square average
                do i=1, nxpoint_k
                  tsdata2d(1:ntpoint_total,-jyw+1:nypoint_k+jyw) = sum(buffer_kplane(1:ntpoint_total,-jyw+1:nypoint_k+jyw,i-ixw:i+ixw,1,m),dim=3)/dble(2*ixw+1)
                  do j=1, nypoint_k
                    tsdata1d(1:ntpoint_total) = sum(tsdata2d(1:ntpoint_total,j-jyw:j+jyw),dim=2)/dble(2*jyw+1)
                    tsdata3d(1:ntpoint_total,j,i) = tsdata1d(1:ntpoint_total)
                  enddo ! end j loop
                enddo ! end i loop
              else
                tsdata3d(1:ntpoint_total,1:nypoint_k,1:nxpoint_k) = buffer_kplane(1:ntpoint_total,1:nypoint_k,1:nxpoint_k,1,m)
              endif

            !print *, 'p = ', buffer_kplane(1,-jywindow+1:1,1,1,1)
            !print *, '', buffer_kplane(1,nypoint_k:nypoint_k+jywindow,1,1,1)
            ! for debug
!            open(12,file='p.dat',status='unknown')
!              write(12,'(a)') 'variables=p'
!              write(12,*) 'zone T = k001, i =', ntpoint_total, ' j =', nypoint_k, ' k =', nxpoint_k
!              do nn=1, nxpoint_k
!                do j=1, nypoint_k
!                  do n=1, ntpoint_total
!                    write(12,*) tsdata3d(n,j,nn)
!                  enddo
!                enddo
!              enddo

              write(unit=fnum4,fmt='(I04.4)') nn
              filename = 'PDF_k'//fnum1//'_i'//fnum2//'-'//fnum3//'_'//trim(varout_k(m))//'_'//fnum4//'.dat'
              open(7,file=trim(filename),status='unknown')
              write(7,*) 'variables =p, PDF'
              if(ilocal_ave.eq.1) write(7,'(a)') '# using local circle average '
              if(ilocal_ave.eq.2) write(7,'(a)') '# using local square average '
              write(7,'(a,E15.8)') '# radius (m) = ', iradius_tmp
              do n=1, num_freq
                ntpoint_tmp2 = 0
                do i = 1, nxpoint_k, 2*ixw+1              !!!!! skip
                  do j = 1, nypoint_k, 2*jyw+1            !!!!! skip
                    call FFTFilter(ntpoint_total,tsdata3d(1:ntpoint_total,j,i),dt_sample,fl(n),fu(n),tsdata1d(1:ntpoint_total))
                    tsdata1d_tmp(ntpoint_tmp2+1:ntpoint_tmp2+ntpoint_total) = tsdata1d(1:ntpoint_total) - sum(tsdata1d)/dble(ntpoint_total)
                    ntpoint_tmp2 = ntpoint_tmp2 + ntpoint_total
                  enddo
                enddo

print *, 'ntpoint_tmp2 = ', ntpoint_tmp2

                call sort(ntpoint_tmp1,tsdata1d_tmp)
                pave = sum(tsdata1d_tmp)/dble(ntpoint_tmp1)
                p2 = sum(tsdata1d_tmp**2)/dble(ntpoint_tmp1)
                prms = sqrt(abs(p2 - pave**2))

                !tsdata1d_tmp = tsdata1d_tmp/prms      ! normalization
                tsdata1d_tmp = tsdata1d_tmp/17345.3
                pave = sum(tsdata1d_tmp)/dble(ntpoint_tmp1)

                std = sqrt( sum((tsdata1d_tmp-pave)**2)/dble(ntpoint_tmp1-1))
                do i=1, ntpoint_tmp1
                  pdf(i) = exp(-0.5*( tsdata1d_tmp(i) - pave )**2/std**2 )/(std*sqrt(8.d0*atan(1.d0)))
                  write(7,*) tsdata1d_tmp(i), pdf(i)
                enddo

              enddo ! end n loop
              close(7)
            enddo ! end nn loop
          enddo ! end m loop
      endif




!      if(ical_skewness.eq.1) then
!          if(allocated(tsdata1d)) deallocate(tsdata1d)
!          allocate(tsdata1d(ntpoint_total))
!          do m = 1, nvar_output
!            filename = 'Skewness_k'//fnum1//'_i'//fnum2//'-'//fnum3//'_'//trim(varout_k(m))//'.dat'
!            open(7,file=trim(filename),status='unknown')
!            write(7,*) 'variables = freq_l, freq_u, pave, prms, Sp, Fp'
!            do n=1, num_freq
!              pave=0.d0; p2=0.d0; p3=0.d0; p4=0.d0;numpt1=0
!              prms = 0.d0; Sp = 0.d0; Fp = 0.d0
!              do i = 1, nxpoint_k
!                do j = 1, nypoint_k
!                  !tsdata1d(1:ntpoint_total) = buffer_kplane(1:ntpoint_total,j,i,1,m)
!                  call FFTFilter(ntpoint_total,buffer_kplane(1:ntpoint_total,j,i,1,m),dt_sample,fl(n),fu(n),tsdata1d(1:ntpoint_total))
!
!                  pave = sum(tsdata1d(1:ntpoint_total))/dble(ntpoint_total)
!                  p2  = sum(tsdata1d(1:ntpoint_total)**2)/dble(ntpoint_total)
!                  p3  = sum(tsdata1d(1:ntpoint_total)**3)/dble(ntpoint_total)
!                  p4  = sum(tsdata1d(1:ntpoint_total)**4)/dble(ntpoint_total)
!
!                  ! print *, 'spect1dave =', real(sum(spect1dave))
!                  prms_tmp = sqrt(abs(p2-pave**2))
!                  Sp_tmp = (p3 - pave**3 - 3.d0*pave*prms_tmp**2 ) /(prms_tmp**3)
!                  Fp_tmp = (p4 - pave**4 - 6.d0*pave**2*prms_tmp**2 - 4.d0*pave*Sp_tmp*prms_tmp**3 )/(prms_tmp**4)
!
!                  prms = prms + prms_tmp
!                  Sp = Sp + Sp_tmp
!                  Fp = Fp + Fp_tmp
!                  numpt1 = numpt1 + 1
!                enddo
!              enddo
!
!              prms = prms/dble(numpt1)
!              Sp = Sp/dble(numpt1)
!              Fp = Fp/dble(numpt1)
!
!              print *, 'pave = ', pave, 'prms = ', prms, 'Sp = ', Sp, 'Fp = ', Fp
!              write(7,'(6E20.12)') fl(n), fu(n), pave, prms, Sp, Fp
!            enddo ! end n loop
!            close(7)
!          enddo ! end m loop
!      endif
!
!      if(ical_skewness.eq.1) then
!          if(allocated(tsdata1d)) deallocate(tsdata1d)
!          allocate(tsdata1d(ntpoint_total))
!          do m = 1, nvar_output
!            filename = 'Skewness_k'//fnum1//'_i'//fnum2//'-'//fnum3//'_'//trim(varout_k(m))//'.dat'
!            open(7,file=trim(filename),status='unknown')
!            write(7,*) 'variables = freq_l, freq_u, pave, prms, Sp, Fp'
!            do n=1, num_freq
!              pave=0.d0; p2=0.d0; p3=0.d0; p4=0.d0;numpt1=0
!              do i = 1, nxpoint_k
!              do j = 1, nypoint_k
!                 !tsdata1d(1:ntpoint_total) = buffer_kplane(1:ntpoint_total,j,i,1,m)
!                 call FFTFilter(ntpoint_total,buffer_kplane(1:ntpoint_total,j,i,1,m),dt_sample,fl(n),fu(n),tsdata1d(1:ntpoint_total))
!
!                 pave = pave + sum(tsdata1d(1:ntpoint_total))
!                 p2  = p2  + sum(tsdata1d(1:ntpoint_total)**2)
!                 p3  = p3  + sum(tsdata1d(1:ntpoint_total)**3)
!                 p4  = p4  + sum(tsdata1d(1:ntpoint_total)**4)
!                 numpt1 = numpt1 + 1
!              enddo
!              enddo
!              pave = pave/dble(ntpoint_total*numpt1)
!              p2 = p2/dble(ntpoint_total*numpt1)
!              p3 = p3/dble(ntpoint_total*numpt1)
!              p4 = p4/dble(ntpoint_total*numpt1)
!              ! print *, 'spect1dave =', real(sum(spect1dave))
!              prms = sqrt(abs(p2-pave**2))
!              Sp = (p3 - pave**3 - 3.d0*pave*prms**2 ) /(prms**3)
!              Fp = (p4 - pave**4 - 6.d0*pave**2*prms**2 - 4.d0*pave*Sp*prms**3 )/(prms**4)
!
!              print *, 'pave = ', pave, 'prms = ', prms, 'Sp = ', Sp, 'Fp = ', Fp
!              write(7,'(6E20.12)') fl(n), fu(n), pave, prms, Sp, Fp
!            enddo ! end n loop
!            close(7)
!          enddo ! end m loop
!      endif


    enddo  ! end looping over k planes
    deallocate(buffer_kplane)
  endif  ! end calculating const-k plane

  ! close file DebugFileName
  close(66)


   ! test FFT
!   call InitFFT1D(5)
!   allocate(spect1dtmp(5))
!   allocate(spect1dave(5))
!   allocate(spect1dxtmp(5))
!   spect1dtmp(1) = 0.
!   spect1dtmp(2) = 1.
!   spect1dtmp(3) = 0.
!   spect1dtmp(4) = 1.
!   spect1dtmp(5) = 0.
!   call FFT1DF(spect1dtmp,spect1dave)
!   print *, 'spect1dtmp = ', spect1dtmp
!   print *, 'spect1dave = ', spect1dave
!   call FFT1DB(spect1dave,spect1dxtmp)
!   print *, 'spect1dxtmp = ', spect1dxtmp

contains

    ! quicksort
    SUBROUTINE sort(n,arr)
      INTEGER n,M,NSTACK
      REAL arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,ir,j,jstack,k,l,istack(NSTACK)
      REAL a,temp

      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do j=l+1,ir
          a=arr(j)
          do i=j-1,l,-1
            if(arr(i).le.a)goto 2
            arr(i+1)=arr(i)
          enddo
          i=l-1
2         arr(i+1)=a
        enddo
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        temp=arr(k)
        arr(k)=arr(l+1)
        arr(l+1)=temp
        if(arr(l).gt.arr(ir))then
          temp=arr(l)
          arr(l)=arr(ir)
          arr(ir)=temp
        endif
        if(arr(l+1).gt.arr(ir))then
          temp=arr(l+1)
          arr(l+1)=arr(ir)
          arr(ir)=temp
        endif
        if(arr(l).gt.arr(l+1))then
          temp=arr(l)
          arr(l)=arr(l+1)
          arr(l+1)=temp
        endif
        i=l+1
        j=ir
        a=arr(l+1)
3       continue
        i=i+1
        if(arr(i).lt.a)goto 3
4       continue
        j=j-1
        if(arr(j).gt.a)goto 4
        if(j.lt.i)goto 5
        temp=arr(i)
        arr(i)=arr(j)
        arr(j)=temp
        goto 3
5       arr(l+1)=arr(j)
        arr(j)=a
        jstack=jstack+2
        if(jstack.gt.NSTACK) pause 'NSTACK too small in sort'
        if(ir-i+1.ge.j-l)then
           istack(jstack)=ir
           istack(jstack-1)=i
           ir=j-1
        else
           istack(jstack)=j-1
           istack(jstack-1)=l
           l=i
        endif
      endif
     goto 1
    END subroutine Sort

    subroutine Input()
      integer, parameter :: nid=5
      integer :: i, j, k, kk, n
      logical :: PlaneNameMatch = .false.

      open(66,file=trim(DebugFileName),status='unknown')
      npath = 1
      read(*,*)
      read(*,*) iFileType
      allocate(fileprop(npath))
      read(*,*)
      do n=1, npath
        read(*,*)
        read(*,*) fileprop(n)%file_be, fileprop(n)%file_end, fileprop(n)%file_skip
        read(*,*)
        read(*,'(a)') fileprop(n)%filepath
      enddo
      read(*,*)
      read(*,*)
      read(*,*) icalPSD, icalcorrxt, icalcorrxy, icalWNSD_ky, icalWNSD_kx, icalphasespeed, icalky_t, icalkx_t, ical3D, ical_skewness, icalkx_ky, icalWNSD_kz
      read(*,*)
      read(*,*) nsection, nsection_kx, ioverlap, iwindow
      read(*,*)
      read(*,*)
      read(*,*) nvar_output
      allocate(varindex_output(nvar_output))
      read(*,*) (varindex_output(i), i = 1, nvar_output)
      read(*,*)
      read(*,*)
      !> timeseries k-plane
      read(*,*)
      read(*,*) ireadkp
      if(ireadkp.gt.0) then
        read(*,*)
        read(*,*) ibe_k, iend_k, jbe_k, jend_k, ixwindowl, ixwindowr !, dx_kp, dy_kp
        read(*,*)
        read(*,*) num_kplane
        allocate(kplane(num_kplane))
        read(*,*) (kplane(i), i = 1, num_kplane)
        if(iFileType.eq.0) then
          print *, 'Reading DNS_index information from input file ... '
          read(*,*)
          read(*,*)
          read(*,*) fileprop(1)%ntpoint, fileprop(1)%nskip, dx_kp, dy_kp, dt_sample
          read(*,*)
          read(*,*) DNSIndex_k%ibe, DNSIndex_k%iend, DNSIndex_k%iskip, &
                    DNSIndex_k%jbe, DNSIndex_k%jend, DNSIndex_k%jskip, DNSIndex_k%ibuffer
        else
          call ReadDNS_index_kplane(fileprop(1)%ntpoint, dx_kp, dy_kp, dt_sample, DNSIndex_k,fileprop(1)%filepath)
          fileprop(1)%nskip = 0
          do n=1, 5
            read(*,*)
          enddo
        endif
        read(*,*)
        read(*,*)
      else
        do n=1, 12
          read(*,*)
        enddo
      endif
      !> timeseries i-plane
      read(*,*)
      read(*,*) ireadip
      if(ireadip.gt.0) then
        read(*,*)
        read(*,*) jbe_i, jend_i, kbe_i, kend_i, ixwindowl, ixwindowr
        read(*,*)
        read(*,*) num_iplane
        allocate(iplane(num_iplane))
        read(*,*) (iplane(i),i=1,num_iplane)
        if(iFileType.eq.0) then
          print *, 'Reading DNS_index information from input file ... '
          read(*,*)
          read(*,*)
          read(*,*) fileprop(1)%ntpoint, fileprop(1)%nskip, dy_ip, dz_ip, dt_sample
          read(*,*)
          read(*,*) DNSIndex_i%jbe, DNSIndex_i%jend, DNSIndex_i%jskip, &
                    DNSIndex_i%kbe, DNSIndex_i%kend, DNSIndex_i%kskip, DNSIndex_i%jbuffer
        else
          call ReadDNS_index_iplane(fileprop(1)%ntpoint, dy_ip, dz_ip, dt_sample, DNSIndex_i,fileprop(1)%filepath)
          fileprop(1)%nskip = 0
          do n=1, 5
            read(*,*)
          enddo
        endif
        read(*,*)
        read(*,*)
      else
        do n=1, 12
          read(*,*)
        enddo
      endif
      !> timeseries j-plane
      read(*,*)
      read(*,*) ireadjp
      if(ireadjp.gt.0) then
        read(*,*)
        read(*,*) ibe_j, iend_j, kbe_j, kend_j, ixwindowl, ixwindowr, kzwindowd, kzwindowu
        read(*,*)
        read(*,*) num_jplane
        allocate(jplane(num_jplane))
        read(*,*) (jplane(i),i=1,num_jplane)
        call ReadDNS_index_jplane(fileprop(1)%ntpoint, dx_jp, dz_jp, dt_sample, DNSIndex_j,fileprop(1)%filepath)
        read(*,*)
        read(*,*) icalcorrxz, icalcorrzt, dx_jp, dz_jp
        read(*,*)
      else
        do n=1, 8
          read(*,*)
        enddo
      endif

      read(*,*)
      read(*,*)
      read(*,*) ireadTSVol
      if(ireadTSVol.gt.0) then
        read(*,*)
        read(*,*) icalcorrzt
        read(*,*)
        read(*,*) ibe_Vol, iend_Vol, jbe_Vol, jend_Vol, kbe_Vol, kend_Vol, ixwindowl, ixwindowr, kzwindowd, kzwindowu
        read(*,*)
        read(*,*) dx_Vol, dy_Vol, dz_Vol, dt_sample, fileprop(1)%ntpoint
        !call ReadDNS_index_Vol(fileprop(1)%ntpoint,dt_sample, DNSIndex_Vol, fileprop(1)%filepath )
        read(*,*)
        read(*,*) DNSIndex_Vol%ibe, DNSIndex_Vol%iend, DNSIndex_Vol%jbe, DNSIndex_Vol%jend, DNSIndex_Vol%kbe, DNSIndex_Vol%kend
      else
        read(*,*)
        read(*,*)
        read(*,*)
        read(*,*)
        read(*,*)
        read(*,*)
        read(*,*)
        read(*,*)
      endif
      read(*,*)
      read(*,*)
      read(*,*)
      read(*,*) ilocal_ave, iradius
      read(*,*)
      read(*,*) num_freq
      if(num_freq.gt.0) then
        allocate(fl(num_freq), fu(num_freq))
        do n=1, num_freq
          read(*,*) fl(n), fu(n)
        enddo
      endif

      if((ireadkp+ireadip+ireadjp+ireadTSVOL).gt.1) then
        print *, 'Only one of ireadkp, ireadip, ireadjp and ireadTSVOL can be 1 ... STOP! '
        stop
      endif

      ! check the input k-plane information
      if(iFileType.ne.0.and.ireadkp.eq.1) then
        if(num_kplane.le.DNSIndex_k%nkplane) then
          do k=1, num_kplane
            do kk=1, DNSIndex_k%nkplane
              if(kplane(k).eq.DNSIndex_k%kplane(kk)) then
                PlaneNameMatch = .true.
              endif
            enddo ! end kk loop
            if(.not.PlaneNameMatch) then
              print *, '#############################################################'
              print *, 'The input kplane name does not match the name in DNSIndex.h5 '
              print *, 'Input kplane name ', kplane(k) ,' does not exist in file DNSIndex.h5 '
              print *, 'Stop ... '
              stop
            else
              PlaneNameMatch = .false.
            endif
          enddo ! end k loop

        else ! num_kplane.gt.DNSIndex_k%nkplane
          print *, 'The input num_kplane does not match the num_kplane in DNSIndex.h5. '
          print *, 'The input num_kplane = ', num_kplane, 'num_kplane in DNSIndex.h5 = ', DNSIndex_k%nkplane
          stop
        endif
      endif

      ! check the input i-plane information
      if(iFileType.ne.0.and.ireadip.eq.1) then
        if(num_iplane.le.DNSIndex_i%niplane) then
          do i=1, num_iplane
            do kk=1, DNSIndex_i%niplane
              if(iplane(i).eq.DNSIndex_i%iplane(kk)) then
                PlaneNameMatch = .true.
              endif
            enddo ! end kk loop
            if(.not.PlaneNameMatch) then
              print *, '#############################################################'
              print *, 'The input iplane name does not match the name in DNSIndex.h5 '
              print *, 'Input iplane name ', iplane(i) ,' does not exist in file DNSIndex.h5 '
              print *, 'Stop ... '
              stop
            else
              PlaneNameMatch = .false.
            endif
          enddo ! end i loop

        else ! num_iplane.gt.DNSIndex_i%niplane
          print *, 'The input num_iplane does not match the num_iplane in DNSIndex.h5. '
          print *, 'The input num_iplane = ', num_iplane, 'num_iplane in DNSIndex.h5 = ', DNSIndex_i%niplane
          stop
        endif
      endif

      ! check the input j-plane information
      if(iFileType.ne.0.and.ireadjp.eq.1) then
        if(num_jplane.le.DNSIndex_j%njplane) then
          do i=1, num_jplane
            do kk=1, DNSIndex_j%njplane
              if(jplane(i).eq.DNSIndex_j%jplane(kk)) then
                PlaneNameMatch = .true.
              endif
            enddo ! end kk loop
            if(.not.PlaneNameMatch) then
              print *, '#############################################################'
              print *, 'The input iplane name does not match the name in DNSIndex.h5 '
              print *, 'Input iplane name ', jplane(i) ,' does not exist in file DNSIndex.h5 '
              print *, 'Stop ... '
              stop
            else
              PlaneNameMatch = .false.
            endif
          enddo ! end i loop
        else ! num_jplane.gt.DNSIndex_j%njplane
          print *, 'The input num_jplane does not match the num_jplane in DNSIndex.h5. '
          print *, 'The input num_jplane = ', num_jplane, 'num_jplane in DNSIndex.h5 = ', DNSIndex_j%njplane
          stop
        endif
        if(icalcorrzt.ne.1) then
          kzwindowd = 0
          kzwindowu = 0
        endif
        print *, 'kzwindowd =', kzwindowd, 'kzwindowu =', kzwindowu
      endif


      ! calculate the total number of time point
      ntpoint_total = 0
      if(iFileType.eq.0) then
        do n=1, npath
          fileprop(n)%ntpoint = fileprop(n)%ntpoint - fileprop(n)%nskip
          ntpoint_total = ntpoint_total + fileprop(n)%ntpoint
        enddo
      elseif(iFileType.ge.1) then
        do n=1, npath
          ntpoint_total = ntpoint_total + ( (fileprop(n)%file_end-fileprop(n)%file_be)/fileprop(n)%file_skip + 1 )*fileprop(n)%ntpoint
        enddo
      endif

      if(icalcorrxt.eq.0.and.icalcorrxy.eq.0.and.icalcorrxz.eq.0) then
         ixwindowl = 0
         ixwindowr = 0
      endif
      print *, 'ixwindowl =', ixwindowl, 'ixwindowr =', ixwindowr
      write(66,*) 'ixwindowl =', ixwindowl, 'ixwindowr =', ixwindowr
      if (ioverlap.ne.0.and.ioverlap.ne.1) then
        print*,'ioverlap can ONLY be 0 or 1'
        stop
      end if
      if(ireadkp.eq.1) then
         print *, 'number of wall-normal locations =', num_kplane
         print *, 'Wall-normal plane indexes k =', kplane
         write(66,*) 'number of wall-normal locations =', num_kplane
         write(66,*) 'Wall-normal plane indexes k =', kplane
      endif
      if(ireadip.eq.1) then
         print *, 'number of streamwise locations =', num_iplane
         print *, 'Streamwise plane indexes i =', iplane
         write(66,*) 'number of streamwise locations =', num_iplane
         write(66,*) 'Streamwise plane indexes i =', iplane
      endif

      if(ireadTSVol.eq.1) then
         if(ibe_Vol.lt.DNSIndex_Vol%ibe.or.jbe_Vol.lt.DNSIndex_Vol%jbe.or.kbe_Vol.lt.DNSIndex_Vol%kbe.or.iend_Vol.gt.DNSIndex_Vol%iend &
            .or.jend_Vol.gt.DNSIndex_Vol%jend.or.kend_Vol.gt.DNSIndex_Vol%kend) then
           print *, 'selected region is out of bound... STOP!'
           print *, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`'
           print *, 'Selected region'
           print *, 'ibe_Vol = ', ibe_Vol, 'iend_Vol = ', iend_Vol
           print *, 'jbe_Vol = ', jbe_Vol, 'jend_Vol = ', jend_Vol
           print *, 'kbe_Vol = ', kbe_Vol, 'kend_Vol = ', kend_Vol
           print *, 'DNS region'
           print *, 'ibe_fromDNS = ', DNSIndex_Vol%ibe, 'iend_fromDNS = ', DNSIndex_Vol%iend
           print *, 'jbe_fromDNS = ', DNSIndex_Vol%jbe, 'jend_fromDNS = ', DNSIndex_Vol%jend
           print *, 'kbe_fromDNS = ', DNSIndex_Vol%kbe, 'kend_fromDNS = ', DNSIndex_Vol%kend
           print *, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`'
           stop
         endif
         print *, 'dx_Vol = ', dx_Vol, 'dy_Vol = ', dy_Vol, 'dz_Vol = ', dz_Vol, 'dt_sample = ', dt_sample
         if(icalcorrzt.ne.1) then
           kzwindowd = 0
           kzwindowu = 0
         endif
         print *, 'kzwindowd =', kzwindowd, 'kzwindowu =', kzwindowu
      endif

      ! close file DebugFileName
      close(66)

      if((ilocal_ave.eq.1.or.ilocal_ave.eq.2).and.ireadkp.eq.1.and.ical_skewness.eq.1) then
        print *, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
        if(ilocal_ave.eq.1) print *, 'using circle local average'
        if(ilocal_ave.eq.2) print *, 'using square local average'
        print *, 'dx_kp = ', dx_kp, 'dy_kp = ', dy_kp, 'iradius = ', iradius
        ixwindowl = int(iradius/dx_kp)
        ixwindowr = ixwindowl
        jywindow  = int(iradius/dy_kp)
        print *, 'ixwindowl = ', ixwindowl, 'ixwindowr = ', ixwindowr, 'jywindow = ', jywindow
        print *, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      else
        jywindow = 0
      endif

      if(num_freq.eq.0) then
        num_freq = 1
        allocate(fl(num_freq), fu(num_freq))
        fl(1) = -1.e3
        fu(1) = 1.e30
      else
        do n=1, num_freq
          if(fl(n).gt.fu(n)) then
            print *, ' lower bound for frequency is greater than the upper bound. STOP! '
            print *, 'freq_lower = ', fl(n), 'freq_upper = ', fu(n)
            stop
          endif
        enddo
      endif

    end subroutine input


end program spectcorr
