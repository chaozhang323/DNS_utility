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
  use modSpect
  use modSpecty
  implicit none
  integer :: icalPSD, icalcorrxt, icalcorrxy, icalphasespeed, iascii
  integer :: ireadkp, ireadip, nodes_k, nodes_i
  integer :: nskip, ixwindowl, ixwindowr
  integer :: nperiod, ioverlap, ntrans, nsection
  integer :: iwindow
  integer :: num_kplane, num_iplane, num_jplane
  integer :: kplane_be, kplane_end, iplane_be, iplane_end, jplane_be, jplane_end
!  integer :: nxpoint, nypoint, nzpoint, imints, imaxts
  integer, dimension(:), allocatable :: kplane, iplane
  real(8), dimension(:), allocatable :: tsdata1d, tsdata1dx

  complex, dimension(:), allocatable :: spect1dtmp, spect1dave
  complex, dimension(:), allocatable :: spect1dxtmp, spect1dxave
  complex, dimension(:,:), allocatable :: spect1dtmp2, spect1dave2
  complex, dimension(:,:), allocatable :: spect1dtmp3, spect1dave3
  complex, dimension(:,:), allocatable :: spect2dtmp, spect2dave
  complex, dimension(:,:), allocatable :: spectspanave, spectxspanave
  real(8), dimension(:), allocatable :: powerspect
  integer :: n, ibe_k, iend_k, jbe_k, jend_k
!  character(200) :: filepath
  type(fprop), dimension(:), allocatable :: fileprop  !!!!!!!!
  character(200) :: fileinput, fileoutput  !!!!!!!!
  character(4) :: fnum1, fnum2, fnum3
  character(4) :: fnum4, fnum5
  integer :: nvar_output
  integer, dimension(:), allocatable :: varindex_output

  integer :: iii, ii, jj, i, j, k, nn, m, numpt1, numpt2, kk
  real(8) :: dt_sample, dx_kp, dy_kp, dx_ip, dy_ip
  real(8) :: tmean, tmeanenergy, tmpenergy, ff, fperiod
  integer, dimension(:), allocatable :: numpt
  integer :: nnn
  real(8),pointer :: zloc_ip(:), dzdk_ip(:)
  logical :: IsZgridRead = .false.
  integer :: klen
  real(8), dimension(:), allocatable :: phasespeed
  real(8) :: zloc_kp, dzdk_kp
  integer :: ntpoint_total                                     !!!!!!!!!!!!!
  real(8), dimension(:,:,:,:,:), allocatable :: buffer_kplane, buffertmp
  type(tp_DNSIndex) :: DNSIndex_k

  !!!!!!!!!!!!!!!!!
  integer :: npath, ireadVOL, iAverage, num_file, ireadMulVol, buffer_MulVol
  integer :: n2, nf, ntpoint_tmp
  character(400) :: filename
  character(8) :: fnum_8

  call Input()
  call InitHDF5()
  allocate(tsdata1d(1:ntpoint_total))
  tsdata1d = 0.
  print *, 'ntpoint_total =', ntpoint_total, 'dt_sample =', dt_sample
  if(icalPSD.eq.1.or.icalcorrxt.eq.1.or.icalPhaseSpeed.eq.1) then
      call InitSpect1d(ntpoint_total, nperiod, iwindow, ioverlap, dt_sample, ntrans, nsection)
  endif

  if(ireadkp.eq.1.and.num_kplane.gt.0) then  ! Using constant-k plane time series

    if(icalPSD.eq.1) then
      if(ireadVOL.eq.0) then
        call InitTSkplane(ntpoint_total, npath, DNSIndex_k, ibe_k, iend_k, jbe_k, jend_k, &
                          ixwindowl, ixwindowr, nvar_output, varindex_output, fileprop)
      else
        call InitVOLTSkplane(ntpoint_total, npath, DNSIndex_k, ibe_k, iend_k, jbe_k, jend_k, &
                             ixwindowl, ixwindowr, nvar_output, varindex_output, fileprop)
      endif
    endif
    if(icalcorrxy.eq.1.or.icalcorrxt.eq.1) then
      call InitVOLTSkplane_NoPrintInfo(ntpoint_total, npath, DNSIndex_k, ibe_k, ibe_k, jbe_k, jend_k, &
                           ixwindowl, ixwindowr, nvar_output, varindex_output, fileprop)
    endif

    if(icalcorrxy.eq.1) then
!      call InitVOLTSkplane(ntpoint_total, npath, DNSIndex_k, ibe_k, ibe_k, jbe_k, jend_k, &
!                           ixwindowl, ixwindowr, nvar_output, varindex_output, fileprop)
      call InitFFT1D(nypoint_k)
      call InitCorry(nypoint_k, dy_kp)
    endif
    allocate( buffer_kplane(ntpoint_total, nypoint_k, imints_k:imaxts_k, 1, nvar_output) )
    allocate( buffertmp(fileprop(1)%ntpoint, nypoint_k, imints_k:imaxts_k, 1, nvar_output) )

    ! Calculate x-y correlation
    if(icalcorrxy.eq.1) then

      num_file = (fileprop(1)%file_end-fileprop(1)%file_be)/fileprop(1)%file_skip + 1

      do k = 1, num_kplane
        write(unit=fnum1,fmt='(I04.4)') kplane(k)
          ! write(unit=fnum2,fmt='(I04.4)') iii
          ! write(unit=fnum3,fmt='(I04.4)') iii
        print *, 'wall-normal location k =', kplane(k)
!        print *, 'streamwise average range i =', iii, iii

        do n=buffer_MulVol+1, num_file-buffer_MulVol

          do iii = ibe_k, iend_k
            if(iii-ibe_k.eq.0) then
              call InitVOLTSkplane(ntpoint_total, npath, DNSIndex_k, iii, iii, jbe_k, jend_k, &
                                   ixwindowl, ixwindowr, nvar_output, varindex_output, fileprop)
            else
              call InitVOLTSkplane_NoPrintInfo(ntpoint_total, npath, DNSIndex_k, iii, iii, jbe_k, jend_k, &
                                               ixwindowl, ixwindowr, nvar_output, varindex_output, fileprop)
            endif
            ntpoint_tmp = 0
            do nn=-buffer_MulVol, buffer_MulVol
              write(unit=fnum_8,fmt='(I08.8)') fileprop(1)%file_be + (n+nn-1)*fileprop(1)%file_skip
              filename = trim(fileprop(1)%filepath)//'timeseries_'//fnum_8//'.h5'
              print *, 'Reading file: ', trim(filename)
              call ReadVOLTSkplane_perFile(filename,k+kplane_be-1,buffertmp) !!!!!!!!!
              buffer_kplane((ntpoint_tmp+1):(ntpoint_tmp+fileprop(1)%ntpoint),1:nypoint_k,imints_k:imaxts_k,1,1:nvar_output) = &
              buffertmp(1:fileprop(1)%ntpoint,1:nypoint_k,imints_k:imaxts_k,1,1:nvar_output)
              ntpoint_tmp = ntpoint_tmp + fileprop(1)%ntpoint
            enddo ! end nn loop

            allocate(spect1dtmp3(-ixwindowl:ixwindowr,nypoint_k), spect1dave3(-ixwindowl:ixwindowr,nypoint_k))
            if(.not.allocated(numpt)) allocate(numpt(-ixwindowl:ixwindowr))
            do m = 1, nvar_output
              ! fileoutput = '_k'//fnum1//'_i'//fnum2//'-'//fnum3//'_'//trim(varout_k(m))//'_'//fnum_8
              fileoutput = '_k'//fnum1//'_'//trim(varout_k(m))//'_'//fnum_8
              numpt = 0
              spect1dtmp3 = (0.,0.); spect1dave3 = (0.,0.)
              do i = 1, nxpoint_k
                do j = 1, ntpoint_total
                  do ii = -min(ixwindowl,i-imints_k), min(ixwindowr,imaxts_k-i)
                    ! print *, 'i = ', i, 'j = ', j, 'ii = ', ii
                    spect1dtmp3(ii,:) = crossspect1dy(buffer_kplane(j,1:nypoint_k,i+ii,1,m),buffer_kplane(j,1:nypoint_k,i,1,m))
                    spect1dave3(ii,:) = spect1dave3(ii,:) + spect1dtmp3(ii,:)
                    numpt(ii) = numpt(ii) + 1
                  enddo
                enddo
              enddo
              forall(ii=-ixwindowl:ixwindowr,numpt(ii).gt.0)
                spect1dave3(ii,:) = spect1dave3(ii,:)/dble(numpt(ii))
              end forall
              call CalCorrxy_HDF5(ixwindowl,ixwindowr,iii,dx_kp,spect1dave3,fileoutput,n,iii-ibe_k)
            enddo ! end m loop
            deallocate(spect1dtmp3, spect1dave3)

          enddo ! end iii loop
        enddo ! end n loop
      enddo ! end k loop
    endif ! end icalcorrxy.eq.1


    ! Calculate x-t correlation
    if(icalcorrxt.eq.1) then
    ! Reading k plane data
      num_file = (fileprop(1)%file_end-fileprop(1)%file_be)/fileprop(1)%file_skip + 1
      do k = 1, num_kplane
        write(unit=fnum1,fmt='(I04.4)') kplane(k)
!        write(unit=fnum2,fmt='(I04.4)') ibe_k
!        write(unit=fnum3,fmt='(I04.4)') iend_k
        print *, 'wall-normal location k =', kplane(k)
!        print *, 'streamwise average range i =', ibe_k, iend_k

        do n=buffer_MulVol+1, num_file-buffer_MulVol

          do iii = ibe_k, iend_k
            if(iii-ibe_k.eq.0) then
              call InitVOLTSkplane(ntpoint_total, npath, DNSIndex_k, iii, iii, jbe_k, jend_k, &
                                 ixwindowl, ixwindowr, nvar_output, varindex_output, fileprop)
            else
              call InitVOLTSkplane_NoPrintInfo(ntpoint_total, npath, DNSIndex_k, iii, iii, jbe_k, jend_k, &
                                               ixwindowl, ixwindowr, nvar_output, varindex_output, fileprop)
            endif
            ntpoint_tmp = 0
            do nn=-buffer_MulVol, buffer_MulVol
              write(unit=fnum_8,fmt='(I08.8)') fileprop(1)%file_be + (n+nn-1)*fileprop(1)%file_skip
              filename = trim(fileprop(1)%filepath)//'timeseries_'//fnum_8//'.h5'
              print *, 'Reading file: ', trim(filename)
              call ReadVOLTSkplane_perFile(filename,k+kplane_be-1,buffertmp) !!!!!!!!!
              buffer_kplane((ntpoint_tmp+1):(ntpoint_tmp+fileprop(1)%ntpoint),1:nypoint_k,imints_k:imaxts_k,1,1:nvar_output) = &
              buffertmp(1:fileprop(1)%ntpoint,1:nypoint_k,imints_k:imaxts_k,1,1:nvar_output)
              ntpoint_tmp = ntpoint_tmp + fileprop(1)%ntpoint
            enddo ! end nn loop

            write(unit=fnum_8,fmt='(I08.8)') fileprop(1)%file_be + (n-1)*fileprop(1)%file_skip
            allocate(spect1dtmp2(ntrans,-ixwindowl:ixwindowr), spect1dave2(ntrans,-ixwindowl:ixwindowr))
            if(.not.allocated(numpt)) allocate(numpt(-ixwindowl:ixwindowr))
            do m = 1, nvar_output
              ! fileoutput = '_k'//fnum1//'_i'//fnum2//'-'//fnum3//'_'//trim(varout_k(m))//'_'//fnum_8
              fileoutput = '_k'//fnum1//'_'//trim(varout_k(m))//'_'//fnum_8
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
              ! call CalCorrxt(ixwindowl,ixwindowr,dx_kp,spect1dave2,fileoutput)
              call CalCorrxt_HDF5(ixwindowl,ixwindowr,iii,dx_kp,spect1dave2,fileoutput,n-buffer_MulVol,iii-ibe_k)
            enddo
            deallocate(spect1dtmp2, spect1dave2)
          enddo ! end iii loop
        enddo ! end n loop
      enddo  ! end looping over k planes
    endif ! end icalcorrxt.eq.1


    ! Calculating PSD
    if(icalPSD.eq.1) then

    ! Reading k plane data
      do k = 1, num_kplane
        write(unit=fnum1,fmt='(I04.4)') kplane(k)
        write(unit=fnum2,fmt='(I04.4)') ibe_k
        write(unit=fnum3,fmt='(I04.4)') iend_k
        if(ireadVOL.eq.0) then
          call ReadTSkplane(kplane(k),buffer_kplane)
        else
          call ReadVOLTSkplane(k+kplane_be-1,kplane(k),buffer_kplane)  !!!!!
        endif
        print *, 'wall-normal location k =', kplane(k)
        print *, 'streamwise average range i =', ibe_k, iend_k

        allocate(spect1dtmp(ntrans), spect1dave(ntrans))
        do m = 1, nvar_output
          spect1dave = 0.d0;  numpt1 = 0
          tmean = 0.d0; tmeanenergy = 0.d0
          fileoutput = '_k'//fnum1//'_i'//fnum2//'-'//fnum3//'_'//trim(varout_k(m))
          print *, 'fileoutput =', trim(fileoutput)
          do i = 1, nxpoint_k
          do j = 1, nypoint_k
            tsdata1d(1:ntpoint_total) = buffer_kplane(1:ntpoint_total,j,i,1,m)
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
     !              print *, 'spect1dave =', real(sum(spect1dave))
          print *, 'tmean =', tmean, 'tmeanenergy =', tmeanenergy
          call CalPowerspect(spect1dave,tmeanenergy,fileoutput)
          call CalCorrt(spect1dave,fileoutput)
        enddo ! end m loop
        deallocate(spect1dtmp, spect1dave)
      enddo ! end k loop
    endif


    deallocate(buffer_kplane,buffertmp)
  endif  ! end calculating const-k plane


! Reading k plane data
!     do k = 1, num_kplane
!       write(unit=fnum1,fmt='(I04.4)') kplane(k)
!       write(unit=fnum2,fmt='(I04.4)') ibe_k
!       write(unit=fnum3,fmt='(I04.4)') iend_k
!       print *, 'wall-normal location k =', kplane(k)
!       print *, 'streamwise average range i =', ibe_k, iend_k
!
!
!      num_file = (fileprop(1)%file_end-fileprop(1)%file_be)/fileprop(1)%file_skip + 1
!
!      do n=buffer_MulVol+1, num_file-buffer_MulVol
!
!        ntpoint_tmp = 0
!        do nn=-buffer_MulVol, buffer_MulVol
!          write(unit=fnum_8,fmt='(I08.8)') fileprop(1)%file_be + (n+nn-1)*fileprop(1)%file_skip
!          filename = trim(fileprop(1)%filepath)//'timeseries_'//fnum_8//'.h5'
!          print *, 'Reading file: ', trim(filename)
!          call ReadVOLTSkplane_perFile(filename,k+kplane_be-1,buffertmp) !!!!!!!!!
!          buffer_kplane((ntpoint_tmp+1):(ntpoint_tmp+fileprop(1)%ntpoint),1:nypoint_k,imints_k:imaxts_k,1,1:nvar_output) = &
!          buffertmp(1:fileprop(1)%ntpoint,1:nypoint_k,imints_k:imaxts_k,1,1:nvar_output)
!          ntpoint_tmp = ntpoint_tmp + fileprop(1)%ntpoint
!        enddo ! end nn loop



!       if(ireadVOL.eq.0) then
!         call ReadTSkplane(kplane(k),buffer_kplane)
!       elseif(ireadVOL.eq.1.and.iSeparate.eq.0) then
!         call ReadVOLTSkplane(k,kplane(k),buffer_kplane)


!      ! Calculate x-y correlation
!      if(icalcorrxy.eq.1) then
!         allocate(spect1dtmp3(-ixwindowl:ixwindowr,nypoint_k), spect1dave3(-ixwindowl:ixwindowr,nypoint_k))
!         if(.not.allocated(numpt)) allocate(numpt(-ixwindowl:ixwindowr))
!         do m = 1, nvar_output
!            fileoutput = '_k'//fnum1//'_i'//fnum2//'-'//fnum3//'_'//trim(varout_k(m))
!            numpt = 0
!            spect1dtmp3 = (0.,0.); spect1dave3 = (0.,0.)
!            do i = 1, nxpoint_k
!            do j = 1, ntpoint_total
!               do ii = -min(ixwindowl,i-imints_k), min(ixwindowr,imaxts_k-i)
!
!         print *, 'i = ', i, 'j = ', j, 'ii = ',ii
!                  spect1dtmp3(ii,:) = crossspect1dy(buffer_kplane(j,1:nypoint_k,i+ii,1,m),buffer_kplane(j,1:nypoint_k,i,1,m))
!                  spect1dave3(ii,:) = spect1dave3(ii,:) + spect1dtmp3(ii,:)
!                  numpt(ii) = numpt(ii) + 1
!               enddo
!            enddo
!            enddo
!            forall(ii=-ixwindowl:ixwindowr,numpt(ii).gt.0)
!               spect1dave3(ii,:) = spect1dave3(ii,:)/dble(numpt(ii))
!            end forall
!            call CalCorrxy(ixwindowl,ixwindowr,dx_kp,spect1dave3,fileoutput)
!         enddo ! end m loop
!         deallocate(spect1dtmp3, spect1dave3)
!
!      endif ! end icalcorrxy.eq.1



!        ! Calculate x-t correlation
!        if(icalcorrxt.eq.1) then
!         allocate(spect1dtmp2(ntrans,-ixwindowl:ixwindowr), spect1dave2(ntrans,-ixwindowl:ixwindowr))
!         if(.not.allocated(numpt)) allocate(numpt(-ixwindowl:ixwindowr))
!         do m = 1, nvar_output
!            fileoutput = '_k'//fnum1//'_i'//fnum2//'-'//fnum3//'_'//trim(varout_k(m))
!            numpt = 0
!            spect1dtmp2 = (0.,0.); spect1dave2 = (0.,0.)
!            do i = 1, nxpoint_k
!            do j = 1, nypoint_k
!               do ii = -min(ixwindowl,i-imints_k), min(ixwindowr,imaxts_k-i)
!                  spect1dtmp2(:,ii) = crossspect1d(buffer_kplane(1:ntpoint_total,j,i+ii,1,m),buffer_kplane(1:ntpoint_total,j,i,1,m))
!                  spect1dave2(:,ii) = spect1dave2(:,ii) + spect1dtmp2(:,ii)
!                  numpt(ii) = numpt(ii) + 1
!               enddo
!            enddo
!            enddo
!            forall(ii=-ixwindowl:ixwindowr,numpt(ii).gt.0)
!               spect1dave2(:,ii) = spect1dave2(:,ii)/dble(numpt(ii))
!            end forall
!            call CalCorrxt(ixwindowl,ixwindowr,dx_kp,spect1dave2,fileoutput)
!         enddo
!         deallocate(spect1dtmp2, spect1dave2)
!      endif






contains
    subroutine Input()
      integer, parameter :: nid=5
      integer :: i, n
!      if (nid.eq.5) then
!        print*,'please input parameters or redirect from a file'
!      else
!        open(nid,file='spectcorr.inp',status='old')
!      end if
      read(*,*)
      read(*,*) npath, ireadVOL, iAverage, buffer_MulVol  !, ireadMulVol
      allocate(fileprop(npath))
      read(*,*)
      do n=1, npath
        read(*,*)
        read(*,*) fileprop(n)%ntpoint, fileprop(n)%nskip
        read(*,*)
        read(*,*) fileprop(n)%file_be, fileprop(n)%file_end, fileprop(n)%file_skip
        read(*,*)
        read(*,'(a)') fileprop(n)%filepath
      enddo
      read(*,*)
      read(*,*)
      read(*,*)icalPSD, icalcorrxt, icalcorrxy, icalphasespeed
      read(*,*)
      read(*,*) nperiod, ioverlap, iwindow, dt_sample
      read(*,*)
      read(*,*)
      read(*,*) nvar_output
      allocate(varindex_output(nvar_output))
      read(*,*) (varindex_output(i), i = 1, nvar_output)
      read(*,*)
      read(*,*)
      read(*,*) ireadkp
      read(*,*)
      read(*,*) ibe_k, iend_k, jbe_k, jend_k, ixwindowl, ixwindowr, dx_kp, dy_kp
      read(*,*)
      read(*,*) DNSIndex_k%ibe, DNSIndex_k%iend, DNSIndex_k%iskip, &
                DNSIndex_k%jbe, DNSIndex_k%jend, DNSIndex_k%jskip, DNSIndex_k%ibuffer
      read(*,*)
      read(*,*) num_kplane, kplane_be, kplane_end
      allocate(kplane(num_kplane))
      read(*,*) (kplane(i), i = 1, num_kplane)
      if (nid.ne.5) close(nid)


      if(ireadVOL.eq.0) buffer_MulVol = 0 !!!!!!!!!!!!!

      ntpoint_total = 0
      if(ireadVOL.eq.0) then
        do n=1, npath
          fileprop(n)%ntpoint = fileprop(n)%ntpoint - fileprop(n)%nskip
          ntpoint_total = ntpoint_total + fileprop(n)%ntpoint
        enddo
      elseif(ireadVOL.eq.1.and.iAverage.eq.0) then
        do n=1, npath
          ntpoint_total = ntpoint_total + ( (fileprop(n)%file_end-fileprop(n)%file_be)/fileprop(n)%file_skip + 1 )*fileprop(n)%ntpoint
        enddo
      endif

      if(ireadVOL.eq.1.and.iAverage.eq.1.and.icalcorrxy.eq.1) then
        ntpoint_total = fileprop(1)%ntpoint
        print *, 'ntpoint_total = ', ntpoint_total
        buffer_MulVol = 0
      endif

      if(ireadVOL.eq.1.and.iAverage.eq.1.and.icalcorrxt.eq.1) then
        if(buffer_MulVol.eq.0) then
          print *, 'please set buffer_MulVol within none 0 value'
          stop
        endif

        ntpoint_total = fileprop(1)%ntpoint*(2*buffer_MulVol+1)
        nperiod = buffer_MulVol*fileprop(1)%ntpoint
        print *, 'ntpoint_total = ', ntpoint_total, 'nperiod = ', nperiod

      endif

      if(icalcorrxt.eq.0.and.icalcorrxy.eq.0) then
         ixwindowl = 0
         ixwindowr = 0
      endif
      print *, 'ixwindowl =', ixwindowl, 'ixwindowr =', ixwindowr
      if (ioverlap.ne.0.and.ioverlap.ne.1) then
        print*,'ioverlap can ONLY be 0 or 1'
        stop
      end if

      if(num_kplane.ne.(kplane_end-kplane_be+1)) then
        print *, 'num_kplane = ', 'kplane_end - kplane_be + 1 = ', kplane_end-kplane_be+1
        print *, 'num_kplane should be equal to kplane_end - kplane_be + 1'
        stop
      endif

      if(ireadkp.eq.1) then
         print *, 'number of wall-normal locations =', num_kplane
         print *, 'Wall-normal plane indexes k =', kplane
      endif



!      if(ireadip.eq.1) then
!         print *, 'number of const-i locations =', num_iplane
!         print *, 'const-i plane indexes i =', iplane
!      endif
!      write(*,*)
    end subroutine input


end program spectcorr
