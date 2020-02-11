! program to compute coherence using timeseries data
program Calcoh
  use MTSKplane
  use modSpect
  use modSpecty
  implicit none
  integer :: nskip, kloc_ref, ixwindowl, ixwindowr, jwindow , icylin
  integer :: nperiod, ioverlap, ntrans, nsection
  integer :: num_kplane, iwindow
  integer :: kplane_be, kplane_end, kplane_skip
  integer, dimension(:), allocatable :: kplane
  integer :: n, ibe_k, iend_k
  integer :: jbe_k, jend_k
  type(fprop), dimension(:), allocatable :: fileprop
  character(400) :: fileinput, fileoutput, fname
  character(4) :: fnum4, fnum4_1, fnum4_2
  integer, parameter :: nvar_output = 1
  integer :: varindex_output(nvar_output)
  character(3) :: varname_output(nvar_output)
  integer :: ii, i, j, k, kk, nn, m, jj, jjj
  real(8) :: dt_sample, dx_kp, dy_kp
  real(8), dimension(:,:,:,:,:), allocatable :: tsdata, tsdata_ref
  real(8), dimension(:), allocatable :: coh2, phase
  complex, dimension(:,:,:), allocatable :: spect12, spect11, spect22
  complex, dimension(:), allocatable ::  corrxtmp
  integer, dimension(:,:), allocatable ::  numpt
  real(8) :: fperiod, corrx
  integer :: numptx
  real(8) :: zloc_ref, dzdk_ref, zloc_kp, dzdk_kp
  integer :: ntpoint_total
  real(8), dimension(:,:,:,:,:), allocatable :: buffer_kplane
  type(tp_DNSIndex) :: DNSIndex_k
  integer :: npath
  integer :: i2D_coh, ix_coh, iy_coh, iAve_data, fdim1, fdim2, fdim3, jpointAve_be, jpointAve_end
  !
  integer :: ilocal_ave, ixw, jyw, iradius_pt, jradius_pt, ixpt, jypt, numpt1, iidx, jidx
  real(8) :: iradius, iradius_tmp, std, dx_tmp
  integer :: nxpoint_k_local, nypoint_k_local, imints_k_local, imaxts_k_local, ixwindowl_local, ixwindowr_local
  real(8), dimension(:,:,:,:,:), allocatable :: tsdata_tmp, tsdata_ref_tmp
  real(8), dimension(:), allocatable :: tsdata1d

  call InitHDF5()
  call Input()

  if(iAve_data.eq.1) then
    call CalAve()
  else
    varindex_output(1)=4 ! '/p/'
    call InitSpect1d(ntpoint_total, nsection, iwindow, ioverlap, dt_sample, ntrans, nperiod)
    fperiod = (ntrans-1)*dt_sample

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

    if(jpointAve_end.gt.nypoint_k) then
      print *, 'STOP! jpointAve_end should be less than or equal to ', nypoint_k
      stop
    endif
    if(jwindow.gt.nypoint_k/2) then
      print *, 'STOP! jwindow should be less than or equal to ', nypoint_k/2
      stop
    endif

    allocate(tsdata_ref(ntpoint_total,nypoint_k,imints_k:imaxts_k,1,1))
    allocate(tsdata(ntpoint_total,nypoint_k,imints_k:imaxts_k,1,nvar_output) )

    write(unit=fnum4_1,fmt='(I04.4)') kloc_ref
    print *, 'reference wall-normal location k =', kloc_ref

    if(i2D_coh.eq.1) then
      !jwindow = nypoint_k/2
      allocate(coh2(ntrans),phase(ntrans))
      allocate(spect12(ntrans,-ixwindowl:ixwindowr,0:jwindow),spect11(ntrans,-ixwindowl:ixwindowr,0:jwindow),spect22(ntrans,-ixwindowl:ixwindowr,0:jwindow))
      allocate(corrxtmp(ntrans),numpt(-ixwindowl:ixwindowr,0:jwindow))
    elseif(ix_coh.eq.1) then
      allocate(coh2(ntrans),phase(ntrans))
      allocate(spect12(ntrans,-ixwindowl:ixwindowr,1),spect11(ntrans,-ixwindowl:ixwindowr,1),spect22(ntrans,-ixwindowl:ixwindowr,1))
      allocate(corrxtmp(ntrans),numpt(-ixwindowl:ixwindowr,1))
    endif

    if(iFileType.eq.0) then
      ! Read reference plane data
      call ReadTSkplane(kloc_ref,tsdata_ref)

      if(ix_coh.eq.1) then
        do k = 1, num_kplane
          write(unit=fnum4,fmt='(I04.4)') kplane(k)
          print *, 'wall-normal location k =', kplane(k)
          call ReadTSkplane(kplane(k),tsdata)

          !do j = 1, nypoint_k
          do j=jpointAve_be, jpointAve_end !nypoint_k
            print *, 'j = ', j
            numpt = 0
            spect11 = (0.d0,0.d0)
            spect22 = (0.d0,0.d0)
            spect12 = (0.d0,0.d0)
            do i = 1, nxpoint_k
              do ii = max(-ixwindowl,imints_k-i), min(ixwindowr,imaxts_k-i)
                spect12(:,ii,1) = spect12(:,ii,1)+crossspect1d(tsdata_ref(1:ntpoint_total,j,i,1,1),tsdata(1:ntpoint_total,j,i+ii,1,1))
                spect11(:,ii,1) = spect11(:,ii,1)+autospect1d(tsdata_ref(1:ntpoint_total,j,i,1,1))
                spect22(:,ii,1) = spect22(:,ii,1)+autospect1d(tsdata(1:ntpoint_total,j,i+ii,1,1))
                numpt(ii,1) = numpt(ii,1) + 1
              enddo
            enddo

            write(unit=fnum4_2,fmt='(I04.4)') j
            fname = 'coherence_dx_kref'//fnum4_1//'_j'//fnum4_2//'.dat'
            open(unit=12,file=trim(fname),status='unknown')
            rewind(12)
            write(12,'(a)') 'variables = f,x,y,spect12_r,spect12_i,spect11,spect22'

            write(12,*) 'zone T = k'//fnum4//', i =', ntrans/2, ' j =', (ixwindowl + ixwindowr + 1), ' k =', 1
            do ii = -ixwindowl,ixwindowr
              if(numpt(ii,1).gt.0) then
                !coh2 = abs(spect12(:,ii,1))**2/(abs(spect11(:,ii,1))*abs(spect22(:,ii,1))+1.e-30)
                !phase = atan2(imag(spect12(:,ii,1)),real(spect12(:,ii,1)))
                do i = 1,ntrans/2
                  write(12,*)  dble(i-1)/fperiod, dble(ii)*dx_kp,0, real(spect12(i,ii,1)),imag(spect12(i,ii,1)), abs(spect11(i,ii,1)),abs(spect22(i,ii,1))
                enddo
              endif
            enddo
            close(12)
          enddo ! end j loop
        enddo  ! end k loop

      elseif(i2D_coh.eq.1) then

        do k = 1, num_kplane
          ! Reading k plane data
          write(unit=fnum4,fmt='(I04.4)') kplane(k)
          print *, 'wall-normal location k =', kplane(k)
          call ReadTSkplane(kplane(k),tsdata)

          !do j = 1, nypoint_k
          do j=jpointAve_be, jpointAve_end !nypoint_k
            numpt = 0
            spect11 = (0.d0,0.d0)
            spect22 = (0.d0,0.d0)
            spect12 = (0.d0,0.d0)
            do i = 1, nxpoint_k
              do jj = 0,jwindow
                jjj = modulo(j+jj,nypoint_k)
                if(jjj.eq.0) jjj = nypoint_k
                do ii = max(-ixwindowl,imints_k-i), min(ixwindowr,imaxts_k-i)
                  spect12(:,ii,jj) = spect12(:,ii,jj)+crossspect1d(tsdata_ref(1:ntpoint_total,j,i,1,1),tsdata(1:ntpoint_total,jjj,i+ii,1,1))
                  spect11(:,ii,jj) = spect11(:,ii,jj)+autospect1d(tsdata_ref(1:ntpoint_total,j,i,1,1))
                  spect22(:,ii,jj) = spect22(:,ii,jj)+autospect1d(tsdata(1:ntpoint_total,jjj,i+ii,1,1))
                  numpt(ii,jj) = numpt(ii,jj) + 1
                enddo
              enddo
            enddo

            write(unit=fnum4_2,fmt='(I04.4)') j
            open(unit=12,file='coherence_kref'//fnum4_1//'_j'//fnum4_2//'.dat')
            rewind(12)
            write(12,'(a)') 'variables = f,x,y,spect12_r,spect12_i,spect11,spect22'

            numptx = count(numpt(:,0).gt.0)
            write(12,*) 'zone T = k'//fnum4//', i =', ntrans/2, ' j =', numptx, ' k =', jwindow+1
            do jj = 0, jwindow
              do ii = -ixwindowl,ixwindowr
                if(numpt(ii,jj).gt.0) then
                  !coh2 = abs(spect12(:,ii,jj))**2/(abs(spect11(:,ii,jj))*abs(spect22(:,ii,jj))+1.e-30)
                  !phase = atan2(imag(spect12(:,ii,jj)),real(spect12(:,ii,jj)))
                  do i = 1,ntrans/2
                    write(12,*)  dble(i-1)/fperiod, dble(ii)*dx_kp, dble(jj)*dy_kp, real(spect12(i,ii,jj)),imag(spect12(i,ii,jj)), abs(spect11(i,ii,jj)),abs(spect22(i,ii,jj))
                  enddo
                 endif
              enddo
            enddo
            close(12)
          enddo ! end j loop
        enddo  ! end k loop

      endif

    else ! iFileType.eq.1

      ! Read reference plane data
      do kk=1, DNSIndex_k%nkplane
        if(kloc_ref.eq.DNSIndex_k%kplane(kk)) then
          kplane_be = kk
        endif
      enddo ! end kk loop
      call ReadVOLTSkplane(kplane_be,kloc_ref,tsdata_ref)

      !print *, 'tsdata_ref = ', tsdata_ref(:,1,1,1,1)
    if(ilocal_ave.eq.0) then

      if(ix_coh.eq.1) then
        do k = 1, num_kplane
          write(unit=fnum4,fmt='(I04.4)') kplane(k)
          print *, 'wall-normal location k =', kplane(k)
          ! Read k-plane data
          do kk=1, DNSIndex_k%nkplane
          if(kplane(k).eq.DNSIndex_k%kplane(kk)) then
            kplane_be = kk
          endif
          enddo ! end kk loop
          call ReadVOLTSkplane(kplane_be,kplane(k),tsdata)

          do j=jpointAve_be, jpointAve_end !nypoint_k
            print *, 'j = ', j
            numpt = 0
            spect11 = (0.d0,0.d0)
            spect22 = (0.d0,0.d0)
            spect12 = (0.d0,0.d0)
            do i = 1, nxpoint_k
              do ii = max(-ixwindowl,imints_k-i), min(ixwindowr,imaxts_k-i)
                spect12(:,ii,1) = spect12(:,ii,1)+crossspect1d(tsdata_ref(1:ntpoint_total,j,i,1,1),tsdata(1:ntpoint_total,j,i+ii,1,1))
                spect11(:,ii,1) = spect11(:,ii,1)+autospect1d(tsdata_ref(1:ntpoint_total,j,i,1,1))
                spect22(:,ii,1) = spect22(:,ii,1)+autospect1d(tsdata(1:ntpoint_total,j,i+ii,1,1))
                numpt(ii,1) = numpt(ii,1) + 1
              enddo
            enddo

            write(unit=fnum4_2,fmt='(I04.4)') j
            fname = 'coherence_dx_kref'//fnum4_1//'_j'//fnum4_2//'.dat'
            open(unit=12,file=trim(fname),status='unknown')
            rewind(12)
            write(12,'(a)') 'variables = f,x,y,spect12_r,spect12_i,spect11,spect22'

            write(12,*) 'zone T = k'//fnum4//', i =', ntrans/2, ' j =', (ixwindowl + ixwindowr + 1), ' k =', 1
            do ii = -ixwindowl,ixwindowr
              if(numpt(ii,1).gt.0) then
                !coh2 = abs(spect12(:,ii,1))**2/(abs(spect11(:,ii,1))*abs(spect22(:,ii,1))+1.e-30)
                !phase = atan2(imag(spect12(:,ii,1)),real(spect12(:,ii,1)))
                do i = 1,ntrans/2
                  write(12,*)  dble(i-1)/fperiod, dble(ii)*dx_kp,0, real(spect12(i,ii,1)),imag(spect12(i,ii,1)), abs(spect11(i,ii,1)),abs(spect22(i,ii,1))
                enddo
              endif
            enddo
            close(12)
          enddo ! end j loop

        enddo ! end k loop

      elseif(i2D_coh.eq.1) then

        do k = 1, num_kplane
          write(unit=fnum4,fmt='(I04.4)') kplane(k)
          print *, 'wall-normal location k =', kplane(k)
          ! Read k-plane data
          do kk=1, DNSIndex_k%nkplane
          if(kplane(k).eq.DNSIndex_k%kplane(kk)) then
            kplane_be = kk
          endif
          enddo ! end kk loop
          call ReadVOLTSkplane(kplane_be,kplane(k),tsdata)

          do j=jpointAve_be, jpointAve_end !nypoint_k
            numpt = 0
            spect11 = (0.d0,0.d0)
            spect22 = (0.d0,0.d0)
            spect12 = (0.d0,0.d0)
            do i = 1, nxpoint_k
              do jj = 0,jwindow
                jjj = modulo(j+jj,nypoint_k)
                if(jjj.eq.0) jjj = nypoint_k
                do ii = max(-ixwindowl,imints_k-i), min(ixwindowr,imaxts_k-i)
                  spect12(:,ii,jj) = spect12(:,ii,jj)+crossspect1d(tsdata_ref(1:ntpoint_total,j,i,1,1),tsdata(1:ntpoint_total,jjj,i+ii,1,1))
                  spect11(:,ii,jj) = spect11(:,ii,jj)+autospect1d(tsdata_ref(1:ntpoint_total,j,i,1,1))
                  spect22(:,ii,jj) = spect22(:,ii,jj)+autospect1d(tsdata(1:ntpoint_total,jjj,i+ii,1,1))
                  numpt(ii,jj) = numpt(ii,jj) + 1
                enddo
              enddo
            enddo

            write(unit=fnum4_2,fmt='(I04.4)') j
            open(unit=12,file='coherence_kref'//fnum4_1//'_j'//fnum4_2//'.dat')
            rewind(12)
            write(12,'(a)') 'variables = f,x,y,spect12_r,spect12_i,spect11,spect22'

            numptx = count(numpt(:,0).gt.0)
            write(12,*) 'zone T = k'//fnum4//', i =', ntrans/2, ' j =', numptx, ' k =', jwindow+1
            do jj = 0, jwindow
              do ii = -ixwindowl,ixwindowr
                if(numpt(ii,jj).gt.0) then
                  !coh2 = abs(spect12(:,ii,jj))**2/(abs(spect11(:,ii,jj))*abs(spect22(:,ii,jj))+1.e-30)
                  !phase = atan2(imag(spect12(:,ii,jj)),real(spect12(:,ii,jj)))
                  do i = 1,ntrans/2
                    write(12,*)  dble(i-1)/fperiod, dble(ii)*dx_kp, dble(jj)*dy_kp, real(spect12(i,ii,jj)),imag(spect12(i,ii,jj)), abs(spect11(i,ii,jj)),abs(spect22(i,ii,jj))
                  enddo
                 endif
              enddo
            enddo
            close(12)
          enddo ! end j loop
        enddo  ! end k loop

      endif

    else ! ilocal_ave.gt.0

      nxpoint_k_local = nxpoint_k/ixpt
      nypoint_k_local = nypoint_k/jypt
      ixwindowl_local = ixwindowl/ixpt
      ixwindowr_local = ixwindowr/ixpt
      imaxts_k_local = nxpoint_k_local + ixwindowr_local
      imints_k_local = -(ixwindowl_local - 1)

      print *, '  '
      print *, 'nxpoint_local = ', nxpoint_k_local
      print *, 'nypoint_local = ', nypoint_k_local
      print *, 'ixwindowl_local = ', ixwindowl_local
      print *, 'ixwindowr_local = ', ixwindowr_local
      print *, 'imints_local = ', imints_k_local
      print *, 'imaxts_local = ', imaxts_k_local

!stop
      !allocate(tsdata_ref(ntpoint_total,nypoint_k,imints_k:imaxts_k,1,1))
      !allocate(tsdata(ntpoint_total,nypoint_k,imints_k:imaxts_k,1,nvar_output) )
      allocate(tsdata_ref_tmp(ntpoint_total,nypoint_k_local,imints_k_local:imaxts_k_local,1,1))
      allocate(tsdata_tmp(ntpoint_total,nypoint_k_local,imints_k_local:imaxts_k_local,1,nvar_output) )
      allocate(tsdata1d(1:ntpoint_total))

!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if(ilocal_ave.eq.1) then
        print *, 'Doing circle local avarage for the reference data ... '

        do i=1, nxpoint_k_local
          do j=1, nypoint_k_local
            numpt1 = 0; tsdata1d = 0;
            iidx = iradius_pt + (i-1)*ixpt + 1
            jidx = jradius_pt + (j-1)*jypt + 1
            do ii=-iradius_pt, iradius_pt
              do jj=-jradius_pt, jradius_pt
                iradius_tmp = sqrt( (dble(ii)*dx_kp)**2 + (dble(jj)*dy_kp)**2 )
                if(iradius_tmp.le.iradius) then
                  !if(i.eq.1.and.j.eq.2) then
                  !  print *, 'iidx+ii = ', iidx+ii, 'jidx+jj = ', jidx+jj
                  !  print *, 'tsdata_ref(1,jidx+jj,iidx+ii,1,1) = ', tsdata_ref(1,jidx+jj,iidx+ii,1,1)
                  !endif
                  tsdata1d(1:ntpoint_total) = tsdata1d(1:ntpoint_total) + tsdata_ref(1:ntpoint_total,jidx+jj,iidx+ii,1,1)
                  numpt1 = numpt1 + 1
                endif
              enddo
            enddo
            tsdata_ref_tmp(1:ntpoint_total,j,i,1,1) = tsdata1d(1:ntpoint_total)/dble(numpt1)
          enddo
        enddo
        print *, 'numpt1 = ', numpt1
      elseif(ilocal_ave.eq.2) then
        print *, 'Doing square local avarage ... '


      endif
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     ! print *, 'p = ', tsdata_ref(1,1,1,1,1)
     ! print *, '', tsdata_ref_tmp(1:ntpoint_total,1,2,1,1)
!stop

      if(ix_coh.eq.1) then
        do k = 1, num_kplane
          write(unit=fnum4,fmt='(I04.4)') kplane(k)
          print *, 'wall-normal location k =', kplane(k)
          ! Read k-plane data
          do kk=1, DNSIndex_k%nkplane
          if(kplane(k).eq.DNSIndex_k%kplane(kk)) then
            kplane_be = kk
          endif
          enddo ! end kk loop
          call ReadVOLTSkplane(kplane_be,kplane(k),tsdata)

!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          if(ilocal_ave.eq.1) then
            print *, 'Doing circle local avarage ... '
            do i=1, nxpoint_k_local
              do j=1, nypoint_k_local
                numpt1 = 0; tsdata1d = 0;
                iidx = iradius_pt + (i-1)*ixpt + 1
                jidx = jradius_pt + (j-1)*jypt + 1
                do ii=-iradius_pt, iradius_pt
                  do jj=-jradius_pt, jradius_pt
                    iradius_tmp = sqrt( (dble(ii)*dx_kp)**2 + (dble(jj)*dy_kp)**2 )
                    if(iradius_tmp.le.iradius) then
                      !if(i.eq.1.and.j.eq.1) then
                      !  print *, 'iidx+ii = ', iidx+ii, 'jidx+jj = ', jidx+jj
                      !  print *, 'tsdata_ref(1,jidx+jj,iidx+ii,1,1) = ', tsdata_ref(1,jidx+jj,iidx+ii,1,1)
                      !endif
                      tsdata1d(1:ntpoint_total) = tsdata1d(1:ntpoint_total) + tsdata(1:ntpoint_total,jidx+jj,iidx+ii,1,1)
                      numpt1 = numpt1 + 1
                    endif
                  enddo
                enddo
                tsdata_tmp(1:ntpoint_total,j,i,1,1) = tsdata1d(1:ntpoint_total)/dble(numpt1)
              enddo
            enddo
            print *, 'numpt1 = ', numpt1
          elseif(ilocal_ave.eq.2) then
            print *, 'Doing square local avarage ... '


          endif
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          do j=jpointAve_be, nypoint_k_local
            print *, 'j = ', j
            numpt = 0
            spect11 = (0.d0,0.d0)
            spect22 = (0.d0,0.d0)
            spect12 = (0.d0,0.d0)
            do i = 1, nxpoint_k_local
              do ii = max(-ixwindowl_local,imints_k_local-i), min(ixwindowr_local,imaxts_k_local-i)
                spect12(:,ii,1) = spect12(:,ii,1)+crossspect1d(tsdata_ref_tmp(1:ntpoint_total,j,i,1,1),tsdata_tmp(1:ntpoint_total,j,i+ii,1,1))
                spect11(:,ii,1) = spect11(:,ii,1)+autospect1d(tsdata_ref_tmp(1:ntpoint_total,j,i,1,1))
                spect22(:,ii,1) = spect22(:,ii,1)+autospect1d(tsdata_tmp(1:ntpoint_total,j,i+ii,1,1))
                numpt(ii,1) = numpt(ii,1) + 1
              enddo
            enddo

            write(unit=fnum4_2,fmt='(I04.4)') j
            fname = 'coherence_dx_kref'//fnum4_1//'_j'//fnum4_2//'.dat'
            open(unit=12,file=trim(fname),status='unknown')
            rewind(12)
            write(12,'(a)') 'variables = f,x,y,spect12_r,spect12_i,spect11,spect22'

            write(12,*) 'zone T = k'//fnum4//', i =', ntrans/2, ' j =', (ixwindowl_local + ixwindowr_local + 1), ' k =', 1
            do ii = -ixwindowl_local, ixwindowr_local
              if(numpt(ii,1).gt.0) then
                !coh2 = abs(spect12(:,ii,1))**2/(abs(spect11(:,ii,1))*abs(spect22(:,ii,1))+1.e-30)
                !phase = atan2(imag(spect12(:,ii,1)),real(spect12(:,ii,1)))
                do i = 1,ntrans/2
                  write(12,*)  dble(i-1)/fperiod, dble(ii*ixpt)*dx_kp,0, real(spect12(i,ii,1)),imag(spect12(i,ii,1)), abs(spect11(i,ii,1)),abs(spect22(i,ii,1))
                enddo
              endif
            enddo
            close(12)
          enddo ! end j loop

        enddo ! end k loop

      elseif(i2D_coh.eq.1) then

        jwindow = nypoint_k_local/2

        do k = 1, num_kplane
          write(unit=fnum4,fmt='(I04.4)') kplane(k)
          print *, 'wall-normal location k =', kplane(k)
          ! Read k-plane data
          do kk=1, DNSIndex_k%nkplane
          if(kplane(k).eq.DNSIndex_k%kplane(kk)) then
            kplane_be = kk
          endif
          enddo ! end kk loop
          call ReadVOLTSkplane(kplane_be,kplane(k),tsdata)

!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          if(ilocal_ave.eq.1) then
            print *, 'Doing circle local avarage ... '
            do i=1, nxpoint_k_local
              do j=1, nypoint_k_local
                numpt1 = 0; tsdata1d = 0;
                iidx = iradius_pt + (i-1)*ixpt + 1
                jidx = jradius_pt + (j-1)*jypt + 1
                do ii=-iradius_pt, iradius_pt
                  do jj=-jradius_pt, jradius_pt
                    iradius_tmp = sqrt( (dble(ii)*dx_kp)**2 + (dble(jj)*dy_kp)**2 )
                    if(iradius_tmp.le.iradius) then
                      !if(i.eq.1.and.j.eq.1) then
                      !  print *, 'iidx+ii = ', iidx+ii, 'jidx+jj = ', jidx+jj
                      !  print *, 'tsdata_ref(1,jidx+jj,iidx+ii,1,1) = ', tsdata_ref(1,jidx+jj,iidx+ii,1,1)
                      !endif
                      tsdata1d(1:ntpoint_total) = tsdata1d(1:ntpoint_total) + tsdata(1:ntpoint_total,jidx+jj,iidx+ii,1,1)
                      numpt1 = numpt1 + 1
                    endif
                  enddo
                enddo
                tsdata_tmp(1:ntpoint_total,j,i,1,1) = tsdata1d(1:ntpoint_total)/dble(numpt1)
              enddo
            enddo
            print *, 'numpt1 = ', numpt1
          elseif(ilocal_ave.eq.2) then
            print *, 'Doing square local avarage ... '


          endif
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          do j=jpointAve_be, nypoint_k_local/2-1
            print *, 'j = ', j
            numpt = 0
            spect11 = (0.d0,0.d0)
            spect22 = (0.d0,0.d0)
            spect12 = (0.d0,0.d0)
            do i = 1, nxpoint_k_local
              do jj = 0,jwindow
                jjj = j+jj
                do ii = max(-ixwindowl_local,imints_k_local-i), min(ixwindowr_local,imaxts_k_local-i)
                  spect12(:,ii,jj) = spect12(:,ii,jj)+crossspect1d(tsdata_ref_tmp(1:ntpoint_total,j,i,1,1),tsdata_tmp(1:ntpoint_total,jjj,i+ii,1,1))
                  spect11(:,ii,jj) = spect11(:,ii,jj)+autospect1d(tsdata_ref_tmp(1:ntpoint_total,j,i,1,1))
                  spect22(:,ii,jj) = spect22(:,ii,jj)+autospect1d(tsdata_tmp(1:ntpoint_total,jjj,i+ii,1,1))
                  numpt(ii,jj) = numpt(ii,jj) + 1
                enddo
              enddo
            enddo

            write(unit=fnum4_2,fmt='(I04.4)') j
            open(unit=12,file='coherence_kref'//fnum4_1//'_j'//fnum4_2//'.dat')
            rewind(12)
            write(12,'(a)') 'variables = f,x,y,spect12_r,spect12_i,spect11,spect22'

            numptx = count(numpt(:,0).gt.0)
            write(12,*) 'zone T = k'//fnum4//', i =', ntrans/2, ' j =', numptx, ' k =', jwindow+1
            do jj = 0, jwindow
              do ii = -ixwindowl_local,ixwindowr_local
                if(numpt(ii,jj).gt.0) then
                  !coh2 = abs(spect12(:,ii,jj))**2/(abs(spect11(:,ii,jj))*abs(spect22(:,ii,jj))+1.e-30)
                  !phase = atan2(imag(spect12(:,ii,jj)),real(spect12(:,ii,jj)))
                  do i = 1,ntrans/2
                    write(12,*)  dble(i-1)/fperiod, dble(ii*ixpt)*dx_kp, dble(jj*jypt)*dy_kp, real(spect12(i,ii,jj)),imag(spect12(i,ii,jj)), abs(spect11(i,ii,jj)),abs(spect22(i,ii,jj))
                  enddo
                 endif
              enddo
            enddo
            close(12)
          enddo ! end j loop
        enddo  ! end k loop

      endif

    endif ! end ilocal_ave.eq.0

    endif ! end iFileType.eq.0

  endif ! end iAve_data.eq.0

  contains

    subroutine CalAve()
      integer :: n, i, j, k, num_file
      real(8), dimension(:,:,:,:), allocatable :: vars, vars_tmp
      character(400) :: fname, variables1, variables2

      allocate(vars(fdim1,fdim2,fdim3,7),vars_tmp(fdim1,fdim2,fdim3,7))
      vars = 0.
      num_file = (fileprop(1)%file_end-fileprop(1)%file_be)/fileprop(1)%file_skip + 1
      do n=1, num_file
        write(unit=fnum4,fmt='(I04.4)') fileprop(1)%file_be + (n-1)*fileprop(1)%file_skip
        fname = trim(fileprop(1)%filepath)//fnum4//'.dat'

        print *, 'open file: ', trim(fname)
        open(7,file=trim(fname),status='old')
          read(7,'(2a)') variables1
          read(7,'(2a)') variables2
          do k=1, fdim3
            do j=1, fdim2
              do i=1, fdim1
                read(7,*) vars_tmp(i,j,k,1:7)
              enddo
            enddo
          enddo
        close(7)
        vars(:,:,:,4:7) = vars(:,:,:,4:7) + vars_tmp(:,:,:,4:7)
      enddo
      vars(:,:,:,1:3) = vars_tmp(:,:,:,1:3)

      write(unit=fnum4_1,fmt='(I04.4)') fileprop(1)%file_be
      write(unit=fnum4_2,fmt='(I04.4)') fileprop(1)%file_end
      fname = trim(fileprop(1)%filepath)//'ave'//fnum4_1//'_'//fnum4_2//'.dat'
      print *, 'writing file: ', trim(fname)
      open(7,file=trim(fname),status='unknown')
        write(7,*) trim(variables1)
        write(7,*) trim(variables2)
        do k=1, fdim3
          do j=1, fdim2
            do i=1, fdim1
              write(7,*) vars(i,j,k,1:7)
            enddo
          enddo
        enddo
      close(7)

    end subroutine CalAve

    subroutine Input()
      integer, parameter :: nid=11
      integer :: i, j, k, kk, n, npt_tmp
      logical :: PlaneNameMatch = .false.
      real(8) :: radius

      npath = 1
      read(*,*)
      read(*,*) iFileType, iAve_data, icylin
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
      read(*,*)
      read(*,*)
      read(*,*) fdim1, fdim2, fdim3
      read(*,*)
      read(*,*)
      read(*,*)
      read(*,*)
      read(*,*) i2D_coh, ix_coh
      read(*,*)
      read(*,*) nsection, ioverlap, iwindow
      read(*,*)
      if(iAve_data.eq.0) then
        read(*,*)
        read(*,*) kloc_ref
        read(*,*)
        read(*,*) ibe_k, iend_k, jbe_k, jend_k, ixwindowl, ixwindowr, jwindow, jpointAve_be, jpointAve_end
        read(*,*)
        read(*,*) num_kplane
        allocate(kplane(num_kplane))
        read(*,*) (kplane(i), i = 1, num_kplane)
        read(*,*)
        read(*,*) ilocal_ave, iradius
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
        endif
        if(num_kplane.ne.1) then
            print*,'num_kplane MUST be 1. STOP!'
            stop
        endif
        if(icylin.ge.1) then
            do kk=1, DNSIndex_k%nkplane
              if(kplane(1).eq.DNSIndex_k%kplane(kk)) then
                kplane_be = kk
              endif
            enddo ! end kk loop
            radius = DNSIndex_k%zloc(kplane_be)
            print*,'cylindrical coordinate, radius at kref',kplane(1),': ',radius
            dy_kp = dy_kp*radius
            print*,'dy_kp: ',dy_kp
        endif
      endif ! end iAve_data.eq.0

      if(iAve_data.eq.0) then
        if(i2D_coh.eq.1.and.ix_coh.eq.1) then
          print *, 'only one of the following can be 1 ... STOP! '
          print *, 'i2D_coh = ', i2D_coh, 'ix_coh = ', ix_coh
          stop
        endif
        ! check kplane location index
        if(iFileType.ne.0) then
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


        if(ioverlap.ne.0.and.ioverlap.ne.1) then
          print*,'ioverlap can ONLY be 0 or 1'
          stop
        end if
        print *, 'number of wall-normal locations =', num_kplane
        print *, 'Wall-normal plane indexes k =', kplane
        write(*,*)
        print *, 'ntpoint_total =', ntpoint_total, 'dt_sample = ', dt_sample

        if((ilocal_ave.eq.1.or.ilocal_ave.eq.2)) then
          print *, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
          if(ilocal_ave.eq.1) print *, 'using circle local average'
          if(ilocal_ave.eq.2) print *, 'using square local average'
          print *, 'dx_kp = ', dx_kp, 'dy_kp = ', dy_kp, 'iradius = ', iradius
          iradius_pt = int(iradius/dx_kp)
          jradius_pt = int(iradius/dy_kp)
          ixpt = iradius_pt*2 + 1
          jypt = jradius_pt*2 + 1
          print *, 'iradius_pt = ', iradius_pt, 'jradius_pt = ', jradius_pt
          print *, 'adjust ibe, iend, jbe, jend, ixwindowl, ixwindowr ... '
          npt_tmp = iend_k - ibe_k + 1
          n = npt_tmp/ixpt
          iend_k = ibe_k + n*ixpt - 1
          npt_tmp = jend_k - jbe_k + 1
          n = npt_tmp/jypt
          jend_k = jbe_k + n*jypt - 1
          print *, 'ibe_ave = ', ibe_k, 'iend_ave = ', iend_k
          print *, 'jbe_ave = ', jbe_k, 'jend_ave = ', jend_k
          print *, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

          n = ixwindowl/ixpt
          ixwindowl = n*ixpt
          n = ixwindowr/ixpt
          ixwindowr = n*ixpt
        endif

        print *, 'ixwindowl =', ixwindowl,  'ixwindowr =', ixwindowr

     !stop

    endif

    end subroutine input


end program Calcoh
