module modSpect
   use MFFT1D
   use MFFT2D
   use MFFT3D
   use MFFTWindow
   implicit none
   integer, private :: nsection, nsection_kx
   integer, private :: iwindow, ioverlap, nperiod, nperiod_kx
   integer, private :: ntpoint, ntrans, nypoint, nxpoint, ntrans_kx
   real(8), private :: dts, dys, dxs
   real(8), private :: pi
   real(8), dimension(:), allocatable, private :: wcoeff_t, wcoeff_kx
   logical :: IsSpectInit = .FALSE.
!                  nsection (# of bins)
!                  iwindow ( windowing type)
!                  ioverlap (whether or not 1/2 overlap between neighboring bins)
!                  ntpoint (total number of points), ntrans (# of points per bin)
!                  ntrans (size of bins used for Fourier transformation)
!                  nsection (number of bins; nsection = 8 for standard Welch method)
contains

     subroutine InitSpect1d(nt, nsect, iwin, ihflap, dt_sample, npbin, npr)
        implicit none
        integer, intent(in) :: nt, nsect, iwin, ihflap
        real(8), intent(in) :: dt_sample  ! sampling time interval
        integer, intent(out) :: npbin, npr

        ntpoint = nt
        nsection = nsect
        iwindow = iwin
        ioverlap = ihflap
        dts = dt_sample
        nperiod = ntpoint/(nsection + ioverlap)
        npr =  nperiod


        pi = 4.*atan(1.)
        ntrans = nperiod*(ioverlap+1)
        npbin = ntrans
        if (ntrans.gt.ntpoint) then
             print*, 'specified nperiod with ioverlap flag greater than total points available'
             print*, 'ntpoint = ', ntpoint, ', nperiod = ', nperiod, 'ioverlap=',ioverlap
             stop
        else if (nperiod.le.0) then
             nperiod = ntpoint
             ntrans = nperiod
             ioverlap = 0
        end if

        print *, '# of data points per segments =', ntrans
        print *, '# of segments =', nsection
        print *, 'Sampling frequency (Herz) =', 1./dts
        print *, 'Interval of sampling (sec) =', dts
        print *, 'Length per segments (sec) =', (ntrans-1)*dts
        if(ioverlap.eq.0) print *, 'Segments with no overlap'
        if(ioverlap.eq.1) print *, 'Segments with half overlap'
        call InitFFT1D(ntrans)
        call InitFFTWindow(ntrans,iwindow)
        IsSpectInit = .TRUE.
    end subroutine InitSpect1d

!     subroutine InitSpect1d(nt, npr, iwin, ihflap, dt_sample, npbin, nsect)
!        implicit none
!        integer, intent(in) :: nt, npr, iwin, ihflap
!        real(8), intent(in) :: dt_sample  ! sampling time interval
!        integer, intent(out) :: npbin, nsect
!
!        ntpoint = nt
!        nperiod = npr
!        iwindow = iwin
!        ioverlap = ihflap
!        dts = dt_sample
!
!        pi = 4.*atan(1.)
!        ntrans = nperiod*(ioverlap+1)
!        npbin = ntrans
!        if (ntrans.gt.ntpoint) then
!             print*, 'specified nperiod with ioverlap flag greater than total points available'
!             print*, 'ntpoint = ', ntpoint, ', nperiod = ', nperiod, 'ioverlap=',ioverlap
!             stop
!        else if (nperiod.le.0) then
!             nperiod = ntpoint
!             ntrans = nperiod
!             ioverlap = 0
!        end if
!        nsection = ntpoint/nperiod - ioverlap
!        nsect = nsection
!        print *, '# of data points per segments =', ntrans
!        print *, '# of segments =', nsection
!        print *, 'Sampling frequency (Herz) =', 1./dts
!        print *, 'Interval of sampling (sec) =', dts
!        print *, 'Length per segments (sec) =', (ntrans-1)*dts
!        if(ioverlap.eq.0) print *, 'Segments with no overlap'
!        if(ioverlap.eq.1) print *, 'Segments with half overlap'
!        call InitFFT1D(ntrans)
!        call InitFFTWindow(ntrans,iwindow)
!        IsSpectInit = .TRUE.
!    end subroutine InitSpect1d

     subroutine InitSpect1dy(ny, nsect, iwin, ihflap, dy_sample, npbin, npr)
        implicit none
        integer, intent(in) :: ny, nsect, iwin, ihflap
        real(8), intent(in) :: dy_sample
        integer, intent(out) :: npbin, npr

        nypoint = ny
        nsection = nsect
        iwindow = iwin
        ioverlap = ihflap
        dys = dy_sample
        nperiod = nypoint/(nsection + ioverlap)
        npr =  nperiod

        pi = 4.*atan(1.)
        ntrans = nperiod*(ioverlap+1)
        npbin = ntrans
        if (ntrans.gt.nypoint) then
             print*, 'specified nperiod with ioverlap flag greater than total points available'
             print*, 'nypoint = ', nypoint, ', nperiod = ', nperiod, 'ioverlap=',ioverlap
             stop
        else if (nperiod.le.0) then
             nperiod = nypoint
             ntrans = nperiod
             ioverlap = 0
        end if

        print *, '# of data points per segments =', ntrans
        print *, '# of segments =', nsection
!        print *, 'Sampling frequency (Herz) =', 1./dts
        print *, 'Interval of sampling (m) =', dys
        print *, 'Length per segments (m) =', (ntrans-1)*dys
        if(ioverlap.eq.0) print *, 'Segments with no overlap'
        if(ioverlap.eq.1) print *, 'Segments with half overlap'
        call InitFFT1D(ntrans)
        call InitFFTWindow(ntrans,iwindow)
        IsSpectInit = .TRUE.
    end subroutine InitSpect1dy


     subroutine InitSpect2d_ty(nt, ny, nsect, iwin, ihflap, dt_sample, dy_sample, npbin, npr)
        implicit none
        integer, intent(in) :: nt, ny, nsect, iwin, ihflap
        real(8), intent(in) :: dt_sample, dy_sample
        integer, intent(out) :: npbin, npr

        ntpoint = nt
        nypoint = ny
        nsection = nsect
        iwindow = iwin
        ioverlap = ihflap
        dts = dt_sample
        dys = dy_sample
        nperiod = ntpoint/(nsection + ioverlap)
!        nperiod = nypoint/(nsection + ioverlap)
        npr =  nperiod

        pi = 4.*atan(1.)
        ntrans = nperiod*(ioverlap+1)
        npbin = ntrans
        if (ntrans.gt.ntpoint) then
             print*, 'specified nperiod with ioverlap flag greater than total points available'
             print*, 'ntpoint = ', ntpoint, ', nperiod = ', nperiod, 'ioverlap=',ioverlap
             stop
        else if (nperiod.le.0) then
             nperiod = ntpoint
             ntrans = nperiod
             ioverlap = 0
        end if
        print *, '###############################################'
        print *, '# of data points per segments in time domain =', ntrans
        print *, '# of segments in time domain =', nsection
        print *, 'Sampling frequency (Herz) =', 1./dts
        print *, 'Interval of sampling (sec) = ', dts
        print *, 'Length per segments (sec) =', (ntrans-1)*dts
        if(ioverlap.eq.0) print *, 'Segments with no overlap'
        if(ioverlap.eq.1) print *, 'Segments with half overlap'

        print *, '###############################################'
        print *, '# of data points per segments in space domain', nypoint
        print *, '# of segments in space domain = 1 '
        print *, 'Interval of sampling (m) =', dys
        print *, 'Length per segments (m) =', (nypoint)*dys
        print *, 'Segments with no overlap'

        call InitFFT2D(ntrans,nypoint)
!        call InitFFT1D(ntrans)
        call InitFFTWindow(ntrans,iwindow)

        IsSpectInit = .TRUE.
    end subroutine InitSpect2d_ty


     subroutine InitSpect3d(nt, ny, nx, nsect, nsect_kx, iwin, ihflap, dt_sample, dy_sample, dx_sample, npbin, npbin_kx)
        implicit none
        integer, intent(in) :: nt, ny, nx, nsect, nsect_kx, iwin, ihflap
        real(8), intent(in) :: dt_sample, dy_sample, dx_sample
        integer, intent(out) :: npbin, npbin_kx

        ntpoint = nt
        nypoint = ny
        nxpoint = nx
        nsection = nsect
        nsection_kx = nsect_kx
        iwindow = iwin
        ioverlap = ihflap
        dts = dt_sample
        dys = dy_sample
        dxs = dx_sample
        nperiod = ntpoint/(nsection + ioverlap)
        nperiod_kx = nxpoint/(nsection_kx + ioverlap)
!        nperiod = nypoint/(nsection + ioverlap)
!        npr =  nperiod

        pi = 4.*atan(1.)
        ntrans = nperiod*(ioverlap+1)
        ntrans_kx = nperiod_kx*(ioverlap+1)
        npbin = ntrans
        npbin_kx = ntrans_kx
        if (ntrans.gt.ntpoint) then
             print*, 'specified nperiod with ioverlap flag greater than total points available'
             print*, 'ntpoint = ', ntpoint, ', nperiod = ', nperiod, 'ioverlap=',ioverlap
             stop
        else if (nperiod.le.0) then
             nperiod = ntpoint
             ntrans = nperiod
             ioverlap = 0
        end if
        if (ntrans_kx.gt.nxpoint) then
             print*, 'specified nperiod with ioverlap flag greater than total points available'
             print*, 'nypoint = ', nypoint, ', nperiod_kx = ', nperiod_kx, 'ioverlap=',ioverlap
             stop
        else if (nperiod_kx.le.0) then
             nperiod_kx = nypoint
             ntrans_kx = nperiod_kx
             ioverlap = 0
        end if

        print *, '###############################################'
        print *, '# of data points per segments in time domain = ', ntrans
        print *, '# of segments in time domain =', nsection
        print *, 'Sampling frequency (Herz) =', 1./dts
        print *, 'Interval of sampling (sec) = ', dts
        print *, 'Length per segments (sec) =', (ntrans-1)*dts
        if(ioverlap.eq.0) print *, 'Segments with no overlap'
        if(ioverlap.eq.1) print *, 'Segments with half overlap'

        print *, '###############################################'
        print *, '# of data points per segments in space domain kx = ', nxpoint
        print *, '# of segments in space domain = ', nsection_kx
        print *, 'Interval of sampling (m) =', dxs
        print *, 'Length per segments (m) =', (ntrans_kx-1)*dxs
        if(ioverlap.eq.0) print *, 'Segments with no overlap'
        if(ioverlap.eq.1) print *, 'Segments with half overlap'

        print *, '###############################################'
        print *, '# of data points per segments in space domain ky = ', nypoint
        print *, '# of segments in space domain = 1 '
        print *, 'Interval of sampling (m) =', dys
        print *, 'Length per segments (m) =', (nypoint-1)*dys
        print *, 'Segments with no overlap'

        call InitFFT3D(ntrans,nypoint,nxpoint)

        allocate(wcoeff_t(ntrans), wcoeff_kx(ntrans_kx))
        call InitFFTWindow_1(ntrans,iwindow,wcoeff_t)
        call InitFFTWindow_1(ntrans_kx,iwindow,wcoeff_kx)

        IsSpectInit = .TRUE.
    end subroutine InitSpect3d

     subroutine InitCoherence(nt, ny, nx, nsect, nsect_kx, iwin, ihflap, dt_sample, dy_sample, dx_sample, npbin, npbin_kx)
        implicit none
        integer, intent(in) :: nt, ny, nx, nsect, nsect_kx, iwin, ihflap
        real(8), intent(in) :: dt_sample, dy_sample, dx_sample
        integer, intent(out) :: npbin, npbin_kx

        ntpoint = nt
        nypoint = ny
        nxpoint = nx
        nsection = nsect
        nsection_kx = nsect_kx
        iwindow = iwin
        ioverlap = ihflap
        dts = dt_sample
        dys = dy_sample
        dxs = dx_sample
        nperiod = ntpoint/(nsection + ioverlap)
        nperiod_kx = nxpoint/(nsection_kx + ioverlap)
!        nperiod = nypoint/(nsection + ioverlap)
!        npr =  nperiod

        pi = 4.*atan(1.)
        ntrans = nperiod*(ioverlap+1)
        ntrans_kx = nperiod_kx*(ioverlap+1)
        npbin = ntrans
        npbin_kx = ntrans_kx
        if (ntrans.gt.ntpoint) then
             print*, 'specified nperiod with ioverlap flag greater than total points available'
             print*, 'ntpoint = ', ntpoint, ', nperiod = ', nperiod, 'ioverlap=',ioverlap
             stop
        else if (nperiod.le.0) then
             nperiod = ntpoint
             ntrans = nperiod
             ioverlap = 0
        end if
        if (ntrans_kx.gt.nxpoint) then
             print*, 'specified nperiod with ioverlap flag greater than total points available'
             print*, 'nypoint = ', nypoint, ', nperiod_kx = ', nperiod_kx, 'ioverlap=',ioverlap
             stop
        else if (nperiod_kx.le.0) then
             nperiod_kx = nxpoint
             ntrans_kx = nperiod_kx
             ioverlap = 0
        end if

        print *, '###############################################'
        print *, '# of data points per segments in time domain = ', ntrans
        print *, '# of segments in time domain =', nsection
        print *, 'Sampling frequency (Herz) =', 1./dts
        print *, 'Interval of sampling (sec) = ', dts
        print *, 'Length per segments (sec) =', (ntrans-1)*dts
        if(ioverlap.eq.0) print *, 'Segments with no overlap'
        if(ioverlap.eq.1) print *, 'Segments with half overlap'

        print *, '###############################################'
        print *, '# of data points per segments in space domain kx = ', ntrans_kx
        print *, '# of segments in space domain = ', nsection_kx
        print *, 'Interval of sampling (m) =', dxs
        print *, 'Length per segments (m) =', (ntrans_kx-1)*dxs
        if(ioverlap.eq.0) print *, 'Segments with no overlap'
        if(ioverlap.eq.1) print *, 'Segments with half overlap'

        print *, '###############################################'
        print *, '# of data points per segments in space domain ky = ', nypoint
        print *, '# of segments in space domain = 1 '
        print *, 'Interval of sampling (m) =', dys
        print *, 'Length per segments (m) =', (nypoint-1)*dys
        print *, 'Segments with no overlap'

        call InitFFT3D(ntrans,nypoint,ntrans_kx)

        allocate(wcoeff_t(ntrans), wcoeff_kx(ntrans_kx))
        call InitFFTWindow_1(ntrans,iwindow,wcoeff_t)
        call InitFFTWindow_1(ntrans_kx,iwindow,wcoeff_kx)

        IsSpectInit = .TRUE.
    end subroutine InitCoherence

    function spect1d(tfun)
      implicit none
      real(8), intent(in) :: tfun(ntpoint)
      complex :: spect1d(ntrans,nsection)
      real(8) :: tmean
      integer :: nn
      complex, dimension(ntrans) :: ctmp

      ! Calculate the mean and energy of time signal
      tmean = sum(tfun(1:ntpoint))/dble(ntpoint)

      ! Compute spectrum for time data length of ntans
      do nn = 1, nsection
         ctmp(1:ntrans) = wcoeff(1:ntrans)*(tfun((nn-1)*nperiod+1:(nn-1)*nperiod+ntrans)-tmean)
         call FFT1DF(ctmp, ctmp)
         spect1d(:,nn) = ctmp
      enddo
    end function spect1d

  function spect1dy(tfun)
      implicit none
      real(8), intent(in) :: tfun(nypoint)
      complex :: spect1dy(ntrans,nsection)
      real(8) :: tmean
      integer :: nn
      complex, dimension(ntrans) :: ctmp

      ! Calculate the mean and energy of time signal
      tmean = sum(tfun(1:nypoint))/dble(nypoint)

      ! Compute spectrum for time data length of ntans
      do nn = 1, nsection
         ctmp(1:ntrans) = wcoeff(1:ntrans)*(tfun((nn-1)*nperiod+1:(nn-1)*nperiod+ntrans)-tmean)
         call FFT1DF(ctmp, ctmp)
         spect1dy(:,nn) = ctmp
      enddo
    end function spect1dy


    function spect2d(tfun)
      implicit none
      real(8), intent(in) :: tfun(ntpoint,nypoint)
      complex :: spect2d(ntrans,nypoint,nsection)
      real(8) :: tmean
      integer :: nn, j
      complex, dimension(ntrans,nypoint) :: ctmp

      tmean = sum(tfun(1:ntpoint,1:nypoint))/dble(ntpoint*nypoint)

      do nn=1, nsection
        do j=1, nypoint
          ctmp(1:ntrans,j) = wcoeff(1:ntrans)*(tfun((nn-1)*nperiod+1:(nn-1)*nperiod+ntrans,j)-tmean)
        enddo
        call FFT2DF(ctmp,ctmp)
        spect2d(1:ntrans,1:nypoint,nn) = ctmp(1:ntrans,1:nypoint)
      enddo ! end nn loop
    end function spect2d

    function spect3d(tfun)
      implicit none
      real(8), intent(in) :: tfun(ntpoint,nypoint,nxpoint)
      real(8) :: tfun_tmp(ntpoint,nypoint,nxpoint)
      complex :: spect3d(ntrans,nypoint,ntrans_kx,nsection,nsection_kx)
      real(8), dimension(nxpoint) :: tmean
      integer :: nn, mm, i, j, jj
      complex, dimension(ntrans,nypoint,ntrans_kx) :: ctmp

      do i=1, nxpoint
        tmean(i) = sum(tfun(1:ntpoint,1:nypoint,i))/dble(ntpoint*nypoint)
        tfun_tmp(1:ntpoint,1:nypoint,i) = tfun(1:ntpoint,1:nypoint,i) - tmean(i)
      enddo

      do mm=1, nsection_kx
        do nn=1, nsection
          do j=1, nypoint
            do jj=1, ntrans_kx
              ctmp(1:ntrans,j,jj) = wcoeff_kx(jj)*wcoeff_t(1:ntrans)*(tfun_tmp( (nn-1)*nperiod+1:(nn-1)*nperiod+ntrans,j, &
                                                                            (mm-1)*nperiod_kx+jj ))
            enddo ! end jj loop
          enddo ! end j loop
          call FFT3DF(ctmp,ctmp)
          spect3d(1:ntrans,1:nypoint,1:ntrans_kx,nn,mm) = ctmp(1:ntrans,1:nypoint,1:ntrans_kx)
        enddo ! end nn loop
      enddo ! end mm loop



    end function spect3d


    function spectra3d_ave(tfun)
      implicit none
      real(8), intent(in) :: tfun(ntpoint,nypoint,nxpoint)
      real(8) :: tfun_tmp(ntpoint,nypoint,nxpoint)
      complex :: spectra3d_ave(ntrans,nypoint,ntrans_kx)
      real(8), dimension(nxpoint) :: tmean
      integer :: nn, mm, i, j, jj
      complex, dimension(ntrans,nypoint,ntrans_kx) :: ctmp

      do i=1, nxpoint
        tmean(i) = sum(tfun(1:ntpoint,1:nypoint,i))/dble(ntpoint*nypoint)
        tfun_tmp(1:ntpoint,1:nypoint,i) = tfun(1:ntpoint,1:nypoint,i) - tmean(i)
      enddo

      spectra3d_ave = 0.d0
      do mm=1, nsection_kx
        do nn=1, nsection
          do j=1, nypoint
            do jj=1, ntrans_kx
              ctmp(1:ntrans,j,jj) = wcoeff_kx(jj)*wcoeff_t(1:ntrans)*(tfun_tmp( (nn-1)*nperiod+1:(nn-1)*nperiod+ntrans,j, &
                                                                            (mm-1)*nperiod_kx+jj ))
            enddo ! end jj loop
          enddo ! end j loop
          call FFT3DF(ctmp,ctmp)
          spectra3d_ave(1:ntrans,1:nypoint,1:ntrans_kx) = spectra3d_ave(1:ntrans,1:nypoint,1:ntrans_kx)+&
                                                          ctmp(1:ntrans,1:nypoint,1:ntrans_kx)*conjg(ctmp(1:ntrans,1:nypoint,1:ntrans_kx))
        enddo ! end nn loop
      enddo ! end mm loop
      spectra3d_ave = spectra3d_ave/dble(nsection*nsection_kx)
    end function spectra3d_ave

    function autospect3d(tfun)
      implicit none
      real(8), intent(in) :: tfun(ntpoint,nypoint,nxpoint)
      complex :: autospect3d(ntrans,nypoint,ntrans_kx)
      complex :: ctmp(ntrans,nypoint,ntrans_kx,nsection,nsection_kx)
      complex :: ctmp_tmp(ntrans,nypoint,ntrans_kx,nsection,nsection_kx)
      integer :: i,j,k

      ctmp = spect3d(tfun)
!      autospect3d = sum(ctmp*conjg(ctmp),dim=3)/dble(nsection)
      ctmp_tmp = ctmp*conjg(ctmp)
      do i=1, ntrans
        do j=1, nypoint
          do k=1, ntrans_kx
            autospect3d(i,j,k) = sum(ctmp_tmp(i,j,k,1:nsection,1:nsection_kx))
          enddo
        enddo
      enddo
!      autospect3d = sum(ctmp_tmp( )
      autospect3d = autospect3d/dble(nsection*nsection_kx)

    end function autospect3d


    function autospect2d(tfun)
      implicit none
      real(8), intent(in) :: tfun(ntpoint,nypoint)
      complex :: autospect2d(ntrans,nypoint)
      complex :: ctmp(ntrans,nypoint,nsection)

      ctmp = spect2d(tfun)
      autospect2d = sum(ctmp*conjg(ctmp),dim=3)/dble(nsection)
    end function autospect2d

    function sumspect1d(tfun)
      implicit none
      real(8), intent(in) :: tfun(ntpoint)
      complex :: sumspect1d(ntrans)
      complex :: ctmp(ntrans,nsection)

      ctmp = spect1d(tfun)
      sumspect1d = sum(ctmp,dim=2)/dble(nsection)
    end function sumspect1d

    function autospect1d(tfun)
      implicit none
      real(8), intent(in) :: tfun(ntpoint)
      complex :: autospect1d(ntrans)
      complex :: ctmp(ntrans,nsection)

      ctmp = spect1d(tfun)
      autospect1d = sum(ctmp*conjg(ctmp),dim=2)/dble(nsection)
    end function autospect1d

    function autospect1dty(tfun)
      implicit none
      real(8), intent(in) :: tfun(nypoint)
      complex :: autospect1dty(ntrans)
      complex :: ctmp(ntrans,nsection)

      ctmp = spect1dy(tfun)
      autospect1dty = sum(ctmp*conjg(ctmp),dim=2)/dble(nsection)
    end function autospect1dty

    function crossspect1d(tfun1,tfun2)
      implicit none
      real(8), intent(in) :: tfun1(ntpoint), tfun2(ntpoint)
      complex :: crossspect1d(ntrans)
      complex, dimension(ntrans,nsection) :: ctmp1, ctmp2

      ctmp1 = spect1d(tfun1)
      ctmp2 = spect1d(tfun2)
      crossspect1d = sum(ctmp1*conjg(ctmp2),dim=2)/dble(nsection)

    end function crossspect1d

    function cross_spectra_fft2db(spectra2d)
        implicit none
        complex, intent(in) :: spectra2d(nypoint,ntrans_kx)
        complex :: cross_spectra_fft2db(nypoint,ntrans_kx)

        call FFT2DB(spectra2d, cross_spectra_fft2db)
    end function cross_spectra_fft2db

subroutine ShiftFFT3D_real(np1,np2,np3,spectin,spectout,wavenum1,wavenum2,wavenum3)
    integer, intent(in) :: np1, np2, np3
    real(8), dimension(np1,np2,np3), intent(in) :: spectin
    real(8), dimension(np1,np2,np3), intent(out) :: spectout
    integer, intent(out) :: wavenum1(np1),wavenum2(np2),wavenum3(np3)
    real(8), dimension(:,:,:),allocatable :: tmp
    integer :: i,whalf

    ! wave 3
    allocate(tmp(np1,np2,np3))
    tmp = spectin
    wavenum3 = (/ (i, i = 0, np3-1) /)
    whalf = wavenum3(np3/2+1)
    wavenum3 = wavenum3 - whalf
    spectout(:,:,np3/2+1:np3) = tmp(:,:,1:np3/2)
    spectout(:,:,1:np3/2) = tmp(:,:,np3/2+1:np3)

    ! wave 2
    wavenum2 = (/ (i, i = 0, np2-1) /)
    whalf = wavenum2(np2/2+1)
    wavenum2 = wavenum2 - whalf
    tmp(:,np2/2+1:np2,:) = spectout(:,1:np2/2,:)
    tmp(:,1:np2/2,:) = spectout(:,np2/2+1:np2,:)

    ! wave 1
    wavenum1 = (/ (i, i = 0, np1-1) /)
    whalf = wavenum1(np1/2+1)
    wavenum1 = wavenum1 - whalf
    spectout(np1/2+1:np1,:,:) = tmp(1:np1/2,:,:)
    spectout(1:np1/2,:,:) = tmp(np1/2+1:np1,:,:)

    deallocate(tmp)
end subroutine ShiftFFT3D_real

subroutine ShiftFFT3D(np1,np2,np3,spectin,spectout,wavenum1,wavenum2,wavenum3)
    integer, intent(in) :: np1, np2, np3
    complex, dimension(np1,np2,np3), intent(in) :: spectin
    complex, dimension(np1,np2,np3), intent(out) :: spectout
    integer, intent(out) :: wavenum1(np1),wavenum2(np2),wavenum3(np3)
    complex, dimension(:,:,:),allocatable :: tmp
    integer :: i,whalf

    ! wave 3
    allocate(tmp(np1,np2,np3))
    tmp = spectin
    wavenum3 = (/ (i, i = 0, np3-1) /)
    whalf = wavenum3(np3/2+1)
    wavenum3 = wavenum3 - whalf
    spectout(:,:,np3/2+1:np3) = tmp(:,:,1:np3/2)
    spectout(:,:,1:np3/2) = tmp(:,:,np3/2+1:np3)

    ! wave 2
    wavenum2 = (/ (i, i = 0, np2-1) /)
    whalf = wavenum2(np2/2+1)
    wavenum2 = wavenum2 - whalf
    tmp(:,np2/2+1:np2,:) = spectout(:,1:np2/2,:)
    tmp(:,1:np2/2,:) = spectout(:,np2/2+1:np2,:)

    ! wave 1
    wavenum1 = (/ (i, i = 0, np1-1) /)
    whalf = wavenum1(np1/2+1)
    wavenum1 = wavenum1 - whalf
    spectout(np1/2+1:np1,:,:) = tmp(1:np1/2,:,:)
    spectout(1:np1/2,:,:) = tmp(np1/2+1:np1,:,:)

    deallocate(tmp)
end subroutine ShiftFFT3D

    subroutine CalPowerspect3D(spect3d,fn)
      implicit none
      complex, intent(in) :: spect3d(ntrans,nypoint,ntrans_kx)
      character(*), intent(in) :: fn
      complex :: powerspect(ntrans,nypoint,ntrans_kx)
      real(8) :: ff, yy, xx, fperiod, kyperiod, kxperiod
      integer :: twave(ntrans), ywave(nypoint), xwave(ntrans_kx)
      integer :: i, j, n

      fperiod = dble(ntrans-1)*dts
      kyperiod = dble(nypoint-1)*dys
      kxperiod = dble(nxpoint-1)*dxs
      call RearrangeFFT3D(ntrans,nypoint,ntrans_kx,spect3d,powerspect,twave,ywave,xwave)

      open(21, file='3D_powerspect'//trim(fn)//'.dat')
      rewind(21)
      write(21,'(a)') 'variables=f,ky, kx,spectrum'
      write(21,'("Zone T=ky-t_PSD, I=", I8, ",J=", I8, ", K=", I8, ", AUXDATA rms2=""", E16.9,"""")') ntrans, nypoint, ntrans_kx, 0.

      do j=1, ntrans_kx
        do i=1, nypoint
          do n=1, ntrans
            ff = dble(twave(n))/fperiod
            yy = 2.d0*pi*dble(ywave(i))/kyperiod
            xx = 2.d0*pi*dble(xwave(j))/kxperiod
            write(21,'(4E20.12)') ff, yy, xx, real(powerspect(n,i,j))
          enddo
        enddo
      enddo
      close(21)

    end subroutine CalPowerspect3D


    subroutine CalPowerspect_ky_t(spect2d,tmeanenergy,fn)
      implicit none
      complex, intent(in) :: spect2d(ntrans,nypoint)
      real(8), intent(in) :: tmeanenergy
      character(*), intent(in) :: fn
      complex :: powerspect(ntrans,nypoint)
      real(8) :: ff, yy, fperiod, kyperiod, tmpenergy
      integer :: twave(ntrans), ywave(nypoint)
      integer :: i, n

      fperiod = dble(ntrans-1)*dts
      kyperiod = dble(nypoint)*dys
      call RearrangeFFT2D(ntrans,nypoint,spect2d,powerspect,twave,ywave)
      tmpenergy = sum(powerspect)/fperiod/kyperiod     

!print *, 'sum(powerspect) = ', real(sum(powerspect))

      open(21, file=trim(fn)//'.dat')
      rewind(21)
      write(21,'(a)') 'variables=f,ky,spectrum'
      write(21,'("Zone T=ky-t_PSD, I=", I8, ",J=", I8, ", AUXDATA rms2=""", E16.9,"""")') ntrans, nypoint, tmeanenergy

      do i=1, nypoint
        do n=1, ntrans
          ff = real(twave(n))/fperiod
          yy = 2.d0*pi*real(ywave(i))/kyperiod
          write(21,'(3E20.12)') ff, yy, real(powerspect(n,i))*tmeanenergy/tmpenergy/(2.d0*pi)
        enddo
      enddo
      close(21)

    end subroutine CalPowerspect_ky_t


    subroutine CalPowerspect_kx_t(spect2d,tmeanenergy,fn)
      implicit none
      complex, intent(in) :: spect2d(ntrans,nypoint)
      real(8), intent(in) :: tmeanenergy
      character(*), intent(in) :: fn
      complex :: powerspect(ntrans,nypoint)
      real(8) :: ff, yy, fperiod, kyperiod, tmpenergy
      integer :: twave(ntrans), ywave(nypoint)
      integer :: i, n

      fperiod = dble(ntrans-1)*dts
      kyperiod = dble(nypoint)*dys

      call RearrangeFFT2D(ntrans,nypoint,spect2d,powerspect,twave,ywave)
      tmpenergy = sum(powerspect)/fperiod/kyperiod

      open(21, file=trim(fn)//'.dat')
      rewind(21)
      write(21,'(a)') 'variables=f,kx,spectrum'
      write(21,'("Zone T=kx-t_PSD, I=", I8, ",J=", I8, ", AUXDATA rms2=""", E16.9,"""")') ntrans, nypoint, tmeanenergy

      do i=1, nypoint
        do n=1, ntrans
          ff = real(twave(n))/fperiod
          yy = 2.d0*pi*real(ywave(i))/kyperiod
          write(21,'(3E20.12)') ff, yy, real(powerspect(n,i))*tmeanenergy/tmpenergy/(2.d0*pi)
        enddo
      enddo
      close(21)

    end subroutine CalPowerspect_kx_t

    subroutine CalPowerspect_test(spect1d,tmeanenergy,fn)
      implicit none
      complex, intent(in) :: spect1d(ntrans)
      real(8), intent(in) :: tmeanenergy
      real(8) :: powerspect(ntrans), tmpenergy
      complex :: specttmp(ntrans)
      character(*), intent(in) :: fn
      real(8) :: ff, fperiod
      integer :: i, n
      character(100) :: vars

      ! Ref: Numerical Recipe, Chapter 13.4
      ! tmeanenergy: mean squre amplitude (1/N)*sum_1^N |c_i|^2
      ! PSD is normalized to preserve time-integral squred amplitude int_0^T |c(t)|^2 dt

      fperiod = dble(ntrans-1)*dts  ! period per segments
      powerspect = spect1d
      !powerspect(2:ntrans/2) = 2.*spect1d(2:ntrans/2)
      tmpenergy = sum(powerspect(1:ntrans/2))/fperiod  ! sum (fi_i * df), with df = 1/fperiod


      call FFT1DB(spect1d, specttmp)
      open(21, file='powerspect'//trim(fn)//'.dat')
      rewind(21)
      write(21,'(a)')'variables=f,spectrum,p'
      write(21,'("Zone T=PSD, I=", I8.8, ", AUXDATA rms2=""", E16.9,"""")') ntrans/2, tmeanenergy
      write(21,'(a,I8.8,a,I3.3)') '# Data points per seg = ', ntrans, '   Num of seg = ', nsection
      write(21,'(a,I1.1,a,I1.1,a,E15.8)') '# Window type = ', iwindow, '   ioverlap = ', ioverlap, '   dt_sample = ', dts
      write(21,'(a,E15.8)') '# meanenergy (mean square amplitude) = ', tmeanenergy
      do i=1,ntrans
        ff = real(i-1)/fperiod
        write(21,'(3E20.12)') ff, powerspect(i),real(specttmp(i)) 
!        write(21,'(2E20.12)') ff, powerspect(i)
      end do
      close(21)
    end subroutine CalPowerspect_test

    subroutine CalPowerspect(spect1d,tmeanenergy,fn,moment)
      implicit none
      complex, intent(in) :: spect1d(ntrans)
      real(8), intent(in) :: tmeanenergy
      real(8) :: powerspect(ntrans/2), tmpenergy
      character(*), intent(in) :: fn
      real(8), intent(in), optional :: moment(4)
      real(8) :: ff, fperiod
      integer :: i, n
      character(100) :: vars

      ! Ref: Numerical Recipe, Chapter 13.4
      ! tmeanenergy: mean squre amplitude (1/N)*sum_1^N |c_i|^2
      ! PSD is normalized to preserve time-integral squred amplitude int_0^T |c(t)|^2 dt

      fperiod = dble(ntrans-1)*dts  ! period per segments
      powerspect(1) = spect1d(1)
      powerspect(2:ntrans/2) = 2.*spect1d(2:ntrans/2)
      tmpenergy = sum(powerspect(1:ntrans/2))/fperiod  ! sum (fi_i * df), with df = 1/fperiod
      open(21, file='powerspect'//trim(fn)//'.dat')
      rewind(21)
      write(21,'(a)')'variables=f,spectrum'
      write(21,'("Zone T=PSD, I=", I8.8, ", AUXDATA rms2=""", E16.9,"""")') ntrans/2, tmeanenergy
      write(21,'(a,I8.8,a,I3.3)') '# Data points per seg = ', ntrans, '   Num of seg = ', nsection
      write(21,'(a,I1.1,a,I1.1,a,E15.8)') '# Window type = ', iwindow, '   ioverlap = ', ioverlap, '   dt_sample = ', dts
      write(21,'(a,E15.8)') '# meanenergy (mean square amplitude) = ', tmeanenergy
      if (present(moment)) then
         write(21,'(a,E15.8)') '# mean     = ', moment(1)
         write(21,'(a,E15.8)') '# rms      = ', moment(2)
         write(21,'(a,E15.8)') '# skewness = ', moment(3)
         write(21,'(a,E15.8)') '# flatness = ', moment(4)
      endif
      do i=1,ntrans/2
        ff = real(i-1)/fperiod
        write(21,'(2E20.12)') ff, powerspect(i)*tmeanenergy/tmpenergy
!        write(21,'(2E20.12)') ff, powerspect(i)
      end do
      close(21)
    end subroutine CalPowerspect

   subroutine CalPowerspecty(spect1d,tmeanenergy,fn)
      implicit none
      complex, intent(in) :: spect1d(ntrans)
      real(8), intent(in) :: tmeanenergy
      real(8) :: wnspect(ntrans/2+1)
      character(*), intent(in) :: fn
      real(8) :: ff, fperiod, tmpenergy
      integer :: i, n
      character(100) :: vars

      ! Ref: Numerical Recipe, Chapter 13.4
      ! tmeanenergy: mean squre amplitude (1/N)*sum_1^N |c_i|^2
      ! PSD is normalized to preserve time-integral squred amplitude int_0^T |c(t)|^2 dt

      fperiod = dble(ntrans)*dys  ! period per segments
      wnspect(1) = spect1d(1)
      wnspect(2:ntrans/2+1) = 2.*spect1d(2:ntrans/2+1)
      tmpenergy = sum(wnspect(1:ntrans/2+1))/fperiod  ! sum (fi_i * df), with df = 1/fperiod

      open(21, file='wnspect'//trim(fn)//'.dat')
      rewind(21)
      write(21,'(a)')'variables=k,spectrum'
      write(21,'("Zone T=PSD, I=", I8.8, ", AUXDATA rms2=""", E16.9,"""")') ntrans/2+1, tmeanenergy
      write(21,'(a,I8.8,a,I3.3)') '# Data points per seg = ', ntrans, '   Num of seg = ', nsection
      write(21,'(a,I1.1,a,I1.1,a,E15.8)') '# Window type = ', iwindow, '   ioverlap = ', ioverlap, '   dy = ', dys
      write(21,'(a,E15.8)') '# meanenergy (mean square amplitude) = ', tmeanenergy
      do i=1,ntrans/2+1
        ff = dble(i-1)/fperiod*8.d0*atan(1.d0)
        write(21,'(2E20.12)') ff, wnspect(i)*tmeanenergy/tmpenergy/(8.d0*atan(1.d0))
      end do
      close(21)
    end subroutine CalPowerspecty

   subroutine CalPowerspectx(spect1d,tmeanenergy,fn)
      implicit none
      complex, intent(in) :: spect1d(ntrans)
      real(8), intent(in) :: tmeanenergy
      real(8) :: wnspect(ntrans/2+1)
      character(*), intent(in) :: fn
      real(8) :: ff, fperiod, tmpenergy
      integer :: i, n

      ! Ref: Numerical Recipe, Chapter 13.4
      ! tmeanenergy: mean squre amplitude (1/N)*sum_1^N |c_i|^2
      ! PSD is normalized to preserve time-integral squred amplitude int_0^T |c(t)|^2 dt

      fperiod = dble(ntrans)*dys  ! period per segments
      wnspect(1) = spect1d(1)
      wnspect(2:ntrans/2+1) = 2.*spect1d(2:ntrans/2+1)
      tmpenergy = sum(wnspect(1:ntrans/2+1))/fperiod  ! sum (fi_i * df), with df = 1/fperiod

      open(21, file='wnspect'//trim(fn)//'.dat')
      rewind(21)
      write(21,'(a)')'variables=k,spectrum'
      write(21,'("Zone T=PSD, I=", I8.8, ", AUXDATA rms2=""", E16.9,"""")') ntrans/2+1, tmeanenergy
      write(21,'(a,I8.8,a,I3.3)') '# Data points per seg = ', ntrans, '   Num of seg = ', nsection
      write(21,'(a,I1.1,a,I1.1,a,E15.8)') '# Window type = ', iwindow, '   ioverlap = ', ioverlap, '   dx = ', dys
      write(21,'(a,E15.8)') '# meanenergy (mean square amplitude) = ', tmeanenergy
      do i=1,ntrans/2+1
        ff = dble(i-1)/fperiod*8.d0*atan(1.d0)
        write(21,'(2E20.12)') ff, wnspect(i)*tmeanenergy/tmpenergy/(8.d0*atan(1.d0))
       ! write(21,'(2E20.12)') ff, wnspect(i) 
      end do
      close(21)
    end subroutine CalPowerspectx

    subroutine CalCoherence(autosp1,autosp2,crosssp,fn)
      implicit none
      complex, intent(in) :: autosp1(ntrans), autosp2(ntrans), crosssp(ntrans)
      character(*), intent(in) :: fn
      complex :: ctmp(ntrans)
      real(8) :: Coh(ntrans)
      real(8) :: ff, fperiod
      integer :: i, n
      character(100) :: vars

      fperiod = dble(ntrans-1)*dts  ! period per segments
      ctmp = abs(crosssp)**2/(abs(autosp1)*abs(autosp2))
      Coh(1) = ctmp(1)
      Coh(2:ntrans/2) = 2.*ctmp(2:ntrans/2)
      open(21, file='powerspect'//trim(fn)//'.dat')
      write(21,'(a)')'variables=f,coherence'
      write(21,'(a,I8.8,a,I3.3)') '# Data points per seg = ', ntrans, '   Num of seg = ', nsection
      write(21,'(a,I1.1,a,I1.1,a,E15.8)') '# Window type = ', iwindow, '   ioverlap = ', ioverlap, '   dt_sample = ', dts
      do i=1,ntrans/2
        ff = real(i-1)/fperiod
        write(21,'(2E20.12)')ff, Coh(i)
      end do
      close(21)
    end subroutine CalCoherence

    subroutine CalCorrt(spect1d,fn)
      implicit none
      complex, intent(in) :: spect1d(ntrans)
      character(*), intent(in) :: fn
      complex(8), dimension(ntrans) :: corrt, corrout
      integer :: twave(ntrans)
      integer ::  n

!      call InitFFT1D(ntrans)
      call FFT1DB(spect1d, corrt)
      call RearrangeFFT1D(ntrans,corrt,corrout,twave)
      open(21, file='corrt'//trim(fn)//'.dat')
      write(21,'(a)')'variables=t,corr'
      write(21,'(a,I8.8,a,I3.3)') '# Data points per seg = ', ntrans, '   Num of seg = ', nsection
      write(21,'(a,I1.1,a,I1.1,a,E15.8)') '# Window type = ', iwindow, '   ioverlap = ', ioverlap, '   dt_sample = ', dts
      do n=1,ntrans
        write(21,'(2E20.12)')real(twave(n))*dts, real(corrout(n))/maxval(real(corrout))
      end do
      close(21)
    end subroutine CalCorrt

    subroutine CalCorrxt(ixhalflenl,ixhalflenr,dxs,spect1d,fn)
      implicit none
      integer, intent(in) :: ixhalflenl,ixhalflenr
      real(8), intent(in) :: dxs
      complex, intent(in) :: spect1d(ntrans,-ixhalflenl:ixhalflenr)
      character(*), intent(in) :: fn
      complex(8), dimension(ntrans) :: ctmp
      real(8), dimension(ntrans,-ixhalflenl:ixhalflenr) :: corrxt
      complex(8), dimension(ntrans,-ixhalflenl:ixhalflenr) :: corrout
      integer :: twave(ntrans)
      integer ::  n,i, iloc(1)
      real(8), dimension(ntrans) :: uc

      if(ixhalflenl.lt.0.or.ixhalflenr.lt.0) then
         print *, 'ixhalflenl, ixhalflenr need to be >= 0. STOP'
      endif
      do i=-ixhalflenl,ixhalflenr
         call FFT1DB(spect1d(:,i), ctmp)
         call RearrangeFFT1D(ntrans,ctmp,corrout(:,i),twave)
     enddo

      open(21, file='corrxt'//trim(fn)//'.dat')
      write(21,'(a)')'variables=t,x,corr'
      write(21,'("Zone T=corrxt, I=", I8, ",J=",I8)') ntrans, ixhalflenl+ixhalflenr+1
      write(21,'(a,I8.8,a,I3.3,a,E15.8)') '# Data points per seg = ', ntrans, '   Num of seg = ', nsection, '   maxval = ', maxval(abs(real(corrout)))
      write(21,'(a,I1.1,a,I1.1,a,E15.8)') '# Window type = ', iwindow, '   ioverlap = ', ioverlap, '   dt_sample = ', dts

!     corrxt = real(corrout)/maxval(abs(real(corrout)))
     print *, 'ixhalflenl =', ixhalflenl, 'ixhalflenr =', ixhalflenr
     print *, 'real(corrout(ntrans/2,0)) =', real(corrout(ntrans/2,0))
     corrxt = real(corrout)/real(corrout(ntrans/2,0))
     do i=-ixhalflenl,ixhalflenr
         do n=1,ntrans
            write(21,'(3E20.12)') real(twave(n))*dts, real(i)*dxs, corrxt(n,i)
         enddo
      enddo
      close(21)

     ! Calculate convection velocity
     ! Ref: Bernardini & Pirozzoli, PhysFluids, 23, 085102 (2011)
      open(21, file='uconv_Pirozzoli_corrxt_'//trim(fn)//'.dat')
      write(21,'(a)')'variables=t,ucorr'
      write(21,'(a,I8.8,a,I3.3)') '# Data points per seg = ', ntrans, '   Num of seg = ', nsection
      write(21,'(a,I1.1,a,I1.1,a,E15.8)') '# Window type = ', iwindow, '   ioverlap = ', ioverlap, '   dt_sample = ', dts
      do n = ntrans/2+1, ntrans
        iloc = maxloc(corrxt(n,:))
        uc(n) = real(iloc(1))*dxs/( real(twave(n))*dts )
        write(21,'(2E20.12)')real(twave(n))*dts, uc(n)
      enddo
      close(21)

    end subroutine CalCorrxt

    subroutine CalCorrzt(kzhalflend,kzhalflenu,dzs,spect1d,fn)
      implicit none
      integer, intent(in) :: kzhalflend,kzhalflenu
      real(8), intent(in) :: dzs
      complex, intent(in) :: spect1d(ntrans,-kzhalflend:kzhalflenu)
      character(*), intent(in) :: fn
      complex(8), dimension(ntrans) :: ctmp
      real(8), dimension(ntrans,-kzhalflend:kzhalflenu) :: corrxt
      complex(8), dimension(ntrans,-kzhalflend:kzhalflenu) :: corrout
      integer :: twave(ntrans)
      integer ::  n,i, kloc(1)
      real(8), dimension(ntrans) :: uc

      if(kzhalflend.lt.0.or.kzhalflenu.lt.0) then
         print *, 'kzhalflend, kzhalflenu need to be >= 0. STOP'
      endif
      do i=-kzhalflend,kzhalflenu
         call FFT1DB(spect1d(:,i), ctmp)
         call RearrangeFFT1D(ntrans,ctmp,corrout(:,i),twave)
     enddo

      open(21, file='corrzt'//trim(fn)//'.dat')
      write(21,'(a)')'variables=t,z,corr'
      write(21,'("Zone T=corrzt, I=", I8, ",J=",I8)') ntrans, kzhalflend+kzhalflenu+1
      write(21,'(a,I8.8,a,I3.3,a,E15.8)') '# Data points per seg = ', ntrans, '   Num of seg = ', nsection, '   maxval = ', maxval(abs(real(corrout)))
      write(21,'(a,I1.1,a,I1.1,a,E15.8)') '# Window type = ', iwindow, '   ioverlap = ', ioverlap, '   dt_sample = ', dts

!     corrxt = real(corrout)/maxval(abs(real(corrout)))
     print *, 'kzhalflend =', kzhalflend, 'kzhalflenu =', kzhalflenu
     print *, 'real(corrout(ntrans/2,0)) =', real(corrout(ntrans/2,0))
     corrxt = real(corrout)/real(corrout(ntrans/2,0))
     do i=-kzhalflend,kzhalflenu
         do n=1,ntrans
            write(21,'(3E20.12)') real(twave(n))*dts, real(i)*dzs, corrxt(n,i)
         enddo
      enddo
      close(21)

     ! Calculate convection velocity
     ! Ref: Bernardini & Pirozzoli, PhysFluids, 23, 085102 (2011)
      open(21, file='uconv_Pirozzoli_corrzt_'//trim(fn)//'.dat')
      write(21,'(a)')'variables=t,ucorr'
      write(21,'(a,I8.8,a,I3.3)') '# Data points per seg = ', ntrans, '   Num of seg = ', nsection
      write(21,'(a,I1.1,a,I1.1,a,E15.8)') '# Window type = ', iwindow, '   ioverlap = ', ioverlap, '   dt_sample = ', dts
      do n = ntrans/2+1, ntrans
        kloc = maxloc(corrxt(n,:))
        uc(n) = real(kloc(1))*dzs/( real(twave(n))*dts )
        write(21,'(2E20.12)')real(twave(n))*dts, uc(n)
      enddo
      close(21)

    end subroutine CalCorrzt

    subroutine CalCorrxtcross(ixhalflenl,ixhalflenr,dxs,spect12,ratio,fn)
      implicit none
      integer, intent(in) :: ixhalflenl,ixhalflenr
      real(8), intent(in) :: dxs
      complex, intent(in) :: spect12(ntrans,-ixhalflenl:ixhalflenr)
      real(8), intent(in) :: ratio(-ixhalflenl:ixhalflenr)
      character(*), intent(in) :: fn
      complex(8), dimension(ntrans) :: ctmp
      real(8), dimension(ntrans,-ixhalflenl:ixhalflenr) :: corrxt
      complex(8), dimension(ntrans,-ixhalflenl:ixhalflenr) :: corrout
      integer :: twave(ntrans)
      integer ::  n,i, iloc(1)
      real(8), dimension(ntrans) :: uc

      open(21, file='corrxtcross'//trim(fn)//'.dat')
      write(21,'(a)')'variables=t,x,corr'
      write(21,'("Zone T=corrxt, I=", I8, ",J=",I8)') ntrans, ixhalflenl+ixhalflenr+1
      write(21,'(a,I8.8,a,I3.3)') '# Data points per seg = ', ntrans, '   Num of seg = ', nsection
      write(21,'(a,I1.1,a,I1.1,a,E15.8)') '# Window type = ', iwindow, '   ioverlap = ', ioverlap, '   dt_sample = ', dts
      do i=-ixhalflenl,ixhalflenr
         call FFT1DB(spect12(:,i), ctmp)
         call RearrangeFFT1D(ntrans,ctmp,corrout(:,i),twave)
     enddo

     corrxt = 0.d0
     do i=-ixhalflenl,ixhalflenr
         if(abs(corrout(ntrans/2,i)).ne.0.d0) corrxt(:,i) = abs(corrout(:,i))/abs(corrout(ntrans/2,i))*ratio(i)
         do n=1,ntrans
            write(21,'(3E20.12)') real(twave(n))*dts, real(i)*dxs, corrxt(n,i)
         enddo
      enddo
      close(21)
!      print *, 'corrxt=',corrxt(ntrans/2,0)

    end subroutine CalCorrxtcross


    subroutine CalPhaseSpeed_VanAtta(spect1d,dxs,phasespeed,phase)
      ! Ref: Stegen & Van Atta, JFM, 42, pp 689-699, 1970
      implicit none
      complex, intent(in) :: spect1d(ntrans)
      real(8), intent(in) :: dxs
      real(8), dimension(ntrans/2), intent(out) :: phasespeed, phase
      integer :: i
      real(8) :: fperiod, ff, ang

      fperiod = dble(ntrans-1)*dts  ! period per segments
      phasespeed=0.; phase=0.
      ang = 0.
      do i=2,ntrans/2
        ff  = (i-1)/fperiod
        !phase(i) = atan2(imag(spect1d(i)),real(spect1d(i))) - 2.d0*pi*min(sign(1.,imag(spect1d(i))),0.)
        phase(i) = atan2(imag(spect1d(i)),real(spect1d(i))) - pi*min(sign(1.,imag(spect1d(i))),0.)
        if (abs(phase(i-1)-phase(i)).ge.1.0) then
          ang = ang+2.d0*pi
          ! print *, 'ang = ', ang
        endif
        if (abs(phase(i)+ang).gt.1.d-14) phasespeed(i)=2.*pi*ff*dxs/abs((phase(i)+ang))

      end do
    end subroutine CalPhaseSpeed_VanAtta


!    subroutine CalPhaseSpeed_VanAtta(spect1d,dxs,fn)
!      ! Ref: Stegen & Van Atta, JFM, 42, pp 689-699, 1970
!      implicit none
!      complex, intent(in) :: spect1d(ntrans)
!      real(8), intent(in) :: dxs
!      character(*), intent(in) :: fn
!      real(8), dimension(ntrans/2) :: phasespeed, phase
!      integer :: i
!      real(8) :: fperiod, ff, ang
!
!      fperiod = dble(ntrans-1)*dts  ! period per segments
!      print *, 'fperiod=',fperiod
!      phasespeed=0.; phase=0.
!      ang = 0.
!      open(21, file='phasespeed_VanAtta_'//trim(fn)//'.dat')
!      write(21,'(a)')'variables=f,phasespeed,phase'
!      write(21,'(a,I8.8,a,I3.3)') '# Data points per seg = ', ntrans, '   Num of seg = ', nsection
!      write(21,'(a,I1.1,a,I1.1,a,E15.8)') '# Window type = ', iwindow, '   ioverlap = ', ioverlap, '   dt_sample = ', dts
!      do i=2,ntrans/2
!        ff  = (i-1)/fperiod
!        phase(i) = atan2(imag(spect1d(i)),real(spect1d(i))) - 2.d0*pi*min(sign(1.,imag(spect1d(i))),0.)
!        if (abs(phase(i-1)-phase(i)).ge.1.0) ang = ang+2.d0*pi
!!        if(imag(spect1d(i)).lt.0.d0) ang = 2.*pi
!        if (abs(phase(i)+ang).gt.1.d-14) phasespeed(i)=2.*pi*ff*dxs/abs((phase(i)+ang))
!!        if (abs(real(spect1d(i))).gt.1.d-14) phase(i) = atan(imag(spect1d(i))/real(spect1d(i)))
!
!!        !phase can be larger than pi/2 due to large distance
!!        if (phase(i-1).ge.0.and.phase(i).lt.0) ang = ang+pi
!!        if (abs(phase(i)+ang).gt.1.d-14) phasespeed(i)=2.*pi*ff*dxs/abs((phase(i)+ang))
!!        if (abs(phase(i)+ang).gt.1.d-14) phasespeed(i)=2.*pi*ff*dxs/abs((phase(i)+ang))
!        write(21,'(3E20.12)')ff, phasespeed(i), phase(i)
!      end do
!      close(21)
!    end subroutine CalPhaseSpeed_VanAtta

    subroutine CalPhaseSpeed_Jimenez(spect1d,spect1dx,fn)
      ! Ref: Juan C. Alamo & Javier Jimenez, JFM, 640, pp 5-26, 2009
      implicit none
      complex, intent(in) :: spect1d(ntrans), spect1dx(ntrans)
      character(*), intent(in) :: fn
      real(8) :: phasespeed
      integer :: i
      real(8) :: fperiod, ff

      fperiod = dble(ntrans-1)*dts  ! period per segments
      print *, 'fperiod=',fperiod
      open(21, file='phasespeed_Jimenez_'//trim(fn)//'.dat')
      write(21,'(a)')'variables=f,phasespeed'
      write(21,'(a,I8.8,a,I3.3)') '# Data points per seg = ', ntrans, '   Num of seg = ', nsection
      write(21,'(a,I1.1,a,I1.1,a,E15.8)') '# Window type = ', iwindow, '   ioverlap = ', ioverlap, '   dt_sample = ', dts
      do i=2,ntrans/2
        ff  = (i-1)/fperiod
        phasespeed = -2.d0*pi*ff*(spect1d(i)*conjg(spect1d(i)))/imag(spect1dx(i)*conjg(spect1d(i)))
        write(21,'(3E20.12)')ff, phasespeed
      end do
      close(21)

    end subroutine CalPhaseSpeed_Jimenez

end module modSpect

