module modSpect
   use MFFT1D
   use MFFT2D
   use MFFT3D
   use MFFTWindow
   implicit none

   integer, private :: iwindow, ioverlap, nperiod_kx, nperiod_ky, nperiod_kz
   integer, private :: nxpoint, nypoint, nzpoint
   real(8), private :: dxs, dys, dzs
   real(8), private :: pi
   real(8), dimension(:), allocatable, private :: wcoeff_kx, wcoeff_ky, wcoeff_kz
   logical :: IsSpectInit = .FALSE.
!                  nsection (# of bins)
!                  iwindow ( windowing type)
!                  ioverlap (whether or not 1/2 overlap between neighboring bins)
!                  ntpoint (total number of points), ntrans (# of points per bin)
!                  ntrans (size of bins used for Fourier transformation)
!                  nsection (number of bins; nsection = 8 for standard Welch method)
contains

     subroutine InitSpect3d(nx, ny, nz, iwin, dx_sample, dy_sample, dz_sample) !, npbin_kx, npbin_ky, npbin_kz)
        implicit none
        integer, intent(in) :: nx, ny, nz, iwin
        real(8), intent(in) :: dx_sample, dy_sample, dz_sample
 
        pi = 4.*atan(1.)
        iwindow = iwin
        nxpoint = nx;         nypoint = ny;         nzpoint = nz 
        dxs = dx_sample;      dys = dy_sample;      dzs = dz_sample
        nperiod_kx = nxpoint; nperiod_ky = nypoint; nperiod_kz = nzpoint
        !ntrans_kx  = nperiod_kx; ntrans_ky = nperiod_ky; ntrans_kz = nperiod_kz;
        !npbin_kx   = ntrans_kx;   npbin_ky = ntrans_ky;   npbin_kz = ntrans_kz

        print *, '###############################################'
        print *, '# of data points per segments in space domain kx = ', nxpoint
        print *, '# of segments in space domain = 1 '
        print *, 'Interval of sampling (m) =', dxs
        print *, 'Length per segments (m) =', dble(nxpoint)*dxs
        !print *, 'Segments with no overlap'
 

        print *, '###############################################'
        print *, '# of data points per segments in space domain ky = ', nypoint
        print *, '# of segments in space domain = 1 '
        print *, 'Interval of sampling (m) =', dys
        print *, 'Length per segments (m) =', dble(nypoint)*dys
        !print *, 'Segments with no overlap'

        print *, '###############################################'
        print *, '# of data points per segments in space domain kz = ', nzpoint
        print *, '# of segments in space domain = 1 '
        print *, 'Interval of sampling (m) =', dzs
        print *, 'Length per segments (m) =', dble(nzpoint)*dzs
        !print *, 'Segments with no overlap'

        !call InitFFT3D(nzpoint,nypoint,nxpoint)
        call InitFFT3D(nxpoint,nypoint,nzpoint)
        allocate(wcoeff_kx(nxpoint), wcoeff_ky(nypoint), wcoeff_kz(nzpoint))
        call InitFFTWindow_1(nxpoint,iwindow,wcoeff_kx)
        call InitFFTWindow_1(nypoint,iwindow,wcoeff_ky)
        call InitFFTWindow_1(nzpoint,iwindow,wcoeff_kz)

        IsSpectInit = .TRUE.
    end subroutine InitSpect3d

subroutine FFTFilter(ns,varin,dts,fl,fu,varout)
    implicit none
    integer, intent(in) :: ns
    real(8), intent(in) :: varin(ns)
    real(8), intent(in) :: dts, fl, fu
    real(8), intent(out) :: varout(ns)
    complex :: ctmp(ns)!, ctmp2(ns)
    real(8) :: ff, fperiod
    integer :: n, i, npoint
    integer :: lensave, lenwork, ierr
    real, dimension(:), allocatable :: wsave, wwork

    npoint = ns
    lensave = 2*npoint+int(log(real(npoint))/log(2.))+4
    lenwork = 2*npoint
    fperiod = dble(npoint)*dts

    if(allocated(wsave)) deallocate(wsave)
    if(allocated(wwork)) deallocate(wwork)
    allocate(wsave(lensave),wwork(lenwork))
    ! initalize fft
    call cfft1i(npoint,wsave,lensave,ierr)
    if(ierr.ne.0) then
      print*, 'FFT1D initialization error!'
      stop
    endif

    ctmp(1:npoint) = varin(1:npoint)
    call cfft1f(npoint,1,ctmp,npoint,wsave,lensave,wwork,lenwork,ierr)
    if(ierr.ne.0) then
      print*, 'Error occured in FFT1DF!'
      stop
    endif

    do i=1, npoint/2
      ff = dble(i)/fperiod
      if(ff.le.fl.or.ff.ge.fu) then
        ctmp(i) = 0.
      endif
    enddo
    do i=npoint/2+1, npoint
      ff =  abs(dble(i-npoint-1)/fperiod)
      if(ff.le.fl.or.ff.ge.fu) then
        ctmp(i) = 0.
      endif
    enddo

    call cfft1b(npoint,1,ctmp,npoint,wsave,lensave,wwork,lenwork,ierr)
    if(ierr.ne.0) then
      print*, 'Error occured in FFT1DB!'
      stop
    endif
    varout = real(ctmp)

end subroutine FFTFilter

    subroutine FFTFilter3d(tfun,varout)
      implicit none
      real(8), intent(in) :: tfun(nxpoint,nypoint,nzpoint)
      real(8), intent(out) :: varout(nxpoint,nypoint,nzpoint)
      integer :: i, j, k
      real(8) :: fp_x, fp_y, fp_z, fl_x, fl_y, fl_z, fu_x, fu_y, fu_z, ff
      complex, dimension(nxpoint,nypoint,nzpoint) :: ctmp, ctmp_2
      real(8), dimension(:), allocatable :: tmp1d
      real(8) :: tmean

      tmean = sum(tfun(:,:,:))/dble(nxpoint*nypoint*nzpoint)
      fp_x = dble(nxpoint)*dxs; fp_y = dble(nypoint)*dys; fp_z = dble(nzpoint)*dzs
      fl_x = 1.d0/fp_x;         fl_y = 1.d0/fp_y;         fl_z = 1.d0/fp_z
      fu_x = 1.d0/(12.d0*dxs);  fu_y = 1.d0/(12.d0*dys);  fu_z = 1.d0/(12.d0*dzs)

      ctmp = tfun - tmean
      call FFT3DF(ctmp,ctmp_2)
      do j=1, nypoint
        do k=1, nzpoint
          do i=1, nxpoint/2
            ff = dble(i)/fp_x
            if(ff.le.fl_x.or.ff.ge.fu_x) then
              ctmp_2(i,j,k) = 0.
            endif
          enddo
          do i=nxpoint/2+1, nxpoint
            ff = abs(dble(i-nxpoint-1))/fp_x
            if(ff.le.fl_x.or.ff.ge.fu_x) then
              ctmp_2(i,j,k) = 0.
            endif
          enddo
        enddo
      enddo

      do i=1, nxpoint
        do k=1, nzpoint
          do j=1, nypoint/2
            ff = dble(j)/fp_y
            if(ff.le.fl_y.or.ff.ge.fu_y) then
              ctmp_2(i,j,k) = 0.
            endif
          enddo
          do j= nypoint/2+1 , nypoint
            ff = abs(dble(j-nypoint-1) )/fp_y
            if(ff.le.fl_y.or.ff.ge.fu_y) then
              ctmp_2(i,j,k) = 0.
            endif
          enddo
        enddo
      enddo

      do i=1, nxpoint
        do j=1, nypoint
          do k=1, nzpoint/2
            ff = dble(k)/fp_z
            if(ff.le.fl_z.or.ff.ge.fu_z) then
              ctmp_2(i,j,k) = 0.
            endif
          enddo
          do k=nzpoint/2+1, nzpoint
            ff = abs(dble(k-nzpoint-1))/fp_z
            if(ff.le.fl_z.or.ff.ge.fu_z) then
              ctmp_2(i,j,k) = 0.
            endif
          enddo
        enddo
      enddo

      call FFT3DB(ctmp_2,ctmp)
      varout = dble(ctmp) + tmean

    end subroutine FFTFilter3d

    subroutine FFTFilter3d_2(tfun,varout)
      implicit none
      real(8), intent(in) :: tfun(nxpoint,nypoint,nzpoint)
      real(8), intent(out) :: varout(nxpoint,nypoint,nzpoint)
      integer :: i, j, k
      real(8) :: fp_x, fp_y, fp_z, fl_x, fl_y, fl_z, fu_x, fu_y, fu_z
      complex, dimension(nxpoint,nypoint,nzpoint) :: ctmp, ctmp_2
      real(8), dimension(:), allocatable :: tmp1d
      real(8) :: tmean

      tmean = sum(tfun(:,:,:))/dble(nxpoint*nypoint*nzpoint)
      fp_x = dble(nxpoint)*dxs; fp_y = dble(nypoint)*dys; fp_z = dble(nzpoint)*dzs
      fl_x = 1.d0/fp_x;         fl_y = 1.d0/fp_y;         fl_z = 1.d0/fp_z
      fu_x = 1.d0/(12.d0*dxs);  fu_y = 1.d0/(12.d0*dys);  fu_z = 1.d0/(12.d0*dzs)

      varout = tfun
      ! filter in x-direction
      allocate(tmp1d(nxpoint))
      do j=1, nypoint
        do k=1, nzpoint
          call FFTFilter(nxpoint,varout(:,j,k),dxs,fl_x,fu_x,tmp1d)
          varout(:,j,k) = tmp1d(:)
        enddo
      enddo
      deallocate(tmp1d)
      ! filter in y-direction
      allocate(tmp1d(nypoint))
      do i=1, nxpoint
        do k=1, nzpoint
          call FFTFilter(nypoint,varout(i,:,k),dys,fl_y,fu_y,tmp1d)
          varout(i,:,k) = tmp1d(:)
        enddo
      enddo
      deallocate(tmp1d)
      ! filter in z-direction
      allocate(tmp1d(nzpoint))
      do i=1, nxpoint
        do j=1, nypoint
          call FFTFilter(nzpoint,varout(i,j,:),dzs,fl_z,fu_z,tmp1d)
          varout(i,j,:) = tmp1d(:)
        enddo
      enddo
      deallocate(tmp1d)

      varout = varout + tmean

    end subroutine FFTFilter3d_2


    function transfer3d(tfun)
      implicit none
      real(8), intent(in) :: tfun(nxpoint,nypoint,nzpoint)
      real(8) :: transfer3d(nxpoint,nypoint,nzpoint)
      integer :: i, j, jj, k
      complex, dimension(nxpoint,nypoint,nzpoint) :: ctmp, ctmp_2
      real(8) :: tmean

      tmean = sum(tfun(:,:,:))/dble(nxpoint*nypoint*nzpoint)
      do k=1, nzpoint
        do i=1, nxpoint
            ctmp(i,1:nypoint,k) = wcoeff_kz(k)*wcoeff_kx(i)*wcoeff_ky(1:nypoint)*(tfun(i,1:nypoint,k)-tmean)
        enddo 
      enddo 
      call FFT3DF(ctmp,ctmp_2)
      call FFT3DB(ctmp_2,ctmp)
      transfer3d = dble(ctmp) + tmean
    end function transfer3d
 
    function forward3d(tfun)
      implicit none
      real(8), intent(in) :: tfun(nxpoint,nypoint,nzpoint)
      real(8) :: tfun_tmp(nxpoint,nypoint,nzpoint)
      complex :: forward3d(nxpoint,nypoint,nzpoint)
      real(8) :: tmean
      integer :: nn, mm, i, j, jj, k
      complex, dimension(nxpoint,nypoint,nzpoint) :: ctmp
 
      tmean = sum(tfun(:,:,:))/dble(nxpoint*nypoint*nzpoint)

      !tfun_tmp(:,:,:) = tfun(:,:,:) - tmean
      tfun_tmp(:,:,:) = tfun(:,:,:)  ! for test

      ! do k=1, nzpoint
      !  do i=1, nxpoint
      !      ctmp(i,1:nypoint,k) = wcoeff_kz(k)*wcoeff_kx(i)*wcoeff_ky(1:nypoint)*(tfun_tmp(i,1:nypoint,k))
      !  enddo 
      !enddo 
      ctmp = tfun_tmp

      call FFT3DF(ctmp,forward3d)

    end function forward3d


    function backward3d(tfun)
      implicit none
      complex, intent(in) :: tfun(nxpoint,nypoint,nzpoint)
      complex :: backward3d(nxpoint,nypoint,nzpoint)
  
      call FFT3DB(tfun,backward3d)

    end function backward3d

    function spect3d(tfun)
      implicit none
      real(8), intent(in) :: tfun(nxpoint,nypoint,nzpoint)
      real(8) :: tfun_tmp(nxpoint,nypoint,nzpoint)
      complex :: spect3d(nxpoint,nypoint,nzpoint)
      real(8) :: tmean
      integer :: nn, mm, i, j, jj, k
      complex, dimension(nxpoint,nypoint,nzpoint) :: ctmp

       ctmp = tfun

      !call FFT3DF(ctmp,ctmp)
      call FFT3DF(ctmp,ctmp)
      !spect3d(:,:,:) = ctmp(:,:,:)

      call FFT3DB(ctmp,spect3d)

    end function spect3d

!    function spect3d(tfun)
!      implicit none
!      real(8), intent(in) :: tfun(nxpoint,nypoint,nzpoint)
!      real(8) :: tfun_tmp(nxpoint,nypoint,nzpoint)
!      complex :: spect3d(nxpoint,nypoint,nzpoint)
!      real(8) :: tmean
!      integer :: nn, mm, i, j, jj, k
!      complex, dimension(nxpoint,nypoint,nzpoint) :: ctmp
! 
!      tmean = sum(tfun(:,:,:))/dble(nxpoint*nypoint*nzpoint)
!
!print *, 'tmean = ', tmean
!      tfun_tmp(:,:,:) = tfun(:,:,:) ! - tmean
! 
!      !do k=1, nzpoint
!      !  do i=1, nxpoint
!      !      ctmp(k,i,1:nypoint) = wcoeff_kz(k)*wcoeff_kx(i)*wcoeff_ky(1:nypoint)*(tfun_tmp( k,i,1:nypoint))
!      !  enddo 
!      !enddo 
!      ctmp = tfun
!
!      call FFT3DF(ctmp,ctmp)
!      !spect3d(:,:,:) = ctmp(:,:,:)
!      call FFT3DB(ctmp,spect3d)
!
!    end function spect3d


!    function autospect1d(tfun)
!      implicit none
!      real(8), intent(in) :: tfun(ntpoint)
!      complex :: autospect1d(ntrans)
!      complex :: ctmp(ntrans,nsection)
!
!      ctmp = spect1d(tfun)
!      autospect1d = sum(ctmp*conjg(ctmp),dim=2)/dble(nsection)
!    end function autospect1d

     subroutine InitSpect2d(nx, ny, iwin, dx_sample, dy_sample)
        implicit none
        integer, intent(in) :: nx, ny,iwin
        real(8), intent(in) :: dx_sample, dy_sample

        nxpoint = nx
        nypoint = ny
        iwindow = iwin
        dxs = dx_sample
        dys = dy_sample
 
        pi = 4.*atan(1.)
        
        print *, '###############################################'
        print *, '# of data points per segments in space domain', nxpoint
        print *, '# of segments in space domain = 1 '
        print *, 'Interval of sampling (m) =', dxs
        print *, 'Length per segments (m) =', (nxpoint)*dxs
        print *, 'Segments with no overlap'
        print *, '###############################################'
        print *, '# of data points per segments in space domain', nypoint
        print *, '# of segments in space domain = 1 '
        print *, 'Interval of sampling (m) =', dys
        print *, 'Length per segments (m) =', (nypoint)*dys
        print *, 'Segments with no overlap'

        call InitFFT2D(nxpoint,nypoint)
        allocate(wcoeff_kx(nxpoint), wcoeff_ky(nypoint))
        call InitFFTWindow_1(nxpoint,iwindow,wcoeff_kx)
        call InitFFTWindow_1(nypoint,iwindow,wcoeff_ky)
        IsSpectInit = .TRUE.
    end subroutine InitSpect2d
 
    function spect2d(tfun)
      implicit none
      real(8), intent(in) :: tfun(nxpoint,nypoint)
      complex :: spect2d(nxpoint,nypoint)
      real(8) :: tmean
      integer :: nn, j
      complex, dimension(nxpoint,nypoint) :: ctmp

      tmean = sum(tfun(1:nxpoint,1:nypoint))/dble(nxpoint*nypoint)

print *, 'tmean = ', tmean
      ctmp(:,:) = tfun(:,:) !-tmean
      call FFT2DF(ctmp,ctmp)
      spect2d(:,:) = ctmp(:,:)
 
    end function spect2d
  
   
 

end module modSpect

