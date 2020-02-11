
!>
module modHelmholtzDecomp
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
   logical :: IsHDecompInit = .FALSE.

   real(8), dimension(:), allocatable, private :: k1, k2, k3
   real(8), dimension(:,:,:), allocatable, private :: ktotsq
   complex, dimension(:,:,:,:), allocatable, private :: buffer_var_c, buffer_sole_c

contains

   subroutine InitHDecomp3d(nx, ny, nz, iwin, dx_sample, dy_sample, dz_sample)
        implicit none
        integer, intent(in) :: nx, ny, nz, iwin
        real(8), intent(in) :: dx_sample, dy_sample, dz_sample
        integer :: i, j, k

        pi = 4.*atan(1.)
        iwindow = iwin
        nxpoint = nx;         nypoint = ny;         nzpoint = nz
        dxs = dx_sample;      dys = dy_sample;      dzs = dz_sample
        nperiod_kx = nxpoint; nperiod_ky = nypoint; nperiod_kz = nzpoint

        print *, '###############################################'
        print *, '# of data points per segments in space domain kx = ', nxpoint
        print *, '# of segments in space domain = 1 '
        print *, 'Interval of sampling (m) =', dxs
        print *, 'Length per segments (m) =', dble(nxpoint)*dxs

        print *, '###############################################'
        print *, '# of data points per segments in space domain ky = ', nypoint
        print *, '# of segments in space domain = 1 '
        print *, 'Interval of sampling (m) =', dys
        print *, 'Length per segments (m) =', dble(nypoint)*dys

        print *, '###############################################'
        print *, '# of data points per segments in space domain kz = ', nzpoint
        print *, '# of segments in space domain = 1 '
        print *, 'Interval of sampling (m) =', dzs
        print *, 'Length per segments (m) =', dble(nzpoint)*dzs

        call InitFFT3D(nxpoint,nypoint,nzpoint)
        allocate(wcoeff_kx(nxpoint), wcoeff_ky(nypoint), wcoeff_kz(nzpoint))
        call InitFFTWindow_1(nxpoint,iwindow,wcoeff_kx)
        call InitFFTWindow_1(nypoint,iwindow,wcoeff_ky)
        call InitFFTWindow_1(nzpoint,iwindow,wcoeff_kz)

        ! calculate wavenumber in i,j,k directions
        allocate(k1(nxpoint),k2(nypoint),k3(nzpoint))
        allocate(ktotsq(nxpoint,nypoint,nzpoint))
        do i=1, nxpoint/2
          !k1(i) = dble(i-1)/dble(nxpoint)
          k1(i) = dble(i-1)/dble(nxpoint*dxs)
        enddo
        do i=nxpoint/2+1,nxpoint
          !k1(i) = dble(i-nxpoint-1)/dble(nxpoint)
          k1(i) = dble(i-nxpoint-1)/dble(nxpoint*dxs)
        enddo
        do j=1, nypoint/2
          !k2(j) = dble(j-1)/dble(nypoint)
          k2(j) = dble(j-1)/dble(nypoint*dys)
        enddo
        do j=nypoint/2+1,nypoint
          !k2(j) = dble(j-nypoint-1)/dble(nypoint)
          k2(j) = dble(j-nypoint-1)/dble(nypoint*dys)
        enddo
        do k=1, nzpoint/2
          !k3(k) = dble(k-1)/dble(nzpoint)
          k3(k) = dble(k-1)/dble(nzpoint*dzs)
        enddo
        do k=nzpoint/2+1, nzpoint
          !k3(k) = dble(k-nzpoint-1)/dble(nzpoint)
          k3(k) = dble(k-nzpoint-1)/dble(nzpoint*dzs)
        enddo
        print *, 'Writing wavenumber files: k1.dat, k2.dat and k3.dat ... '
        open(77,file='k1.dat',status='unknown')
          write(77,*) 'variables=i,k1'
          do i=1, nxpoint
            write(77,*) i, k1(i)
          enddo
        close(77)
        open(77,file='k2.dat',status='unknown')
          write(77,*) 'variables=i,k2'
          do j=1, nypoint
            write(77,*) j, k2(j)
          enddo
        close(77)
        open(77,file='k3.dat',status='unknown')
          write(77,*) 'variables=k,k3'
          do k=1, nzpoint
            write(77,*) k, k3(k)
          enddo
        close(77)
        do i=1, nxpoint
        do j=1, nypoint
        do k=1, nzpoint
          ktotsq(i,j,k) = k1(i)**2 + k2(j)**2 + k3(k)**2 + 1.e-30
        enddo
        enddo
        enddo

        allocate(buffer_var_c(nxpoint,nypoint,nzpoint,3), buffer_sole_c(nxpoint,nypoint,nzpoint,3))

        IsHDecompInit = .TRUE.
    end subroutine InitHDecomp3d

subroutine FFTFilter2(ns,varin,dts,fl,fu,varout)
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
    real(8) :: varave, var2, varrms
    real(8) :: varave_out, var2_out, varrms_out

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

!    ! original data
!    varave = sum(varin)/dble(npoint)
!    var2 = sum(varin**2)/dble(npoint)
!    varrms = sqrt(abs(var2-varave**2))
!
!    varave_out = sum(varout)/dble(npoint)
!    var2_out = sum(varout**2)/dble(npoint)
!    varrms_out = sqrt(abs(var2_out - varave_out**2))
!
!    varout = varout/varrms_out*varrms

end subroutine FFTFilter2

    subroutine FFTFilter3d_2(tfun,varout,fl_x,fu_x,fl_y,fu_y,fl_z,fu_z)
      implicit none
      real(8), intent(in) :: tfun(nxpoint,nypoint,nzpoint)
      real(8), intent(out) :: varout(nxpoint,nypoint,nzpoint)
      real(8), intent(in) :: fl_x, fu_x, fl_y, fu_y, fl_z, fu_z
      integer :: i, j, k
      complex, dimension(nxpoint,nypoint,nzpoint) :: ctmp, ctmp_2
      real(8), dimension(:), allocatable :: tmp1d
      real(8) :: tmean

      tmean = sum(tfun(:,:,:))/dble(nxpoint*nypoint*nzpoint)
      varout = tfun
      ! filter in x-direction
!      allocate(tmp1d(nxpoint))
!      do j=1, nypoint
!        do k=1, nzpoint
!          call FFTFilter2(nxpoint,varout(:,j,k),dxs,fl_x,fu_x,tmp1d)
!          varout(:,j,k) = tmp1d(:)
!        enddo
!      enddo
!      deallocate(tmp1d)
      ! filter in y-direction
      allocate(tmp1d(nypoint))
      do i=1, nxpoint
        do k=1, nzpoint
          call FFTFilter2(nypoint,varout(i,:,k),dys,fl_y,fu_y,tmp1d)
          varout(i,:,k) = tmp1d(:)
        enddo
      enddo
      deallocate(tmp1d)
      ! filter in z-direction
      allocate(tmp1d(nzpoint))
      do i=1, nxpoint
        do j=1, nypoint
          call FFTFilter2(nzpoint,varout(i,j,:),dzs,fl_z,fu_z,tmp1d)
          varout(i,j,:) = tmp1d(:)
        enddo
      enddo
      deallocate(tmp1d)

      varout = varout + tmean

    end subroutine FFTFilter3d_2

    subroutine FFTFilter3d_3(tfun,varout)
      implicit none
      real(8), intent(in) :: tfun(nxpoint,nypoint,nzpoint)
      real(8), intent(out) :: varout(nxpoint,nypoint,nzpoint)
      complex, dimension(:,:,:), allocatable :: buffer_c
      integer :: i, j, k

      allocate(buffer_c(nxpoint,nypoint,nzpoint))

      buffer_c = forward3d(tfun)

      i = nxpoint/2
      j = nypoint/2
      k = nzpoint/2
      print *, 'i = ', i, 'j = ', j, 'k = ', k
      buffer_c(i-105:i+105,j-108:j+108,k-108:k+108 ) = cmplx(0.d0,0.d0)
!      !buffer_c = 0.

      varout = dble(backward3d(buffer_c))

    end subroutine FFTFilter3d_3


    subroutine DoHDecomp3d(buffer_var, buffer_sole, buffer_dila)
      implicit none
      real(8), intent(in) ::  buffer_var(nxpoint,nypoint,nzpoint,3)
      real(8), intent(out):: buffer_sole(nxpoint,nypoint,nzpoint,3), buffer_dila(nxpoint,nypoint,nzpoint,3)
      integer :: i, j, k, n

      do n=1, 3
        buffer_var_c(:,:,:,n) = forward3d(buffer_var(:,:,:,n))
      enddo
      do i=1, nxpoint
        do j=1, nypoint
          do k=1, nzpoint
            buffer_sole_c(i,j,k,1) = buffer_var_c(i,j,k,1)*( 1.d0 - k1(i)*k1(i)/ktotsq(i,j,k) ) +  buffer_var_c(i,j,k,2)*(      - k1(i)*k2(j)/ktotsq(i,j,k) ) + buffer_var_c(i,j,k,3)*(      - k1(i)*k3(k)/ktotsq(i,j,k) )
            buffer_sole_c(i,j,k,2) = buffer_var_c(i,j,k,1)*(      - k1(i)*k2(j)/ktotsq(i,j,k) ) +  buffer_var_c(i,j,k,2)*( 1.d0 - k2(j)*k2(j)/ktotsq(i,j,k) ) + buffer_var_c(i,j,k,3)*(      - k2(j)*k3(k)/ktotsq(i,j,k) )
            buffer_sole_c(i,j,k,3) = buffer_var_c(i,j,k,1)*(      - k1(i)*k3(k)/ktotsq(i,j,k) ) +  buffer_var_c(i,j,k,2)*(      - k2(j)*k3(k)/ktotsq(i,j,k) ) + buffer_var_c(i,j,k,3)*( 1.d0 - k3(k)*k3(k)/ktotsq(i,j,k) )
          enddo
        enddo
      enddo
      do n=1, 3
        buffer_sole(:,:,:,n) = dble( backward3d(buffer_sole_c(:,:,:,n)))
      enddo
      do n=1, 3
        buffer_dila(:,:,:,n) = buffer_var(:,:,:,n) - buffer_sole(:,:,:,n)
      enddo

    end subroutine DoHDecomp3d

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
      tfun_tmp(:,:,:) = tfun(:,:,:)            ! for test

      ctmp = tfun_tmp
      call FFT3DF(ctmp,forward3d)

    end function forward3d


    function backward3d(tfun)
      implicit none
      complex, intent(in) :: tfun(nxpoint,nypoint,nzpoint)
      complex :: backward3d(nxpoint,nypoint,nzpoint)

      call FFT3DB(tfun,backward3d)

    end function backward3d

    function backward3d2(tfun)
      implicit none
      complex, intent(inout) :: tfun(nxpoint,nypoint,nzpoint)
      complex :: backward3d2(nxpoint,nypoint,nzpoint)
      integer :: i, j, k

      i = nxpoint/2
      j = nypoint/2
      k = nzpoint/2

      tfun(i-100:i+100,j-100:j+100,k-100:k+100 ) = cmplx(0.d0,0.d0)
      !buffer_c = 0.

      call FFT3DB(tfun,backward3d2)

    end function backward3d2


end module modHelmholtzDecomp
