module modSpecty
   use MFFT1D
   use MFFT2D
   use MFFT3D
   use MFFTWindow
   implicit none
   integer, private :: nypoint, nxpoint, iwindow
   real(8), private :: dys, dxs
   real(8), private :: pi
   real(8), dimension(:), allocatable, private :: wcoeff_kx
   logical, save :: IsCorryInit =  .false.
!                  iwindow ( windowing type)

contains
     subroutine InitCorry(ny, dy_sample)
        implicit none
        integer, intent(in) :: ny
        real(8), intent(in) :: dy_sample

        nypoint = ny
        dys = dy_sample
        pi = 4.*atan(1.)
!        call InitFFT1D(nypoint)
        IsCorryInit = .true.

     end subroutine InitCorry

     function spect1dy(tfun)
      implicit none
      real(8), intent(in) :: tfun(nypoint)
      complex :: spect1dy(nypoint)
      real(8) :: tmean
      integer :: i, j, n, nn, nt
      complex, dimension(nypoint) :: ctmp
      real(8) :: tbin(nypoint)

      ! Calculate the mean and energy of time signal
      tmean = sum(tfun(1:nypoint))/dble(nypoint)

      ! Compute spectrum for time data length of ntans
      tbin(1:nypoint) = tfun(1:nypoint) - tmean
      ctmp(1:nypoint) = tbin(1:nypoint)
      call FFT1DF(ctmp, ctmp)
      spect1dy =  ctmp

    end function spect1dy

     subroutine InitSpect2d_ky_kx(ny, nx, iwin, dy_sample, dx_sample)
        implicit none
        integer, intent(in) :: ny, nx, iwin
        real(8), intent(in) :: dy_sample, dx_sample

        nypoint = ny; nxpoint = nx; iwindow = iwin
        dys = dy_sample; dxs = dx_sample
        pi = 4.*atan(1.)

        print *, '###############################################'
        print *, '# of data points per segments in ky-domain =', nypoint
        print *, '# of segments in space domain = 1 '
        print *, 'Interval of sampling (m) =', dys
        print *, 'Length per segments (m) =', (nypoint)*dys
        print *, 'Segments with no overlap'

        print *, '###############################################'
        print *, '# of data points per segments in kx-domain', nxpoint
        print *, '# of segments in space domain = 1 '
        print *, 'Interval of sampling (m) =', dxs
        print *, 'Length per segments (m) =', (nxpoint)*dxs
        print *, 'Segments with no overlap'
        print *, '###############################################'

        call InitFFT2D(nypoint,nxpoint)
        call InitFFTWindow(nxpoint,iwindow)

        !IsSpectInit = .TRUE.
    end subroutine InitSpect2d_ky_kx


   function spect2d(tfun)
      implicit none
      real(8), intent(in) :: tfun(nypoint,nxpoint)
      complex :: spect2d(nypoint,nxpoint)
      real(8) :: tmean
      integer :: j
      complex, dimension(nypoint,nxpoint) :: ctmp

      tmean = sum(tfun(1:nypoint,1:nxpoint))/dble(nypoint*nxpoint)

      do j=1, nypoint
        ctmp(j,1:nxpoint) = wcoeff(1:nxpoint)*(tfun(j,1:nxpoint)-tmean)
      enddo
      call FFT2DF(ctmp,ctmp)
      spect2d(1:nypoint,1:nxpoint) = ctmp(1:nypoint,1:nxpoint)

    end function spect2d

    function autospect2d_ky_kx(tfun)
      implicit none
      real(8), intent(in) :: tfun(nypoint,nxpoint)
      complex :: autospect2d_ky_kx(nypoint,nxpoint)
      complex :: ctmp(nypoint,nxpoint)

      ctmp = spect2d(tfun)
      autospect2d_ky_kx = ctmp*conjg(ctmp)
    end function autospect2d_ky_kx


    function autospect1dy(tfun)
      implicit none
      real(8), intent(in) :: tfun(nypoint)
      complex :: autospect1dy(nypoint)
      complex :: ctmp(nypoint)

      ctmp = spect1dy(tfun)
      autospect1dy = ctmp*conjg(ctmp)
    end function autospect1dy

    function crossspect1dy(tfun1,tfun2)
      implicit none
      real(8), intent(in) :: tfun1(nypoint), tfun2(nypoint)
      complex :: crossspect1dy(nypoint)
      complex, dimension(nypoint) :: ctmp1, ctmp2

      ctmp1 = spect1dy(tfun1)
      ctmp2 = spect1dy(tfun2)
      crossspect1dy = ctmp1*conjg(ctmp2)

    end function crossspect1dy

    subroutine CalPowerspect_kx_ky(spect2d,tmeanenergy,fn)
      implicit none
      complex, intent(in) :: spect2d(nypoint,nxpoint)
      real(8), intent(in) :: tmeanenergy
      character(*), intent(in) :: fn
      complex :: powerspect(nypoint,nxpoint)
      real(8) :: yy, xx, kyperiod, kxperiod, tmpenergy
      integer :: ywave(nypoint), xwave(nxpoint)
      integer :: i, n

      kyperiod = dble(nypoint)*dys
      kxperiod = dble(nxpoint)*dxs
      call RearrangeFFT2D(nypoint,nxpoint,spect2d,powerspect,ywave,xwave)
      tmpenergy = sum(powerspect)/kyperiod/kxperiod

print *, 'sum(powerspect) = ', real(sum(powerspect))

      open(21, file=trim(fn)//'.dat')
      rewind(21)
      write(21,'(a)') 'variables=ky,kx,spectrum'
      write(21,'("Zone T=ky-t_PSD, I=", I8, ",J=", I8, ", AUXDATA rms2=""", E16.9,"""")') nypoint, nxpoint, tmeanenergy

      do i=1, nxpoint
        do n=1, nypoint
          yy = 2.d0*pi*real(ywave(n))/kyperiod
          xx = 2.d0*pi*real(xwave(i))/kxperiod
          write(21,'(3E20.12)') yy, xx, real(powerspect(n,i))*tmeanenergy/tmpenergy/(2.d0*pi)/(2.d0*pi)
        enddo
      enddo
      close(21)

    end subroutine CalPowerspect_kx_ky


    subroutine CalCorrxy(ixhalflenl,ixhalflenr,dxs,spect1d,fn)
      implicit none
      integer, intent(in) :: ixhalflenl, ixhalflenr
      real(8), intent(in) :: dxs
      complex, intent(in) :: spect1d(-ixhalflenl:ixhalflenr,nypoint)
      character(*), intent(in) :: fn
      complex(8), dimension(nypoint) :: ctmp
      real(8), dimension(-ixhalflenl:ixhalflenr,nypoint) :: corrxy
      complex(8), dimension(-ixhalflenl:ixhalflenr,nypoint) :: corrout
      integer :: ywave(nypoint)
      integer ::  i, j

      if(.not.IsCorryInit) then
         print *, 'call InitCorry first'
         stop
      endif
      if(ixhalflenl.lt.0.or.ixhalflenr.lt.0) then
         print *, 'ixhalflenl, ixhalflenr need to be > 0. STOP'
      endif
      open(21, file='corrxy'//trim(fn)//'.dat')
      write(21,'(a)')'variables=x,y,corr'
      write(21,'("Zone T=corrxy, I=", I8, ",J=",I8)') ixhalflenl+ixhalflenr+1,nypoint
      do i=-ixhalflenl,ixhalflenr
         call FFT1DB(spect1d(i,:), ctmp)
         call RearrangeFFT1D(nypoint,ctmp,corrout(i,:),ywave)
     enddo
!     corrxy = real(corrout)/maxval(real(corrout))
     corrxy = real(corrout)/abs(real(corrout(0,nypoint/2))+1.e-30)

     do j=1,nypoint
        do i=-ixhalflenl,ixhalflenr
           write(21,'(3E20.12)') real(i)*dxs, real(ywave(j))*dys, corrxy(i,j)
        enddo
     enddo
     close(21)

    end subroutine CalCorrxy

    subroutine CalCorrxz(ixhalflenl,ixhalflenr,dxs,spect1d,fn)
      implicit none
      integer, intent(in) :: ixhalflenl, ixhalflenr
      real(8), intent(in) :: dxs
      complex, intent(in) :: spect1d(-ixhalflenl:ixhalflenr,nypoint)
      character(*), intent(in) :: fn
      complex(8), dimension(nypoint) :: ctmp
      real(8), dimension(-ixhalflenl:ixhalflenr,nypoint) :: corrxy
      complex(8), dimension(-ixhalflenl:ixhalflenr,nypoint) :: corrout
      integer :: ywave(nypoint)
      integer ::  i, j

      if(.not.IsCorryInit) then
         print *, 'call InitCorry first'
         stop
      endif
      if(ixhalflenl.lt.0.or.ixhalflenr.lt.0) then
         print *, 'ixhalflenl, ixhalflenr need to be > 0. STOP'
      endif
      open(21, file='corrxz'//trim(fn)//'.dat')
      write(21,'(a)')'variables=x,z,corr'
      write(21,'("Zone T=corrxy, I=", I8, ",J=",I8)') ixhalflenl+ixhalflenr+1,nypoint
      do i=-ixhalflenl,ixhalflenr
         call FFT1DB(spect1d(i,:), ctmp)
         call RearrangeFFT1D(nypoint,ctmp,corrout(i,:),ywave)
     enddo
!     corrxy = real(corrout)/maxval(real(corrout))
     corrxy = real(corrout)/abs(real(corrout(0,nypoint/2))+1.e-30)

     do j=1,nypoint
        do i=-ixhalflenl,ixhalflenr
           write(21,'(3E20.12)') real(i)*dxs, real(ywave(j))*dys, corrxy(i,j)
        enddo
     enddo
     close(21)

    end subroutine CalCorrxz

end module modSpecty
