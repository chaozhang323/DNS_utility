!contains wrapper of fftpack5 and fft window module
!for details about the FFT, please refer to the documents in fftpack5
!website: http://www.cisl.ucar.edu/css/software/fftpack5/

!1D FFT wrapper
module MFFT1D
  implicit none
  integer, private :: npoint
  integer, private :: lensave, lenwork, ierr
  real, dimension(:), private, allocatable :: wsave, wwork
  contains
    !initialize FFT given number of points (np)
    subroutine InitFFT1D(np)
      integer, intent(in) :: np
      npoint = np
      lensave = 2*npoint+int(log(real(npoint))/log(2.))+4
      lenwork = 2*npoint

      if (allocated(wsave)) deallocate(wsave)
      if (allocated(wwork)) deallocate(wwork)
      allocate(wsave(lensave))
      allocate(wwork(lenwork))
      ! initalize fft
      call cfft1i(npoint,wsave,lensave,ierr)
      if (ierr.ne.0) then
        print*, 'FFT1D initialization error!'
        stop
      end if
    end subroutine InitFFT1D

    !forward FFT
    subroutine FFT1DF(varin, varout)
      complex, dimension(npoint), intent(in) :: varin
      complex, dimension(npoint), intent(out) :: varout
      varout(:) = varin(:)
      call cfft1f(npoint,1,varout,npoint,wsave,lensave,wwork,lenwork,ierr)
      if (ierr.ne.0) then
        print*, 'Error occured in FFT1DF!'
        stop
      end if
    end subroutine FFT1DF

    !backward FFT
    subroutine FFT1DB(varin, varout)
      complex, dimension(npoint), intent(in) :: varin
      complex, dimension(npoint), intent(out) :: varout
      varout(:) = varin(:)
      call cfft1b(npoint,1,varout,npoint,wsave,lensave,wwork,lenwork,ierr)
      if (ierr.ne.0) then
        print*, 'Error occured in FFT1DB!'
        stop
      end if
    end subroutine FFT1DB

end module MFFT1D

!2D FFT wrapper
module MFFT2D
  implicit none
  integer, private :: npoint1, npoint2
  integer, private :: lensave, lenwork, ierr
  real, dimension(:), private, allocatable :: wsave, wwork
  contains
    !initialize 2D FFT given the two dimesions (np1 and np2)
    subroutine InitFFT2D(np1, np2)
      integer, intent(in) :: np1, np2
      npoint1 = np1
      npoint2 = np2
!      lensave = 2*(npoint1+npoint2)+int(log(real(npoint1)))&
!              + int(log(real(npoint2)))+8
      lensave = 2*(npoint1+npoint2)+int(log(real(npoint1))/log(2.))&
              + int(log(real(npoint2))/log(2.))+8
      lenwork = 2*npoint1*npoint2

      if (allocated(wsave)) deallocate(wsave)
      if (allocated(wwork)) deallocate(wwork)
      allocate(wsave(lensave))
      allocate(wwork(lenwork))
      ! initalize fft
      call cfft2i(npoint1,npoint2,wsave,lensave,ierr)
      if (ierr.ne.0) then
        print*, 'FFT2D initialization error!'
        stop
      end if
    end subroutine InitFFT2D

    !forward FFT
    subroutine FFT2DF(varin, varout)
      complex, dimension(npoint1,npoint2), intent(in) :: varin
      complex, dimension(npoint1,npoint2), intent(out) :: varout
      varout = varin
      call cfft2f(npoint1,npoint1,npoint2,varout,wsave,lensave,wwork,lenwork,ierr)
      if (ierr.ne.0) then
        print*, 'Error occured in FFT2DF!'
        stop
      end if
    end subroutine FFT2DF

    !backward FFT
    subroutine FFT2DB(varin, varout)
      complex, dimension(npoint1,npoint2), intent(in) :: varin
      complex, dimension(npoint1,npoint2), intent(out) :: varout
      varout = varin
      call cfft2b(npoint1,npoint1,npoint2,varout,wsave,lensave,wwork,lenwork,ierr)
      if (ierr.ne.0) then
        print*, 'Error occured in FFT2DB!'
        stop
      end if
    end subroutine FFT2DB
end module MFFT2D

!3D FFT wrapper
module MFFT3D
  implicit none
  integer, private :: npoint1, npoint2, npoint3
  integer, private :: lensave, lenwork, ierr
  integer, private :: lensave1d, lenwork1d
  real, dimension(:), private, allocatable :: wsave, wwork
  real, dimension(:), private, allocatable :: wsave1d, wwork1d
  contains
    !initialize 3D FFT given the 3 dimensions (np1, np2, np3)
    subroutine InitFFT3D(np1, np2, np3)
      integer, intent(in) :: np1, np2, np3
      npoint1 = np1
      npoint2 = np2
      npoint3 = np3
      ! initalize fft
!      lensave = 2*(npoint1+npoint2)+int(log(real(npoint1)))&
!              + int(log(real(npoint2)))+8
      lensave = 2*(npoint1+npoint2)+int(log(real(npoint1))/log(2.))&
              + int(log(real(npoint2))/log(2.))+8
      lenwork = 2*npoint1*npoint2
      !lensave1d = 2*npoint3+int(log(real(npoint3)))+4
      lensave1d = 2*npoint3+int(log(real(npoint3))/log(2.))+4
      lenwork1d = 2*npoint3

      if (allocated(wsave)) deallocate(wsave)
      if (allocated(wwork)) deallocate(wwork)
      if (allocated(wsave1d)) deallocate(wsave1d)
      if (allocated(wwork1d)) deallocate(wwork1d)
      allocate(wsave(lensave))
      allocate(wwork(lenwork))
      allocate(wsave1d(lensave1d))
      allocate(wwork1d(lenwork1d))
      call cfft2i(npoint1,npoint2,wsave,lensave,ierr)
      if (ierr.ne.0) then
        print*, 'FFT2D initialization in FFT3D error!'
        stop
      end if

      call cfft1i(npoint3,wsave1d,lensave1d,ierr)

      if (ierr.ne.0) then
        print*, 'FFT1D initialization in FFT3D error!'
        stop
      end if
    end subroutine InitFFT3D

    !forward FFT
    subroutine FFT3DF(varin, varout)
      complex, dimension(npoint1,npoint2,npoint3), intent(in) :: varin
      complex, dimension(npoint1,npoint2,npoint3), intent(out) :: varout
      integer :: i,j,k, lenc
      varout = varin
      !fft in first 2 dimension
      do k=1,npoint3
        call cfft2f(npoint1,npoint1,npoint2,varout(1,1,k)&
             ,wsave,lensave,wwork,lenwork,ierr)
        if (ierr.ne.0) then
          print*, 'Error occured in FFT2DF of FFT3D!'
          stop
        end if
      end do
      !fft in 3rd dimension
      lenc = (npoint1*npoint2)*(npoint3-1)+1
      do j=1,npoint2
        do i=1,npoint1
          call cfft1f(npoint3,npoint1*npoint2,varout(i,j,1),lenc&
               ,wsave1d,lensave1d,wwork1d, lenwork1d,ierr)
          if (ierr.ne.0) then
            print*, 'Error occured in FFT1DF of FFT3D!'
            stop
          end if
        end do
      end do
    end subroutine FFT3DF

    !backward FFT
    subroutine FFT3DB(varin, varout)
      complex, dimension(npoint1,npoint2,npoint3), intent(in) :: varin
      complex, dimension(npoint1,npoint2,npoint3), intent(out) :: varout
      integer :: i,j,k, lenc
      varout = varin
      !fft in first 2 dimension
      do k=1,npoint3
        call cfft2b(npoint1,npoint1,npoint2,varout(1,1,k)&
             ,wsave,lensave,wwork,lenwork,ierr)
!print *, 'npoint1 = ', npoint1, 'k = ', k

        if (ierr.ne.0) then
          print*, 'Error occured in FFT2DF of FFT3D!'
          stop
        end if
      end do
      !fft in 3rd dimension
      lenc = (npoint1*npoint2)*(npoint3-1)+1
      do j=1,npoint2
        do i=1,npoint1
          call cfft1b(npoint3,npoint1*npoint2,varout(i,j,1),lenc&
               ,wsave1d,lensave1d,wwork1d, lenwork1d,ierr)
          if (ierr.ne.0) then
            print*, 'Error occured in FFT1DF of FFT3D!'
            stop
          end if
        end do
      end do
    end subroutine FFT3DB
end module MFFT3D

!module for computing window functions
module MFFTWindow
  integer, private :: npoint, ntype
  real, dimension(:), allocatable :: wcoeff !window function values
  real :: energy_ratio  !=1/window energy
  contains
    !initialize window function computing
    !np (input): number of points used
    !nt (input): type of the window to be computed
    !            0: flat top window
    !            1: Hanning  window
    !            2: modified Hanning  window
    subroutine InitFFTWindow(np, nt)
      implicit none
      integer, intent(in) :: np, nt
      npoint = np
      ntype = nt
      if (allocated(wcoeff)) deallocate(wcoeff)
      allocate(wcoeff(npoint))
      call CalWindow()
      if (.true.) call WriteWindow('window.dat')
    end subroutine InitFFTWindow

    subroutine InitFFTWindow_1(np, nt, wcoeff_output)
      implicit none
      integer, intent(in) :: np, nt
      real(8), intent(out) :: wcoeff_output(np)
      npoint = np
      ntype = nt
      if (allocated(wcoeff)) deallocate(wcoeff)
      allocate(wcoeff(npoint))
      call CalWindow()
      wcoeff_output = wcoeff
      if (.true.) call WriteWindow('window.dat')
    end subroutine InitFFTWindow_1

    subroutine InitFFTWindow_noprint(np, nt)
      implicit none
      integer, intent(in) :: np, nt
      npoint = np
      ntype = nt
      if (allocated(wcoeff)) deallocate(wcoeff)
      allocate(wcoeff(npoint))
      call CalWindow_noprint()
      if (.true.) call WriteWindow('window.dat')
    end subroutine InitFFTWindow_noprint
    
    !subroutine to compute the window functions
    subroutine CalWindow()
      integer :: i
      real :: pi
      pi  = 4.*atan(1.)
      select case(ntype)
      case(0)  !flat top window
        print*,'calculating flat top window'
        wcoeff = 1.
      case(1)  !Hanning
        print*,'calculating hanning window'
        do i=1, npoint
          wcoeff(i) = 0.5*(1.-cos(2.*pi*real(i-1)/real(npoint)))
        end do
      case(2)  !modified Hanning with flat top
        print*,'calculating modified hanning window'
        wcoeff = 1.
        do i=1, npoint/8-1
          wcoeff(i) = 0.5*(1.-cos(8.*pi*real(i-1)/real(npoint)))
        end do
        do i=7*npoint/8+1, npoint
          wcoeff(i) = 0.5*(1.-cos(8.*pi*real(i-1)/real(npoint)))
        end do
      case(3) ! Hamming
        print*,'calculating hamming window'
        do i=1, npoint
          wcoeff(i) = 0.54 - 0.46*cos(2.*pi*real(i-1)/real(npoint))
        end do
      case default
        print*,'unknow window type with ntype = ', ntype
        stop
      end select
      energy_ratio = 1./(sum(wcoeff*wcoeff)/npoint)
      wcoeff = wcoeff*sqrt(energy_ratio)
    end subroutine CalWindow

    !subroutine to compute the window functions
    subroutine CalWindow_noprint()
      integer :: i
      real :: pi
      pi  = 4.*atan(1.)
      select case(ntype)
      case(0)  !flat top window
   !     print*,'calculating flat top window'
        wcoeff = 1.
      case(1)  !Hanning
   !     print*,'calculating hanning window'
        do i=1, npoint
          wcoeff(i) = 0.5*(1.-cos(2.*pi*real(i-1)/real(npoint)))
        end do
      case(2)  !modified Hanning with flat top
   !     print*,'calculating modified hanning window'
        wcoeff = 1.
        do i=1, npoint/8-1
          wcoeff(i) = 0.5*(1.-cos(8.*pi*real(i-1)/real(npoint)))
        end do
        do i=7*npoint/8+1, npoint
          wcoeff(i) = 0.5*(1.-cos(8.*pi*real(i-1)/real(npoint)))
        end do
      case(3) ! Hamming
   !     print*,'calculating hamming window'
        do i=1, npoint
          wcoeff(i) = 0.54 - 0.46*cos(2.*pi*real(i-1)/real(npoint))
        end do
      case default
        print*,'unknow window type with ntype = ', ntype
        stop
      end select
      energy_ratio = 1./(sum(wcoeff*wcoeff)/npoint)
      wcoeff = wcoeff*sqrt(energy_ratio)
    end subroutine CalWindow_noprint

    !output window function to a file
    subroutine WriteWindow(fn)
      character(*),intent(in) :: fn
      integer :: i
      open(11,file=fn)
      write(11,'(a)')'variables=i,w'
      do i=1,npoint
        write(11,'(I8, E20.11)')i,wcoeff(i)
      end do
      close(11)
    end subroutine WriteWindow
end module MFFTWindow
