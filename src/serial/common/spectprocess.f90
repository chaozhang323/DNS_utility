!function to do spectral interpolation, grid should be uniformly distributed
!input data should have no extra point for periodicity
!and neither does the interpolated output.
!input:
!    varin: array contains the data to be interpolated
!    nold:  dimension of array varin
!    nnew:  dimension of array varout
!output:
!    varout: array to store the interpolated data
    subroutine FFTInterp(varin,nold,varout,nnew)
      use MFFT1D
      implicit none
      integer, intent(in) :: nold, nnew
      real, dimension(nold), intent(in) :: varin
      real, dimension(nnew), intent(out) :: varout
      complex, dimension(nnew) :: tmp
      integer :: ns, n
      if (nnew.lt.nold) then
        print*,'output array must have greater size than the input array!'
        stop
      else if (nnew.eq.nold) then
        varout = varin
        return
      end if
      call InitFFT1D(nold)
      call FFT1DF(cmplx(varin,0.),tmp(1:nold))
      ns = nold/2+1
      do n=nold, ns+1, -1
        tmp(nnew-nold+n)=tmp(n)
      end do
      do n=ns+1,nnew-nold+ns
        tmp(n) = 0.
      end do
      call InitFFT1D(nnew)
      call FFT1DB(tmp,tmp)
      varout = real(tmp)
    end subroutine
