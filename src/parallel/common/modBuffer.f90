!module used to compute buffer zone values
!Subroutines:
!  InitBuffer, initialize the module
!  DoBuffer,  compute values in the buffer zone using scalar input
!  DoBufferBlend,  compute values in the buffer zone using vector input
module MBuffer
  implicit none
  integer, private :: ibufferlen, itype
  real(8), private :: bufferstrength
  real(8), private, dimension(:), allocatable :: coeff
  contains
    !subroutine to initialize the module and compute buffer coefficient
    !ibl (input): length of the buffer zone in terms of number of points
    !bcs (input): strength coefficient if tanh profile is used
    !it (input): type of the transient function
    !            1: tanh profile, coeff(i) = 0.5*(1+tanh(bcs*(i-ibl/2)))
    !            2: half cosine profile, coef(i) = 0.5*(1-cos((i-1)pi/(ibl-1)))
    subroutine InitBuffer(ibl,bcs,it)
      integer, intent(in) :: ibl
      integer, intent(in), optional :: it
      real(8), intent(in) :: bcs
      real(8) :: pi
      integer :: i
      ibufferlen = ibl
      bufferstrength = bcs
      if (present(it)) then
        itype = it
      else
        itype = 1
      end if
      if (allocated(coeff)) deallocate(coeff)
      allocate(coeff(ibufferlen))
      if (itype.eq.1) then
        do i = 1, ibufferlen
          coeff(i) = 0.5*(1.d0+tanh(bufferstrength*((2.*i-float(ibufferlen))/2.)))
        end do
      else if (itype.eq.2) then
        pi = 4.*atan(1.)
        do i = 1, ibufferlen
          coeff(i) = 0.5*(1.d0-cos(float(i-1)*pi/float(ibufferlen-1)))
        end do
      else
        print*,'unknown buffer type itype = ', itype
        stop
      end if
      open(11,file='buffer.dat',status='unknown')
      write(11,'(a)')'varialbles=idx,coeff'
      do i = 1, ibufferlen
        write(11,*)i, coeff(i)
      end do
      close(11)
    end subroutine InitBuffer

    !compute buffer zone values using scalar input
    !uin1, uin2 (input): scalar values at the beginning and end of the buffer zone
    !uout (output): values in between
    !               uout(i) = uin1*(1-coeff(i))+uin2*coeff(i)
    subroutine DoBuffer(uin1, uin2, uout)
      real(8), intent(in) :: uin1, uin2
      real(8), intent(out), dimension(ibufferlen) :: uout
      integer :: i
      do i = 1, ibufferlen
        uout(i) = uin1*(1.d0-coeff(i))+uin2*coeff(i)
      end do
    end subroutine DoBuffer

    !compute buffer zone values using vector input
    !uin1, uin2 (input): vector values at the beginning and end of the buffer zone
    !uout (output): values in between
    !               uout(i) = uin1(i)*(1-coeff(i))+uin2(i)*coeff(i)
    subroutine DoBufferBlend(uin1, uin2, uout)
      real(8), intent(in), dimension(ibufferlen) :: uin1
      real(8), intent(in) :: uin2
      real(8), intent(out), dimension(ibufferlen) :: uout
      integer :: i
      do i = 1, ibufferlen
        uout(i) = uin1(i)*(1.d0-coeff(i))+uin2*coeff(i)
      end do
    end subroutine DoBufferBlend
end module MBuffer
