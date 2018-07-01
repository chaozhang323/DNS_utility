!module used to compute buffer zone values
!Subroutines:
!  InitBuffer, initialize the module
!  DoBuffer,  compute values in the buffer zone using scalar input
!  DoBufferBlend,  compute values in the buffer zone using vector input
module MBuffer
  implicit none
  integer, private :: itype
  integer, private :: ibufferlen
  real(8), private :: bufferstrength
  real(8), private, dimension(:), allocatable :: coeff
  contains
  
      !subroutine to initialize the module and compute buffer coefficient
    !ibl (input): length of the buffer zone in terms of number of points
    !xloc(ibl) (input): x coordinate (in meters) of the buffer zone of dimension ibl
    !bcs (input): strength coefficient (in 1/meters) if tanh profile is used
    !it (input): type of the transient function (default: = 1)
    !            1: tanh profile, coeff(x) = 0.5*(1+tanh(bcs*(x-xmid)/(xend-xst))/tanh(bcs/2.0))
    !            2: half cosine profile, coef(x) = 0.5*(1-cos((x-xst)pi/(xend-xst)))
    subroutine InitBuffer(ibl,xloc,bcs,it)
      integer, intent(in) :: ibl
      real(8), intent(in) :: xloc(ibl)
      integer, intent(in), optional :: it
      real(8), intent(in) :: bcs
      real(8) :: xst, xend, xmid
      real(8) :: pi
      integer :: i
     
      ibufferlen = ibl
      if(size(xloc).ne.ibufferlen) then
         print *, 'Dimension of buffer index number NOT consistent with coordinate dimension! STOP'
         stop
      endif   
      
      xst = xloc(1)
      xend = xloc(ibufferlen)      
      if(xend.lt.xst) then
         print *, 'coordinate of buffer end should be larger than that of buffer beg! STOP!!!'
         stop
      endif   

!      xbufferlen = xend - xst ! buffer length in meters
      xmid = (xend + xst)/2.d0  ! buffer half point coord
      
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
          coeff(i) = 0.5*(1.d0+tanh(bufferstrength*((xloc(i)-xmid)/(xend-xst)))/tanh(bufferstrength/2.d0))
        end do
      else if (itype.eq.2) then
        pi = 4.*atan(1.)
        do i = 1, ibufferlen
          coeff(i) = 0.5*(1.d0-cos((xloc(i)-xst)*pi/(xend-xst)))
        end do
      else
        print*,'unknown buffer type itype = ', itype
        stop
      end if
      open(11,file='buffer.dat',status='unknown')
      write(11,'(a)')'variables=idx,xloc,coeff'
      write(11,*) 'Zone T=buffer, i =', ibufferlen
      do i = 1, ibufferlen
        write(11,*)i, xloc(i), coeff(i)
      end do
      close(11)
    end subroutine InitBuffer
    
    !subroutine to initialize the module and compute buffer coefficient
    !ibl (input): length of the buffer zone in terms of number of points
    !bcs (input): strength coefficient if tanh profile is used
    !it (input): type of the transient function
    !            1: tanh profile, coeff(i) = 0.5*(1+tanh(bcs*(i-ibl/2)))
    !            2: half cosine profile, coef(i) = 0.5*(1-cos((i-1)pi/(ibl-1)))
    subroutine InitBuffer_index(ibl,bcs,it)
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
      open(11,file='buffer_index.dat',status='unknown')
      write(11,'(a)')'varialbles=idx,coeff'
      do i = 1, ibufferlen
        write(11,*)i, coeff(i)
      end do
      close(11)
    end subroutine InitBuffer_index

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
