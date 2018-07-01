!subroutine to generate uniform random numbers
!time is used to seed the random_number intrinsic call
!input:
!   ndim: number of random numbers to output
!output:
!   x: array of dimension ndim to store the numbers
  subroutine Rand_uniform(x,ndim)
    integer, intent(in) :: ndim
    real, dimension(ndim) :: x
    integer, dimension(:), allocatable :: seed
    integer :: n, clock
    call random_seed(size = n)
    allocate(seed(n))
    call system_clock(count=clock)
    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    call random_seed(put = seed)
    call random_number(x)
    deallocate(seed)
  end subroutine Rand_uniform

!subroutine to generate gaussian random numbers
!input:
!   ndim: number of random numbers to output
!output:
!   x: array of dimension ndim to store the numbers
  subroutine GaussRandom(x,ndim)
    integer, intent(in) :: ndim
    real, dimension(ndim), intent(out) :: x
    integer :: nbuffer
    real, dimension(:), allocatable :: buf
    real :: u, v, s, tt
    integer :: i, nfilled
          
    nbuffer = min(ndim*2,2000)
    allocate(buf(nbuffer))
    nfilled=0
    do while(nfilled.lt.ndim)
      call rand_uniform(buf,nbuffer)
      buf = 2.*buf-1  !-1 to 1
      i=1
      do while (i.lt.nbuffer)
        u=buf(i)
        v=buf(i+1)
        s=u*u+v*v
        if (s.lt.1.and.s.gt.0) then
          tt = sqrt(-2.*log(s)/s)
          if (nfilled+1.le.ndim) x(nfilled+1) = u*tt
          if (nfilled+2.le.ndim) x(nfilled+2)=v*tt
          nfilled = nfilled+2
          if (nfilled.ge.ndim)  exit
        end if
        i=i+2
      end do
    end do
    deallocate(buf)
  end subroutine GaussRandom
