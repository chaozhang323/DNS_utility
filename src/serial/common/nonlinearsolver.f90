!Newton Rapson solver
!input fuction f(x) and its derivative f'(x)
!with initial guess x0
!output the xout value where f(xout)=0
  subroutine NewtonRapson(xout,func,funcd,x0)
    real, intent(in) :: x0
    real, intent(out) :: xout
    real, parameter :: eps = 1.e-10
    real :: fv, dx
    interface
      real function func(x)
        implicit none
        real, intent(in) :: x
      end function func
    end interface
    interface
      real function funcd(x)
        implicit none
        real, intent(in) :: x
      end function funcd
    end interface
    xout = x0
    fv = func(xout)
    do while(abs(fv).gt.eps)
      dx = -fv/funcd(xout)
      xout = xout+dx
      fv = func(xout)
    end do
  end subroutine NewtonRapson
