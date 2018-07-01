!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!     Search using Bisection method                                !!!
!!!     input:                                                       !!!
!!!        xarray:  array contains bunch of points in increasing     !!!
!!!                 order                                            !!!
!!!        xpt:     the point to be searched                         !!!
!!!     output:                                                      !!!
!!!        index:   the index of the 1st point that xpt lines        !!!
!!!                 between, ie x(index+1)>xpt>=x(index)             !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine BisectionSearch(xarray, ndim, xpt, index)
      implicit none
      integer, intent(in) :: ndim
      real(8), dimension(ndim), intent(in) :: xarray
      real(8), intent(in)     :: xpt
      integer, intent(out) :: index
      real(8), parameter :: err=1.e-4
      integer :: ib, ie, ihalf
      if (xpt-err.gt.xarray(ndim)) then
        index = ndim + 10
        return
      else if (xpt+err.lt.xarray(1)) then
        index = 1- 10
        return
      endif
      ihalf   = (1+ndim)/2
      if (xpt.ge.xarray(ihalf)) then
        ib = ihalf
        ie = ndim
      else
        ib = 1
        ie = ihalf
      end if
      do while (ie-ib.gt.1)
        ihalf = (ib+ie)/2
        if (xpt.ge.xarray(ihalf)) then
          ib = ihalf
        else
          ie = ihalf
        end if
      end do
      index = ib
      return
    end subroutine BisectionSearch

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!     Search using Bisection method                                !!!
!!!     input:                                                       !!!
!!!        xarray:  array contains bunch of points in increasing     !!!
!!!                 order                                            !!!
!!!        xinput:     the point to be searched
!!!        ndim:    dimension for input data
!!!        nnew:    dimension for ouput data
!!!     output:                                                      !!!
!!!        id:      the index of the 1st point that xpt lines        !!!
!!!                 between, ie x(idx+1)>xpt>=x(idx)
!!!        xnormal:      normalized xinput with
!!!                       xnormal in [0,1]            !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine BisectionSearch_updated(xarray, ndim, xinput, nnew, id, xnormal)
      implicit none
      integer, intent(in) :: ndim, nnew
      real(8), dimension(ndim), intent(in) :: xarray
      real(8), dimension(nnew), intent(in) :: xinput
      integer, dimension(nnew), intent(out) :: id
      real(8), dimension(nnew), intent(out) :: xnormal
      integer :: i
      real(8), parameter :: err=1.e-4
      integer :: ib, ie, ihalf
      integer :: index
      real(8) :: xpt

      do i = 1, nnew
         xpt = xinput(i)
         if (xpt-err.gt.xarray(ndim)) then
           index = ndim + 10
         else if (xpt+err.lt.xarray(1)) then
           index = 1- 10
         endif
         ihalf   = (1+ndim)/2
         if (xpt.ge.xarray(ihalf)) then
           ib = ihalf
           ie = ndim
         else
           ib = 1
           ie = ihalf
         end if
         do while (ie-ib.gt.1)
           ihalf = (ib+ie)/2
           if (xpt.ge.xarray(ihalf)) then
             ib = ihalf
           else
             ie = ihalf
           end if
         end do
         index = ib
         id(i) = index
         xnormal(i) = (xpt - xarray(index))/(xarray(index+1)-xarray(index))
      enddo
      return
    end subroutine BisectionSearch_updated




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!     Linear interploation using 2 points                          !!!
!!!     input:                                                       !!!
!!!        x1, x2:  coordinates of the 2 points                      !!!
!!!        Q1, Q2:  values at these two points                       !!!
!!!        x     :  coordinate of the  point to be interpolated      !!!
!!!     output:                                                      !!!
!!!        Qinterp: output interpolated value                        !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine LinearInterp2P(x1, x2, Q1, Q2, x, Qinterp)
      implicit none
      real(8), intent(in)  :: x1, x2, x
      real(8), intent(in)  :: Q1, Q2
      real(8), intent(out) :: Qinterp
      real(8), parameter :: ep=1.e-5
      if (x2-x1.le.ep) then
        Qinterp = 0.5*(Q1+Q2)
      else
        Qinterp = (Q2-Q1)*(x-x1)/(x2-x1) + Q1
      end if
    end subroutine LinearInterp2P

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!     Linear interploation with array inputs                       !!!
!!!     input:                                                       !!!
!!!        x, y: x points and function values y=f(x)                 !!!
!!!        ndim: dimension of x and y array                          !!!
!!!        xnew: x location whose y value is to be interpolated      !!!
!!!     output:                                                      !!!
!!!        ynew: the interpolated value at x=xnew                    !!!
!!!        iextrap:  what to do if xnew lies outside x array         !!!
!!!                 <0 no extrapolation, exit with error             !!!
!!!                 =0 zeroth order extrapolation                    !!!
!!!                 =1 first order extrapolation                     !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine LinearInterp(x,y,ndim,xnew,ynew,iextrap)
      implicit none
      integer, intent(in) :: ndim, iextrap
      real(8), dimension(ndim), intent(in) :: x, y
      real(8), intent(in) :: xnew
      real(8), intent(out) :: ynew
      integer :: idx
      call BisectionSearch(x,ndim,xnew,idx)
      if (idx.lt.1.or.idx.ge.ndim) then
        if (iextrap.lt.0) then
          print*,'point outside the orignal data range!'
          stop
        else if (iextrap.eq.0) then
          if (idx.lt.1) then
            ynew = y(1)
          else if (idx.ge.ndim) then
            ynew = y(ndim)
          end if
        else
          if (idx.lt.1) then
            call LinearInterp2P(x(1),x(2),y(1),y(2),xnew,ynew)
          else if (idx.ge.ndim) then
            call LinearInterp2P(x(ndim-1),x(ndim),y(ndim-1),y(ndim),xnew,ynew)
          end if
        end if
      end if
      call LinearInterp2P(x(idx),x(idx+1),y(idx),y(idx+1),xnew,ynew)
    end subroutine LinearInterp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Subroutine to compute spline interpolation coefficient            !!!
!! For each segment, map x to t so that t in [0,1]                   !!!
!! Thus we have for each segment y=a+bt+ct^2+dt^3, where t [0,1]     !!!
!! input:                                                            !!!
!!    x: array contains x locations                                  !!!
!!    y: array contains function value                               !!!
!!    ndim: size of the x and y arrays                               !!!
!!    y2beg: 2nd order derivative at the beginning of the spline     !!!
!!    y2end: 2nd order derivative at the end of the spline           !!!
!! output:                                                           !!!
!!    coef: output coefficients, 1st dim in the order of a, b, c, d  !!!
!!         dimension is 4x(ndim-1)                                   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine SplineCoefEx(x,y,ndim,coef,y2beg, y2end)
      implicit none
      integer, intent(in) :: ndim
      real(8), dimension(ndim), intent(in) :: x,y
      real(8), dimension(4,ndim-1), intent(out) :: coef
      real(8), intent(in) :: y2beg, y2end
      real(8), dimension(ndim) :: aa, bb, cc, rr, tmp
      integer :: i
      !tri-diagonal system
      aa = 1.
      bb = 4.; bb(1) = 2.; bb(ndim)=2.
      cc = 1.
      do i=2,ndim-1
        rr(i) = 3.*(y(i+1)-y(i-1))
      end do
      rr(1) = 3.*(y(2)-y(1))-0.5*y2beg*(x(2)-x(1))**2
      rr(ndim)=3.*(y(ndim)-y(ndim-1))+0.5*y2end*(x(ndim)-x(ndim-1))**2
      call tridiagonal_solve(ndim,aa,bb,cc,rr,tmp)
      !test
      !do i=2,ndim-1
      !  coef(1,i) = rr(i)-(aa(i)*tmp(i-1)+bb(i)*tmp(i)+cc(i)*tmp(i+1))
      !end do
      !print*,maxval(abs(coef(1,2:ndim-1)))

      do i=1,ndim-1
        coef(1,i) = y(i)
        coef(2,i) = tmp(i)
        coef(3,i) = 3.*(y(i+1)-y(i))-2.*tmp(i)-tmp(i+1)
        coef(4,i) = 2.*(y(i)-y(i+1))+tmp(i)+tmp(i+1)
      end do
      return
    end subroutine SplineCoefEx

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Subroutine to compute NATURAL spline interpolation coefficient    !!!
!! Second order derivate at the beginning and end are set to zero    !!!
!! For each segment, map x to t so that t in [0,1]                   !!!
!! Thus we have for each segment y=a+bt+ct^2+dt^3, where t [0,1]     !!!
!! input:                                                            !!!
!!    y: array contains function value                               !!!
!!    ndim: size of the x and y arrays                               !!!
!! output:                                                           !!!
!!    coef: output coefficients, 1st dim in the order of a, b, c, d  !!!
!!         dimension is 4x(ndim-1)                                   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine SplineCoef_natural(y,ndim,coef)
      implicit none
      integer, intent(in) :: ndim
      real(8), dimension(ndim), intent(in) :: y
      real(8), dimension(4,ndim-1), intent(out) :: coef
      real(8), dimension(ndim) :: aa, bb, cc, rr, tmp
      integer :: i
      !tri-diagonal system
      aa = 1.
      bb = 4.; bb(1) = 2.; bb(ndim)=2.
      cc = 1.
      do i=2,ndim-1
        rr(i) = 3.*(y(i+1)-y(i-1))
      end do
      rr(1) = 3.*(y(2)-y(1)); rr(ndim)=3.*(y(ndim)-y(ndim-1))
      call tridiagonal_solve(ndim,aa,bb,cc,rr,tmp)
      do i=1,ndim-1
        coef(1,i) = y(i)
        coef(2,i) = tmp(i)
        coef(3,i) = 3.*(y(i+1)-y(i))-2.*tmp(i)-tmp(i+1)
        coef(4,i) = 2.*(y(i)-y(i+1))+tmp(i)+tmp(i+1)
      end do
      return
    end subroutine SplineCoef_natural

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Subroutine to compute PERIODIC spline interpolation coefficient   !!!
!! Second order derivate at the beginning and end are set to zero    !!!
!! For each segment, map x to t so that t in [0,1]                   !!!
!! Thus we have for each segment y=a+bt+ct^2+dt^3, where t [0,1]     !!!
!! input:                                                            !!!
!!    y: array contains function value                               !!!
!!    ndim: size of the x and y arrays                               !!!
!! output:                                                           !!!
!!    coef: output coefficients, 1st dim in the order of a, b, c, d  !!!
!!         dimension is 4x(ndim-1)                                   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine SplineCoef_periodic(y,ndim,coef)
      implicit none
      integer, intent(in) :: ndim
      real(8), dimension(ndim), intent(in) :: y
      real(8), dimension(4,ndim-1), intent(out) :: coef
      real(8), dimension(ndim) :: aa, bb, cc, rr, tmp
      integer :: i
      !tri-diagonal system
      aa = 1.
      bb = 4.
      cc = 1.
      do i=2,ndim-1
        rr(i) = 3.*(y(i+1)-y(i-1))
      end do
      rr(1) = 3.*(y(2)-y(ndim)); rr(ndim)=3.*(y(1)-y(ndim-1))
      call cyclictridiagonal_solve(ndim,aa,bb,cc,rr,tmp)
      do i=1,ndim-1
        coef(1,i) = y(i)
        coef(2,i) = tmp(i)
        coef(3,i) = 3.*(y(i+1)-y(i))-2.*tmp(i)-tmp(i+1)
        coef(4,i) = 2.*(y(i)-y(i+1))+tmp(i)+tmp(i+1)
      end do
      return
    end subroutine SplineCoef_periodic

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! pointwise spline interpolation after coefficients are computed    !!!
!! input:                                                            !!!
!!   x: array of the input x locations                               !!!
!!   coef: spline coefficient matrix computed by SplineCoef(EX)      !!!
!!         dimension is 4xndim-1                                     !!!
!!   ndim: number of dimension for the input data                    !!!
!!   xnew: x location where y value is to be interoploated           !!!
!! output:                                                           !!!
!!   ynew: interpolated value at x=xnew                              !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine SplineInterpPT(x,coef,ndim,xnew,ynew)
      implicit none
      integer, intent(in) :: ndim
      real(8), dimension(ndim), intent(in) :: x
      real(8), dimension(4,ndim-1), intent(in) :: coef
      real(8), intent(in) :: xnew
      real(8), intent(out) :: ynew
      real(8) :: xt, xt2, xt3
      integer :: idx
      real(8),parameter :: eps = 1.d-30
      call BisectionSearch(x,ndim,xnew,idx)
      if (idx.lt.1.or.idx.gt.ndim) then
!        print*,'warning: point outside the original region in SplineInterpPT'
        if (idx.lt.1) then
          xt = (xnew-x(1))/(x(2)-x(1)+eps)
          idx = 1
        else
          xt = (xnew-x(ndim-1))/(x(ndim)-x(ndim-1)+eps)
          idx = ndim-1
        end if
      else
!        print *, 'idx=',idx, 'x(idx+1)=',x(idx+1),'x(idx)=',x(idx)
        xt = (xnew-x(idx))/(x(idx+1)-x(idx)+eps)
      end if
      xt2 = xt*xt; xt3 = xt2*xt
      ynew = coef(1,idx)+coef(2,idx)*xt+coef(3,idx)*xt2+coef(4,idx)*xt3
      return
    end subroutine SplineInterpPT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! arraywise spline interpolation                                    !!!
!! input:                                                            !!!
!!   x: array of the input x locations                               !!!
!!   y: array of the function values at x                            !!!
!!   ndim1: number of dimension for the input data (x and y)         !!!
!!   xnew: array of x locations to be interoploated                  !!!
!!   ndim2: number of dimension for the interpolated data (xnew ynew)!!!
!!   iextrap: what to do for extraploation                           !!!
!!            <0 no extrapolation, exit with error                   !!!
!!            =0 zeroth order extrapolation                          !!!
!!            =1 first  order extrapolation                          !!!
!!            =2 extrapolation using the spline                      !!!
!! output:                                                           !!!
!!   ynew: interpolated values at x=xnew                             !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine SplineInterp(x,y,ndim1,xnew,ynew,ndim2,iextrap)
      implicit none
      integer, intent(in) :: ndim1, ndim2, iextrap
      real(8), dimension(ndim1), intent(in) :: x, y
      real(8), dimension(ndim2), intent(in) :: xnew
      real(8), dimension(ndim2), intent(out) :: ynew
      real(8), dimension(4,ndim1-1) :: coef
      integer :: i
      call SplineCoef_natural(y,ndim1,coef)
      do i=1,ndim2
        if (xnew(i).lt.x(1)) then
          if (iextrap.lt.0) then
            print*,'point outside original region in SplineInterp'
            stop
          else if (iextrap.eq.0) then
            ynew(i) = y(1)
          else if (iextrap.eq.1) then
            call LinearInterp2P(x(1),x(2),y(1),y(2),xnew(i),ynew(i))
          else
             call SplineInterpPT(x,coef,ndim1,xnew(i),ynew(i))
          end if
        else if (xnew(i).gt.x(ndim1)) then
          if (iextrap.lt.0) then
            print*,'point outside original region in SplineInterp'
            stop
          else if (iextrap.eq.0) then
            ynew(i) = y(ndim1)
          else if (iextrap.eq.1) then
            call LinearInterp2P(x(ndim1-1),x(ndim1),y(ndim1-1),y(ndim1),xnew(i),ynew(i))
          else
            call SplineInterpPT(x,coef,ndim1,xnew(i),ynew(i))
          end if
        else
          call SplineInterpPT(x,coef,ndim1,xnew(i),ynew(i))
        end if
      end do
    end subroutine SplineInterp
     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! compute derivatives of spline based on coefficien                 !!!
!! input:                                                            !!!
!!   x: array of the input x locations                               !!!
!!   coef: spline coefficients computed by SplineCoeff(Ex)           !!!
!!         dimension is 4xndim-1                                     !!!
!!   ndim: number of dimension for the input data x                  !!!
!!   xnew: array of x locations to be interoploated                  !!!
!!   nderiv: order of derivative to be computed (>0 and <4)          !!!
!! output:                                                           !!!
!!   dynew: derivative computed at xnew using spline                 !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine SplineDerivative(x,coef,ndim,xnew,dynew,nderiv)
      implicit none
      integer, intent(in) :: ndim, nderiv
      real(8), dimension(ndim), intent(in) :: x
      real(8), dimension(4,ndim-1), intent(in) :: coef
      real(8), intent(in) :: xnew
      real(8), intent(out) :: dynew
      real :: xt, xt2, xl
      integer :: idx
      if (nderiv.lt.0) then
        print*,'Error in SplineDerivative: derivative order less than zero!'
        stop
      else if (nderiv.eq.0) then
        call SplineInterpPT(x, coef, ndim, xnew, dynew)
        return
      else if (nderiv.ge.4) then
        dynew = 0.
        return
      end if
      call BisectionSearch(x,ndim,xnew,idx)
      if (idx.lt.1.or.idx.gt.ndim) then
!        print*,'warning: point outside the original region in SplineDerivative'
        if (idx.lt.1) then
          xl = (x(2)-x(1))
          xt = (xnew-x(1))/xl
          idx = 1
        else
          xl = (x(ndim)-x(ndim-1))
          xt = (xnew-x(ndim-1))/xl
          idx = ndim-1
        end if
      else
        xl = (x(idx+1)-x(idx))
        xt = (xnew-x(idx))/xl
      end if
      if (nderiv.eq.1) then
        xt2 = xt*xt
        dynew = (coef(2,idx)+2.*coef(3,idx)*xt+3.*coef(4,idx)*xt2)/xl
      else if (nderiv.eq.2) then
        dynew = (2.*coef(3,idx)+6.*coef(4,idx)*xt)/(xl*xl)
      else if (nderiv.eq.3) then
        dynew = 6.*coef(4,idx)/(xl*xl*xl)
      end if
      return
    end subroutine SplineDerivative

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Compute integration of spline at each node point                  !!!
!! The integration is zero at the first input x point                !!!
!! input:                                                            !!!
!!   x: array of the input x locations                               !!!
!!   coef: spline coefficients computed by SplineCoeff(Ex)           !!!
!!         dimension is 4xndim-1                                     !!!
!!   ndim: number of dimension for the input data x                  !!!
!! output:                                                           !!!
!!   nodeval: integral values at each x locations                    !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine SplineIntNode(x,coef,ndim,nodeval)
      implicit none
      integer, intent(in) :: ndim
      real(8), dimension(ndim), intent(in) :: x
      real(8), dimension(4,ndim-1), intent(in) :: coef
      real(8), intent(out), dimension(ndim) :: nodeval
      real :: xl
      integer :: idx
      nodeval(1) = 0.
      do idx=1,ndim-1
        xl = (x(idx+1)-x(idx))
        nodeval(idx+1) = nodeval(idx)+xl*(coef(1,idx)+0.5*coef(2,idx)+coef(3,idx)/3.+0.25*coef(4,idx))
      end do
      return
    end subroutine SplineIntNode

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Compute integration of spline at a given x point                  !!!
!! The integration is zero at the first input x point                !!!
!! input:                                                            !!!
!!   x: array of the input x locations                               !!!
!!   coef: spline coefficients computed by SplineCoeff(Ex)           !!!
!!         dimension is 4xndim-1                                     !!!
!!   ndim: number of dimension for the input data x                  !!!
!!   nodeval: integral values at each x locations                    !!!
!!            computed by subroutine SplineIntNode                   !!!
!!   xnew: the end point of the integration                          !!!
!! output:                                                           !!!
!!   yint: the computed integration from x(1) to xnew                !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine SplineIntegral(x,coef,nodeval,ndim,xnew,yint)
      implicit none
      integer, intent(in) :: ndim
      real(8), dimension(ndim), intent(in) :: x, nodeval
      real(8), dimension(4,ndim-1), intent(in) :: coef
      real(8), intent(in) :: xnew
      real(8), intent(out) :: yint
      real :: xt, xt2, xt3, xt4, xl
      integer :: idx
      call BisectionSearch(x,ndim,xnew,idx)
      if (idx.lt.1.or.idx.gt.ndim) then
!        print*,'warning: point outside the original region in SplineIntegral'
        if (idx.lt.1) then
          xl = (x(2)-x(1))
          xt = (xnew-x(1))/xl
          idx = 1
        else
          xl = (x(ndim)-x(ndim-1))
          xt = (xnew-x(ndim-1))/xl
          idx = ndim-1
        end if
      else
        xl = (x(idx+1)-x(idx))
        xt = (xnew-x(idx))/xl
      end if
      xt2 = xt*xt; xt3 = xt2*xt; xt4 = xt3*xt
      yint = xl*(coef(1,idx)*xt+0.5*coef(2,idx)*xt2+coef(3,idx)*xt3/3.+0.25*coef(4,idx)*xt4)+nodeval(idx)
      return
    end subroutine SplineIntegral

    
    !Jacobi iterative method for solving linear system Ax=b
    !initially used for spline interpolation
    !not used anymore after tri-diagonal solver was implemented
    subroutine JacobiIteration(A, b, x, ndim, err)
      implicit none
      integer, intent(in) :: ndim
      real(8), intent(in) :: A(ndim,ndim), b(ndim)
      real(8), intent(inout) :: x(ndim)
      real(8), intent(in) :: err
      real :: LU(ndim,ndim), x1(ndim), errmax, error
      integer i, j
    
      if (err.le.0) then
        errmax=1.e-5 !default error torrelance
      else
        errmax=err
      end if
      LU = A
      error = 1.e10
      do while (error.ge.errmax)
        do i=1,ndim
          LU(i,i)=0.
        end do
        x1 = matmul(LU,x)
        do i=1,ndim
          x1(i) = (b(i)-x1(i))/A(i,i)
        end do
        error = 0.
        do i=1,ndim
          error=error+(x1(i)-x(i))**2
        end do
        error=sqrt(error)
        x  = x1
      end do
      return
    end subroutine JacobiIteration
