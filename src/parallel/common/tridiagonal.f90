   !subroutine to solve tridiagnoal system Ax=r
   !ref: numerical recipes in C chapter 2.4
   !array a, b and c store values along the diagonal
   !a(1) and c(nmax) are not used
   !output:
   !   xout: the solution vector
   subroutine tridiagonal_solve(nmax,a,b,c,r,xout)
     implicit none
     integer, intent(in) :: nmax
     real(8), dimension(nmax), intent(in) :: a, b, c
     real(8), dimension(nmax), intent(in) :: r
     real(8), dimension(nmax), intent(out) :: xout
     real(8) :: gam(nmax), bet
     integer :: n
     if (b(1).eq.0) then
       print*,'error in solving tridiagonal system, b(1) is zero'
       stop
     end if
     bet = b(1)
     xout(1) = r(1)/bet
     do n=2,nmax
       gam(n)=c(n-1)/bet
       bet=b(n)-a(n)*gam(n)
       if (bet.eq.0) then
         print*,'error in solving tridiagonal system, b(1) is zero'
         stop
       end if
       xout(n)=(r(n)-a(n)*xout(n-1))/bet
     end do
     !back substitution
     do n=nmax-1,1,-1
       xout(n) = xout(n) - gam(n+1)*xout(n+1)
     end do
   end subroutine tridiagonal_solve

   !subroutine to solve cyclic tridiagnoal system
   !ref: numerical recipes in C chapter 2.7
   !alpha is stored at a(1) and gamma is stored at c(nmax)
   !output:
   !   xout: the solution vector
   !  |b1     c1     0    ...            beta|
   !  |a2     b2     c2   ...                |
   !  |                   ...                |
   !  |                        aN-1 bN-1 cN-1|
   !  |alpha                    0   aN   bN  |
   subroutine cyclictridiagonal_solve(nmax,a,b,c,r,xout)
     implicit none
     integer, intent(in) :: nmax
     real(8), dimension(nmax), intent(in) :: a,b,c
     real(8), dimension(nmax), intent(in) :: r
     real(8), dimension(nmax), intent(out) :: xout
     real(8) :: gamma, factor, bb(nmax), uu(nmax), zz(nmax)
     if (nmax.le.2) then
       print*,'nmax must be greater than 2 for cyclic system'
       stop
     end if
     gamma = -b(1) !choose gamma to be -b(1)
     bb(1) = b(1)-gamma
     bb(nmax) = b(nmax)-a(1)*c(nmax)/gamma
     bb(2:nmax-1) = b(2:nmax-1)
     call tridiagonal_solve(nmax,a,bb,c,r,xout)
     uu(1) = gamma
     uu(nmax) = a(1) !alpha
     uu(2:nmax-1) = 0.d0
     call tridiagonal_solve(nmax,a,bb,c,uu,zz)
     factor = (xout(1)+c(nmax)*xout(nmax)/gamma)/(1.+zz(1)+c(nmax)*zz(nmax)/gamma)
     xout = xout-zz*factor
   end subroutine cyclictridiagonal_solve

   !subroutine to inverse tridiagonal matrix
   !xout is the inversed matrix with 1st dimension being row
   subroutine tridiagonal_inverse(nmax,a,b,c,xout)
     implicit none
     integer, intent(in) :: nmax
     real(8), dimension(nmax), intent(in) :: a,b,c
     real(8), dimension(nmax,nmax), intent(out) :: xout
     real(8) :: r(nmax)
     integer :: n
     do n=1,nmax
       r = 0.d0
       r(n) = 1.d0
       call tridiagonal_solve(nmax,a,b,c,r,xout(:,n))
     end do
   end subroutine tridiagonal_inverse

   !subroutine to inverse cyclic tridiagonal matrix
   !xout is the inversed matrix with 1st dimension being row
   subroutine cyclictridiagonal_inverse(nmax,a,b,c,xout)
     implicit none
     integer, intent(in) :: nmax
     real(8), dimension(nmax), intent(in) :: a,b,c
     real(8), dimension(nmax,nmax), intent(out) :: xout
     real(8), dimension(nmax) :: aa, bb, cc, uu, vv, ww, zz
     real(8) :: gamma, lambda
     integer :: i,j
     if (nmax.le.2) then
       print*,'nmax must be greater than 2 for cyclic system'
       stop
     end if
     gamma = -b(1) !choose gamma to be -b(1)
     bb(1) = b(1)-gamma
     bb(nmax) = b(nmax)-a(1)*c(nmax)/gamma
     bb(2:nmax-1) = b(2:nmax-1)
     uu(1) = gamma
     uu(nmax) = a(1) !alpha
     uu(2:nmax-1) = 0.d0
     call tridiagonal_solve(nmax,a,bb,c,uu,zz) !solve Az=u
     vv(1) = 1.d0
     vv(nmax) = c(nmax)/gamma
     vv(2:nmax-1) = 0.d0
     aa(2:nmax)=c(1:nmax-1)
     cc(1:nmax-1)=a(2:nmax)
     call tridiagonal_solve(nmax,aa,bb,cc,vv,ww) !solve A^Tw=v
     call tridiagonal_inverse(nmax,a,bb,c,xout) !inverse the tridiagonal part
     lambda = zz(1)+c(nmax)*zz(nmax)/gamma
     do j=1,nmax
       do i=1,nmax
         xout(i,j) = xout(i,j)-(zz(i)*ww(j))/(1.+lambda)
       end do
     end do
   end subroutine cyclictridiagonal_inverse
