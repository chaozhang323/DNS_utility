      module MGrid1d
      implicit none
contains

!****************************************************************************
! Generate uniform distribution
! Inputs:
!      N:  Number of grid points 
! Outputs:
!      s: 1D uniform distribution s=[0,1]

      pure function dist_uniform(N)
         integer, intent(in) :: N
         real(8) :: dist_uniform(N)
         integer :: i
      
         dist_uniform = 0.0
         forall(i=1:N, N.gt.1)
           dist_uniform(i) = dble(i-1)/dble(N-1)
         end forall
      end function dist_uniform

!****************************************************************************
! Generate geometric distribution 
! Inputs:      
!      N:  Number of grid points
!    ratio: grid stretching resolution (ratio = zL/z2)

      function dist_geometric(N,ratio)
        integer, intent(in) :: N
        real(8), intent(in) :: ratio
        real(8) :: dist_geometric(N)
        real(8) :: a1, a2, f1, f2, fdf, ds1, alpha
        integer, parameter :: maxiter = 100000
        real(8), parameter :: tol = 1.d-5
        integer :: nn, i

        if(N.lt.1) then
           print *, 'N too small for geometric dist. STOP'
           stop
        endif
        ds1 = 1.d0/ratio
        a1 = 1.001; a2 = 1.002
        f1 = (a1**(N-1) - 1.0)/(a1 - 1.0) - ratio
        f2 = (a2**(N-1) - 1.0)/(a2 - 1.0) - ratio
        do nn = 1, maxiter
           fdf = f2*(a2 - a1)/(f2 - f1)
           if (abs(fdf).gt.0.1) fdf = fdf/abs(fdf)*0.1
           a1 = a2
           f1 = f2
           a2 = a2 - fdf
           f2 = (a2**(N-1) - 1.0)/(a2 - 1.0) - ratio
           if (abs(f2).lt.tol) exit
        enddo
        if(nn.ge.maxiter) then
           write(*,*) ' convergence failure in calculating grid ratio '
           stop
        endif
        alpha = a2
        dist_geometric = dist_uniform(N)
        forall(i=2:N, alpha.ne.1.0)
           dist_geometric(i) = ds1*(alpha**(i-1)-1.d0)/(alpha-1.d0)
        end forall
        dist_geometric = dist_geometric/dist_geometric(N)
      end function dist_geometric

!*****************************************************************************

! Generate distribution with grid clustering at both ends
! Inputs:
!      N:  Number of grid points 
!      ds1, ds2: grid resolution at both ends
! Outputs:
!      s: 1D distribution s=[0,1]
! Stretching function (Ref: Marcel Vinokur, JCP, 50, 215-234, 1983)

! u(eta) = 1/2*( 1+tanh(d*(eta-1/2))/tanh(d/2) )
! d can be calculated by solving sinh(d)/d = B, with B = 1/(N*sqrt(ds1*ds2)) if B>1
! d can be calculated by solving sin(d)/d = B, with B = 1/(N*sqrt(ds1*ds2)) if B<1

! s(eta) = u(eta)/( A+(1-A)*u(eta) ), with A = sqrt(ds2/ds1)
! the actual grid can be generated using z = z1 + (z2-z1)*s(eta), s = [0,1]

      pure function dist_double_stretch(N,ds1,ds2)
         integer, intent(in) :: N
         real(8), intent(in) :: ds1, ds2
         real(8) :: dist_double_stretch(N)
         real(8) :: A, B, B1, B2, B3, B4, B5, B6
         real(8) :: d, v, w1, w2, w3, w4
         real(8) :: u(N), s(N)
         real(8), parameter :: pi = 3.14159265
         real(8), parameter :: eps = 1.d-3
         integer :: i

         A = sqrt(ds2/ds1)
         B = 1.d0/(dble(N-1)*sqrt(ds1*ds2))
         if(B.gt.1.d0+eps) then
           ! solve sinh(d)/d = B
           if(B.lt.2.7829681) then
              B1 = B-1.0; B2 = B1*B1; B3 = B1*B2; B4 = B2*B2; B5 = B2*B3
              d = sqrt(6.0*B1)*(1.0-0.15*B1+0.057321429*B2 &
                 -0.024907295*B3+0.0077424461*B4-0.0010794123*B5)
            else
              v = log(B); w1 = 1.d0/B-0.028527431
              w2 = w1*w1; w3 = w1*w2; w4 = w2*w2
              d = v + (1.d0+1.d0/v)*log(2.d0*v)-0.02041793+0.24902722*w1 &
                  +1.9496443*w2-2.6294547*w3+8.56795911*w4
            endif
         endif

         if(B.lt.1.d0-eps) then
             if(B.lt.0.26938972) then
                B1 = B; B2 = B1*B1; B3 = B1*B2; B4 = B2*B2; B5 = B2*B3; B6 = B3*B3
                d = pi*(1.0-B1+B2-(1.0+pi**2/6.0)*B3+6.794732*B4 &
                    -13.205501*B5+11.726095*B6)

             else
                B1=1.0-B; B2 = B1*B1; B3 = B1*B2; B4 = B2*B2; B5 = B2*B3
                d = sqrt(6.0*B1)*(1.0+0.15*B1+0.057321429*B2+0.048774238*B3 &
                     -0.053337753*B4+0.075845134*B5)
             endif
         endif
     
!      print *, 'B=',B, 'd=',d
         forall(i=1:N)
            u(i) = 1.0/2.0*( 1.0 + &
                   tanh(d*(dble(i-1)/(N-1)-1.0/2.0))/tanh(d/2.0) )
         end forall

         if(abs(B-1.0).le.eps) then
            forall(i=1:N)
               u(i) = dble(i-1)/(N-1)*( 1.0+2.0*(B-1.0)*(dble(i-1)/(N-1)-0.5) &
                            *(1.0-dble(i-1)/(N-1)) )
            end forall
         endif
         s = u/( A + (1.d0-A)*u )
         dist_double_stretch = s
      end function dist_double_stretch
      end module MGrid1d
