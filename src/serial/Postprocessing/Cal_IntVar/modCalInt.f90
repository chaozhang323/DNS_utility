


    module modCalInt
      implicit none
      real(8), parameter :: R=8314.3D0
      real(8) :: Tr, Tinf, Ma, Tw, Rm, rbar, cp, uinf, rhoinf
      integer :: iCalMu

    contains

      !subroutine Cal_Int(buffer_grid, kmax, kinf, nvarin, buffer, nvarout, buffer_output, Tinf, Ma)
      subroutine Cal_Int(buffer_grid, kmax, kinf, nvarin, buffer, nvarout, buffer_output)
        implicit none
        integer, intent(in) :: kmax, kinf, nvarin, nvarout
        real(8), intent(in) :: buffer_grid(kmax,2), buffer(kmax,nvarin)
        real(8), intent(out) :: buffer_output(2:nvarout)

        real(8) :: kdim(kmax)
        real(8) :: uave(kmax), ruave(kmax), rhoave(kmax), dudk(kmax), tave(kmax), dtdk(kmax), divave(kmax), dkdz(kmax)
        real(8) :: delta, theta, dstar, Cf, qw, Ch, Cf_inf, Ch_inf
        integer :: i, j, k, kbl
        integer :: k_lower, k_upper
        real(8) :: muw, muinf, kappa, Pr = 0.71

        integer :: temp(1)
        integer :: ks, ke

        kdim(1:kmax) = buffer_grid(1:kmax,1)
        dkdz(1:kmax) = buffer_grid(1:kmax,2)
        uave(1:kmax) = buffer(1:kmax,1)
        ruave(1:kmax) = buffer(1:kmax,2)
        rhoave(1:kmax) = buffer(1:kmax,3)
        dudk(1:kmax) = buffer(1:kmax,4)
        tave(1:kmax) = buffer(1:kmax,5)
        dtdk(1:kmax) = buffer(1:kmax,6)
        divave(1:kmax) = buffer(1:kmax,7)

        temp = minloc(divave(5:kmax-1))
        ks = temp(1) + 4
        !kinf = int(0.9*ks)            ! setup kinf
        !kinf = ks - 20

        !print *, 'ks = ', ks, 'kinf = ', kinf

        do k=1, kmax
          if((0.99*uave(kinf)).ge.uave(k)) then
            k_lower = k
          elseif((0.99*uave(kinf)).le.uave(k)) then
            k_upper = k
            exit
          endif
        enddo
        ! print *, 'k_lower = ', k_lower , 'k_upper = ', k_upper
        if(k_lower.ge.k_upper) then
          print *, 'cannot find the range for the velocity'
          stop
        endif
        if((k_upper-k_lower).gt.1) then
          print *, 'k_upper and k_lower is wrong'
          stop
        endif
        ke = k_upper
        delta = (0.99*uave(kinf)-uave(k_lower))/(uave(k_upper)-uave(k_lower))*(kdim(k_upper)-kdim(k_lower)) + kdim(k_lower)

        theta = 0.
        dstar = 0.
        do k=1, kinf-1
          theta = ( ruave(k)*(1.0-uave(k)/uave(kinf))/(ruave(kinf))+ruave(k+1)*(1.0-uave(k+1)/uave(kinf))/(ruave(kinf)) )*(kdim(k+1)-kdim(k))*0.5 + theta

          dstar = ((1.0-ruave(k)/ruave(kinf))+(1.0-ruave(k+1)/ruave(kinf)))*(kdim(k+1)-kdim(k))*0.5 + dstar
        enddo

!!******************************************
!      print *, 'Calculate dudk, dtdk, dkdz using input data '
!      dudk(1) = -25./12.*uave(1) + 4.*uave(2) - 3.*uave(3) + 4./3.*uave(4) - 0.25*uave(5)
!      dtdk(1) = -25./12.*tave(1) + 4.*tave(2) - 3.*tave(3) + 4./3.*tave(4) - 0.25*tave(5)
!      dkdz(1) = 20576.86
!!******************************************
        Tinf = tave(kinf)
        Tw = tave(1)
        Ma = uave(kinf)/sqrt(1.4*287*Tinf)

        call Cal_Mu(iCalMu,tave(1),muw)
        call Cal_Mu(iCalMu,Tinf,muinf)

        !muw = (1.458E-6)*((tave(1))**1.5)/(tave(1)+110.4)
        !muinf = (1.458E-6)*((Tinf)**1.5)/(Tinf+110.4)

        Cf_inf = muw*dudk(1)*dkdz(1)/(0.5*rhoave(kinf)*uave(kinf)**2)
        Cf = muw*dudk(1)*dkdz(1)/(0.5*rhoave(k_upper)*uave(k_upper)**2)

        kappa = cp*muw/Pr
        qw = kappa*dtdk(1)*dkdz(1)
        Tr = Tinf*(1.d0+0.89*0.2*Ma**2)
        Ch = qw/(rhoave(k_upper)*uave(k_upper)*cp*(Tr-tave(1)))
        Ch_inf = qw/(cp*(Tr-tave(1)))/rhoave(kinf)/uave(kinf)

        buffer_output(2) = delta
        buffer_output(3) = theta
        buffer_output(4) = dstar
        buffer_output(5) = Cf
        buffer_output(6) = qw
        buffer_output(7) = Ch
        buffer_output(8) = Cf_inf
        buffer_output(9) = Ch_inf

        buffer_output(10) = sqrt(abs(muw*dudk(1)*dkdz(1)/rhoave(1))) ! utau
        buffer_output(11) = muw/rhoave(1)/buffer_output(10)          ! ztau
        buffer_output(12) = delta/buffer_output(11)                 ! Retau

        buffer_output(13) = rhoave(kinf)*uave(kinf)*theta/muinf     ! Retheta
        buffer_output(14) = rhoave(kinf)*uave(kinf)*theta/muw       ! Redelta2
        buffer_output(15) = dstar/theta                             ! H
        !buffer_output(16) = cf_VanDriestII(buffer_output(13))       ! Cf_VanDriest
        call VanDriestII(buffer_output(13),buffer_output(19),buffer_output(16)) ! Cf_VanDriest
        buffer_output(17) = ks                                      ! kindex for shock
        buffer_output(18) = ke                                      ! kindex for boundary layer edge

!!************************************************
!        Retau = rhoave(1)*utau(1)*delta/muw
!        temp = minloc(abs( uave(1:kmax)/uave(kinf) - 0.99 ))
!        kbl = temp(1)
!        delta = kdim(kbl)
!        theta = sum(ruave(1:kinf)/ruave(kinf)*(1.0-uave(1:kinf)/uave(kinf))*1.0/dkdz(1:kinf)  )
!        dstar = sum((1.0-ruave(1:kinf)/ruave(kinf) )*1.0/dkdz(1:kinf)  )
!!************************************************

      end subroutine Cal_Int

      subroutine Cal_Mu(iCalMu,Tin,Muout)
        integer, intent(in) :: iCalMu
        real(8), intent(in) :: Tin
        real(8), intent(out) :: Muout

        if(iCalMu.eq.0) then ! Sutheland law for air
          Muout = (1.458d-6)*((Tin)**1.5)/(Tin+110.4)
        else  ! Keyes' Law for N2
          Muout = 1.418d-6*(Tin**1.5)/(Tin+116.4*10.**(-5./Tin))
        endif

      end subroutine Cal_Mu

      real function mu_Sutherland(T)
        real, intent(in) :: T
        mu_Sutherland=(1.458e-6)*T**(1.5)/(T+110.4)
      end function mu_Sutherland

      real function Cf_Karman(Reth)
        real, intent(in) :: Reth
        Cf_Karman=1.0/(17.08*(log10(Reth))**2+25.11*(log10(Reth))+6.012)
      end function Cf_Karman

      subroutine VanDriestII(Re_theta,cf_bar,cf_VanDriestII)
        real(8), intent(in) :: Re_theta
        real(8), intent(out) :: cf_bar, cf_VanDriestII
        real(8) :: FRe0,Fc,alpha,beta,A,B
        real(8) :: r,Taw, muw, muinf

        r=0.89
        Taw=Tinf*(1.d0+0.5*r*(1.4-1.d0)*(Ma**2))
        alpha=sqrt(0.5*r*(1.4-1.d0)*(Ma**2)*Tinf/Tw)
        beta=Taw/Tw-1.d0
        A=(2*alpha**2-beta)/(sqrt(beta**2+4.d0*alpha**2))
        B=beta/(sqrt(beta**2+4.d0*alpha**2))
        !FRe0=((Tinf/Tw)**(1.5))*(Tw+110.4)/(Tinf+110.4)
        call Cal_Mu(iCalMu,Tw,muw)
        call Cal_Mu(iCalMu,Tinf,muinf)
        FRe0 = muinf/muw
        Fc=((Taw/Tinf)-1.d0)/((asin(A)+asin(B))**2)
        cf_bar=Cf_Karman(FRe0*Re_theta)
        cf_VanDriestII=cf_bar/Fc
      end subroutine VanDriestII

!      real function cf_VanDriestII(Re_theta)
!        real, intent(in) :: Re_theta
!        real :: FRe0,Fc,alpha,beta,A,B
!        real :: r,Taw,cf_bar
!        r=0.89
!        Taw=Tinf*(1.d0+0.5*r*(1.4-1.d0)*(Ma**2))
!        alpha=sqrt(0.5*r*(1.4-1.d0)*(Ma**2)*Tinf/Tw)
!        beta=Taw/Tw-1.d0
!        A=(2*alpha**2-beta)/(sqrt(beta**2+4.d0*alpha**2))
!        B=beta/(sqrt(beta**2+4.d0*alpha**2))
!        FRe0=((Tinf/Tw)**(1.5))*(Tw+110.4)/(Tinf+110.4)
!        Fc=((Taw/Tinf)-1.d0)/((asin(A)+asin(B))**2)
!        cf_bar=Cf_Karman(FRe0*Re_theta)
!        cf_VanDriestII=cf_bar/Fc
!      end function cf_VanDriestII

    end module modCalInt














