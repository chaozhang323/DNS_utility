


    module modCalInt
      implicit none
      real(8), parameter :: gamma = 1.4
      real(8) :: Tr, Tinf, Ma, Tw
      real(8) :: uinf, rhoinf, rbar, Rm
      integer :: iCalMu
      real(8) :: Cp, Cv, Pr
      real(8) :: zplus_ref

    contains

      subroutine Cal_Int_dat(kmax, kinf, nvarin, buffer, nvarout, buffer_output)
        implicit none
        integer, intent(in) :: kmax, kinf, nvarin, nvarout
        real(8), intent(in) :: buffer(kmax,nvarin)
        real(8), intent(out) :: buffer_output(1:kmax,1:nvarout)

        real(8) :: zloc(kmax)
        real(8), dimension(kmax) :: uave, ru, rhoave, dudk, tave, dtdk, divave, dkdz, wave, ruw, rw, rwt, pave, uw, w2
        real(8), dimension(kmax) :: rupwp, rwptp, dTdu
        real(8) :: delta, dstar_k
        integer :: i, j, k, kbl
        integer :: k_lower, k_upper
        real(8) :: muw, muinf, Kappa_w, rg, qw, tauw, Trg
        real(8), dimension(kmax) :: mu, dmudz, dmudk, drhodz, drhodk, vand1, vand2, vand3, muT, dzdk
        real(8), dimension(kmax) :: util, dutildk, dutildz, utau_star, ytau_star, dutildz_star, ZII_plus
        real(8) :: Prt_infty, d_T, Bplus
        real(8) :: ysub_plus, ybuf_plus, zplus

        !Prt_infty = 0.75
        !Bplus = 12.d0

        zloc(1:kmax)   = buffer(1:kmax,1)
        rhoave(1:kmax) = buffer(1:kmax,2)
        uave(1:kmax)   = buffer(1:kmax,3)
        wave(1:kmax)   = buffer(1:kmax,5)
        tave(1:kmax)   = buffer(1:kmax,6)
        pave(1:kmax)   = buffer(1:kmax,7)

        !dkdz(1:kmax) = buffer_grid(1:kmax,2)
        !ru(1:kmax) = buffer(1:kmax,2)
        !dudk(1:kmax) = buffer(1:kmax,4)
        !dtdk(1:kmax) = buffer(1:kmax,6)
        !divave(1:kmax) = buffer(1:kmax,7)
        !ruw(1:kmax) = buffer(1:kmax,9)
        !rw(1:kmax) = buffer(1:kmax,10)
        !rwt(1:kmax) = buffer(1:kmax,11)
        !uw(1:kmax) = buffer(1:kmax,13)
        !w2(1:kmax) = buffer(1:kmax,14)
        !util(1:kmax) = ru(1:kmax)/rhoave(1:kmax)

        vand1 = 0.d0; vand2 = 0.d0; vand3 = 0.d0
        buffer_output = 0.d0

        do k=1, kmax
          if((0.99*uave(kinf)).ge.uave(k)) then
            k_lower = k
          elseif((0.99*uave(kinf)).le.uave(k)) then
            k_upper = k
            exit
          endif
        enddo
        if(k_lower.ge.k_upper) then
          print *, 'cannot find the range for the velocity'
          stop
        endif
        if((k_upper-k_lower).gt.1) then
          print *, 'k_upper and k_lower is wrong'
          stop
        endif

        delta = (0.99*uave(kinf)-uave(k_lower))/(uave(k_upper)-uave(k_lower))*(zloc(k_upper)-zloc(k_lower)) + zloc(k_lower)

        Tinf = tave(kinf)
        Tw = tave(1)
        Ma = uave(kinf)/sqrt(1.4*rbar*Tinf)
        Tr = Tinf*(1.d0+0.89*0.2*Ma**2)

        call Derivative(tave,dtdk)
        call Derivative(uave,dudk)
        call Derivative(zloc,dzdk)
        dkdz = 1.d0/dzdk
        do k=1, kmax
          call Cal_Mu(iCalMu,tave(k),mu(k))
          dTdu(k) = dtdk(k)/(dudk(k)+1.e-30)
        enddo

        Kappa_w = mu(1)*Cp/Pr
        qw = - Kappa_w*dtdk(1)*dkdz(1)
        tauw = mu(1)*dudk(1)*dkdz(1)
        rg = (tave(1) - tave(k_upper))/( uave(k_upper)**2/(2.d0*Cp) ) - 2.d0*Pr/uave(k_upper)*qw/tauw
        Trg = tave(k_upper) + rg*uave(k_upper)**2/(2.d0*Cp)

        call Derivative(mu,dmudk)
        call Derivative(rhoave,drhodk)

        dmudz(1:kmax)  = dmudk(1:kmax)*dkdz(1:kmax)
        drhodz(1:kmax) = drhodk(1:kmax)*dkdz(1:kmax)

        buffer_output(:,1) = zloc(:)
        buffer_output(:,2) = sqrt(abs(tauw/rhoave(1)))                   ! utau
        buffer_output(:,3) = mu(1)/rhoave(1)/buffer_output(1,2)          ! ztau
        buffer_output(:,4) = delta/buffer_output(1,3)                    ! Retau
        buffer_output(:,5) = (abs(tauw)*rhoave(k_upper))**0.5*delta/mu(k_upper)      ! Retau_star

        do k=2, kmax
          vand1(k) = sqrt(rhoave(k)/rhoave(1))*(uave(k)-uave(k-1))
          vand2(k) = sqrt(rhoave(k)/rhoave(1))*(uave(k)-uave(k-1))*(1.d0+0.5/rhoave(k)*drhodz(k)*zloc(k)-1.d0/mu(k)*dmudz(k)*zloc(k))
        enddo
        do k=2, kmax
          buffer_output(k,6) = sum(vand1(1:k))/buffer_output(1,2)                        ! Uvd
          buffer_output(k,7) = sum(vand2(1:k))/buffer_output(1,2)                        ! Uvd_Trettel
          buffer_output(k,8) = (abs(mu(1)*dudk(1)*dkdz(1))*rhoave(k))**0.5*zloc(k)/mu(k) ! zplus_Trettel


      !    rupwp(k) = ruw(k)-ru(k)*wave(k)-rw(k)*uave(k)+rhoave(k)*uave(k)*wave(k)
      !    rwptp(k) = rwt(k)-pave(k)/rbar*wave(k)-rw(k)*tave(k)+rhoave(k)*wave(k)*tave(k)
      !  buffer_output(k,11) = buffer_output(k,10)*( (1.d0+wave(k)*(ru(k)-rhoave(k)*uave(k))/(rupwp(k)+1.e-30 ))/(1.d0+wave(k)*(pave(k)/rbar-rhoave(k)*tave(k))/(rwptp(k)+1.e-30 ) )  ) !!

          buffer_output(k,9) = uave(k)/2.d0*dTdu(k)/( tave(k) - uave(k)/2.d0*dTdu(1)-tave(1) )  ! Pre  eq. 4.9
        enddo

        buffer_output(:,10) = delta
        !buffer_output(:,11) = tave(1)/tave(k_upper) + (Tr-tave(1))/tave(k_upper)*uave(:)/uave(k_upper) + (tave(k_upper)-Tr)/tave(k_upper)*(uave(:)/uave(k_upper) )**2    ! Walz
        !buffer_output(:,12) = tave(1)/tave(k_upper) + (Trg-tave(1))/tave(k_upper)*uave(:)/uave(k_upper) + (tave(k_upper)-Trg)/tave(k_upper)*(uave(:)/uave(k_upper) )**2  ! Walz_g
        buffer_output(:,11) = tave(1)/tave(kinf) + (Tr-tave(1))/tave(kinf)*uave(:)/uave(kinf) + (tave(kinf)-Tr)/tave(kinf)*(uave(:)/uave(kinf) )**2    ! Walz
        buffer_output(:,12) = tave(1)/tave(kinf) + (Trg-tave(1))/tave(kinf)*uave(:)/uave(kinf)+ (tave(kinf)-Trg)/tave(kinf)*(uave(:)/uave(kinf) )**2  ! Walz_g
        buffer_output(:,13) = uave(:)/uave(kinf)   ! u_uinf
        buffer_output(:,14) = tave(:)/tave(kinf)   ! T_Tinf


        !buffer_output(:,27) = dTdk(:)*dkdz(:)
        !buffer_output(:,28) = dudk(:)*dkdz(:)


      end subroutine Cal_Int_dat

      subroutine Cal_Int(buffer_grid, kmax, kinf, nvarin, buffer, nvarout, buffer_output)
        implicit none
        integer, intent(in) :: kmax, kinf, nvarin, nvarout
        real(8), intent(in) :: buffer_grid(kmax,2), buffer(kmax,nvarin)
        real(8), intent(out) :: buffer_output(1:kmax,3:nvarout)

        real(8) :: zloc(kmax)
        real(8), dimension(kmax) :: uave, ru, rhoave, dudk, tave, dtdk, divave, dkdz, wave, ruw, rw, rwt, pave, uw, w2, u2, t2, t0ave
        real(8), dimension(kmax) :: rupwp, rwptp, dTdu
        real(8) :: delta, theta, dstar, dstar_k
        integer :: i, j, k, kbl
        integer :: k_lower, k_upper
        real(8) :: muw, muinf, Kappa_w, rg, qw, tauw, Trg, urms, trms, Ma_loc
        real(8), dimension(kmax) :: mu, dmudz, dmudk, drhodz, drhodk, vand1, vand2, vand3, vand4, muT, dt0dk, dRetaustdk
        real(8), dimension(kmax) :: util, dutildk, dutildz, utau_star, ytau_star, dutildz_star, ZII_plus
        real(8) :: Prt_infty, d_T, Bplus
        real(8) :: ysub_plus, ybuf_plus, zplus

        Prt_infty = 0.75
        Bplus = 12.d0

        zloc(1:kmax) = buffer_grid(1:kmax,1)
        dkdz(1:kmax) = buffer_grid(1:kmax,2)
        uave(1:kmax) = buffer(1:kmax,1)
        ru(1:kmax) = buffer(1:kmax,2)
        rhoave(1:kmax) = buffer(1:kmax,3)
        dudk(1:kmax) = buffer(1:kmax,4)
        tave(1:kmax) = buffer(1:kmax,5)
        dtdk(1:kmax) = buffer(1:kmax,6)
        divave(1:kmax) = buffer(1:kmax,7)
        wave(1:kmax) = buffer(1:kmax,8)
        ruw(1:kmax) = buffer(1:kmax,9)
        rw(1:kmax) = buffer(1:kmax,10)
        rwt(1:kmax) = buffer(1:kmax,11)
        pave(1:kmax) = buffer(1:kmax,12)
        uw(1:kmax) = buffer(1:kmax,13)
        w2(1:kmax) = buffer(1:kmax,14)
        u2(1:kmax) = buffer(1:kmax,15)
        t2(1:kmax) = buffer(1:kmax,16)
        t0ave(1:kmax) = buffer(1:kmax,17)

        util(1:kmax) = ru(1:kmax)/rhoave(1:kmax)

        vand1 = 0.d0; vand2 = 0.d0; vand3 = 0.d0; vand4 = 0.d0
        buffer_output = 0.d0

        ! find the boundary layer edge
        do k=1, kmax
          if((0.99*uave(kinf)).ge.uave(k)) then
            k_lower = k
          elseif((0.99*uave(kinf)).le.uave(k)) then
            k_upper = k
            exit
          endif
        enddo
        if(k_lower.ge.k_upper) then
          print *, 'cannot find the range for the velocity'
          stop
        endif
        if((k_upper-k_lower).gt.1) then
          print *, 'k_upper and k_lower is wrong'
          stop
        endif

        delta = (0.99*uave(kinf)-uave(k_lower))/(uave(k_upper)-uave(k_lower))*(zloc(k_upper)-zloc(k_lower)) + zloc(k_lower)
        theta = 0.
        dstar = 0.
        do k=1, kinf-1
          theta = ( ru(k)*(1.0-uave(k)/uave(kinf))/(ru(kinf))+ru(k+1)*(1.0-uave(k+1)/uave(kinf))/(ru(kinf)) )*(zloc(k+1)-zloc(k))*0.5 + theta

          dstar = ((1.0-ru(k)/ru(kinf))+(1.0-ru(k+1)/ru(kinf)))*(zloc(k+1)-zloc(k))*0.5 + dstar
        enddo

        dstar_k = 0.d0
        do k=1, kinf-1
         dstar_k = ( (1.d0-uave(k+1)/uave(kinf)) + (1.d0-uave(k)/uave(kinf)) )*(zloc(k+1)-zloc(k))/2.d0 + dstar_k
        enddo
        Tinf = tave(kinf)
        Tw = tave(1)
        Ma = uave(kinf)/sqrt(1.4*rbar*Tinf)
        !Tr = Tinf*(1.d0+0.89*0.2*Ma**2)
        Tr = Tinf*(1.d0+0.9*0.2*Ma**2)

        do k=1, kmax
          call Cal_Mu(iCalMu,tave(k),mu(k))
          dTdu(k) = dtdk(k)/(dudk(k)+1.e-30)
        enddo
        Kappa_w = mu(1)*Cp/Pr
        qw = - Kappa_w*dtdk(1)*dkdz(1)
        tauw = mu(1)*dudk(1)*dkdz(1)
        rg = (tave(1) - tave(kinf))/( uave(kinf)**2/(2.d0*Cp) ) - 2.d0*Pr/uave(kinf)*qw/tauw
        Trg = tave(kinf) + rg*uave(kinf)**2/(2.d0*Cp)

        call Derivative(mu,dmudk)
        call Derivative(rhoave,drhodk)
        call Derivative(util,dutildk)
        call Derivative(t0ave,dt0dk)

        dmudz(1:kmax)  = dmudk(1:kmax)*dkdz(1:kmax)
        drhodz(1:kmax) = drhodk(1:kmax)*dkdz(1:kmax)
        dutildz(1:kmax) = dutildk(1:kmax)*dkdz(1:kmax)

        do k=1, kmax
          utau_star(k) = sqrt(tauw/rhoave(k))
          ytau_star(k) = mu(k)/(rhoave(k)*utau_star(k))
          dutildz_star(k) = dutildz(k)/(utau_star(k)/ytau_star(k))
        enddo

!        buffer_output(:,32) = utau_star(:) ! ustar
!        buffer_output(:,33) = zloc(:)*utau_star(:)/(mu(:)/rhoave(:))     ! zstar

        buffer_output(:,3) = sqrt(abs(tauw/rhoave(1)))                   ! utau
        buffer_output(:,4) = mu(1)/rhoave(1)/buffer_output(1,3)          ! ztau
        buffer_output(:,5) = delta/buffer_output(1,4)                    ! Retau

        buffer_output(:,6) = (abs(tauw)*rhoave(k_upper))**0.5*delta/mu(k_upper)      ! Retau_star

        do k=2, kmax
          vand1(k) = sqrt(rhoave(k)/rhoave(1))*(uave(k)-uave(k-1))
          vand2(k) = sqrt(rhoave(k)/rhoave(1))*(uave(k)-uave(k-1))*(1.d0+0.5/rhoave(k)*drhodz(k)*zloc(k)-1.d0/mu(k)*dmudz(k)*zloc(k))
        enddo
        do k=2, kmax
          buffer_output(k,7) = sum(vand1(1:k))/buffer_output(1,3)                        ! Uvd
          buffer_output(k,8) = sum(vand2(1:k))/buffer_output(1,3)                        ! Uvd_Trettel
          buffer_output(k,9) = (abs(mu(1)*dudk(1)*dkdz(1))*rhoave(k))**0.5*zloc(k)/mu(k) ! zplus_Trettel

          rupwp(k) = ruw(k)-ru(k)*wave(k)-rw(k)*uave(k)+rhoave(k)*uave(k)*wave(k)
          rwptp(k) = rwt(k)-pave(k)/rbar*wave(k)-rw(k)*tave(k)+rhoave(k)*wave(k)*tave(k)

          buffer_output(k,10) = rupwp(k)/(rwptp(k)+1.e-30)*dtdk(k)/(dudk(k)+1.e-30)                       ! Prt
          buffer_output(k,11) = buffer_output(k,10)*( (1.d0+wave(k)*(ru(k)-rhoave(k)*uave(k))/(rupwp(k)+1.e-30 ))/(1.d0+wave(k)*(pave(k)/rbar-rhoave(k)*tave(k))/(rwptp(k)+1.e-30 ) )  )
          buffer_output(k,12) = uave(k)/2.d0*dTdu(k)/( tave(k) - uave(k)/2.d0*dTdu(1)-tave(1) )  ! Pre  eq. 4.9
          ! d_T = 1.d0- exp( -zloc(k)/buffer_output(1,4)/Bplus )
          ! buffer_output(k,13) = Prt_infty/d_T
        enddo

        do k=1, kmax
          buffer_output(k,14) = -rhoave(k)*( uw(k)-uave(k)*wave(k) ) ! tau_xz    -<r><u'w'>
          buffer_output(k,15) = -rhoave(k)*( w2(k)-wave(k)**2 )      ! tau_zz    -<r><w'w'>
        enddo

        ! theta_11, theta_22, theta_12 for Poggie 2015
        call Cal_energyFlux(kmax,buffer_output(:,14),buffer_output(:,15),dtdk(:)*dkdz(:),dudk(:)*dkdz(:),buffer_output(:,16),buffer_output(:,17),buffer_output(:,18))
        ! x2  for Poggie 2015
        buffer_output(:,19) = (0.68/0.41)*uave(k_upper)/buffer_output(1,3)*dstar_k !x2_smooth

        buffer_output(:,20) = delta  ! delta
        buffer_output(:,21) = tave(1)/tave(kinf) + (Tr-tave(1))/tave(kinf)*uave(:)/uave(kinf) + (tave(kinf)-Tr)/tave(kinf)*(uave(:)/uave(kinf) )**2    ! Walz
        buffer_output(:,22) = tave(1)/tave(kinf) +(Trg-tave(1))/tave(kinf)*uave(:)/uave(kinf) +(tave(kinf)-Trg)/tave(kinf)*(uave(:)/uave(kinf) )**2  ! Walz_g
        buffer_output(:,23) = uave(:)/uave(kinf)   ! u_uinf
        buffer_output(:,24) = tave(:)/tave(kinf)   ! T_Tinf

        ! calculate zplus using Wu's method
        ysub_plus =  9.8*(tave(1)/Tr)**(-0.5)
        ybuf_plus = 43.8*(tave(1)/Tr)**(-0.5)
        do k=1, kmax
          zplus = zloc(k)/buffer_output(1,4)
          if(zplus.lt.ysub_plus) then
            ZII_plus(k) = ( (sqrt(rhoave(k)/rhoave(1))/(mu(k)/mu(1)) )**(2.d0/3.d0) )*((ysub_plus/9.8)**(1.d0/3.d0) )/((ybuf_plus/43.8)**(2.d0/3.d0))*zplus
          elseif(ysub_plus.le.zplus.and.zplus.lt.ybuf_plus) then
            ZII_plus(k) = ( (sqrt(rhoave(k)/rhoave(1))/(mu(k)/mu(1)) )**0.5 )/((ybuf_plus/43.8)**0.5)*zplus
          else
            ZII_plus(k) = (sqrt(rhoave(k)/rhoave(1)))/(mu(k)/mu(1))*zplus
          endif
        enddo
        do k=2, kmax
          vand3(k) = dutildz_star(k)*(ZII_plus(k)-ZII_plus(k-1))
        enddo
        buffer_output(:,25) = ZII_plus(:)             ! zplus_Wu
        do k=2, kmax
          buffer_output(k,26) = sum(vand3(1:k))       ! Uvd_Wu
        enddo

        ! calculate SRA using Huang's and Zhang's method
        do k=2, kmax
         urms = sqrt(abs(u2(k)-uave(k)**2))
         trms = sqrt(abs(t2(k)-tave(k)**2))
         Ma_loc = uave(k)/sqrt(gamma*rbar*tave(k))
         buffer_output(k,27) = trms/tave(k)*(1.d0-dt0dk(k)/(dtdk(k)+1.e-30))/( (gamma-1)*Ma_loc**2*urms/uave(k)+1.e-30 )*buffer_output(k,10)  ! SRA_Huang
         buffer_output(k,28) = trms/tave(k)*(1.d0-dt0dk(k)/(dtdk(k)+1.e-30))/( (gamma-1)*Ma_loc**2*urms/uave(k)+1.e-30 )*buffer_output(k,11)  ! SRA_Zhang
        enddo

        do k=1, kmax
          buffer_output(k,29) = mu(k)*dudk(k)*dkdz(k)                                                       ! mu*dudz
          buffer_output(k,30) = - ( ruw(k) - ru(k)*wave(k) - rw(k)*uave(k) + rhoave(k)*uave(k)*wave(k) )    ! - <ru'w'>
          buffer_output(k,31) = buffer_output(k,29) + buffer_output(k,30)                                   ! mu*dudz + (- <ru'w'>)
        enddo
        buffer_output(:,32) = dkdz(:)  ! dkdz

        do k=1, kmax
          buffer_output(k,33) =  buffer_output(k,5)*sqrt(rhoave(k)/rhoave(1))/(mu(k)/mu(1))   ! Retau_star_Patel
          buffer_output(k,34) = zloc(k)*buffer_output(k,33)/buffer_output(k,20)             ! zstar_Patel
        enddo

        call Derivative(buffer_output(1:kmax,33),dRetaustdk)
        do k=2, kmax
          vand4(k) = ( 1.d0+zloc(k)/buffer_output(k,33)*dRetaustdk(k)*dkdz(k))*(buffer_output(k,7)-buffer_output(k-1,7))
        enddo
        do k=2, kmax
          buffer_output(k,35) = sum(vand4(1:k))  ! ustar_Patel
        enddo

        do k=1, kmax  ! Davide Modesti & Sergio Pirozzoli, <Reynolds and Mach number effects in compressible turbulent channel flow>
          buffer_output(k,36) = ( (buffer_output(k,3)*delta)**(0.5d0) )/sqrt(abs(dudk(k)*dkdz(k)+1.e-30))  ! l12
          buffer_output(k,37) = ( (sqrt(abs(tauw/rhoave(k)))*delta)**(0.5d0) )/sqrt(abs(dutildk(k)*dkdz(k)+1.e-30))  ! l12_star
          buffer_output(k,38) = buffer_output(k,3)/(dudk(k)*dkdz(k)+1.e-30) ! lm
          buffer_output(k,39) = sqrt(abs(tauw/rhoave(k)))/(dutildk(k)*dkdz(k)+1.e-30) ! lm_star
        enddo
   
        buffer_output(:,40) = theta
        buffer_output(:,41) = dstar
        buffer_output(:,42) = dstar/theta
        call Cal_Mu(iCalMu,Tinf,muinf)
        buffer_output(:,43) = rhoave(kinf)*uave(kinf)*theta/muinf
        call Cal_Mu(iCalMu,tave(1),muw)
        buffer_output(:,44) = rhoave(kinf)*uave(kinf)*theta/muw

!print *, 'buffer_output(k,3)*delta ', buffer_output(2,3)*delta, 'dudk(k)*dkdz(k) = ', dudk(2)*dkdz(2)

      end subroutine Cal_Int

      ! inner and outer boundary layer region  z^+ = 30
      subroutine Cal_Mu_turbulent_viscosity(kdim,ztau,dstar_k,Ue,delta,rho,x2,du1dk,dkdz,muout)
        integer, intent(in) :: kdim
        real(8), intent(in) :: ztau, dstar_k, Ue, delta
        real(8), dimension(kdim), intent(in) :: rho, x2, du1dk, dkdz
        real(8), dimension(kdim), intent(out) :: muout
        real(8) :: kp, Aplus, C, Cc, Ck, nK
        real(8) :: zplus!, zplus_ref
        integer :: k

        kp = 0.41; Aplus=26.d0
        Cc = 0.018; Ck = 1.2
        C  = 0.68
        nK = 6.d0

        do k=1, kdim
          zplus = x2(k)/ztau

          if(zplus.le.zplus_ref) then !! inner region
            muout(k) = rho(k)*(( kp*x2(k)*(1.d0-exp(-(x2(k)/ztau)/Aplus) )  )**2.d0 )*abs(du1dk(k)*dkdz(k))
          else         !! outer region
            muout(k) = Cc*rho(k)*Ue*dstar_k/(1.d0+Ck*(x2(k)/delta)**nK)
          endif

        enddo

      end subroutine Cal_Mu_turbulent_viscosity

      subroutine Cal_energyFlux(kdim,tau12,tau22,dTdx2,du1dx2,theta1,theta2,theta3)
        integer, intent(in) :: kdim
        real(8), dimension(kdim), intent(in) :: tau12, tau22, dTdx2, du1dx2
        real(8), dimension(kdim), intent(out) :: theta1, theta2, theta3
        real(8) :: tau_theta, sigma_theta, sigma_e, tau_u, tau_e, a1, Cmu
        integer :: k

        sigma_theta = 0.28/gamma
        sigma_e = 0.72/gamma
        a1 = 0.28; Cmu=0.09

        !Cp = gamma*rbar/(gamma-1.d0)
        !Cv = 2.5*rbar

        do k=1, kdim
          tau_u = a1/Cmu/(du1dx2(k)+1.e-30)
          tau_theta = sigma_theta*tau_u
          tau_e = sigma_e*tau_u

          theta1(k) = tau12(k)*Cp*dTdx2(k)*tau_theta - tau22(k)*Cp*dTdx2(k)*du1dx2(k)*tau_theta**2
          theta2(k) = tau22(k)*Cp*dTdx2(k)*tau_theta
          !theta3(k) = -2.d0*tau22(k)*(Cp*dTdx2(k))**2*tau_theta*tau_e
          theta3(k) = -2.d0*tau22(k)*(Cp*dTdx2(k))*(Cp*dTdx2(k))*tau_theta*tau_e

        enddo

        theta1 = theta1/Cv  !!!!!!!!!!!!!!!!!!!
        theta2 = theta2/Cv
        theta3 = theta3/Cv/Cv

      end subroutine Cal_energyFlux

      subroutine Derivative(var,dvardx)
        real(8), intent(in) :: var(:)
        real(8), intent(out) :: dvardx(size(var))
        integer :: i, ilen

        ilen = size(var)
        do i=3, ilen-2
          dvardx(i) = (-var(i+2)+8.0*var(i+1)-8.0*var(i-1)+var(i-2))/(12.d0)
        enddo
        dvardx(1) = -25./12.*var(1)+4.*var(2)-3.*var(3)+4./3.*var(4)-0.25*var(5)
        dvardx(2) = -0.25*var(1)-5./6.*var(2)+1.5*var(3)-0.5*var(4)+1./12.*var(5)
        dvardx(ilen-1) = -1./12.*var(ilen-4)+0.5*var(ilen-3)-1.5*var(ilen-2) + 5./6.*var(ilen-1)+0.25*var(ilen)
        dvardx(ilen) = 0.25*var(ilen-4)-4./3.*var(ilen-3)+3.*var(ilen-2)- 4.*var(ilen-1)+25./12.*var(ilen)

      end subroutine Derivative



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

      real function cf_VanDriestII(Re_theta)
        real, intent(in) :: Re_theta
        real :: FRe0,Fc,alpha,beta,A,B
        real :: r,Taw,cf_bar
        r=0.89
        Taw=Tinf*(1.d0+0.5*r*(1.4-1.d0)*(Ma**2))
        alpha=sqrt(0.5*r*(1.4-1.d0)*(Ma**2)*Tinf/Tw)
        beta=Taw/Tw-1.d0
        A=(2*alpha**2-beta)/(sqrt(beta**2+4.d0*alpha**2))
        B=beta/(sqrt(beta**2+4.d0*alpha**2))
        FRe0=((Tinf/Tw)**(1.5))*(Tw+110.4)/(Tinf+110.4)
        Fc=((Taw/Tinf)-1.d0)/((asin(A)+asin(B))**2)
        cf_bar=Cf_Karman(FRe0*Re_theta)
        cf_VanDriestII=cf_bar/Fc
      end function cf_VanDriestII

    end module modCalInt














