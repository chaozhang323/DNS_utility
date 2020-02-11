


    module modCalInt
      implicit none

    contains

!      subroutine Cal_Int(kdim, kmax, kinf, uave, ruave, rhoave, dudk, tave, dkdz, delta, theta, dstar, Cf)
      subroutine Cal_Int(buffer_grid, kmax, kinf, nvarin, buffer, nvarout, buffer_output, Tinf, Ma)
        implicit none
        integer, intent(in) :: kmax, kinf, nvarin, nvarout
        real(8), intent(in) :: buffer_grid(kmax,2), buffer(kmax,nvarin)
        real(8), intent(in) :: Tinf, Ma
        real(8), intent(out) :: buffer_output(nvarout-1)


        ! integer :: kmax, kinf
        real(8) :: kdim(kmax)
        real(8) :: uave(kmax), ruave(kmax),rhoave(kmax), dudk(kmax), tave(kmax), dtdk(kmax),dkdz(kmax)
        real(8) :: delta, theta, dstar, Cf, qw, Ch, Cf_inf, Ch_inf
        integer :: i, j, k, kbl
        integer :: k_lower, k_upper
        real(8) :: muw, kappa, cp = 1004.5, Pr = 0.71
        real(8) :: Tr
        integer :: temp(1)

        kdim(1:kmax) = buffer_grid(1:kmax,1)
        dkdz(1:kmax) = buffer_grid(1:kmax,2)
        uave(1:kmax) = buffer(1:kmax,1)
        ruave(1:kmax) = buffer(1:kmax,2)
        rhoave(1:kmax) = buffer(1:kmax,3)
        dudk(1:kmax) = buffer(1:kmax,4)
        tave(1:kmax) = buffer(1:kmax,5)
        dtdk(1:kmax) = buffer(1:kmax,6)
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


        muw = (1.458E-6)*((tave(1))**1.5)/(tave(1)+110.4)
        ! Cf = muw*dudk(1)*dkdz(1)
        Cf_inf = muw*dudk(1)*dkdz(1)
        Cf = muw*dudk(1)*dkdz(1)/(0.5*rhoave(k_upper)*uave(k_upper)**2)

        kappa = cp*muw/Pr
        qw = kappa*dtdk(1)*dkdz(1)
        Tr = Tinf*(1.d0+0.89*0.2*Ma**2)
        Ch = qw/(rhoave(k_upper)*uave(k_upper)*cp*(Tr-tave(1)))
        Ch_inf = qw/(cp*(Tr-tave(1)))

        buffer_output(1) = delta
        buffer_output(2) = theta
        buffer_output(3) = dstar
        buffer_output(4) = Cf
        buffer_output(5) = qw
        buffer_output(6) = Ch
        buffer_output(7) = Cf_inf
        buffer_output(8) = Ch_inf

!!************************************************
!        Retau = rhoave(1)*utau(1)*delta/muw
!        temp = minloc(abs( uave(1:kmax)/uave(kinf) - 0.99 ))
!        kbl = temp(1)
!        delta = kdim(kbl)
!        theta = sum(ruave(1:kinf)/ruave(kinf)*(1.0-uave(1:kinf)/uave(kinf))*1.0/dkdz(1:kinf)  )
!        dstar = sum((1.0-ruave(1:kinf)/ruave(kinf) )*1.0/dkdz(1:kinf)  )
!!************************************************

      end subroutine Cal_Int

    end module modCalInt














