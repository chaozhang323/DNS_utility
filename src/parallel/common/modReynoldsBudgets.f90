!----------------------------------------------------------------------------------------------
! MODULE: MReynoldsBudgets
!
!> @author
!> Chao Zhang
!
! DESCRIPTION:
!> Calcluate the Reynolds stress Budgets variables.
!
! REVISION HISTORY:
! 08/21/2014 - Initial Version
!----------------------------------------------------------------------------------------------

    module MReynoldsBudgets
      use MDerivative
      use MAverage
      implicit none

      integer, private :: idim, jdim, kdim
      real(8),dimension(:,:,:),allocatable, private :: x,y,z
      real(8),dimension(:,:,:),allocatable, private :: uave,vave,wave,rhoave,pave
      real(8),dimension(:,:,:),allocatable, private :: uave2,vave2,wave2
      real(8),dimension(:,:,:),allocatable, private :: uflu,vflu,wflu,pflu,rhoflu,mu
      real(8),dimension(:,:,:),allocatable, private :: u2flu,v2flu,w2flu
      real(8),dimension(:,:,:),allocatable, private :: u2fluave,v2fluave,w2fluave
      real(8),dimension(:,:,:),allocatable, private :: S11,S22,S33,S12,S13,S23
      real(8),dimension(:,:,:),allocatable, private :: Sig11,Sig22,Sig33,Sig12,Sig13,Sig23
      real(8),dimension(:,:,:),allocatable, private :: Sig11ave,Sig22ave,Sig33ave,Sig12ave,Sig13ave,Sig23ave
      real(8),dimension(:,:,:),allocatable, private :: Sig11ave2,Sig22ave2,Sig33ave2,Sig12ave2,Sig13ave2,Sig23ave2
      real(8),dimension(:,:,:),allocatable, private :: Sig11flu,Sig22flu,Sig33flu,Sig12flu,Sig13flu,Sig23flu
      real(8),dimension(:,:,:),allocatable, private :: Sig112flu,Sig222flu,Sig332flu,Sig122flu,Sig132flu,Sig232flu
      real(8),dimension(:,:,:),allocatable, private :: tau11,tau22,tau33,tau12,tau13,tau23
      real(8),dimension(:,:,:,:,:),allocatable,private :: met

contains

!> @brief
!> Calculate the average and fluctuation value of each variables.
!> @details
!> <u> = ensembel avg. u'  = u - <u>. u'' = u - <rho.u>/<rho>. uave = <u>, uave2 = <rho.u>/<rho>.
!> uflu=u', u2flu=u''.
!> \f$ \overline{u}, \overline{v}, \overline{w}, \overline{p}, \overline{\rho}, \widetilde{u} \f$
!> \f$ \widetilde{v}, \widetilde{w}, u', v', w', p', \rho', \mu, u'', v'', w'' \f$
!> \f$ \overline{u''}, \overline{v''}, \overline{w''} \f$
!>

!> @param[in] n
!> @param[in] idimin i dimension
!> @param[in] jdimin j dimension
!> @param[in] kdimin k dimension
!> @param[in] xin x coordinate
!> @param[in] xin y coordinate
!> @param[in] xin z coordinate
!> @param[in] u stream-wise velocity
!> @param[in] v spanwise velocity
!> @param[in] w wall-normal velocity
!> @param[in] p pressure
!> @param[in] t temperature
!> @param[in] rho density



    Subroutine Cal_Average(n,idimin,jdimin,kdimin,xin,yin,zin,u,v,w,p,t,rho)
      Integer, intent(in) :: n,idimin,jdimin,kdimin
      Real(8),dimension(:,:,:), intent(in) :: xin,yin,zin,u,v,w,p,t,rho
      Integer :: i,j,k

      idim=idimin
      jdim=jdimin
      kdim=kdimin

      If(n==1) then
        Call Init(idim,jdim,kdim)
      End If

      x=xin
      y=yin
      z=zin

      met = planeCalmet(idim,jdim,kdim,planeCalmm(idim,jdim,kdim,x,y,z))

      Do k=1,kdim
        Do i=1,idim
           uave(i,:,k) = ave1D(jdim,u(i,:,k))
           vave(i,:,k) = ave1D(jdim,v(i,:,k))
           wave(i,:,k) = ave1D(jdim,w(i,:,k))
           pave(i,:,k) = ave1D(jdim,p(i,:,k))
           rhoave(i,:,k) = ave1D(jdim,rho(i,:,k))

           uave2(i,:,k) = ave1D(jdim,u(i,:,k)*rho(i,:,k))/rhoave(i,1,k)
           vave2(i,:,k) = ave1D(jdim,v(i,:,k)*rho(i,:,k))/rhoave(i,1,k)
           wave2(i,:,k) = ave1D(jdim,w(i,:,k)*rho(i,:,k))/rhoave(i,1,k)
           Do j=1,jdim
             uflu(i,j,k)=u(i,j,k)-uave(i,j,k)
             vflu(i,j,k)=v(i,j,k)-vave(i,j,k)
             wflu(i,j,k)=w(i,j,k)-wave(i,j,k)
             pflu(i,j,k)=p(i,j,k)-pave(i,j,k)
             rhoflu(i,j,k)=rho(i,j,k)-rhoave(i,j,k)
             mu(i,j,k)=(1.458e-6)*t(i,j,k)**(1.5)/(t(i,j,k)+110.4)

             u2flu(i,j,k)=u(i,j,k)-uave2(i,j,k)
             v2flu(i,j,k)=v(i,j,k)-vave2(i,j,k)
             w2flu(i,j,k)=w(i,j,k)-wave2(i,j,k)
           End Do
           u2fluave(i,:,k)=ave1D(jdim,u2flu(i,:,k))
           v2fluave(i,:,k)=ave1D(jdim,v2flu(i,:,k))
           w2fluave(i,:,k)=ave1D(jdim,w2flu(i,:,k))
         End Do
       End Do

       Do k=1,kdim
          Do i=1,idim
            tau11(i,:,k)=ave1D(jdim,u2flu(i,:,k)*u2flu(i,:,k)*rho(i,:,k))
            tau22(i,:,k)=ave1D(jdim,v2flu(i,:,k)*v2flu(i,:,k)*rho(i,:,k))
            tau33(i,:,k)=ave1D(jdim,w2flu(i,:,k)*w2flu(i,:,k)*rho(i,:,k))
            tau12(i,:,k)=ave1D(jdim,u2flu(i,:,k)*v2flu(i,:,k)*rho(i,:,k))
            tau13(i,:,k)=ave1D(jdim,u2flu(i,:,k)*w2flu(i,:,k)*rho(i,:,k))
            tau23(i,:,k)=ave1D(jdim,v2flu(i,:,k)*w2flu(i,:,k)*rho(i,:,k))
            Do j=1,jdim
              S11(i,j,k) = ddgrid2(i,j,k,x,y,z,u,met,1)
              S22(i,j,k) = ddgrid2(i,j,k,x,y,z,v,met,2)
              S33(i,j,k) = ddgrid2(i,j,k,x,y,z,w,met,3)
              S12(i,j,k) = 0.5*(ddgrid2(i,j,k,x,y,z,u,met,2)+ddgrid2(i,j,k,x,y,z,v,met,1))
              S13(i,j,k) = 0.5*(ddgrid2(i,j,k,x,y,z,u,met,3)+ddgrid2(i,j,k,x,y,z,w,met,1))
              S23(i,j,k) = 0.5*(ddgrid2(i,j,k,x,y,z,v,met,3)+ddgrid2(i,j,k,x,y,z,w,met,2))

              Sig11(i,j,k)=2*mu(i,j,k)*S11(i,j,k)-2.0/3.0*mu(i,j,k)*S11(i,j,k)&
                           -2.0/3.0*mu(i,j,k)*S22(i,j,k)-2.0/3.0*mu(i,j,k)*S33(i,j,k)
              Sig22(i,j,k)=2*mu(i,j,k)*S22(i,j,k)-2.0/3.0*mu(i,j,k)*S11(i,j,k)&
                           -2.0/3.0*mu(i,j,k)*S22(i,j,k)-2.0/3.0*mu(i,j,k)*S33(i,j,k)
              Sig33(i,j,k)=2*mu(i,j,k)*S33(i,j,k)-2.0/3.0*mu(i,j,k)*S11(i,j,k)&
                           -2.0/3.0*mu(i,j,k)*S22(i,j,k)-2.0/3.0*mu(i,j,k)*S33(i,j,k)
              Sig12(i,j,k)=2*mu(i,j,k)*S12(i,j,k)
              Sig13(i,j,k)=2*mu(i,j,k)*S13(i,j,k)
              Sig23(i,j,k)=2*mu(i,j,k)*S23(i,j,k)
            End Do ! End j loop
            Sig11ave(i,:,k)=ave1D(jdim,Sig11(i,:,k))
            Sig22ave(i,:,k)=ave1D(jdim,Sig22(i,:,k))
            Sig33ave(i,:,k)=ave1D(jdim,Sig33(i,:,k))
            Sig12ave(i,:,k)=ave1D(jdim,Sig12(i,:,k))
            Sig13ave(i,:,k)=ave1D(jdim,Sig13(i,:,k))
            Sig23ave(i,:,k)=ave1D(jdim,Sig23(i,:,k))

            Sig11ave2(i,:,k)=ave1D(jdim,Sig11(i,:,k)*rho(i,:,k))/rhoave(i,1,k)
            Sig22ave2(i,:,k)=ave1D(jdim,Sig22(i,:,k)*rho(i,:,k))/rhoave(i,1,k)
            Sig33ave2(i,:,k)=ave1D(jdim,Sig33(i,:,k)*rho(i,:,k))/rhoave(i,1,k)
            Sig12ave2(i,:,k)=ave1D(jdim,Sig12(i,:,k)*rho(i,:,k))/rhoave(i,1,k)
            Sig13ave2(i,:,k)=ave1D(jdim,Sig13(i,:,k)*rho(i,:,k))/rhoave(i,1,k)
            Sig23ave2(i,:,k)=ave1D(jdim,Sig23(i,:,k)*rho(i,:,k))/rhoave(i,1,k)
            Do j=1,jdim
              Sig11flu(i,j,k)=Sig11(i,j,k)-Sig11ave(i,1,k)
              Sig22flu(i,j,k)=Sig22(i,j,k)-Sig22ave(i,1,k)
              Sig33flu(i,j,k)=Sig33(i,j,k)-Sig33ave(i,1,k)
              Sig12flu(i,j,k)=Sig12(i,j,k)-Sig12ave(i,1,k)
              Sig13flu(i,j,k)=Sig13(i,j,k)-Sig13ave(i,1,k)
              Sig23flu(i,j,k)=Sig23(i,j,k)-Sig23ave(i,1,k)

              Sig112flu(i,j,k)=Sig11(i,j,k)-Sig11ave2(i,1,k)
              Sig222flu(i,j,k)=Sig22(i,j,k)-Sig22ave2(i,1,k)
              Sig332flu(i,j,k)=Sig33(i,j,k)-Sig33ave2(i,1,k)
              Sig122flu(i,j,k)=Sig12(i,j,k)-Sig12ave2(i,1,k)
              Sig132flu(i,j,k)=Sig13(i,j,k)-Sig13ave2(i,1,k)
              Sig232flu(i,j,k)=Sig23(i,j,k)-Sig23ave2(i,1,k)
            End Do ! End j loop
          End Do ! End  i loop
        End Do ! End k loop
    End Subroutine Cal_Average






!> @brief
!> Allocate space for each variables.
!> @param[in] idim i dimension
!> @param[in] jdim j dimension
!> @param[in] kdim k dimension

    Subroutine Init(idim,jdim,kdim)
      Integer, intent(in) :: idim,jdim,kdim

      allocate(x(idim,jdim,kdim),y(idim,jdim,kdim),z(idim,jdim,kdim))
      allocate( uave(idim,jdim,kdim) )
      allocate( vave(idim,jdim,kdim) )
      allocate( wave(idim,jdim,kdim) )
      allocate( pave(idim,jdim,kdim) )
      allocate( rhoave(idim,jdim,kdim) )
      allocate( mu(idim,jdim,kdim) )

      allocate( uave2(idim,jdim,kdim) )
      allocate( vave2(idim,jdim,kdim) )
      allocate( wave2(idim,jdim,kdim) )

      allocate( uflu(idim,jdim,kdim) )
      allocate( vflu(idim,jdim,kdim) )
      allocate( wflu(idim,jdim,kdim) )
      allocate( pflu(idim,jdim,kdim) )
      allocate( rhoflu(idim,jdim,kdim) )

      allocate( u2flu(idim,jdim,kdim) )
      allocate( v2flu(idim,jdim,kdim) )
      allocate( w2flu(idim,jdim,kdim) )

      allocate( u2fluave(idim,jdim,kdim) )
      allocate( v2fluave(idim,jdim,kdim) )
      allocate( w2fluave(idim,jdim,kdim) )

      allocate( S11(idim,jdim,kdim) )
      allocate( S22(idim,jdim,kdim) )
      allocate( S33(idim,jdim,kdim) )
      allocate( S12(idim,jdim,kdim) )
      allocate( S13(idim,jdim,kdim) )
      allocate( S23(idim,jdim,kdim) )

      allocate( Sig11(idim,jdim,kdim) )
      allocate( Sig22(idim,jdim,kdim) )
      allocate( Sig33(idim,jdim,kdim) )
      allocate( Sig12(idim,jdim,kdim) )
      allocate( Sig13(idim,jdim,kdim) )
      allocate( Sig23(idim,jdim,kdim) )

      allocate( Sig11ave(idim,jdim,kdim) )
      allocate( Sig22ave(idim,jdim,kdim) )
      allocate( Sig33ave(idim,jdim,kdim) )
      allocate( Sig12ave(idim,jdim,kdim) )
      allocate( Sig13ave(idim,jdim,kdim) )
      allocate( Sig23ave(idim,jdim,kdim) )

      allocate( Sig11ave2(idim,jdim,kdim) )
      allocate( Sig22ave2(idim,jdim,kdim) )
      allocate( Sig33ave2(idim,jdim,kdim) )
      allocate( Sig12ave2(idim,jdim,kdim) )
      allocate( Sig13ave2(idim,jdim,kdim) )
      allocate( Sig23ave2(idim,jdim,kdim) )

      allocate( Sig11flu(idim,jdim,kdim) )
      allocate( Sig22flu(idim,jdim,kdim) )
      allocate( Sig33flu(idim,jdim,kdim) )
      allocate( Sig12flu(idim,jdim,kdim) )
      allocate( Sig13flu(idim,jdim,kdim) )
      allocate( Sig23flu(idim,jdim,kdim) )

      allocate( Sig112flu(idim,jdim,kdim) )
      allocate( Sig222flu(idim,jdim,kdim) )
      allocate( Sig332flu(idim,jdim,kdim) )
      allocate( Sig122flu(idim,jdim,kdim) )
      allocate( Sig132flu(idim,jdim,kdim) )
      allocate( Sig232flu(idim,jdim,kdim) )

      allocate( tau11(idim,jdim,kdim) )
      allocate( tau22(idim,jdim,kdim) )
      allocate( tau33(idim,jdim,kdim) )
      allocate( tau12(idim,jdim,kdim) )
      allocate( tau13(idim,jdim,kdim) )
      allocate( tau23(idim,jdim,kdim) )

      allocate( met(3,3,idim,jdim,kdim) )

    End Subroutine Init

!> @brief
!> Calculate the Production term \f$ P_{ij} \f$.
!> @details
!> \f[
!> P_{ij}=-(\tau_{ik}\frac{\partial{\widetilde{u}_j}}{\partial{x_k}}+\frac{\partial{\widetilde{u}_i}}{\partial{x_k}}\tau_{kj})
!> \f]
!> @param[out] P11 Production term
!> @param[out] P22 Production term
!> @param[out] P33 Production term
!> @param[out] P12 Production term
!> @param[out] P13 Production term
!> @param[out] P23 Production term


    Subroutine Cal_Production(P11,P22,P33,P12,P13,P23)
      Real(8), dimension(idim,jdim,kdim), intent(out) :: P11,P22,P33,P12,P13,P23
      Integer :: i,j,k

      Do k=1,kdim
        Do i=1,idim
          Do j=1,jdim
              P11(i,j,k)=-2*(tau11(i,j,k)*ddgrid2(i,j,k,x,y,z,uave2,met,1)+tau12(i,j,k)*ddgrid2(i,j,k,x,y,z,uave2,met,2)&
                          +tau13(i,j,k)*ddgrid2(i,j,k,x,y,z,uave2,met,3))
              P22(i,j,k)=-2*(tau12(i,j,k)*ddgrid2(i,j,k,x,y,z,vave2,met,1)+tau22(i,j,k)*ddgrid2(i,j,k,x,y,z,vave2,met,2)&
                          +tau23(i,j,k)*ddgrid2(i,j,k,x,y,z,vave2,met,3))
              P33(i,j,k)=-2*(tau13(i,j,k)*ddgrid2(i,j,k,x,y,z,wave2,met,1)+tau23(i,j,k)*ddgrid2(i,j,k,x,y,z,wave2,met,2)&
                          +tau33(i,j,k)*ddgrid2(i,j,k,x,y,z,wave2,met,3))
              P12(i,j,k)=-(tau12(i,j,k)*ddgrid2(i,j,k,x,y,z,uave2,met,1)+tau22(i,j,k)*ddgrid2(i,j,k,x,y,z,uave2,met,2)&
                          +tau23(i,j,k)*ddgrid2(i,j,k,x,y,z,uave2,met,3))-(tau11(i,j,k)*ddgrid2(i,j,k,x,y,z,vave2,met,1)&
                          +tau12(i,j,k)*ddgrid2(i,j,k,x,y,z,vave2,met,2)+tau13(i,j,k)*ddgrid2(i,j,k,x,y,z,vave2,met,3))
              P13(i,j,k)=-(tau13(i,j,k)*ddgrid2(i,j,k,x,y,z,uave2,met,1)+tau23(i,j,k)*ddgrid2(i,j,k,x,y,z,uave2,met,2)&
                          +tau33(i,j,k)*ddgrid2(i,j,k,x,y,z,uave2,met,3))-(tau11(i,j,k)*ddgrid2(i,j,k,x,y,z,wave2,met,1)&
                          +tau12(i,j,k)*ddgrid2(i,j,k,x,y,z,wave2,met,2)+tau13(i,j,k)*ddgrid2(i,j,k,x,y,z,wave2,met,3))
              P23(i,j,k)=-(tau13(i,j,k)*ddgrid2(i,j,k,x,y,z,vave2,met,1)+tau23(i,j,k)*ddgrid2(i,j,k,x,y,z,vave2,met,2)&
                          +tau33(i,j,k)*ddgrid2(i,j,k,x,y,z,vave2,met,3))-(tau12(i,j,k)*ddgrid2(i,j,k,x,y,z,wave2,met,1)&
                          +tau22(i,j,k)*ddgrid2(i,j,k,x,y,z,wave2,met,2)+tau23(i,j,k)*ddgrid2(i,j,k,x,y,z,wave2,met,3))
            End Do ! End j loop
          End Do ! End i loop
        End Do ! End K loop
    End Subroutine Cal_Production


!> @brief
!> Calculate the Turbulent Transport term \f$ T_{ij} \f$.
!> @details
!> \f[
!> T_{ij} = -\frac{\partial}{\partial x_k}(\overline{\rho u_i''u_j''u_k''})
!> \f]

!> @param[in] rho density
!> @param[out] TT11 Turbulent Transport term
!> @param[out] TT22 Turbulent Transport term
!> @param[out] TT33 Turbulent Transport term
!> @param[out] TT12 Turbulent Transport term
!> @param[out] TT13 Turbulent Transport term
!> @param[out] TT23 Turbulent Transport term

   Subroutine Cal_TurTransport(rho,TT11,TT22,TT33,TT12,TT13,TT23)
      Real(8), dimension(idim,jdim,kdim), intent(in) :: rho
      Real(8), dimension(idim,jdim,kdim), intent(out) :: TT11,TT22,TT33,TT12,TT13,TT23
      real(8),dimension(idim,jdim,kdim) :: TT11temp1,TT11temp2,TT11temp3
      real(8),dimension(idim,jdim,kdim) :: TT22temp1,TT22temp2,TT22temp3
      real(8),dimension(idim,jdim,kdim) :: TT33temp1,TT33temp2,TT33temp3
      real(8),dimension(idim,jdim,kdim) :: TT12temp1,TT12temp2,TT12temp3
      real(8),dimension(idim,jdim,kdim) :: TT13temp1,TT13temp2,TT13temp3
      real(8),dimension(idim,jdim,kdim) :: TT23temp1,TT23temp2,TT23temp3
      Integer :: i,j,k

        Do k=1,kdim
          Do i=1,idim

            TT11temp1(i,:,k) = ave1D(jdim,rho(i,:,k)*u2flu(i,:,k)*u2flu(i,:,k)*u2flu(i,:,k))
            TT11temp2(i,:,k) = ave1D(jdim,rho(i,:,k)*u2flu(i,:,k)*u2flu(i,:,k)*v2flu(i,:,k))
            TT11temp3(i,:,k) = ave1D(jdim,rho(i,:,k)*u2flu(i,:,k)*u2flu(i,:,k)*w2flu(i,:,k))
            TT22temp1(i,:,k) = ave1D(jdim,rho(i,:,k)*v2flu(i,:,k)*v2flu(i,:,k)*u2flu(i,:,k))
            TT22temp2(i,:,k) = ave1D(jdim,rho(i,:,k)*v2flu(i,:,k)*v2flu(i,:,k)*v2flu(i,:,k))
            TT22temp3(i,:,k) = ave1D(jdim,rho(i,:,k)*v2flu(i,:,k)*v2flu(i,:,k)*w2flu(i,:,k))
            TT33temp1(i,:,k) = ave1D(jdim,rho(i,:,k)*w2flu(i,:,k)*w2flu(i,:,k)*u2flu(i,:,k))
            TT33temp2(i,:,k) = ave1D(jdim,rho(i,:,k)*w2flu(i,:,k)*w2flu(i,:,k)*v2flu(i,:,k))
            TT33temp3(i,:,k) = ave1D(jdim,rho(i,:,k)*w2flu(i,:,k)*w2flu(i,:,k)*w2flu(i,:,k))
            TT12temp1(i,:,k) = ave1D(jdim,rho(i,:,k)*u2flu(i,:,k)*v2flu(i,:,k)*u2flu(i,:,k))
            TT12temp2(i,:,k) = ave1D(jdim,rho(i,:,k)*u2flu(i,:,k)*v2flu(i,:,k)*v2flu(i,:,k))
            TT12temp3(i,:,k) = ave1D(jdim,rho(i,:,k)*u2flu(i,:,k)*v2flu(i,:,k)*w2flu(i,:,k))
            TT13temp1(i,:,k) = ave1D(jdim,rho(i,:,k)*u2flu(i,:,k)*w2flu(i,:,k)*u2flu(i,:,k))
            TT13temp2(i,:,k) = ave1D(jdim,rho(i,:,k)*u2flu(i,:,k)*w2flu(i,:,k)*v2flu(i,:,k))
            TT13temp3(i,:,k) = ave1D(jdim,rho(i,:,k)*u2flu(i,:,k)*w2flu(i,:,k)*w2flu(i,:,k))
            TT23temp1(i,:,k) = ave1D(jdim,rho(i,:,k)*v2flu(i,:,k)*w2flu(i,:,k)*u2flu(i,:,k))
            TT23temp2(i,:,k) = ave1D(jdim,rho(i,:,k)*v2flu(i,:,k)*w2flu(i,:,k)*v2flu(i,:,k))
            TT23temp3(i,:,k) = ave1D(jdim,rho(i,:,k)*v2flu(i,:,k)*w2flu(i,:,k)*w2flu(i,:,k))
          End Do
        End Do

        Do k=1,kdim
          Do i=1,idim
            Do j=1,jdim

              TT11(i,j,k)=-ddgrid2(i,j,k,x,y,z,TT11temp1,met,1)-ddgrid2(i,j,k,x,y,z,TT11temp2,met,2)&
                          -ddgrid2(i,j,k,x,y,z,TT11temp3,met,3)
              TT22(i,j,k)=-ddgrid2(i,j,k,x,y,z,TT22temp1,met,1)-ddgrid2(i,j,k,x,y,z,TT22temp2,met,2)&
                          -ddgrid2(i,j,k,x,y,z,TT22temp3,met,3)
              TT33(i,j,k)=-ddgrid2(i,j,k,x,y,z,TT33temp1,met,1)-ddgrid2(i,j,k,x,y,z,TT33temp2,met,2)&
                          -ddgrid2(i,j,k,x,y,z,TT33temp3,met,3)
              TT12(i,j,k)=-ddgrid2(i,j,k,x,y,z,TT12temp1,met,1)-ddgrid2(i,j,k,x,y,z,TT12temp2,met,2)&
                          -ddgrid2(i,j,k,x,y,z,TT12temp3,met,3)
              TT13(i,j,k)=-ddgrid2(i,j,k,x,y,z,TT13temp1,met,1)-ddgrid2(i,j,k,x,y,z,TT13temp2,met,2)&
                          -ddgrid2(i,j,k,x,y,z,TT13temp3,met,3)
              TT23(i,j,k)=-ddgrid2(i,j,k,x,y,z,TT23temp1,met,1)-ddgrid2(i,j,k,x,y,z,TT23temp2,met,2)&
                          -ddgrid2(i,j,k,x,y,z,TT23temp3,met,3)
            End Do ! End j loop
          End Do ! End i loop
        End Do ! End k loop
     End Subroutine Cal_TurTransport


!> @brief
!> Calculate Pressure-Dilatation term \f$ \Pi^d_{ij} \f$.
!> @details
!> \f[
!> \Pi_{ij}^d = \overline{p'(\frac{\partial{u_i''}}{\partial{x_j}}+\frac{\partial{u_j''}}{\partial{x_i}})}
!> =\overline{p'(\frac{\partial{u_i'}}{\partial{x_j}}+\frac{\partial{u_j'}}{\partial{x_i}})}
!> \f]

!> @param[out] Pid11 Pressure-Dilatation term
!> @param[out] Pid22 Pressure-Dilatation term
!> @param[out] Pid33 Pressure-Dilatation term
!> @param[out] Pid12 Pressure-Dilatation term
!> @param[out] Pid13 Pressure-Dilatation term
!> @param[out] Pid23 Pressure-Dilatation term

     Subroutine Cal_PreDilatation(Pid11,Pid22,Pid33,Pid12,Pid13,Pid23)
      Real(8), dimension(idim,jdim,kdim), intent(out) :: Pid11,Pid22,Pid33,Pid12,Pid13,Pid23
      real(8) :: Pid11t,Pid22t,Pid33t,Pid12t,Pid13t,Pid23t
      Integer :: i,j,k

        Do k=1,kdim
          Do i=1,idim
            Pid11t=0
            Pid22t=0
            Pid33t=0
            Pid12t=0
            Pid13t=0
            Pid23t=0
            Do j=1,jdim
              Pid11t=Pid11t+2*pflu(i,j,k)*ddgrid2(i,j,k,x,y,z,uflu,met,1)
              Pid22t=Pid22t+2*pflu(i,j,k)*ddgrid2(i,j,k,x,y,z,vflu,met,2)
              Pid33t=Pid33t+2*pflu(i,j,k)*ddgrid2(i,j,k,x,y,z,wflu,met,3)
              Pid12t=Pid12t+pflu(i,j,k)*(ddgrid2(i,j,k,x,y,z,uflu,met,2)+ddgrid2(i,j,k,x,y,z,vflu,met,1))
              Pid13t=Pid13t+pflu(i,j,k)*(ddgrid2(i,j,k,x,y,z,uflu,met,3)+ddgrid2(i,j,k,x,y,z,wflu,met,1))
              Pid23t=Pid23t+pflu(i,j,k)*(ddgrid2(i,j,k,x,y,z,vflu,met,3)+ddgrid2(i,j,k,x,y,z,wflu,met,2))
            End Do ! End j loop
            Pid11(i,:,k)=Pid11t/dble(jdim)
            Pid22(i,:,k)=Pid22t/dble(jdim)
            Pid33(i,:,k)=Pid33t/dble(jdim)
            Pid12(i,:,k)=Pid12t/dble(jdim)
            Pid13(i,:,k)=Pid13t/dble(jdim)
            Pid23(i,:,k)=Pid23t/dble(jdim)

          End Do
        End Do
     End Subroutine Cal_PreDilatation
 

!> @brief
!> Calculate Pressure Transport term \f$ \Pi^t_{ij} \f$.
!> @details
!> \f[
!> \Pi_{ij}^t = -\frac{\partial}{\partial x_k}(\delta_{ik}\overline{p'u_j''}+\overline{p'u_i''}\delta_{jk})
!> \f]

!> @param[out] Pit11 Pressure Transport term.
!> @param[out] Pit22 Pressure Transport term.
!> @param[out] Pit33 Pressure Transport term.
!> @param[out] Pit12 Pressure Transport term.
!> @param[out] Pit13 Pressure Transport term.
!> @param[out] Pit23 Pressure Transport term.



   Subroutine Cal_PreTransport(Pit11,Pit22,Pit33,Pit12,Pit13,Pit23)
      Real(8), dimension(idim,jdim,kdim), intent(out) :: Pit11,Pit22,Pit33,Pit12,Pit13,Pit23
      Real(8), dimension(idim,jdim,kdim) :: Pit11temp1,Pit22temp1,Pit33temp1
      Real(8), dimension(idim,jdim,kdim) :: Pit12temp1,Pit12temp2,Pit13temp1
      Real(8), dimension(idim,jdim,kdim) :: Pit13temp2,Pit23temp1,Pit23temp2
      Integer :: i,j,k

       Do k=1,kdim
          Do i=1,idim
             Pit11temp1(i,:,k) = 2*ave1D(jdim,pflu(i,:,k)*u2flu(i,:,k))
             Pit22temp1(i,:,k) = 2*ave1D(jdim,pflu(i,:,k)*v2flu(i,:,k))
             Pit33temp1(i,:,k) = 2*ave1D(jdim,pflu(i,:,k)*w2flu(i,:,k))
             Pit12temp1(i,:,k) = ave1D(jdim,pflu(i,:,k)*v2flu(i,:,k))
             Pit12temp2(i,:,k) = ave1D(jdim,pflu(i,:,k)*u2flu(i,:,k))
             Pit13temp1(i,:,k) = ave1D(jdim,pflu(i,:,k)*w2flu(i,:,k))
             Pit13temp2(i,:,k) = ave1D(jdim,pflu(i,:,k)*u2flu(i,:,k))
             Pit23temp1(i,:,k) = ave1D(jdim,pflu(i,:,k)*w2flu(i,:,k))
             Pit23temp2(i,:,k) = ave1D(jdim,pflu(i,:,k)*v2flu(i,:,k))

          End Do
        End Do

        Do k=1,kdim
          Do i=1,idim
            Do j=1,jdim
              Pit11(i,j,k) = -ddgrid2(i,j,k,x,y,z,Pit11temp1,met,1)
              Pit22(i,j,k) = -ddgrid2(i,j,k,x,y,z,Pit22temp1,met,2)
              Pit33(i,j,k) = -ddgrid2(i,j,k,x,y,z,Pit33temp1,met,3)
              Pit12(i,j,k) = -ddgrid2(i,j,k,x,y,z,Pit12temp1,met,1)-ddgrid2(i,j,k,x,y,z,Pit12temp2,met,2)
              Pit13(i,j,k) = -ddgrid2(i,j,k,x,y,z,Pit13temp1,met,1)-ddgrid2(i,j,k,x,y,z,Pit13temp2,met,3)
              Pit23(i,j,k) = -ddgrid2(i,j,k,x,y,z,Pit23temp1,met,2)-ddgrid2(i,j,k,x,y,z,Pit23temp2,met,3)
            End Do ! End j loop
          End Do ! End i loop
        End Do ! End k loop
     End Subroutine Cal_PreTransport

!> @brief
!> Calculate Viscous Dissipation term \f$ \epsilon_{ij} \f$.
!> @details
!> \f[
!> \epsilon_{ij}=\overline{\sigma_{ik}'\frac{\partial{u_j''}}{\partial{x_k}}}+\overline{\sigma_{jk}'\frac{\partial{u_i''}}{\partial{x_k}}} 
!> =\overline{\sigma_{ik}'\frac{\partial{u_j'}}{\partial{x_k}}}+\overline{\sigma_{jk}'\frac{\partial{u_i'}}{\partial{x_k}}}
!> \f]

!> @param[out] Delta11 Viscous Dissipation term.
!> @param[out] Delta22 Viscous Dissipation term.
!> @param[out] Delta33 Viscous Dissipation term.
!> @param[out] Delta12 Viscous Dissipation term.
!> @param[out] Delta13 Viscous Dissipation term.
!> @param[out] Delta23 Viscous Dissipation term.



     Subroutine Cal_VisDissipation(Delta11, Delta22, Delta33, Delta12, Delta13, Delta23)
       Real(8), dimension(idim,jdim,kdim), intent(out) :: Delta11, Delta22, Delta33, Delta12, Delta13, Delta23
       Real(8) :: Delta11t1,Delta11t2,Delta11t3,Delta22t1,Delta22t2,Delta22t3,Delta33t1,Delta33t2,Delta33t3
       Real(8) :: Delta12t1,Delta12t2,Delta12t3,Delta12t4,Delta12t5,Delta12t6
       Real(8) :: Delta13t1,Delta13t2,Delta13t3,Delta13t4,Delta13t5,Delta13t6
       Real(8) :: Delta23t1, Delta23t2, Delta23t3,Delta23t4,Delta23t5,Delta23t6
       Integer :: i,j,k

            Delta11t1=0
            Delta11t2=0
            Delta11t3=0
            Delta22t1=0
            Delta22t2=0
            Delta22t3=0
            Delta33t1=0
            Delta33t2=0
            Delta33t3=0
            Delta12t1=0
            Delta12t2=0
            Delta12t3=0
            Delta12t4=0
            Delta12t5=0
            Delta12t6=0
            Delta13t1=0
            Delta13t2=0
            Delta13t3=0
            Delta13t4=0
            Delta13t5=0
            Delta13t6=0
            Delta23t1=0
            Delta23t2=0
            Delta23t3=0
            Delta23t4=0
            Delta23t5=0
            Delta23t6=0

      Do k=1,kdim
        Do i=1,idim
          Do j=1,jdim
            Delta11t1=Delta11t1+Sig11flu(i,j,k)*ddgrid2(i,j,k,x,y,z,uflu,met,1)
            Delta11t2=Delta11t2+Sig12flu(i,j,k)*ddgrid2(i,j,k,x,y,z,uflu,met,2)
            Delta11t3=Delta11t3+Sig13flu(i,j,k)*ddgrid2(i,j,k,x,y,z,uflu,met,3)

            Delta22t1=Delta22t1+Sig12flu(i,j,k)*ddgrid2(i,j,k,x,y,z,vflu,met,1)
            Delta22t2=Delta22t2+Sig22flu(i,j,k)*ddgrid2(i,j,k,x,y,z,vflu,met,2)
            Delta22t3=Delta22t3+Sig23flu(i,j,k)*ddgrid2(i,j,k,x,y,z,vflu,met,3)

            Delta33t1=Delta33t1+Sig13flu(i,j,k)*ddgrid2(i,j,k,x,y,z,wflu,met,1)
            Delta33t2=Delta33t2+Sig23flu(i,j,k)*ddgrid2(i,j,k,x,y,z,wflu,met,2)
            Delta33t3=Delta33t3+Sig33flu(i,j,k)*ddgrid2(i,j,k,x,y,z,wflu,met,3)

            Delta12t1=Delta12t1+Sig12flu(i,j,k)*ddgrid2(i,j,k,x,y,z,uflu,met,1)
            Delta12t2=Delta12t2+Sig22flu(i,j,k)*ddgrid2(i,j,k,x,y,z,uflu,met,2)
            Delta12t3=Delta12t3+Sig23flu(i,j,k)*ddgrid2(i,j,k,x,y,z,uflu,met,3)
            Delta12t4=Delta12t4+Sig11flu(i,j,k)*ddgrid2(i,j,k,x,y,z,vflu,met,1)
            Delta12t5=Delta12t5+Sig12flu(i,j,k)*ddgrid2(i,j,k,x,y,z,vflu,met,2)
            Delta12t6=Delta12t6+Sig13flu(i,j,k)*ddgrid2(i,j,k,x,y,z,vflu,met,3)

            Delta13t1=Delta13t1+Sig13flu(i,j,k)*ddgrid2(i,j,k,x,y,z,uflu,met,1)
            Delta13t2=Delta13t2+Sig23flu(i,j,k)*ddgrid2(i,j,k,x,y,z,uflu,met,2)
            Delta13t3=Delta13t3+Sig33flu(i,j,k)*ddgrid2(i,j,k,x,y,z,uflu,met,3)
            Delta13t4=Delta13t4+Sig11flu(i,j,k)*ddgrid2(i,j,k,x,y,z,wflu,met,1)
            Delta13t5=Delta13t5+Sig12flu(i,j,k)*ddgrid2(i,j,k,x,y,z,wflu,met,2)
            Delta13t6=Delta13t6+Sig13flu(i,j,k)*ddgrid2(i,j,k,x,y,z,wflu,met,3)

            Delta23t1=Delta23t1+Sig13flu(i,j,k)*ddgrid2(i,j,k,x,y,z,vflu,met,1)
            Delta23t2=Delta23t2+Sig23flu(i,j,k)*ddgrid2(i,j,k,x,y,z,vflu,met,2)
            Delta23t3=Delta23t3+Sig33flu(i,j,k)*ddgrid2(i,j,k,x,y,z,vflu,met,3)
            Delta23t4=Delta23t4+Sig12flu(i,j,k)*ddgrid2(i,j,k,x,y,z,wflu,met,1)
            Delta23t5=Delta23t5+Sig22flu(i,j,k)*ddgrid2(i,j,k,x,y,z,wflu,met,2)
            Delta23t6=Delta23t6+Sig23flu(i,j,k)*ddgrid2(i,j,k,x,y,z,wflu,met,3)
          End Do ! End j loop

          Delta11(i,:,k)=2*(Delta11t1/dble(jdim)+Delta11t2/dble(jdim)+Delta11t3/dble(jdim))
          Delta22(i,:,k)=2*(Delta22t1/dble(jdim)+Delta22t2/dble(jdim)+Delta22t3/dble(jdim))
          Delta33(i,:,k)=2*(Delta33t1/dble(jdim)+Delta33t2/dble(jdim)+Delta33t3/dble(jdim))
          Delta12(i,:,k)=Delta12t1/dble(jdim)+Delta12t2/dble(jdim)+Delta12t3/dble(jdim)&
                         +Delta12t4/dble(jdim)+Delta12t5/dble(jdim)+Delta12t6/dble(jdim)
          Delta13(i,:,k)=Delta13t1/dble(jdim)+Delta13t2/dble(jdim)+Delta13t3/dble(jdim)&
                         +Delta13t4/dble(jdim)+Delta13t5/dble(jdim)+Delta13t6/dble(jdim)
          Delta23(i,:,k)=Delta23t1/dble(jdim)+Delta23t2/dble(jdim)+Delta23t3/dble(jdim)&
                         +Delta23t4/dble(jdim)+Delta23t5/dble(jdim)+Delta23t6/dble(jdim)

        End Do
      End Do

     End Subroutine Cal_VisDissipation


!> @brief
!> Calculate Density Fluctuation term \f$ M_{ij} \f$.
!> @details
!> \f[
!> M_{ij}=\overline{u_i''}(\frac{\partial{\overline{\sigma}_{jk}}}{\partial{x_k}}-\frac{\partial{\overline{p}}}{\partial{x_j}})
!> +\overline{u_j''}(\frac{\partial{\overline{\sigma}_{ik}}}{\partial{x_k}}-\frac{\partial{\overline{p}}}{\partial{x_i}})
!> =\frac{\overline{\rho'u_i'}}{\overline{\rho}}(\frac{\partial{\overline{p}}}{\partial{x_j}}
!> -\frac{\partial{\overline{\sigma}_{jk}}}{\partial{x_k}})
!> +\frac{\overline{\rho'u_j'}}{\overline{\rho}}(\frac{\partial{\overline{p}}}{\partial{x_i}}
!> -\frac{\partial{\overline{\sigma}_{ik}}}{\partial{x_k}})
!> \f]

!> @param[out] M11 Density Fluctuation term.
!> @param[out] M22 Density Fluctuation term.
!> @param[out] M33 Density Fluctuation term.
!> @param[out] M12 Density Fluctuation term.
!> @param[out] M13 Density Fluctuation term.
!> @param[out] M23 Density Fluctuation term.



     Subroutine Cal_DenFluctuation(M11,M22,M33,M12,M13,M23)
       Real(8), dimension(idim,jdim,kdim), intent(out) :: M11,M22,M33,M12,M13,M23
       Real(8), dimension(idim,jdim,kdim) :: M11t,M22t,M33t
       Integer :: i,j,k
  
     Do k=1,kdim
          Do i=1,idim
            M11t(i,:,k)=ave1D(jdim,uflu(i,:,k)*rhoflu(i,:,k))/rhoave(i,1,k)
            M22t(i,:,k)=ave1D(jdim,vflu(i,:,k)*rhoflu(i,:,k))/rhoave(i,1,k)
            M33t(i,:,k)=ave1D(jdim,wflu(i,:,k)*rhoflu(i,:,k))/rhoave(i,1,k)
            Do j=1,jdim
              M11(i,j,k)=2*M11t(i,j,k)*(3*ddgrid2(i,j,k,x,y,z,pave,met,1)-ddgrid2(i,j,k,x,y,z,Sig11ave,met,1)&
                          -ddgrid2(i,j,k,x,y,z,Sig12ave,met,2)-ddgrid2(i,j,k,x,y,z,Sig13ave,met,3))
              M22(i,j,k)=2*M22t(i,j,k)*(3*ddgrid2(i,j,k,x,y,z,pave,met,2)-ddgrid2(i,j,k,x,y,z,Sig12ave,met,1)&
                          -ddgrid2(i,j,k,x,y,z,Sig22ave,met,2)-ddgrid2(i,j,k,x,y,z,Sig23ave,met,3))
              M33(i,j,k)=2*M33t(i,j,k)*(3*ddgrid2(i,j,k,x,y,z,pave,met,3)-ddgrid2(i,j,k,x,y,z,Sig13ave,met,1)&
                          -ddgrid2(i,j,k,x,y,z,Sig23ave,met,2)-ddgrid2(i,j,k,x,y,z,Sig33ave,met,3))
              M12(i,j,k)=M11t(i,j,k)*(3*ddgrid2(i,j,k,x,y,z,pave,met,2)-ddgrid2(i,j,k,x,y,z,Sig12ave,met,1)&
                          -ddgrid2(i,j,k,x,y,z,Sig22ave,met,2)-ddgrid2(i,j,k,x,y,z,Sig23ave,met,3))&
                          +M22t(i,j,k)*(3*ddgrid2(i,j,k,x,y,z,pave,met,1)-ddgrid2(i,j,k,x,y,z,Sig11ave,met,1)&
                          -ddgrid2(i,j,k,x,y,z,Sig12ave,met,2)-ddgrid2(i,j,k,x,y,z,Sig13ave,met,3))
              M13(i,j,k)=M11t(i,j,k)*(3*ddgrid2(i,j,k,x,y,z,pave,met,3)-ddgrid2(i,j,k,x,y,z,Sig13ave,met,1)&
                          -ddgrid2(i,j,k,x,y,z,Sig23ave,met,2)-ddgrid2(i,j,k,x,y,z,Sig33ave,met,3))&
                          +M33t(i,j,k)*(3*ddgrid2(i,j,k,x,y,z,pave,met,1)-ddgrid2(i,j,k,x,y,z,Sig11ave,met,1)&
                          -ddgrid2(i,j,k,x,y,z,Sig12ave,met,2)-ddgrid2(i,j,k,x,y,z,Sig13ave,met,3))
              M23(i,j,k)=M22t(i,j,k)*(3*ddgrid2(i,j,k,x,y,z,pave,met,3)-ddgrid2(i,j,k,x,y,z,Sig13ave,met,1)&
                          -ddgrid2(i,j,k,x,y,z,Sig23ave,met,2)-ddgrid2(i,j,k,x,y,z,Sig33ave,met,3))&
                          +M33t(i,j,k)*(3*ddgrid2(i,j,k,x,y,z,pave,met,2)-ddgrid2(i,j,k,x,y,z,Sig12ave,met,1)&
                          -ddgrid2(i,j,k,x,y,z,Sig22ave,met,2)-ddgrid2(i,j,k,x,y,z,Sig23ave,met,3))
            End Do ! End j loop
          End Do ! End i loop
        End Do ! End k loop

     End Subroutine Cal_DenFluctuation


!> @brief
!> Calculate Viscous Diffusion term \f$ D_{ij} \f$.
!> @details
!> \f[
!> D_{ij} = \frac{\partial}{\partial{x_k}}(\overline{(\sigma_{ik}''u_j''+\sigma_{jk}''u_i'')})
!> = \frac{\partial}{\partial{x_k}}\overline{(\sigma_{ik}'u_j'+\sigma_{jk}'u_i')}
!> \f]
!> @param[out] D11 Viscous Diffusion term.
!> @param[out] D22 Viscous Diffusion term.
!> @param[out] D33 Viscous Diffusion term.
!> @param[out] D12 Viscous Diffusion term.
!> @param[out] D13 Viscous Diffusion term.
!> @param[out] D23 Viscous Diffusion term.


     Subroutine Cal_VisDiffusion(D11,D22,D33,D12,D13,D23)
       Real(8), dimension(idim,jdim,kdim), intent(out) :: D11,D22,D33,D12,D13,D23
       Real(8) :: Dt1,Dt2,Dt3,Dt4,Dt5,Dt6,Dt7,Dt8,Dt9
      real(8),dimension(idim,jdim,kdim) :: D11temp1,D11temp2,D11temp3,D22temp1,D22temp2,D22temp3
      real(8),dimension(idim,jdim,kdim) :: D33temp1,D33temp2,D33temp3
      real(8),dimension(idim,jdim,kdim) :: D12temp1,D12temp2,D12temp3,D13temp1,D13temp2,D13temp3
      real(8),dimension(idim,jdim,kdim) :: D23temp1,D23temp2,D23temp3
       Integer :: i,j,k

        Do k=1,kdim
          Do i=1,idim
              D11temp1(i,:,k) = 2*ave1D(jdim,Sig112flu(i,:,k)*u2flu(i,:,k))
              D11temp2(i,:,k) = 2*ave1D(jdim,Sig122flu(i,:,k)*u2flu(i,:,k))
              D11temp3(i,:,k) = 2*ave1D(jdim,Sig132flu(i,:,k)*u2flu(i,:,k))

              D22temp1(i,:,k) = 2*ave1D(jdim,Sig122flu(i,:,k)*v2flu(i,:,k))
              D22temp2(i,:,k) = 2*ave1D(jdim,Sig222flu(i,:,k)*v2flu(i,:,k))
              D22temp3(i,:,k) = 2*ave1D(jdim,Sig232flu(i,:,k)*v2flu(i,:,k))

              D33temp1(i,:,k) = 2*ave1D(jdim,Sig132flu(i,:,k)*w2flu(i,:,k))
              D33temp2(i,:,k) = 2*ave1D(jdim,Sig232flu(i,:,k)*w2flu(i,:,k))
              D33temp3(i,:,k) = 2*ave1D(jdim,Sig332flu(i,:,k)*w2flu(i,:,k))

              D12temp1(i,:,k) = ave1D(jdim,Sig112flu(i,:,k)*v2flu(i,:,k)+Sig122flu(i,:,k)*u2flu(i,:,k))
              D12temp2(i,:,k) = ave1D(jdim,Sig122flu(i,:,k)*v2flu(i,:,k)+Sig222flu(i,:,k)*u2flu(i,:,k))
              D12temp3(i,:,k) = ave1D(jdim,Sig132flu(i,:,k)*v2flu(i,:,k)+Sig232flu(i,:,k)*u2flu(i,:,k))

              D13temp1(i,:,k) = ave1D(jdim,Sig112flu(i,:,k)*w2flu(i,:,k)+Sig132flu(i,:,k)*u2flu(i,:,k))
              D13temp2(i,:,k) = ave1D(jdim,Sig122flu(i,:,k)*w2flu(i,:,k)+Sig232flu(i,:,k)*u2flu(i,:,k))
              D13temp3(i,:,k) = ave1D(jdim,Sig132flu(i,:,k)*w2flu(i,:,k)+Sig332flu(i,:,k)*u2flu(i,:,k))

              D23temp1(i,:,k) = ave1D(jdim,Sig122flu(i,:,k)*w2flu(i,:,k)+Sig132flu(i,:,k)*v2flu(i,:,k))
              D23temp2(i,:,k) = ave1D(jdim,Sig222flu(i,:,k)*w2flu(i,:,k)+Sig232flu(i,:,k)*v2flu(i,:,k))
              D23temp3(i,:,k) = ave1D(jdim,Sig232flu(i,:,k)*w2flu(i,:,k)+Sig332flu(i,:,k)*v2flu(i,:,k))
          End Do
        End Do


        Do k=1,kdim
          Do i=1,idim
            Do j=1,jdim

              D11(i,j,k) = ddgrid2(i,j,k,x,y,z,D11temp1,met,1)+ddgrid2(i,j,k,x,y,z,D11temp2,met,2)&
                         +ddgrid2(i,j,k,x,y,z,D11temp3,met,3)
              D22(i,j,k) = ddgrid2(i,j,k,x,y,z,D22temp1,met,1)+ddgrid2(i,j,k,x,y,z,D22temp2,met,2)&
                         +ddgrid2(i,j,k,x,y,z,D22temp3,met,3)
              D33(i,j,k) = ddgrid2(i,j,k,x,y,z,D33temp1,met,1)+ddgrid2(i,j,k,x,y,z,D33temp2,met,2)&
                         +ddgrid2(i,j,k,x,y,z,D33temp3,met,3)
              D12(i,j,k) = ddgrid2(i,j,k,x,y,z,D12temp1,met,1)+ddgrid2(i,j,k,x,y,z,D12temp2,met,2)&
                         +ddgrid2(i,j,k,x,y,z,D12temp3,met,3)
              D13(i,j,k) = ddgrid2(i,j,k,x,y,z,D13temp1,met,1)+ddgrid2(i,j,k,x,y,z,D13temp2,met,2)&
                         +ddgrid2(i,j,k,x,y,z,D13temp3,met,3)
              D23(i,j,k) = ddgrid2(i,j,k,x,y,z,D23temp1,met,1)+ddgrid2(i,j,k,x,y,z,D23temp2,met,2)&
                         +ddgrid2(i,j,k,x,y,z,D23temp3,met,3)

            End Do ! End j loop
          End Do ! End i loop
        End Do ! End k loop

     End Subroutine Cal_VisDiffusion

   End module MReynoldsBudgets



