!----------------------------------------------------------------------------------------------
! MODULE: MTKEBudgets
!
!> @author
!> Chao Zhang
!
! DESCRIPTION:
!> Calcluate the Turbulent Kinetic Energy Budgets variables.
!
! REVISION HISTORY:
! 08/21/2014 - Initial Version
!----------------------------------------------------------------------------------------------

    module MTKEBudgets
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
              S11(i,j,k) = ddgrid(i,j,k,x,y,z,u,1)
              S22(i,j,k) = ddgrid(i,j,k,x,y,z,v,2)
              S33(i,j,k) = ddgrid(i,j,k,x,y,z,w,3)
              S12(i,j,k) = 0.5*(ddgrid(i,j,k,x,y,z,u,2)+ddgrid(i,j,k,x,y,z,v,1))
              S13(i,j,k) = 0.5*(ddgrid(i,j,k,x,y,z,u,3)+ddgrid(i,j,k,x,y,z,w,1))
              S23(i,j,k) = 0.5*(ddgrid(i,j,k,x,y,z,v,3)+ddgrid(i,j,k,x,y,z,w,2))

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


    End Subroutine Init




!> @brief
!> Calculate the Production term \f$ P_{ij} \f$.
!> @details
!> \f[
!> P=\frac{P_{ii}}{2}=-\tau_{ik}\frac{\partial \widetilde{u_i}}{\partial x_k}
!> = -\tau_{ik}\widetilde{S_{ki}}\
!> \f]
!> @param[out] PP Production term

    
    Subroutine Cal_Production(PP)
      Real(8), dimension(idim,jdim,kdim), intent(out) :: PP
      Integer :: i,j,k

      Do k=1,kdim
        Do i=1,idim
          Do j=1,jdim
              PP(i,j,k)=-tau11(i,j,k)*ddgrid(i,j,k,x,y,z,uave2,1)-tau12(i,j,k)*ddgrid(i,j,k,x,y,z,uave2,2)&
                        -tau13(i,j,k)*ddgrid(i,j,k,x,y,z,uave2,3)-tau12(i,j,k)*ddgrid(i,j,k,x,y,z,vave2,1)&
                        -tau22(i,j,k)*ddgrid(i,j,k,x,y,z,vave2,2)-tau23(i,j,k)*ddgrid(i,j,k,x,y,z,vave2,3)&
                        -tau13(i,j,k)*ddgrid(i,j,k,x,y,z,wave2,1)-tau23(i,j,k)*ddgrid(i,j,k,x,y,z,wave2,2)&
                        -tau33(i,j,k)*ddgrid(i,j,k,x,y,z,wave2,3)

          End Do
        End Do
      End Do
    End Subroutine Cal_Production

!> @brief
!> Calculate the Turbulent Transport term \f$ T_{ij} \f$.
!> @details
!> \f[
!> T = \frac{T_{ii}}{2}
!>   = -\frac{\partial }{\partial x_k}(\frac{1}{2}\overline{\rho}\widetilde{u_i''u_i''u_k''})
!> \f]

!> @param[in] rho density
!> @param[out] TT Turbulent Transport term

   Subroutine Cal_TurTransport(rho,TT)
      Real(8), dimension(idim,jdim,kdim), intent(in) :: rho
      Real(8), dimension(idim,jdim,kdim), intent(out) :: TT
      Real(8) :: TT1,TT2,TT3,TT4,TT5,TT6,TT7,TT8,TT9
      Real(8), dimension(idim,jdim,kdim) :: TTtemp1,TTtemp2,TTtemp3
      Integer :: i,j,k

      Do k=1,kdim
        Do i=1,idim
            TT1 = ave1D(jdim,rho(i,1:jdim,k)*u2flu(i,1:jdim,k)*u2flu(i,1:jdim,k)*u2flu(i,1:jdim,k))/2.0
            TT2 = ave1D(jdim,rho(i,1:jdim,k)*v2flu(i,1:jdim,k)*v2flu(i,1:jdim,k)*u2flu(i,1:jdim,k))/2.0
            TT3 = ave1D(jdim,rho(i,1:jdim,k)*w2flu(i,1:jdim,k)*w2flu(i,1:jdim,k)*u2flu(i,1:jdim,k))/2.0
            TT4 = ave1D(jdim,rho(i,1:jdim,k)*u2flu(i,1:jdim,k)*u2flu(i,1:jdim,k)*v2flu(i,1:jdim,k))/2.0
            TT5 = ave1D(jdim,rho(i,1:jdim,k)*v2flu(i,1:jdim,k)*v2flu(i,1:jdim,k)*v2flu(i,1:jdim,k))/2.0
            TT6 = ave1D(jdim,rho(i,1:jdim,k)*w2flu(i,1:jdim,k)*w2flu(i,1:jdim,k)*v2flu(i,1:jdim,k))/2.0
            TT7 = ave1D(jdim,rho(i,1:jdim,k)*u2flu(i,1:jdim,k)*u2flu(i,1:jdim,k)*w2flu(i,1:jdim,k))/2.0
            TT8 = ave1D(jdim,rho(i,1:jdim,k)*v2flu(i,1:jdim,k)*v2flu(i,1:jdim,k)*w2flu(i,1:jdim,k))/2.0
            TT9 = ave1D(jdim,rho(i,1:jdim,k)*w2flu(i,1:jdim,k)*w2flu(i,1:jdim,k)*w2flu(i,1:jdim,k))/2.0

            TTtemp1(i,:,k) = TT1+TT2+TT3
            TTtemp2(i,:,k) = TT4+TT5+TT6
            TTtemp3(i,:,k) = TT7+TT8+TT9
         End Do
       End Do

       Do k=1,kdim
         Do i=1,idim
           Do j=1,jdim
       TT(i,j,k) = -ddgrid(i,j,k,x,y,z,TTtemp1,1)-ddgrid(i,j,k,x,y,z,TTtemp2,2)-ddgrid(i,j,k,x,y,z,TTtemp3,3)

           End Do
         End Do
       End Do
     End Subroutine Cal_TurTransport

!> @brief
!> Calculate Pressure-Dilatation term \f$ \Pi^d_{ij} \f$.
!> @details
!> \f[
!> \Pi^d = \frac{\Pi_{ii}^d}{2} = \overline{p'\frac{\partial u_i'}{\partial x_i}}
!> \f]

!> @param[out] Pid Pressure-Dilatation term

     Subroutine Cal_PreDilatation(Pid)
      Real(8), dimension(idim,jdim,kdim), intent(out) :: Pid
      Real(8) :: Pid1,Pid2,Pid3
      Integer :: i,j,k

      Do k=1,kdim
        Do i=1,idim
          Pid1=0
          Pid2=0
          Pid3=0
          Do j=1,jdim
            Pid1 = Pid1+pflu(i,j,k)*ddgrid(i,j,k,x,y,z,uflu,1)
            Pid2 = Pid2+pflu(i,j,k)*ddgrid(i,j,k,x,y,z,vflu,2)
            Pid3 = Pid3+pflu(i,j,k)*ddgrid(i,j,k,x,y,z,wflu,3)
          End Do
            Pid(i,:,k)=Pid1/dble(jdim)+Pid2/dble(jdim)+Pid3/dble(jdim)
        End Do
      End Do
     End Subroutine Cal_PreDilatation


!> @brief
!> Calculate Pressure Transport term \f$ \Pi^t_{ij} \f$.
!> @details
!> \f[
!> \Pi^t = \frac{\Pi_{ii}^t}{2} = -\frac{\partial}{\partial x_k} (\overline{p'u_i'}\delta_{ik})
!> \f]

!> @param[out] Pit Pressure Transport term.

    Subroutine Cal_PreTransport(Pit)
      real(8),dimension(idim,jdim,kdim) :: Pittemp1,Pittemp2,Pittemp3
      Real(8), dimension(idim,jdim,kdim), intent(out) :: Pit
      Integer :: i,j,k

       Do k=1,kdim
         Do i=1,idim
           Pittemp1(i,:,k) = ave1D(jdim,pflu(i,:,k)*uflu(i,:,k))
           Pittemp2(i,:,k) = ave1D(jdim,pflu(i,:,k)*vflu(i,:,k))
           Pittemp3(i,:,k) = ave1D(jdim,pflu(i,:,k)*wflu(i,:,k))
         End Do
       End Do

       Do k=1,kdim
         Do i=1,idim
           Do j=1,jdim
        Pit(i,j,k) = -ddgrid(i,j,k,x,y,z,Pittemp1,1)-ddgrid(i,j,k,x,y,z,Pittemp2,2)-ddgrid(i,j,k,x,y,z,Pittemp3,3)
           End Do
         End Do
       End Do
     End Subroutine Cal_PreTransport


!> @brief
!> Calculate Viscous Dissipation term \f$ \epsilon_{ij} \f$.
!> @details
!> \f$ \epsilon=\frac{\epsilon_{ii}}{2}=\overline{\sigma_{ik}'\frac{\partial u_i'}{\partial x_k}} \f$

!> @param[out] Delta Viscous Dissipation term.

     Subroutine Cal_VisDissipation(Delta)
       Real(8), dimension(idim,jdim,kdim), intent(out) :: Delta
       Real(8) :: Deltat1,Deltat2,Deltat3,Deltat4,Deltat5,Deltat6,Deltat7,Deltat8,Deltat9
       Integer :: i,j,k

       Do k=1,kdim
         Do i=1,idim
           Deltat1=0
           Deltat2=0
           Deltat3=0
           Deltat4=0
           Deltat5=0
           Deltat6=0
           Deltat7=0
           Deltat8=0
           Deltat9=0

           Do j=1,jdim
             Deltat1=Deltat1+Sig11flu(i,j,k)*ddgrid(i,j,k,x,y,z,uflu,1)
             Deltat2=Deltat2+Sig12flu(i,j,k)*ddgrid(i,j,k,x,y,z,uflu,2)
             Deltat3=Deltat3+Sig13flu(i,j,k)*ddgrid(i,j,k,x,y,z,uflu,3)
             Deltat4=Deltat4+Sig12flu(i,j,k)*ddgrid(i,j,k,x,y,z,vflu,1)
             Deltat5=Deltat5+Sig22flu(i,j,k)*ddgrid(i,j,k,x,y,z,vflu,2)
             Deltat6=Deltat6+Sig23flu(i,j,k)*ddgrid(i,j,k,x,y,z,vflu,3)
             Deltat7=Deltat7+Sig13flu(i,j,k)*ddgrid(i,j,k,x,y,z,wflu,1)
             Deltat8=Deltat8+Sig23flu(i,j,k)*ddgrid(i,j,k,x,y,z,wflu,2)
             Deltat9=Deltat9+Sig33flu(i,j,k)*ddgrid(i,j,k,x,y,z,wflu,3)
           End Do
           Delta(i,:,k)=Deltat1/dble(jdim)+Deltat2/dble(jdim)+Deltat3/dble(jdim)+Deltat4/dble(jdim)&
                        +Deltat5/dble(jdim)+Deltat6/dble(jdim)+Deltat7/dble(jdim)+Deltat8/dble(jdim)&
                        +Deltat9/dble(jdim)
         End Do
       End Do
     End Subroutine Cal_VisDissipation




!> @brief
!> Calculate Density Fluctuation term \f$ M_{ij} \f$.
!> @details
!> \f[
!> M=\frac{M_{ii}}{2}
!> =\overline{u_i''}(\frac{\partial \overline{\sigma}_{ik}}{\partial x_k}
!> -\frac{\partial \overline{p}}{\partial x_i})
!> \f]

!> @param[out] M Density Fluctuation term.

     Subroutine Cal_DenFluctuation(M)
       Real(8), dimension(idim,jdim,kdim), intent(out) :: M
       Integer :: i,j,k

       Do k=1,kdim
         Do i=1,idim
           Do j=1,jdim
              M(i,j,k)=u2fluave(i,j,k)*(ddgrid(i,j,k,x,y,z,Sig11ave,1)-ddgrid(i,j,k,x,y,z,pave,1))&
                      +u2fluave(i,j,k)*(ddgrid(i,j,k,x,y,z,Sig12ave,2)-ddgrid(i,j,k,x,y,z,pave,1))&
                      +u2fluave(i,j,k)*(ddgrid(i,j,k,x,y,z,Sig13ave,3)-ddgrid(i,j,k,x,y,z,pave,1))&
                      +v2fluave(i,j,k)*(ddgrid(i,j,k,x,y,z,Sig12ave,1)-ddgrid(i,j,k,x,y,z,pave,2))&
                      +v2fluave(i,j,k)*(ddgrid(i,j,k,x,y,z,Sig22ave,2)-ddgrid(i,j,k,x,y,z,pave,2))&
                      +v2fluave(i,j,k)*(ddgrid(i,j,k,x,y,z,Sig23ave,3)-ddgrid(i,j,k,x,y,z,pave,2))&
                      +w2fluave(i,j,k)*(ddgrid(i,j,k,x,y,z,Sig13ave,1)-ddgrid(i,j,k,x,y,z,pave,3))&
                      +w2fluave(i,j,k)*(ddgrid(i,j,k,x,y,z,Sig23ave,2)-ddgrid(i,j,k,x,y,z,pave,3))&
                      +w2fluave(i,j,k)*(ddgrid(i,j,k,x,y,z,Sig33ave,3)-ddgrid(i,j,k,x,y,z,pave,3))
           End Do
         End Do
       End Do
     End Subroutine Cal_DenFluctuation


!> @brief
!> Calculate Viscous Diffusion term \f$ D_{ij} \f$.
!> @details
!> \f$ D=\frac{D_{ii}}{2} = \frac{\partial}{\partial x_k}(\overline{\sigma_{ik}'u_i'}) \f$
!> @param[out] Viscous Diffusion term.

     Subroutine Cal_VisDiffusion(D)
       Real(8), dimension(idim,jdim,kdim), intent(out) :: D
       Real(8) :: Dt1,Dt2,Dt3,Dt4,Dt5,Dt6,Dt7,Dt8,Dt9
       Real(8),dimension(idim,jdim,kdim) :: Dtemp1,Dtemp2,Dtemp3
       Integer :: i,j,k

       Do k=1,kdim
         Do i=1,idim
           Dt1=ave1D(jdim,Sig11flu(i,:,k)*uflu(i,:,k))
           Dt2=ave1D(jdim,Sig12flu(i,:,k)*vflu(i,:,k))
           Dt3=ave1D(jdim,Sig13flu(i,:,k)*wflu(i,:,k))
           Dt4=ave1D(jdim,Sig12flu(i,:,k)*uflu(i,:,k))
           Dt5=ave1D(jdim,Sig22flu(i,:,k)*vflu(i,:,k))
           Dt6=ave1D(jdim,Sig23flu(i,:,k)*wflu(i,:,k))
           Dt7=ave1D(jdim,Sig13flu(i,:,k)*uflu(i,:,k))
           Dt8=ave1D(jdim,Sig23flu(i,:,k)*vflu(i,:,k))
           Dt9=ave1D(jdim,Sig33flu(i,:,k)*wflu(i,:,k))

           Dtemp1(i,:,k) = Dt1+Dt2+Dt3
           Dtemp2(i,:,k) = Dt4+Dt5+Dt6
           Dtemp3(i,:,k) = Dt7+Dt8+Dt9
         End Do
       End Do

       Do k=1,kdim
         Do i=1,idim
           Do j=1,jdim
        D(i,j,k) = ddgrid(i,j,k,x,y,z,Dtemp1,1)+ddgrid(i,j,k,x,y,z,Dtemp2,2)+ddgrid(i,j,k,x,y,z,Dtemp3,3)
           End Do
         End Do
       End Do
     End Subroutine Cal_VisDiffusion


  End module MTKEBudgets
