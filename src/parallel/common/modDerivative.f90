      module MDerivative
      use modMetrics
      implicit none
      real(8), parameter, private :: c12i = 1.d0/12.d0
contains
      pure function dfdi(idim,f)
      integer, intent(in) :: idim
      real(8), intent(in) :: f(idim)
      real(8) :: dfdi(idim)
      integer :: i

      if(idim.eq.1) then
         dfdi(1) = 1.d0
      elseif(idim.eq.2) then
         dfdi(1) = f(2)-f(1) 
         dfdi(2) = f(2)-f(1) 
      elseif(idim.le.4) then
         dfdi(1) = -1.5*f(1) + 2*f(2) - 0.5*f(3)
         forall(i=2:idim-1)
            dfdi(i) = -0.5*f(i-1) + 0.5*f(i+1)
         end forall
         dfdi(idim) = 0.5*f(idim-2) - 2.0*f(idim-1) + 1.5*f(idim)
      else 
         dfdi(1)=-25./12.*f(1)+4.*f(2)-3.*f(3)+4./3.*f(4)-0.25*f(5)
         dfdi(2)=-0.25*f(1)-5./6.*f(2)+1.5*f(3)-0.5*f(4)+1./12.*f(5)
         forall(i=3:idim-2)
           dfdi(i) = ( f(i-2)-8.d0*(f(i-1)-f(i+1))-f(i+2) )/12.d0
         end forall
         dfdi(idim-1)=-1./12.*f(idim-4)+0.5*f(idim-3)-1.5*f(idim-2) &
                          + 5./6.*f(idim-1)+0.25*f(idim)
         dfdi(idim)=0.25*f(idim-4)-4./3.*f(idim-3)+3.*f(idim-2) &
                         - 4.*f(idim-1)+25./12.*f(idim)
      endif 

      end function dfdi

      function dfdipt(iin,idim,f)
! INPUT
!    iin: point of derivative evaluation
!    idim: dimension of input array
!    f: input array with size idim
! OUTPUT
!    dfdipt: deriviative evaluated at point iin

      integer, intent(in) :: iin, idim
      real(8), intent(in) :: f(idim)
      real(8) :: dfdipt
      integer :: i

      if(iin.gt.idim) then
         print *, 'iin (evaluation point) should be less than or equal to idim (size of the array). STOP ...'
         stop
      endif
      if(idim.eq.1) then
         dfdipt = 1.d0
      elseif(idim.eq.2) then
         dfdipt = f(2)-f(1) 
      elseif(idim.le.4) then
         if(iin.eq.1) dfdipt = -1.5*f(1) + 2*f(2) - 0.5*f(3)
         if(iin.gt.1.and.iin.lt.idim) dfdipt = -0.5*f(i-1) + 0.5*f(i+1)
         if(iin.eq.idim) dfdipt = 0.5*f(idim-2) - 2.0*f(idim-1) + 1.5*f(idim)
      else 
         if(iin.eq.1) dfdipt =-25./12.*f(1)+4.*f(2)-3.*f(3)+4./3.*f(4)-0.25*f(5)
         if(iin.eq.2) dfdipt =-0.25*f(1)-5./6.*f(2)+1.5*f(3)-0.5*f(4)+1./12.*f(5)
         if(iin.gt.2.and.iin.lt.idim-1) dfdipt = ( f(i-2)-8.d0*(f(i-1)-f(i+1))-f(i+2) )/12.d0
         if(iin.eq.idim-1) dfdipt=-1./12.*f(idim-4)+0.5*f(idim-3)-1.5*f(idim-2) &
                                         + 5./6.*f(idim-1)+0.25*f(idim)
         if(iin.eq.idim) dfdipt=0.25*f(idim-4)-4./3.*f(idim-3)+3.*f(idim-2) &
                                     - 4.*f(idim-1)+25./12.*f(idim)
      endif 

      end function dfdipt

    real(8) function ddgrid(iin,jin,kin,xin,yin,zin,vin,n)
      integer, intent(in) :: iin,jin,kin,n
      real(8), intent(in) :: xin(:,:,:), yin(:,:,:), zin(:,:,:)
      real(8), intent(in) :: vin(:,:,:)
      integer :: idim, jdim, kdim
      real(8) :: ddi, ddj, ddk
      real(8) :: met(3,3)

      idim = size(xin,dim=1)
      jdim = size(xin,dim=2)
      kdim = size(xin,dim=3)
      if(size(xin).ne.size(yin).or.size(xin).ne.size(zin).or.size(xin).ne.size(vin)) then
          print *, 'The size for input arrays (xin,yin,zin,vin) should be the same. STOP'
          stop
      endif

      if(idim.ge.5) then
        if(iin.ge.3.and.iin.le.idim-2) then
           ddi = (vin(iin-2,jin,kin)-8.*(vin(iin-1,jin,kin)-vin(iin+1,jin,kin))-vin(iin+2,jin,kin))*c12i
        elseif(iin.eq.1) then
           ddi = -25./12.*vin(1,jin,kin)+4.*vin(2,jin,kin)-3.*vin(3,jin,kin)+4./3.*vin(4,jin,kin)-0.25*vin(5,jin,kin)
        elseif(iin.eq.2) then
           ddi = -0.25*vin(1,jin,kin)-5./6.*vin(2,jin,kin)+1.5*vin(3,jin,kin)-0.5*vin(4,jin,kin)+1./12.*vin(5,jin,kin)
        elseif(iin.eq.idim-1) then
           ddi = -1./12.*vin(idim-4,jin,kin)+0.5*vin(idim-3,jin,kin)-1.5*vin(idim-2,jin,kin)&
                          + 5./6.*vin(idim-1,jin,kin)+0.25*vin(idim,jin,kin)
        elseif(iin.eq.idim) then
           ddi = 0.25*vin(idim-4,jin,kin)-4./3.*vin(idim-3,jin,kin)+3.*vin(idim-2,jin,kin)&
                         - 4.*vin(idim-1,jin,kin)+25./12.*vin(idim,jin,kin)
        endif ! end if(iin)
     elseif(idim.eq.3.or.idim.eq.4) then
        if(iin.ge.2.and.iin.le.3) then
           ddi = -0.5*vin(iin-1,jin,kin) + 0.5*vin(iin+1,jin,kin)
        elseif(iin.eq.1) then
           ddi = -1.5*vin(1,jin,kin) + 2.0*vin(2,jin,kin) - 0.5*vin(3,jin,kin)
        elseif(iin.eq.idim) then
           ddi = 0.5*vin(idim-2,jin,kin) - 2.0*vin(idim-1,jin,kin) + 1.5*vin(idim,jin,kin)
       endif
     elseif(idim.eq.2) then
           ddi = vin(idim,jin,kin) - vin(idim-1,jin,kin)
     elseif(idim.eq.1) then
           ddi = 0.0
     endif

    if(jdim.ge.5) then
       if(jin.ge.3.and.jin.le.jdim-2) then
           ddj = (vin(iin,jin-2,kin)-8.*(vin(iin,jin-1,kin)-vin(iin,jin+1,kin))-vin(iin,jin+2,kin))*c12i
       elseif(jin.eq.1) then
           ddj = -25./12.*vin(iin,1,kin)+4.*vin(iin,2,kin)-3.*vin(iin,3,kin)+4./3.*vin(iin,4,kin)-0.25*vin(iin,5,kin)
       elseif(jin.eq.2) then
           ddj = -0.25*vin(iin,1,kin)-5./6.*vin(iin,2,kin)+1.5*vin(iin,3,kin)-0.5*vin(iin,4,kin)+1./12.*vin(iin,5,kin)
       elseif(jin.eq.jdim-1) then
           ddj = -1./12.*vin(iin,jdim-4,kin)+0.5*vin(iin,jdim-3,kin)-1.5*vin(iin,jdim-2,kin)&
                          + 5./6.*vin(iin,jdim-1,kin)+0.25*vin(iin,jdim,kin)
       elseif(jin.eq.jdim) then
           ddj = 0.25*vin(iin,jdim-4,kin)-4./3.*vin(iin,jdim-3,kin)+3.*vin(iin,jdim-2,kin)&
                        - 4.*vin(iin,jdim-1,kin)+25./12.*vin(iin,jdim,kin)
       endif
    elseif(jdim.eq.3.or.jdim.eq.4) then
       if(jin.ge.2.and.jin.le.3) then
           ddj = -0.5*vin(iin,jin-1,kin) + 0.5*vin(iin,jin+1,kin)
       elseif(jin.eq.1) then
           ddj = -1.5*vin(iin,1,kin) + 2.0*vin(iin,2,kin) - 0.5*vin(iin,3,kin)
       elseif(jin.eq.jdim) then
           ddj = 0.5*vin(iin,jdim-2,kin) - 2.0*vin(iin,jdim-1,kin) + 1.5*vin(iin,jdim,kin)
       endif
    elseif(jdim.eq.2) then
           ddj = vin(iin,jdim,kin) - vin(iin,jdim-1,kin)
    elseif(jdim.eq.1) then
           ddj = 0.0
    endif

    if(kdim.ge.5) then
       if(kin.ge.3.and.kin.le.kdim-2) then
           ddk = (vin(iin,jin,kin-2)-8.*(vin(iin,jin,kin-1)-vin(iin,jin,kin+1))-vin(iin,jin,kin+2))*c12i
       elseif(kin.eq.1) then
           ddk = -25./12.*vin(iin,jin,1)+4.*vin(iin,jin,2)-3.*vin(iin,jin,3)+4./3.*vin(iin,jin,4)-0.25*vin(iin,jin,5)
       elseif(kin.eq.2) then
           ddk = -0.25*vin(iin,jin,1)-5./6.*vin(iin,jin,2)+1.5*vin(iin,jin,3)-0.5*vin(iin,jin,4)+1./12.*vin(iin,jin,5)
       elseif(kin.eq.kdim-1) then
           ddk = -1./12.*vin(iin,jin,kdim-4)+0.5*vin(iin,jin,kdim-3)-1.5*vin(iin,jin,kdim-2)&
                          + 5./6.*vin(iin,jin,kdim-1)+0.25*vin(iin,jin,kdim)
       elseif(kin.eq.kdim) then
           ddk = 0.25*vin(iin,jin,kdim-4)-4./3.*vin(iin,jin,kdim-3)+3.*vin(iin,jin,kdim-2)&
                        - 4.*vin(iin,jin,kdim-1)+25./12.*vin(iin,jin,kdim)
       endif
    elseif(kdim.eq.3.or.kdim.eq.4) then
       if(kin.ge.2.and.kin.le.3) then
           ddk = -0.5*vin(iin,jin,kin-1) + 0.5*vin(iin,jin,kin+1)
       elseif(kin.eq.1) then
           ddk = -1.5*vin(iin,jin,1) + 2.0*vin(iin,jin,2) - 0.5*vin(iin,jin,3)
       elseif(kin.eq.kdim) then
           ddk = 0.5*vin(iin,jin,kdim-2) - 2.0*vin(iin,jin,kdim-1) + 1.5*vin(iin,jin,kdim)
       endif
    elseif(kdim.eq.2) then
           ddk = vin(iin,jin,kdim) - vin(iin,jin,kdim-1)
    elseif(kdim.eq.1) then
           ddk = 0.0
    endif

 !   met = Calmet(Calmm(iin,jin,kin,xin,yin,zin))
    select case(n)
      case(1) !ddx
        ddgrid = ddi*met(1,1)+ddj*met(2,1)+ddk*met(3,1)
      case(2) !ddy
        ddgrid = ddi*met(1,2)+ddj*met(2,2)+ddk*met(3,2)
      case(3) !ddz
        ddgrid = ddi*met(1,3)+ddj*met(2,3)+ddk*met(3,3)
      case default
        print*,'n parameter out of range'
        stop
      end select
    end function ddgrid





    real(8) function ddgrid_P(iin,jin,kin,idim,jdim,kdim,xin,yin,zin,vin,n)
      integer, intent(in) :: iin,jin,kin,n
      integer, intent(in) :: idim, jdim, kdim
      real(8), intent(in) :: xin(1:kdim,-1:idim-2,-1:jdim-2), yin(1:kdim,-1:idim-2,-1:jdim-2), zin(1:kdim,-1:idim-2,-1:jdim-2)
      real(8), intent(in) :: vin(1:kdim,-1:idim-2,-1:jdim-2)
      real(8) :: ddi, ddj, ddk
      real(8) :: met(3,3)


      if(size(xin).ne.size(yin).or.size(xin).ne.size(zin).or.size(xin).ne.size(vin)) then
          print *, 'The size for input arrays (xin,yin,zin,vin) should be the same. STOP'
          stop
      endif

      ddi = (vin(kin,iin-2,jin)-8.*(vin(kin,iin-1,jin)-vin(kin,iin+1,jin))-vin(kin,iin+2,jin))*c12i
      ddj = (vin(kin,iin,jin-2)-8.*(vin(kin,iin,jin-1)-vin(kin,iin,jin+1))-vin(kin,iin,jin+2))*c12i


       if(kin.ge.3.and.kin.le.kdim-2) then
           ddk = (vin(kin-2,iin,jin)-8.*(vin(kin-1,iin,jin)-vin(kin+1,iin,jin))-vin(kin+2,iin,jin))*c12i
       elseif(kin.eq.1) then
           ddk = -25./12.*vin(1,iin,jin)+4.*vin(2,iin,jin)-3.*vin(3,iin,jin)+4./3.*vin(4,iin,jin)-0.25*vin(5,iin,jin)
       elseif(kin.eq.2) then
           ddk = -0.25*vin(1,iin,jin)-5./6.*vin(2,iin,jin)+1.5*vin(3,iin,jin)-0.5*vin(4,iin,jin)+1./12.*vin(5,iin,jin)
       elseif(kin.eq.kdim-1) then
           ddk = -1./12.*vin(kdim-4,iin,jin)+0.5*vin(kdim-3,iin,jin)-1.5*vin(kdim-2,iin,jin)&
                          + 5./6.*vin(kdim-1,iin,jin)+0.25*vin(kdim,iin,jin)
       elseif(kin.eq.kdim) then
           ddk = 0.25*vin(kdim-4,iin,jin)-4./3.*vin(kdim-3,iin,jin)+3.*vin(kdim-2,iin,jin)&
                        - 4.*vin(kdim-1,iin,jin)+25./12.*vin(kdim,iin,jin)
       endif


    met = Calmet_P(Calmm_P(iin,jin,kin,idim,jdim,kdim,xin,yin,zin))
    select case(n)
      case(1) !ddx
        ddgrid_P = ddi*met(1,1)+ddj*met(2,1)+ddk*met(3,1)
      case(2) !ddy
        ddgrid_P = ddi*met(1,2)+ddj*met(2,2)+ddk*met(3,2)
      case(3) !ddz
       ddgrid_P = ddi*met(1,3)+ddj*met(2,3)+ddk*met(3,3)
      case default
        print*,'n parameter out of range'
        stop
      end select


    end function ddgrid_P


















    real(8) function ddgrid_tsiplane(iin,jin,kin,jdim,kdim,buffer_grid,vin,dvindi,n)
      integer, intent(in) :: iin,jin,kin,n
      real(8), intent(in) :: buffer_grid(:)
      real(8), intent(in) :: vin(:,:)
      real(8), intent(in) :: dvindi
      integer, intent(in) :: jdim, kdim
!      real(8) :: ddi, ddj, ddk
      real(8) :: ddj, ddk
      real(8) :: met(3,3)

!      idim = size(xin,dim=1)
!      jdim = size(xin,dim=2)
!      kdim = size(xin,dim=3)
!      if(size(xin).ne.size(yin).or.size(xin).ne.size(zin).or.size(xin).ne.size(vin)) then
!          print *, 'The size for input arrays (xin,yin,zin,vin) should be the same. STOP'
!          stop
!      endif

    if(jdim.ge.5) then
       if(jin.ge.3.and.jin.le.jdim-2) then
           ddj = (vin(jin-2,kin)-8.*(vin(jin-1,kin)-vin(jin+1,kin))-vin(jin+2,kin))*c12i
       elseif(jin.eq.1) then
           ddj = -25./12.*vin(1,kin)+4.*vin(2,kin)-3.*vin(3,kin)+4./3.*vin(4,kin)-0.25*vin(5,kin)
       elseif(jin.eq.2) then
           ddj = -0.25*vin(1,kin)-5./6.*vin(2,kin)+1.5*vin(3,kin)-0.5*vin(4,kin)+1./12.*vin(5,kin)
       elseif(jin.eq.jdim-1) then
           ddj = -1./12.*vin(jdim-4,kin)+0.5*vin(jdim-3,kin)-1.5*vin(jdim-2,kin)&
                          + 5./6.*vin(jdim-1,kin)+0.25*vin(jdim,kin)
       elseif(jin.eq.jdim) then
           ddj = 0.25*vin(jdim-4,kin)-4./3.*vin(jdim-3,kin)+3.*vin(jdim-2,kin)&
                        - 4.*vin(jdim-1,kin)+25./12.*vin(jdim,kin)
       endif
    elseif(jdim.eq.3.or.jdim.eq.4) then
       if(jin.ge.2.and.jin.le.3) then
           ddj = -0.5*vin(jin-1,kin) + 0.5*vin(jin+1,kin)
       elseif(jin.eq.1) then
           ddj = -1.5*vin(1,kin) + 2.0*vin(2,kin) - 0.5*vin(3,kin)
       elseif(jin.eq.jdim) then
           ddj = 0.5*vin(jdim-2,kin) - 2.0*vin(jdim-1,kin) + 1.5*vin(jdim,kin)
       endif
    elseif(jdim.eq.2) then
           ddj = vin(jdim,kin) - vin(jdim-1,kin)
    elseif(jdim.eq.1) then
           ddj = 0.0
    endif

    if(kdim.ge.5) then
       if(kin.ge.3.and.kin.le.kdim-2) then
           ddk = (vin(jin,kin-2)-8.*(vin(jin,kin-1)-vin(jin,kin+1))-vin(jin,kin+2))*c12i
       elseif(kin.eq.1) then
           ddk = -25./12.*vin(jin,1)+4.*vin(jin,2)-3.*vin(jin,3)+4./3.*vin(jin,4)-0.25*vin(jin,5)
       elseif(kin.eq.2) then
           ddk = -0.25*vin(jin,1)-5./6.*vin(jin,2)+1.5*vin(jin,3)-0.5*vin(jin,4)+1./12.*vin(jin,5)
       elseif(kin.eq.kdim-1) then
           ddk = -1./12.*vin(jin,kdim-4)+0.5*vin(jin,kdim-3)-1.5*vin(jin,kdim-2)&
                          + 5./6.*vin(jin,kdim-1)+0.25*vin(jin,kdim)
       elseif(kin.eq.kdim) then
           ddk = 0.25*vin(jin,kdim-4)-4./3.*vin(jin,kdim-3)+3.*vin(jin,kdim-2)&
                        - 4.*vin(jin,kdim-1)+25./12.*vin(jin,kdim)
       endif
    elseif(kdim.eq.3.or.kdim.eq.4) then
       if(kin.ge.2.and.kin.le.3) then
           ddk = -0.5*vin(jin,kin-1) + 0.5*vin(jin,kin+1)
       elseif(kin.eq.1) then
           ddk = -1.5*vin(jin,1) + 2.0*vin(jin,2) - 0.5*vin(jin,3)
       elseif(kin.eq.kdim) then
           ddk = 0.5*vin(jin,kdim-2) - 2.0*vin(jin,kdim-1) + 1.5*vin(jin,kdim)
       endif
    elseif(kdim.eq.2) then
           ddk = vin(jin,kdim) - vin(jin,kdim-1)
    elseif(kdim.eq.1) then
           ddk = 0.0
    endif

!    met = Calmet(Calmm(iin,jin,kin,xin,yin,zin))
    select case(n)
      case(1) !ddx
        ddgrid_tsiplane = dvindi*buffer_grid(4)+ddj*buffer_grid(5)+ddk*buffer_grid(6)
      case(2) !ddy
        ddgrid_tsiplane = dvindi*buffer_grid(7)+ddj*buffer_grid(8)+ddk*buffer_grid(9)
      case(3) !ddz
        ddgrid_tsiplane = dvindi*buffer_grid(10)+ddj*buffer_grid(11)+ddk*buffer_grid(12)
      case default
        print*,'n parameter out of range'
        stop
      end select
    end function ddgrid_tsiplane










!    real(8) function ddgrid_tsjplane(iin,jin,kin,xin,yin,zin,didx,djdx,dkdx,didy,djdy,dkdy,didz,djdz,dkdz,vin,dvindj,n)
    real(8) function ddgrid_tsjplane(iin,jin,kin,idim,kdim,buffer_grid,vin,dvindj,n)
      integer, intent(in) :: iin,jin,kin,n
      real(8), intent(in) :: buffer_grid(:)
      real(8), intent(in) :: vin(:,:)
      real(8), intent(in) :: dvindj
      integer, intent(in) :: idim, kdim
!      real(8) :: ddi, ddj, ddk
      real(8) :: ddi, ddk
      real(8) :: met(3,3)

!      idim = size(xin,dim=1)
!      jdim = size(xin,dim=2)
!      kdim = size(xin,dim=3)
!      if(size(xin).ne.size(yin).or.size(xin).ne.size(zin).or.size(xin).ne.size(vin)) then
!          print *, 'The size for input arrays (xin,yin,zin,vin) should be the same. STOP'
!          stop
!      endif

      if(idim.ge.5) then
        if(iin.ge.3.and.iin.le.idim-2) then
           ddi = (vin(iin-2,kin)-8.*(vin(iin-1,kin)-vin(iin+1,kin))-vin(iin+2,kin))*c12i
        elseif(iin.eq.1) then
           ddi = -25./12.*vin(1,kin)+4.*vin(2,kin)-3.*vin(3,kin)+4./3.*vin(4,kin)-0.25*vin(5,kin)
        elseif(iin.eq.2) then
           ddi = -0.25*vin(1,kin)-5./6.*vin(2,kin)+1.5*vin(3,kin)-0.5*vin(4,kin)+1./12.*vin(5,kin)
        elseif(iin.eq.idim-1) then
           ddi = -1./12.*vin(idim-4,kin)+0.5*vin(idim-3,kin)-1.5*vin(idim-2,kin)&
                          + 5./6.*vin(idim-1,kin)+0.25*vin(idim,kin)
        elseif(iin.eq.idim) then
           ddi = 0.25*vin(idim-4,kin)-4./3.*vin(idim-3,kin)+3.*vin(idim-2,kin)&
                         - 4.*vin(idim-1,kin)+25./12.*vin(idim,kin)
        endif ! end if(iin)
     elseif(idim.eq.3.or.idim.eq.4) then
        if(iin.ge.2.and.iin.le.3) then
           ddi = -0.5*vin(iin-1,kin) + 0.5*vin(iin+1,kin)
        elseif(iin.eq.1) then
           ddi = -1.5*vin(1,kin) + 2.0*vin(2,kin) - 0.5*vin(3,kin)
        elseif(iin.eq.idim) then
           ddi = 0.5*vin(idim-2,kin) - 2.0*vin(idim-1,kin) + 1.5*vin(idim,kin)
       endif
     elseif(idim.eq.2) then
           ddi = vin(idim,kin) - vin(idim-1,kin)
     elseif(idim.eq.1) then
           ddi = 0.0
     endif

    if(kdim.ge.5) then
       if(kin.ge.3.and.kin.le.kdim-2) then
           ddk = (vin(iin,kin-2)-8.*(vin(iin,kin-1)-vin(iin,kin+1))-vin(iin,kin+2))*c12i
       elseif(kin.eq.1) then
           ddk = -25./12.*vin(iin,1)+4.*vin(iin,2)-3.*vin(iin,3)+4./3.*vin(iin,4)-0.25*vin(iin,5)
       elseif(kin.eq.2) then
           ddk = -0.25*vin(iin,1)-5./6.*vin(iin,2)+1.5*vin(iin,3)-0.5*vin(iin,4)+1./12.*vin(iin,5)
       elseif(kin.eq.kdim-1) then
           ddk = -1./12.*vin(iin,kdim-4)+0.5*vin(iin,kdim-3)-1.5*vin(iin,kdim-2)&
                          + 5./6.*vin(iin,kdim-1)+0.25*vin(iin,kdim)
       elseif(kin.eq.kdim) then
           ddk = 0.25*vin(iin,kdim-4)-4./3.*vin(iin,kdim-3)+3.*vin(iin,kdim-2)&
                        - 4.*vin(iin,kdim-1)+25./12.*vin(iin,kdim)
       endif
    elseif(kdim.eq.3.or.kdim.eq.4) then
       if(kin.ge.2.and.kin.le.3) then
           ddk = -0.5*vin(iin,kin-1) + 0.5*vin(iin,kin+1)
       elseif(kin.eq.1) then
           ddk = -1.5*vin(iin,1) + 2.0*vin(iin,2) - 0.5*vin(iin,3)
       elseif(kin.eq.kdim) then
           ddk = 0.5*vin(iin,kdim-2) - 2.0*vin(iin,kdim-1) + 1.5*vin(iin,kdim)
       endif
    elseif(kdim.eq.2) then
           ddk = vin(iin,kdim) - vin(iin,kdim-1)
    elseif(kdim.eq.1) then
           ddk = 0.0
    endif

!    met = Calmet(Calmm(iin,jin,kin,xin,yin,zin))
    select case(n)
      case(1) !ddx
        ddgrid_tsjplane = ddi*buffer_grid(4)+dvindj*buffer_grid(5)+ddk*buffer_grid(6)
      case(2) !ddy
        ddgrid_tsjplane = ddi*buffer_grid(7)+dvindj*buffer_grid(8)+ddk*buffer_grid(9)
      case(3) !ddz
        ddgrid_tsjplane = ddi*buffer_grid(10)+dvindj*buffer_grid(11)+ddk*buffer_grid(12)
      case default
        print*,'n parameter out of range'
        stop
      end select
    end function ddgrid_tsjplane




!input: iin,jin,kin,n
!       xin(:,:,:), yin(:,:,:), zin(:,:,:)
!       vin(:,:,:),met(:,:,:,:,:)



    real(8) function ddgrid2(iin,jin,kin,xin,yin,zin,vin,met,n)
      integer, intent(in) :: iin,jin,kin,n
      real(8), intent(in) :: xin(:,:,:), yin(:,:,:), zin(:,:,:)
      real(8), intent(in) :: vin(:,:,:),met(:,:,:,:,:)
      integer :: idim, jdim, kdim
      real(8) :: ddi, ddj, ddk
!      real(8) :: met(3,3)

      idim = size(xin,dim=1)
      jdim = size(xin,dim=2)
      kdim = size(xin,dim=3)
      if(size(xin).ne.size(yin).or.size(xin).ne.size(zin).or.size(xin).ne.size(vin)) then
          print *, 'The size for input arrays (xin,yin,zin,vin) should be the same. STOP'
          stop
      endif

      if(idim.ge.5) then
        if(iin.ge.3.and.iin.le.idim-2) then
           ddi = (vin(iin-2,jin,kin)-8.*(vin(iin-1,jin,kin)-vin(iin+1,jin,kin))-vin(iin+2,jin,kin))*c12i
        elseif(iin.eq.1) then
           ddi = -25./12.*vin(1,jin,kin)+4.*vin(2,jin,kin)-3.*vin(3,jin,kin)+4./3.*vin(4,jin,kin)-0.25*vin(5,jin,kin)
        elseif(iin.eq.2) then
           ddi = -0.25*vin(1,jin,kin)-5./6.*vin(2,jin,kin)+1.5*vin(3,jin,kin)-0.5*vin(4,jin,kin)+1./12.*vin(5,jin,kin)
        elseif(iin.eq.idim-1) then
           ddi = -1./12.*vin(idim-4,jin,kin)+0.5*vin(idim-3,jin,kin)-1.5*vin(idim-2,jin,kin)&
                          + 5./6.*vin(idim-1,jin,kin)+0.25*vin(idim,jin,kin)
        elseif(iin.eq.idim) then
           ddi = 0.25*vin(idim-4,jin,kin)-4./3.*vin(idim-3,jin,kin)+3.*vin(idim-2,jin,kin)&
                         - 4.*vin(idim-1,jin,kin)+25./12.*vin(idim,jin,kin)
        endif ! end if(iin)
     elseif(idim.eq.3.or.idim.eq.4) then
        if(iin.ge.2.and.iin.le.3) then
           ddi = -0.5*vin(iin-1,jin,kin) + 0.5*vin(iin+1,jin,kin)
        elseif(iin.eq.1) then
           ddi = -1.5*vin(1,jin,kin) + 2.0*vin(2,jin,kin) - 0.5*vin(3,jin,kin)
        elseif(iin.eq.idim) then
           ddi = 0.5*vin(idim-2,jin,kin) - 2.0*vin(idim-1,jin,kin) + 1.5*vin(idim,jin,kin)
       endif
     elseif(idim.eq.2) then
           ddi = vin(idim,jin,kin) - vin(idim-1,jin,kin)
     elseif(idim.eq.1) then
           ddi = 0.0
     endif

    if(jdim.ge.5) then
       if(jin.ge.3.and.jin.le.jdim-2) then
           ddj = (vin(iin,jin-2,kin)-8.*(vin(iin,jin-1,kin)-vin(iin,jin+1,kin))-vin(iin,jin+2,kin))*c12i
       elseif(jin.eq.1) then
           ddj = -25./12.*vin(iin,1,kin)+4.*vin(iin,2,kin)-3.*vin(iin,3,kin)+4./3.*vin(iin,4,kin)-0.25*vin(iin,5,kin)
       elseif(jin.eq.2) then
           ddj = -0.25*vin(iin,1,kin)-5./6.*vin(iin,2,kin)+1.5*vin(iin,3,kin)-0.5*vin(iin,4,kin)+1./12.*vin(iin,5,kin)
       elseif(jin.eq.jdim-1) then
           ddj = -1./12.*vin(iin,jdim-4,kin)+0.5*vin(iin,jdim-3,kin)-1.5*vin(iin,jdim-2,kin)&
                          + 5./6.*vin(iin,jdim-1,kin)+0.25*vin(iin,jdim,kin)
       elseif(jin.eq.jdim) then
           ddj = 0.25*vin(iin,jdim-4,kin)-4./3.*vin(iin,jdim-3,kin)+3.*vin(iin,jdim-2,kin)&
                        - 4.*vin(iin,jdim-1,kin)+25./12.*vin(iin,jdim,kin)
       endif
    elseif(jdim.eq.3.or.jdim.eq.4) then
       if(jin.ge.2.and.jin.le.3) then
           ddj = -0.5*vin(iin,jin-1,kin) + 0.5*vin(iin,jin+1,kin)
       elseif(jin.eq.1) then
           ddj = -1.5*vin(iin,1,kin) + 2.0*vin(iin,2,kin) - 0.5*vin(iin,3,kin)
       elseif(jin.eq.jdim) then
           ddj = 0.5*vin(iin,jdim-2,kin) - 2.0*vin(iin,jdim-1,kin) + 1.5*vin(iin,jdim,kin)
       endif
    elseif(jdim.eq.2) then
           ddj = vin(iin,jdim,kin) - vin(iin,jdim-1,kin)
    elseif(jdim.eq.1) then
           ddj = 0.0
    endif

    if(kdim.ge.5) then
       if(kin.ge.3.and.kin.le.kdim-2) then
           ddk = (vin(iin,jin,kin-2)-8.*(vin(iin,jin,kin-1)-vin(iin,jin,kin+1))-vin(iin,jin,kin+2))*c12i
       elseif(kin.eq.1) then
           ddk = -25./12.*vin(iin,jin,1)+4.*vin(iin,jin,2)-3.*vin(iin,jin,3)+4./3.*vin(iin,jin,4)-0.25*vin(iin,jin,5)
       elseif(kin.eq.2) then
           ddk = -0.25*vin(iin,jin,1)-5./6.*vin(iin,jin,2)+1.5*vin(iin,jin,3)-0.5*vin(iin,jin,4)+1./12.*vin(iin,jin,5)
       elseif(kin.eq.kdim-1) then
           ddk = -1./12.*vin(iin,jin,kdim-4)+0.5*vin(iin,jin,kdim-3)-1.5*vin(iin,jin,kdim-2)&
                          + 5./6.*vin(iin,jin,kdim-1)+0.25*vin(iin,jin,kdim)
       elseif(kin.eq.kdim) then
           ddk = 0.25*vin(iin,jin,kdim-4)-4./3.*vin(iin,jin,kdim-3)+3.*vin(iin,jin,kdim-2)&
                        - 4.*vin(iin,jin,kdim-1)+25./12.*vin(iin,jin,kdim)
       endif
    elseif(kdim.eq.3.or.kdim.eq.4) then
       if(kin.ge.2.and.kin.le.3) then
           ddk = -0.5*vin(iin,jin,kin-1) + 0.5*vin(iin,jin,kin+1)
       elseif(kin.eq.1) then
           ddk = -1.5*vin(iin,jin,1) + 2.0*vin(iin,jin,2) - 0.5*vin(iin,jin,3)
       elseif(kin.eq.kdim) then
           ddk = 0.5*vin(iin,jin,kdim-2) - 2.0*vin(iin,jin,kdim-1) + 1.5*vin(iin,jin,kdim)
       endif
    elseif(kdim.eq.2) then
           ddk = vin(iin,jin,kdim) - vin(iin,jin,kdim-1)
    elseif(kdim.eq.1) then
           ddk = 0.0
    endif

!    met = Calmet(Calmm(iin,jin,kin,xin,yin,zin))
    select case(n)
      case(1) !ddx
        ddgrid2 = ddi*met(1,1,iin,jin,kin)+ddj*met(2,1,iin,jin,kin)+ddk*met(3,1,iin,jin,kin)
      case(2) !ddy
        ddgrid2 = ddi*met(1,2,iin,jin,kin)+ddj*met(2,2,iin,jin,kin)+ddk*met(3,2,iin,jin,kin)
      case(3) !ddz
        ddgrid2 = ddi*met(1,3,iin,jin,kin)+ddj*met(2,3,iin,jin,kin)+ddk*met(3,3,iin,jin,kin)
      case default
        print*,'n parameter out of range'
        stop
      end select
    end function ddgrid2










    function planeddgrid(idim,jdim,kdim,xin,yin,zin,vin,met,n)
      integer, intent(in) :: idim,jdim,kdim,n
      real(8), intent(in) :: xin(:,:,:), yin(:,:,:), zin(:,:,:)
      real(8), intent(in) :: vin(:,:,:), met(:,:,:,:,:)
      integer :: i, j, k
      real(8), dimension(:,:,:),allocatable :: ddi, ddj, ddk
 !     real(8) :: met(3,3,idim,jdim,kdim)
      real(8) :: planeddgrid(idim,jdim,kdim)

      allocate(ddi(idim,jdim,kdim))
      allocate(ddj(idim,jdim,kdim))
      allocate(ddk(idim,jdim,kdim))

      if(size(xin).ne.size(yin).or.size(xin).ne.size(zin).or.size(xin).ne.size(vin)) then
          print *, 'The size for input arrays (xin,yin,zin,vin) should be the same. STOP'
          stop
      endif

      if(idim.ge.5) then

        Do k=1,kdim
          Do j=1,jdim
            ddi(1,j,k)= -25./12.*vin(1,j,k)+4.*vin(2,j,k)-3.*vin(3,j,k)+4./3.*vin(4,j,k)-0.25*vin(5,j,k)
            ddi(2,j,k)= -0.25*vin(1,j,k)-5./6.*vin(2,j,k)+1.5*vin(3,j,k)-0.5*vin(4,j,k)+1./12.*vin(5,j,k)
            ddi(idim-1,j,k)= -1./12.*vin(idim-4,j,k)+0.5*vin(idim-3,j,k)-1.5*vin(idim-2,j,k)&
                          + 5./6.*vin(idim-1,j,k)+0.25*vin(idim,j,k)
            ddi(idim,j,k)= 0.25*vin(idim-4,j,k)-4./3.*vin(idim-3,j,k)+3.*vin(idim-2,j,k)&
                         - 4.*vin(idim-1,j,k)+25./12.*vin(idim,j,k)
            Do i=3,idim-2
              ddi(i,j,k)= (vin(i-2,j,k)-8.*(vin(i-1,j,k)-vin(i+1,j,k))-vin(i+2,j,k))*c12i
            End Do
          End Do
        End Do

      End If



!        if(iin.ge.3.and.iin.le.idim-2) then
!           ddi = (vin(iin-2,jin,kin)-8.*(vin(iin-1,jin,kin)-vin(iin+1,jin,kin))-vin(iin+2,jin,kin))*c12i
!        elseif(iin.eq.1) then
!           ddi = -25./12.*vin(1,jin,kin)+4.*vin(2,jin,kin)-3.*vin(3,jin,kin)+4./3.*vin(4,jin,kin)-0.25*vin(5,jin,kin)
!        elseif(iin.eq.2) then
!           ddi = -0.25*vin(1,jin,kin)-5./6.*vin(2,jin,kin)+1.5*vin(3,jin,kin)-0.5*vin(4,jin,kin)+1./12.*vin(5,jin,kin)
!        elseif(iin.eq.idim-1) then
!           ddi = -1./12.*vin(idim-4,jin,kin)+0.5*vin(idim-3,jin,kin)-1.5*vin(idim-2,jin,kin)&
!                          + 5./6.*vin(idim-1,jin,kin)+0.25*vin(idim,jin,kin)
!        elseif(iin.eq.idim) then
!           ddi = 0.25*vin(idim-4,jin,kin)-4./3.*vin(idim-3,jin,kin)+3.*vin(idim-2,jin,kin)&
!                         - 4.*vin(idim-1,jin,kin)+25./12.*vin(idim,jin,kin)
!        endif ! end if(iin)
!     elseif(idim.eq.3.or.idim.eq.4) then
!        if(iin.ge.2.and.iin.le.3) then
!           ddi = -0.5*vin(iin-1,jin,kin) + 0.5*vin(iin+1,jin,kin)
!        elseif(iin.eq.1) then
!           ddi = -1.5*vin(1,jin,kin) + 2.0*vin(2,jin,kin) - 0.5*vin(3,jin,kin)
!        elseif(iin.eq.idim) then
!           ddi = 0.5*vin(idim-2,jin,kin) - 2.0*vin(idim-1,jin,kin) + 1.5*vin(idim,jin,kin)
!       endif
!     elseif(idim.eq.2) then
!           ddi = vin(idim,jin,kin) - vin(idim-1,jin,kin)
!     elseif(idim.eq.1) then
!           ddi = 0.0
!     endif

      if(jdim.ge.5) then

        Do k=1,kdim
          Do i=1,idim
            ddj(i,1,k)= -25./12.*vin(i,1,k)+4.*vin(i,2,k)-3.*vin(i,3,k)+4./3.*vin(i,4,k)-0.25*vin(i,5,k)
            ddj(i,2,k)= -0.25*vin(i,1,k)-5./6.*vin(i,2,k)+1.5*vin(i,3,k)-0.5*vin(i,4,k)+1./12.*vin(i,5,k)
            ddj(i,jdim-1,k)= -1./12.*vin(i,jdim-4,k)+0.5*vin(i,jdim-3,k)-1.5*vin(i,jdim-2,k)&
                          + 5./6.*vin(i,jdim-1,k)+0.25*vin(i,jdim,k)
            ddj(i,jdim,k)= 0.25*vin(i,jdim-4,k)-4./3.*vin(i,jdim-3,k)+3.*vin(i,jdim-2,k)&
                        - 4.*vin(i,jdim-1,k)+25./12.*vin(i,jdim,k)
            Do j=3,jdim-2
              ddj(i,j,k)= (vin(i,j-2,k)-8.*(vin(i,j-1,k)-vin(i,j+1,k))-vin(i,j+2,k))*c12i
            End Do
          End Do
        End Do

      End If

!    if(jdim.ge.5) then
!       if(jin.ge.3.and.jin.le.jdim-2) then
!           ddj = (vin(iin,jin-2,kin)-8.*(vin(iin,jin-1,kin)-vin(iin,jin+1,kin))-vin(iin,jin+2,kin))*c12i
!       elseif(jin.eq.1) then
!           ddj = -25./12.*vin(iin,1,kin)+4.*vin(iin,2,kin)-3.*vin(iin,3,kin)+4./3.*vin(iin,4,kin)-0.25*vin(iin,5,kin)
!       elseif(jin.eq.2) then
!           ddj = -0.25*vin(iin,1,kin)-5./6.*vin(iin,2,kin)+1.5*vin(iin,3,kin)-0.5*vin(iin,4,kin)+1./12.*vin(iin,5,kin)
!       elseif(jin.eq.jdim-1) then
!           ddj = -1./12.*vin(iin,jdim-4,kin)+0.5*vin(iin,jdim-3,kin)-1.5*vin(iin,jdim-2,kin)&
!                          + 5./6.*vin(iin,jdim-1,kin)+0.25*vin(iin,jdim,kin)
!       elseif(jin.eq.jdim) then
!           ddj = 0.25*vin(iin,jdim-4,kin)-4./3.*vin(iin,jdim-3,kin)+3.*vin(iin,jdim-2,kin)&
!                        - 4.*vin(iin,jdim-1,kin)+25./12.*vin(iin,jdim,kin)
!       endif
!    elseif(jdim.eq.3.or.jdim.eq.4) then
!       if(jin.ge.2.and.jin.le.3) then
!           ddj = -0.5*vin(iin,jin-1,kin) + 0.5*vin(iin,jin+1,kin)
!       elseif(jin.eq.1) then
!           ddj = -1.5*vin(iin,1,kin) + 2.0*vin(iin,2,kin) - 0.5*vin(iin,3,kin)
!       elseif(jin.eq.jdim) then
!           ddj = 0.5*vin(iin,jdim-2,kin) - 2.0*vin(iin,jdim-1,kin) + 1.5*vin(iin,jdim,kin)
!       endif
!    elseif(jdim.eq.2) then
!           ddj = vin(iin,jdim,kin) - vin(iin,jdim-1,kin)
!    elseif(jdim.eq.1) then
!           ddj = 0.0
!    endif

      if(kdim.ge.5) then

        Do i=1,idim
          Do j=1,jdim
            ddk(i,j,1)= -25./12.*vin(i,j,1)+4.*vin(i,j,2)-3.*vin(i,j,3)+4./3.*vin(i,j,4)-0.25*vin(i,j,5)
            ddk(i,j,2)= -0.25*vin(i,j,1)-5./6.*vin(i,j,2)+1.5*vin(i,j,3)-0.5*vin(i,j,4)+1./12.*vin(i,j,5)
            ddk(i,j,kdim-1)= -1./12.*vin(i,j,kdim-4)+0.5*vin(i,j,kdim-3)-1.5*vin(i,j,kdim-2)&
                          + 5./6.*vin(i,j,kdim-1)+0.25*vin(i,j,kdim)
            ddk(i,j,kdim)= 0.25*vin(i,j,kdim-4)-4./3.*vin(i,j,kdim-3)+3.*vin(i,j,kdim-2)&
                        - 4.*vin(i,j,kdim-1)+25./12.*vin(i,j,kdim)
            Do k=3,kdim-2
              ddk(i,j,k)= (vin(i,j,k-2)-8.*(vin(i,j,k-1)-vin(i,j,k+1))-vin(i,j,k+2))*c12i
            End Do
          End Do
        End Do

      End If



!    if(kdim.ge.5) then
!       if(kin.ge.3.and.kin.le.kdim-2) then
!           ddk = (vin(iin,jin,kin-2)-8.*(vin(iin,jin,kin-1)-vin(iin,jin,kin+1))-vin(iin,jin,kin+2))*c12i
!       elseif(kin.eq.1) then
!           ddk = -25./12.*vin(iin,jin,1)+4.*vin(iin,jin,2)-3.*vin(iin,jin,3)+4./3.*vin(iin,jin,4)-0.25*vin(iin,jin,5)
!       elseif(kin.eq.2) then
!           ddk = -0.25*vin(iin,jin,1)-5./6.*vin(iin,jin,2)+1.5*vin(iin,jin,3)-0.5*vin(iin,jin,4)+1./12.*vin(iin,jin,5)
!       elseif(kin.eq.kdim-1) then
!           ddk = -1./12.*vin(iin,jin,kdim-4)+0.5*vin(iin,jin,kdim-3)-1.5*vin(iin,jin,kdim-2)&
!                          + 5./6.*vin(iin,jin,kdim-1)+0.25*vin(iin,jin,kdim)
!       elseif(kin.eq.kdim) then
!           ddk = 0.25*vin(iin,jin,kdim-4)-4./3.*vin(iin,jin,kdim-3)+3.*vin(iin,jin,kdim-2)&
!                        - 4.*vin(iin,jin,kdim-1)+25./12.*vin(iin,jin,kdim)
!       endif
!    elseif(kdim.eq.3.or.kdim.eq.4) then
!       if(kin.ge.2.and.kin.le.3) then
!           ddk = -0.5*vin(iin,jin,kin-1) + 0.5*vin(iin,jin,kin+1)
!       elseif(kin.eq.1) then
!           ddk = -1.5*vin(iin,jin,1) + 2.0*vin(iin,jin,2) - 0.5*vin(iin,jin,3)
!       elseif(kin.eq.kdim) then
!           ddk = 0.5*vin(iin,jin,kdim-2) - 2.0*vin(iin,jin,kdim-1) + 1.5*vin(iin,jin,kdim)
!       endif
!    elseif(kdim.eq.2) then
!           ddk = vin(iin,jin,kdim) - vin(iin,jin,kdim-1)
!    elseif(kdim.eq.1) then
!           ddk = 0.0
!    endif

!    met = planeCalmet(idim,jdim,kdim,planeCalmm(idim,jdim,kdim,xin,yin,zin))

    Do k=1,kdim
      Do j=1,jdim
        Do i=1,idim
    select case(n)
      case(1) !ddx
        planeddgrid(i,j,k) = ddi(i,j,k)*met(1,1,i,j,k)+ddj(i,j,k)*met(2,1,i,j,k)+ddk(i,j,k)*met(3,1,i,j,k)
      case(2) !ddy
        planeddgrid(i,j,k) = ddi(i,j,k)*met(1,2,i,j,k)+ddj(i,j,k)*met(2,2,i,j,k)+ddk(i,j,k)*met(3,2,i,j,k)
      case(3) !ddz
        planeddgrid(i,j,k) = ddi(i,j,k)*met(1,3,i,j,k)+ddj(i,j,k)*met(2,3,i,j,k)+ddk(i,j,k)*met(3,3,i,j,k)
      case default
        print*,'n parameter out of range'
        stop
      end select
          End Do
        End Do
     End Do
    end function planeddgrid












!      subroutine Deriv_i(idim,jdim,kdim,x,y,z,dxdi,dydi,dzdi)
!      integer, intent(in) :: idim jdim, kdim
!      real(8), dimension(idim,jdim,kdim),intent(in) :: x, y, z
!      real(8), dimension(idim,jdim,kdim),intent(out) :: dxdi,dydi,dzdi
!
!      if(idim.eq.1) then
!         dxdi = 1.0
!         dydi = 1.0
!         dzdi = 1.0
!      else      
!         forall(j=1:jdim, k=1:kdim)
!            dxdi(:,j,k) = dfdi(idim,x(:,j,k))
!            dydi(:,j,k) = dfdi(idim,y(:,j,k))
!            dzdi(:,j,k) = dfdi(idim,z(:,j,k))
!         end forall
!      endif
!    
!      end subroutine Deriv_i
!
!      subroutine Deriv_i(idim,jdim,kdim,x,y,z,dxdi,dydi,dzdi)
!      integer, intent(in) :: idim jdim, kdim
!      real(8), dimension(idim,jdim,kdim),intent(in) :: x, y, z
!      real(8), dimension(idim,jdim,kdim),intent(out) :: dxdi,dydi,dzdi
!
!      forall(i=1:idim, k=1:kdim)
!         dxdj(i,:,k) = dfdi(jdim,x(i,:,k))
!         dydj(i,:,k) = dfdi(jdim,y(i,:,k))
!         dzdj(i,:,k) = dfdi(jdim,z(i,:,k))
!      end forall
!
!      forall(i=1:idim, j=1:jdim)
!         dxdk(i,j,:) = dfdi(kdim,x(i,j,:))
!         dydk(i,j,:) = dfdi(kdim,y(i,j,:))
!         dzdk(i,j,:) = dfdi(kdim,z(i,j,:))
!      end forall
!
!      end subroutine CalDerivative
!
!
!      subroutine CalDerivative(idim,jdim,kdim,x,y,z,dxdi,dxdj,dxdk,&
!                               dydi,dydj,dydk,dzdi,dzdj,dzdk)
!      integer, intent(in) :: idim jdim, kdim
!      real(8), dimension(idim,jdim,kdim),intent(in) :: x, y, z
!      real(8), dimension(idim,jdim,kdim),intent(out) :: dxdi,dxdj,dxdk, &
!                                                        dydi,dydj,dydk, &
!                                                        dzdi,dzdj,dzdk
!      
!      forall(j=1:jdim, k=1:kdim)
!         dxdi(:,j,k) = dfdi(idim,x(:,j,k))
!         dydi(:,j,k) = dfdi(idim,y(:,j,k))
!         dzdi(:,j,k) = dfdi(idim,z(:,j,k))
!      end forall
!
!      forall(i=1:idim, k=1:kdim)
!         dxdj(i,:,k) = dfdi(jdim,x(i,:,k))
!         dydj(i,:,k) = dfdi(jdim,y(i,:,k))
!         dzdj(i,:,k) = dfdi(jdim,z(i,:,k))
!      end forall
!
!      forall(i=1:idim, j=1:jdim)
!         dxdk(i,j,:) = dfdi(kdim,x(i,j,:))
!         dydk(i,j,:) = dfdi(kdim,y(i,j,:))
!         dzdk(i,j,:) = dfdi(kdim,z(i,j,:))
!      end forall
!
!      end subroutine CalDerivative
!

!      subroutine CalJacobian(idim,jdim,kdim,x,y,z,kapi,kapj,kapk,vol)
!     kapi(1) = didx, kapi(2) = didy, kapi(3) = didz, kapi(4) = sqrt(didx**2+didy**2+didz**2)
!     kapj(1) = djdx, kapj(2) = djdy, kapj(3) = djdz, kapj(4) = sqrt(djdx**2+djdy**2+djdz**2)
!     kapk(1) = dkdx, kapk(2) = dkdy, kapk(3) = dkdz, kapk(4) = sqrt(dkdx**2+dkdy**2+dkdz**2)
!      integer, intent(in) :: idim jdim, kdim
!      real(8), dimension(idim,jdim,kdim),intent(in) :: x, y, z
!      real(8), dimension(idim,jdim,kdim),intent(out) :: x, y, z
!      real :: xdi, xdj, xdk, ydi, ydj, ydk, zdi, zdj, zdk, jacobi
 
!      do i = 1, idim
!      do j = 1, jdim
!      do k = 1, kdim
!      xdi  = dxdi(i,j,k)
!      xdj  = dxdj(iin,jin,kin)
!      xdk  = dxdk(iin,jin,kin)
!      ydi  = dydi(iin,jin,kin)
!      ydj  = dydj(iin,jin,kin)
!      ydk  = dydk(iin,jin,kin)
!      zdi  = dzdi(iin,jin,kin)
!      zdj  = dzdj(iin,jin,kin)
!      zdk  = dzdk(iin,jin,kin)
!      jacobi = xdi*(ydj*zdk-ydk*zdj)-ydi*(xdj*zdk-xdk*zdj)&
!              +zdi*(xdj*ydk-xdk*ydj)
!      kapi(i,j,k,1) =  (ydj*zdk-ydk-zdj)/jacobi
!      kapi(i,j,k,2) = -(xdj*zdk-xdk*zdj)/jacobi
!      kapi(i,j,k,3) =  (xdj*ydk-xdk*ydj)/jacobi
!      kapi(i,j,k,4) = sqrt( kapi(i,j,k,1)**2 + kapi(i,j,k,2)**2 + kapi(i,j,k,3)**2 )
!      kapj(i,j,k,1) = -(ydi*zdk-ydk*zdi)/jacobi
!      kapj(i,j,k,2) =  (xdi*zdk-xdk*zdi)/jacobi
!      kapj(i,j,k,3) = -(xdi*ydk-xdk*ydi)/jacobi
!      kapj(i,j,k,4) = sqrt( kapj(i,j,k,1)**2 + kapj(i,j,k,2)**2 + kapj(i,j,k,3)**2 )
!      dkdx =  (ydi*zdj-ydj*zdi)/jacobi
!      dkdy = -(xdi*zdj-xdj*zdi)/jacobi   ! Check
!      dkdz =  (xdi*ydj-xdj*ydi)/jacobi
!      vol = jacobi
!      ddj = sqrt( djdx**2 + djdy**2 + djdz**2 )
!      ddk = sqrt( dkdx**2 + dkdy**2 + dkdz**2 )
!      end subroutine CalJacobian
!
!      subroutine CalJacobian(iin,jin,kin)
!      integer, intent(in) :: iin, jin, kin
!      real :: xdi, xdj, xdk, ydi, ydj, ydk, zdi, zdj, zdk, jacobi
!
!      xdi  = dxdi(iin,jin,kin)
!      xdj  = dxdj(iin,jin,kin)
!      xdk  = dxdk(iin,jin,kin)
!      ydi  = dydi(iin,jin,kin)
!      ydj  = dydj(iin,jin,kin)
!      ydk  = dydk(iin,jin,kin)
!      zdi  = dzdi(iin,jin,kin)
!      zdj  = dzdj(iin,jin,kin)
!      zdk  = dzdk(iin,jin,kin)
!      jacobi = xdi*(ydj*zdk-ydk*zdj)-ydi*(xdj*zdk-xdk*zdj)&
!              +zdi*(xdj*ydk-xdk*ydj)
!      didx =  (ydj*zdk-ydk-zdj)/jacobi
!      didy = -(xdj*zdk-xdk*zdj)/jacobi
!      didz =  (xdj*ydk-xdk*ydj)/jacobi
!      djdx = -(ydi*zdk-ydk*zdi)/jacobi
!      djdy =  (xdi*zdk-xdk*zdi)/jacobi
!      djdz = -(xdi*ydk-xdk*ydi)/jacobi
!      dkdx =  (ydi*zdj-ydj*zdi)/jacobi
!      dkdy = -(xdi*zdj-xdj*zdi)/jacobi   ! Check
!      dkdz =  (xdi*ydj-xdj*ydi)/jacobi
!      vol = jacobi
!      ddi = sqrt( didx**2 + didy**2 + didz**2 )
!      ddj = sqrt( djdx**2 + djdy**2 + djdz**2 )
!      ddk = sqrt( dkdx**2 + dkdy**2 + dkdz**2 )
!      end subroutine CalJacobian
!
!      pure real(8) function diffcenter2(val)
!      real(8), dimension(-2:2), intent(in) :: val
!      diffcenter2 = (val(-2)-8.0*(val(-1)-val(1))-val(2))/12.
!      end function diffcenter2

!      pure real(8) function diff0p4(val)
!      real(8), dimension(0:4), intent(in) :: val
!      diff0p4 = -25./12.*val(0)+4.*val(1)-3.*val(2)+4./3.*val(3)-0.25*val(4)
!      end function diff0p4

!      pure real(8) function diffm1p3(val)
!      real(8), dimension(-1:3), intent(in) :: val
!      diffm1p3 = -0.25*val(-1)-5./6.*val(0)+1.5*val(1)-0.5*val(2)+1./12.*val(3)
!      end function diffm1p3

!      pure real(8) function diffm3p1(val)
!      real(8), dimension(-3:1), intent(in) :: val
!      diffm3p1 = -1./12.*val(-3)+0.5*val(-2)-1.5*val(-1) &
!                          + 5./6.*val(0)+0.25*val(1)
!      end function diffm3p1

!      pure real(8) function diffm40(val)
!      real(8), dimension(-4:0), intent(in) :: val
!      diffm40 = 0.25*val(-4)-4./3.*val(-3)+3.*val(-2) &
!                        - 4.*val(-1)+25./12.*val(0)
!      end function diffm40

!      pure real(8) function ddx(iin,jin,kin,val)
!      integer, intent(in) :: iin, jin, kin
!      real(8), dimension(:,:,:), target, intent(in) :: val
!      real(8), dimension(:), pointer :: vx, vy, vz
!      vx => val(iin-2:iin+2,jin,kin)
!      vy => val(iin,jin-2:jin+2,kin)
!      vz => val(iin,jin,kin-2:kin+2)
!      ddx = didx*diffcenter2(vx)+djdx*diffcenter2(vy)+dkdx*diffcenter2(vz)
!      end function ddx

!      pure, real(8) function ddy(iin, jin, kin,val)
!      integer, intent(in) :: iin, jin, kin
!      real(8), dimension(:,:,:), target, intent(in) :: val
!!      real(8), dimension(:), pointer :: vx, vy, vz
!      vx => val(iin-2:iin+2,jin,kin)
!      vy => val(iin,jin-2:jin+2,kin)
!      vz => val(iin,jin,kin-2:kin+2)
!      ddy = didy*diffcenter2(vx)+djdy*diffcenter2(vy)+dkdy*diffcenter2(vz)
!      end function ddy

!      pure, real(8) function ddz(iin, jin, kin,val)
!      integer, intent(in) :: iin, jin, kin
!      real(8), dimension(:,:,:), target, intent(in) :: val
!      real(8), dimension(:), pointer :: vx, vy, vz
!      vx => val(iin-2:iin+2,jin,kin)
!      vy => val(iin,jin-2:jin+2,kin)
!      vz => val(iin,jin,kin-2:kin+2)
!      ddz = didz*diffcenter2(vx)+djdz*diffcenter2(vy)+dkdz*diffcenter2(vz)
!      end function ddz
!




      end module MDerivative

