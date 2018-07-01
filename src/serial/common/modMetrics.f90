module modMetrics
implicit none
real(8), parameter, private :: c12i = 1.d0/12.d0

contains





    function Calmm_2D(iin,jin, x, y)

    !  mm(2,2)
    !  mm(1,1) = dxdi, mm(1,2) = dxdj
    !  mm(2,1) = dydi, mm(2,2) = dydj

    integer, intent(in) :: iin, jin
    real(8), intent(in) :: x(:,:), y(:,:)
    real(8) :: Calmm_2D(2,2)
    integer :: idim, jdim
    real(8) :: mm(2,2)

    idim = size(x,dim=1)
    jdim = size(x,dim=2)

    if(idim.ge.5) then
       if(iin.ge.3.and.iin.le.idim-2) then
          mm(1,1) = (x(iin-2,jin)-8.*(x(iin-1,jin)-x(iin+1,jin))-x(iin+2,jin))*c12i
          mm(2,1) = (y(iin-2,jin)-8.*(y(iin-1,jin)-y(iin+1,jin))-y(iin+2,jin))*c12i
       elseif(iin.eq.1) then
          mm(1,1) = -25./12.*x(1,jin)+4.*x(2,jin)-3.*x(3,jin)+4./3.*x(4,jin)-0.25*x(5,jin)
          mm(2,1) = -25./12.*y(1,jin)+4.*y(2,jin)-3.*y(3,jin)+4./3.*y(4,jin)-0.25*y(5,jin)
       elseif(iin.eq.2) then
          mm(1,1) = -0.25*x(1,jin)-5./6.*x(2,jin)+1.5*x(3,jin)-0.5*x(4,jin)+1./12.*x(5,jin)
          mm(2,1) = -0.25*y(1,jin)-5./6.*y(2,jin)+1.5*y(3,jin)-0.5*y(4,jin)+1./12.*y(5,jin)
       elseif(iin.eq.idim-1) then
          mm(1,1) = -1./12.*x(idim-4,jin)+0.5*x(idim-3,jin)-1.5*x(idim-2,jin)&
                          + 5./6.*x(idim-1,jin)+0.25*x(idim,jin)
          mm(2,1) = -1./12.*y(idim-4,jin)+0.5*y(idim-3,jin)-1.5*y(idim-2,jin)&
                          + 5./6.*y(idim-1,jin)+0.25*y(idim,jin)
       elseif(iin.eq.idim) then
          mm(1,1) = 0.25*x(idim-4,jin)-4./3.*x(idim-3,jin)+3.*x(idim-2,jin)&
                         - 4.*x(idim-1,jin)+25./12.*x(idim,jin)
          mm(2,1) = 0.25*y(idim-4,jin)-4./3.*y(idim-3,jin)+3.*y(idim-2,jin)&
                         - 4.*y(idim-1,jin)+25./12.*y(idim,jin)
       endif ! end if(iin)
    elseif(idim.eq.3.or.idim.eq.4) then
       if(iin.ge.2.and.iin.le.3) then
          mm(1,1) = -0.5*x(iin-1,jin) + 0.5*x(iin+1,jin)
          mm(2,1) = -0.5*y(iin-1,jin) + 0.5*y(iin+1,jin)
       elseif(iin.eq.1) then
          mm(1,1) = -1.5*x(1,jin) + 2.0*x(2,jin) - 0.5*x(3,jin)
          mm(2,1) = -1.5*y(1,jin) + 2.0*y(2,jin) - 0.5*y(3,jin)
       elseif(iin.eq.idim) then
          mm(1,1) = 0.5*x(idim-2,jin) - 2.0*x(idim-1,jin) + 1.5*x(idim,jin)
          mm(2,1) = 0.5*y(idim-2,jin) - 2.0*y(idim-1,jin) + 1.5*y(idim,jin)
       endif
    elseif(idim.eq.2) then
          mm(1,1) = x(idim,jin) - x(idim-1,jin)
          mm(2,1) = y(idim,jin) - y(idim-1,jin)
    elseif(idim.eq.1) then
          mm(1,1) = 0.0
          mm(2,1) = 0.0
    endif

    if(jdim.ge.5) then
       if(jin.ge.3.and.jin.le.jdim-2) then
          mm(1,2) = (x(iin,jin-2)-8.*(x(iin,jin-1)-x(iin,jin+1))-x(iin,jin+2))*c12i
          mm(2,2) = (y(iin,jin-2)-8.*(y(iin,jin-1)-y(iin,jin+1))-y(iin,jin+2))*c12i
       elseif(jin.eq.1) then
          mm(1,2) = -25./12.*x(iin,1)+4.*x(iin,2)-3.*x(iin,3)+4./3.*x(iin,4)-0.25*x(iin,5)
          mm(2,2) = -25./12.*y(iin,1)+4.*y(iin,2)-3.*y(iin,3)+4./3.*y(iin,4)-0.25*y(iin,5)
       elseif(jin.eq.2) then
          mm(1,2) = -0.25*x(iin,1)-5./6.*x(iin,2)+1.5*x(iin,3)-0.5*x(iin,4)+1./12.*x(iin,5)
          mm(2,2) = -0.25*y(iin,1)-5./6.*y(iin,2)+1.5*y(iin,3)-0.5*y(iin,4)+1./12.*y(iin,5)
       elseif(jin.eq.jdim-1) then
          mm(1,2) = -1./12.*x(iin,jdim-4)+0.5*x(iin,jdim-3)-1.5*x(iin,jdim-2)&
                          + 5./6.*x(iin,jdim-1)+0.25*x(iin,jdim)
          mm(2,2) = -1./12.*y(iin,jdim-4)+0.5*y(iin,jdim-3)-1.5*y(iin,jdim-2)&
                          + 5./6.*y(iin,jdim-1)+0.25*y(iin,jdim)
       elseif(jin.eq.jdim) then
          mm(1,2) = 0.25*x(iin,jdim-4)-4./3.*x(iin,jdim-3)+3.*x(iin,jdim-2)&
                        - 4.*x(iin,jdim-1)+25./12.*x(iin,jdim)
          mm(2,2) = 0.25*y(iin,jdim-4)-4./3.*y(iin,jdim-3)+3.*y(iin,jdim-2)&
                        - 4.*y(iin,jdim-1)+25./12.*y(iin,jdim)
       endif
    elseif(jdim.eq.3.or.jdim.eq.4) then
       if(jin.ge.2.and.jin.le.3) then
          mm(1,2) = -0.5*x(iin,jin-1) + 0.5*x(iin,jin+1)
          mm(2,2) = -0.5*y(iin,jin-1) + 0.5*y(iin,jin+1)
       elseif(jin.eq.1) then
          mm(1,2) = -1.5*x(iin,1) + 2.0*x(iin,2) - 0.5*x(iin,3)
          mm(2,2) = -1.5*y(iin,1) + 2.0*y(iin,2) - 0.5*y(iin,3)
       elseif(jin.eq.jdim) then
          mm(1,2) = 0.5*x(iin,jdim-2) - 2.0*x(iin,jdim-1) + 1.5*x(iin,jdim)
          mm(2,2) = 0.5*y(iin,jdim-2) - 2.0*y(iin,jdim-1) + 1.5*y(iin,jdim)
       endif
    elseif(jdim.eq.2) then
          mm(1,2) = x(iin,jdim) - x(iin,jdim-1)
          mm(2,2) = y(iin,jdim) - y(iin,jdim-1)
    elseif(jdim.eq.1) then
          mm(1,2) = 0.0
          mm(2,2) = 0.0
    endif

    Calmm_2D = mm
    end function Calmm_2D


    function Calmm(iin,jin,kin, x, y, z)

    !  mm(3,3)
    !  mm(1,1) = dxdi, mm(1,2) = dxdj, mm(1,3) = dxdk
    !  mm(2,1) = dydi, mm(2,2) = dydj, mm(2,3) = dydk
    !  mm(3,1) = dzdi, mm(3,2) = dzdj, mm(3,3) = dzdk

    integer, intent(in) :: iin, jin, kin
    real(8), intent(in) :: x(:,:,:), y(:,:,:), z(:,:,:)
    real(8) :: Calmm(3,3)
    integer :: idim, jdim, kdim
    real(8) :: mm(3,3)

    idim = size(x,dim=1)
    jdim = size(x,dim=2)
    kdim = size(x,dim=3)

    if(idim.ge.5) then
       if(iin.ge.3.and.iin.le.idim-2) then
          mm(1,1) = (x(iin-2,jin,kin)-8.*(x(iin-1,jin,kin)-x(iin+1,jin,kin))-x(iin+2,jin,kin))*c12i
          mm(2,1) = (y(iin-2,jin,kin)-8.*(y(iin-1,jin,kin)-y(iin+1,jin,kin))-y(iin+2,jin,kin))*c12i
          mm(3,1) = (z(iin-2,jin,kin)-8.*(z(iin-1,jin,kin)-z(iin+1,jin,kin))-z(iin+2,jin,kin))*c12i
       elseif(iin.eq.1) then
          mm(1,1) = -25./12.*x(1,jin,kin)+4.*x(2,jin,kin)-3.*x(3,jin,kin)+4./3.*x(4,jin,kin)-0.25*x(5,jin,kin)
          mm(2,1) = -25./12.*y(1,jin,kin)+4.*y(2,jin,kin)-3.*y(3,jin,kin)+4./3.*y(4,jin,kin)-0.25*y(5,jin,kin)
          mm(3,1) = -25./12.*z(1,jin,kin)+4.*z(2,jin,kin)-3.*z(3,jin,kin)+4./3.*z(4,jin,kin)-0.25*z(5,jin,kin)
       elseif(iin.eq.2) then
          mm(1,1) = -0.25*x(1,jin,kin)-5./6.*x(2,jin,kin)+1.5*x(3,jin,kin)-0.5*x(4,jin,kin)+1./12.*x(5,jin,kin)
          mm(2,1) = -0.25*y(1,jin,kin)-5./6.*y(2,jin,kin)+1.5*y(3,jin,kin)-0.5*y(4,jin,kin)+1./12.*y(5,jin,kin)
          mm(3,1) = -0.25*z(1,jin,kin)-5./6.*z(2,jin,kin)+1.5*z(3,jin,kin)-0.5*z(4,jin,kin)+1./12.*z(5,jin,kin)
       elseif(iin.eq.idim-1) then
          mm(1,1) = -1./12.*x(idim-4,jin,kin)+0.5*x(idim-3,jin,kin)-1.5*x(idim-2,jin,kin)&
                          + 5./6.*x(idim-1,jin,kin)+0.25*x(idim,jin,kin)
          mm(2,1) = -1./12.*y(idim-4,jin,kin)+0.5*y(idim-3,jin,kin)-1.5*y(idim-2,jin,kin)&
                          + 5./6.*y(idim-1,jin,kin)+0.25*y(idim,jin,kin)
          mm(3,1) = -1./12.*z(idim-4,jin,kin)+0.5*z(idim-3,jin,kin)-1.5*z(idim-2,jin,kin)&
                          + 5./6.*z(idim-1,jin,kin)+0.25*z(idim,jin,kin)
       elseif(iin.eq.idim) then
          mm(1,1) = 0.25*x(idim-4,jin,kin)-4./3.*x(idim-3,jin,kin)+3.*x(idim-2,jin,kin)&
                         - 4.*x(idim-1,jin,kin)+25./12.*x(idim,jin,kin)
          mm(2,1) = 0.25*y(idim-4,jin,kin)-4./3.*y(idim-3,jin,kin)+3.*y(idim-2,jin,kin)&
                         - 4.*y(idim-1,jin,kin)+25./12.*y(idim,jin,kin)
          mm(3,1) = 0.25*z(idim-4,jin,kin)-4./3.*z(idim-3,jin,kin)+3.*z(idim-2,jin,kin)&
                         - 4.*z(idim-1,jin,kin)+25./12.*z(idim,jin,kin)
       endif ! end if(iin)
    elseif(idim.eq.3.or.idim.eq.4) then
       if(iin.ge.2.and.iin.le.3) then
          mm(1,1) = -0.5*x(iin-1,jin,kin) + 0.5*x(iin+1,jin,kin)
          mm(2,1) = -0.5*y(iin-1,jin,kin) + 0.5*y(iin+1,jin,kin)
          mm(3,1) = -0.5*z(iin-1,jin,kin) + 0.5*z(iin+1,jin,kin)
       elseif(iin.eq.1) then        
          mm(1,1) = -1.5*x(1,jin,kin) + 2.0*x(2,jin,kin) - 0.5*x(3,jin,kin)
          mm(2,1) = -1.5*y(1,jin,kin) + 2.0*y(2,jin,kin) - 0.5*y(3,jin,kin)
          mm(3,1) = -1.5*z(1,jin,kin) + 2.0*z(2,jin,kin) - 0.5*z(3,jin,kin)
       elseif(iin.eq.idim) then
          mm(1,1) = 0.5*x(idim-2,jin,kin) - 2.0*x(idim-1,jin,kin) + 1.5*x(idim,jin,kin)
          mm(2,1) = 0.5*y(idim-2,jin,kin) - 2.0*y(idim-1,jin,kin) + 1.5*y(idim,jin,kin)
          mm(3,1) = 0.5*z(idim-2,jin,kin) - 2.0*z(idim-1,jin,kin) + 1.5*z(idim,jin,kin)
       endif
    elseif(idim.eq.2) then
          mm(1,1) = x(idim,jin,kin) - x(idim-1,jin,kin)
          mm(2,1) = y(idim,jin,kin) - y(idim-1,jin,kin)
          mm(3,1) = z(idim,jin,kin) - z(idim-1,jin,kin)
    elseif(idim.eq.1) then
          mm(1,1) = 0.0
          mm(2,1) = 0.0
          mm(3,1) = 0.0
    endif

    if(jdim.ge.5) then
       if(jin.ge.3.and.jin.le.jdim-2) then
          mm(1,2) = (x(iin,jin-2,kin)-8.*(x(iin,jin-1,kin)-x(iin,jin+1,kin))-x(iin,jin+2,kin))*c12i
          mm(2,2) = (y(iin,jin-2,kin)-8.*(y(iin,jin-1,kin)-y(iin,jin+1,kin))-y(iin,jin+2,kin))*c12i
          mm(3,2) = (z(iin,jin-2,kin)-8.*(z(iin,jin-1,kin)-z(iin,jin+1,kin))-z(iin,jin+2,kin))*c12i
       elseif(jin.eq.1) then
          mm(1,2) = -25./12.*x(iin,1,kin)+4.*x(iin,2,kin)-3.*x(iin,3,kin)+4./3.*x(iin,4,kin)-0.25*x(iin,5,kin)
          mm(2,2) = -25./12.*y(iin,1,kin)+4.*y(iin,2,kin)-3.*y(iin,3,kin)+4./3.*y(iin,4,kin)-0.25*y(iin,5,kin)
          mm(3,2) = -25./12.*z(iin,1,kin)+4.*z(iin,2,kin)-3.*z(iin,3,kin)+4./3.*z(iin,4,kin)-0.25*z(iin,5,kin)
       elseif(jin.eq.2) then
          mm(1,2) = -0.25*x(iin,1,kin)-5./6.*x(iin,2,kin)+1.5*x(iin,3,kin)-0.5*x(iin,4,kin)+1./12.*x(iin,5,kin)
          mm(2,2) = -0.25*y(iin,1,kin)-5./6.*y(iin,2,kin)+1.5*y(iin,3,kin)-0.5*y(iin,4,kin)+1./12.*y(iin,5,kin)
          mm(3,2) = -0.25*z(iin,1,kin)-5./6.*z(iin,2,kin)+1.5*z(iin,3,kin)-0.5*z(iin,4,kin)+1./12.*z(iin,5,kin)
       elseif(jin.eq.jdim-1) then
          mm(1,2) = -1./12.*x(iin,jdim-4,kin)+0.5*x(iin,jdim-3,kin)-1.5*x(iin,jdim-2,kin)&
                          + 5./6.*x(iin,jdim-1,kin)+0.25*x(iin,jdim,kin)
          mm(2,2) = -1./12.*y(iin,jdim-4,kin)+0.5*y(iin,jdim-3,kin)-1.5*y(iin,jdim-2,kin)&
                          + 5./6.*y(iin,jdim-1,kin)+0.25*y(iin,jdim,kin)
          mm(3,2) = -1./12.*z(iin,jdim-4,kin)+0.5*z(iin,jdim-3,kin)-1.5*z(iin,jdim-2,kin)&
                          + 5./6.*z(iin,jdim-1,kin)+0.25*z(iin,jdim,kin)
       elseif(jin.eq.jdim) then
          mm(1,2) = 0.25*x(iin,jdim-4,kin)-4./3.*x(iin,jdim-3,kin)+3.*x(iin,jdim-2,kin)&
                        - 4.*x(iin,jdim-1,kin)+25./12.*x(iin,jdim,kin)
          mm(2,2) = 0.25*y(iin,jdim-4,kin)-4./3.*y(iin,jdim-3,kin)+3.*y(iin,jdim-2,kin)&
                        - 4.*y(iin,jdim-1,kin)+25./12.*y(iin,jdim,kin)
          mm(3,2) = 0.25*z(iin,jdim-4,kin)-4./3.*z(iin,jdim-3,kin)+3.*z(iin,jdim-2,kin)&
                        - 4.*z(iin,jdim-1,kin)+25./12.*z(iin,jdim,kin)
       endif
    elseif(jdim.eq.3.or.jdim.eq.4) then
       if(jin.ge.2.and.jin.le.3) then
          mm(1,2) = -0.5*x(iin,jin-1,kin) + 0.5*x(iin,jin+1,kin)
          mm(2,2) = -0.5*y(iin,jin-1,kin) + 0.5*y(iin,jin+1,kin)
          mm(3,2) = -0.5*z(iin,jin-1,kin) + 0.5*z(iin,jin+1,kin)
       elseif(jin.eq.1) then        
          mm(1,2) = -1.5*x(iin,1,kin) + 2.0*x(iin,2,kin) - 0.5*x(iin,3,kin)
          mm(2,2) = -1.5*y(iin,1,kin) + 2.0*y(iin,2,kin) - 0.5*y(iin,3,kin)
          mm(3,2) = -1.5*z(iin,1,kin) + 2.0*z(iin,2,kin) - 0.5*z(iin,3,kin)
       elseif(jin.eq.jdim) then
          mm(1,2) = 0.5*x(iin,jdim-2,kin) - 2.0*x(iin,jdim-1,kin) + 1.5*x(iin,jdim,kin)
          mm(2,2) = 0.5*y(iin,jdim-2,kin) - 2.0*y(iin,jdim-1,kin) + 1.5*y(iin,jdim,kin)
          mm(3,2) = 0.5*z(iin,jdim-2,kin) - 2.0*z(iin,jdim-1,kin) + 1.5*z(iin,jdim,kin)
       endif
    elseif(jdim.eq.2) then
          mm(1,2) = x(iin,jdim,kin) - x(iin,jdim-1,kin)
          mm(2,2) = y(iin,jdim,kin) - y(iin,jdim-1,kin)
          mm(3,2) = z(iin,jdim,kin) - z(iin,jdim-1,kin)
    elseif(jdim.eq.1) then
          mm(1,2) = 0.0
          mm(2,2) = 0.0
          mm(3,2) = 0.0
    endif



    if(kdim.ge.5) then
       if(kin.ge.3.and.kin.le.kdim-2) then
          mm(1,3) = (x(iin,jin,kin-2)-8.*(x(iin,jin,kin-1)-x(iin,jin,kin+1))-x(iin,jin,kin+2))*c12i
          mm(2,3) = (y(iin,jin,kin-2)-8.*(y(iin,jin,kin-1)-y(iin,jin,kin+1))-y(iin,jin,kin+2))*c12i
          mm(3,3) = (z(iin,jin,kin-2)-8.*(z(iin,jin,kin-1)-z(iin,jin,kin+1))-z(iin,jin,kin+2))*c12i
       elseif(kin.eq.1) then
          mm(1,3) = -25./12.*x(iin,jin,1)+4.*x(iin,jin,2)-3.*x(iin,jin,3)+4./3.*x(iin,jin,4)-0.25*x(iin,jin,5)
          mm(2,3) = -25./12.*y(iin,jin,1)+4.*y(iin,jin,2)-3.*y(iin,jin,3)+4./3.*y(iin,jin,4)-0.25*y(iin,jin,5)
          mm(3,3) = -25./12.*z(iin,jin,1)+4.*z(iin,jin,2)-3.*z(iin,jin,3)+4./3.*z(iin,jin,4)-0.25*z(iin,jin,5)
       elseif(kin.eq.2) then
          mm(1,3) = -0.25*x(iin,jin,1)-5./6.*x(iin,jin,2)+1.5*x(iin,jin,3)-0.5*x(iin,jin,4)+1./12.*x(iin,jin,5)
          mm(2,3) = -0.25*y(iin,jin,1)-5./6.*y(iin,jin,2)+1.5*y(iin,jin,3)-0.5*y(iin,jin,4)+1./12.*y(iin,jin,5)
          mm(3,3) = -0.25*z(iin,jin,1)-5./6.*z(iin,jin,2)+1.5*z(iin,jin,3)-0.5*z(iin,jin,4)+1./12.*z(iin,jin,5)
       elseif(kin.eq.kdim-1) then
          mm(1,3) = -1./12.*x(iin,jin,kdim-4)+0.5*x(iin,jin,kdim-3)-1.5*x(iin,jin,kdim-2)&
                          + 5./6.*x(iin,jin,kdim-1)+0.25*x(iin,jin,kdim)
          mm(2,3) = -1./12.*y(iin,jin,kdim-4)+0.5*y(iin,jin,kdim-3)-1.5*y(iin,jin,kdim-2)&
                          + 5./6.*y(iin,jin,kdim-1)+0.25*y(iin,jin,kdim)
          mm(3,3) = -1./12.*z(iin,jin,kdim-4)+0.5*z(iin,jin,kdim-3)-1.5*z(iin,jin,kdim-2)&
                          + 5./6.*z(iin,jin,kdim-1)+0.25*z(iin,jin,kdim)
       elseif(kin.eq.kdim) then
          mm(1,3) = 0.25*x(iin,jin,kdim-4)-4./3.*x(iin,jin,kdim-3)+3.*x(iin,jin,kdim-2)&
                        - 4.*x(iin,jin,kdim-1)+25./12.*x(iin,jin,kdim)
          mm(2,3) = 0.25*y(iin,jin,kdim-4)-4./3.*y(iin,jin,kdim-3)+3.*y(iin,jin,kdim-2)&
                        - 4.*y(iin,jin,kdim-1)+25./12.*y(iin,jin,kdim)
          mm(3,3) = 0.25*z(iin,jin,kdim-4)-4./3.*z(iin,jin,kdim-3)+3.*z(iin,jin,kdim-2)&
                        - 4.*z(iin,jin,kdim-1)+25./12.*z(iin,jin,kdim)
       endif
    elseif(kdim.eq.3.or.kdim.eq.4) then
       if(kin.ge.2.and.kin.le.3) then
          mm(1,3) = -0.5*x(iin,jin,kin-1) + 0.5*x(iin,jin,kin+1)
          mm(2,3) = -0.5*y(iin,jin,kin-1) + 0.5*y(iin,jin,kin+1)
          mm(3,3) = -0.5*z(iin,jin,kin-1) + 0.5*z(iin,jin,kin+1)
       elseif(kin.eq.1) then        
          mm(1,3) = -1.5*x(iin,jin,1) + 2.0*x(iin,jin,2) - 0.5*x(iin,jin,3)
          mm(2,3) = -1.5*y(iin,jin,1) + 2.0*y(iin,jin,2) - 0.5*y(iin,jin,3)
          mm(3,3) = -1.5*z(iin,jin,1) + 2.0*z(iin,jin,2) - 0.5*z(iin,jin,3)
       elseif(kin.eq.kdim) then
          mm(1,3) = 0.5*x(iin,jin,kdim-2) - 2.0*x(iin,jin,kdim-1) + 1.5*x(iin,jin,kdim)
          mm(2,3) = 0.5*y(iin,jin,kdim-2) - 2.0*y(iin,jin,kdim-1) + 1.5*y(iin,jin,kdim)
          mm(3,3) = 0.5*z(iin,jin,kdim-2) - 2.0*z(iin,jin,kdim-1) + 1.5*z(iin,jin,kdim)
       endif
    elseif(kdim.eq.2) then
          mm(1,3) = x(iin,jin,kdim) - x(iin,jin,kdim-1)
          mm(2,3) = y(iin,jin,kdim) - y(iin,jin,kdim-1)
          mm(3,3) = z(iin,jin,kdim) - z(iin,jin,kdim-1)
    elseif(kdim.eq.1) then
          mm(1,3) = 0.0
          mm(2,3) = 0.0
          mm(3,3) = 0.0
    endif

    Calmm = mm
    end function Calmm

    ! x(kmax,imax,jmax)
    function CalmmShift(iin,jin,kin, x, y, z)

    !  mm(3,3)
    !  mm(1,1) = dxdi, mm(1,2) = dxdj, mm(1,3) = dxdk
    !  mm(2,1) = dydi, mm(2,2) = dydj, mm(2,3) = dydk
    !  mm(3,1) = dzdi, mm(3,2) = dzdj, mm(3,3) = dzdk

    integer, intent(in) :: iin, jin, kin
    real(8), intent(in) :: x(:,:,:), y(:,:,:), z(:,:,:)
    real(8) :: CalmmShift(3,3)
    integer :: idim, jdim, kdim
    real(8) :: mm(3,3)

    idim = size(x,dim=2)
    jdim = size(x,dim=3)
    kdim = size(x,dim=1)

    if(idim.ge.5) then
       if(iin.ge.3.and.iin.le.idim-2) then
          mm(1,1) = (x(kin,iin-2,jin)-8.*(x(kin,iin-1,jin)-x(kin,iin+1,jin))-x(kin,iin+2,jin))*c12i
          mm(2,1) = (y(kin,iin-2,jin)-8.*(y(kin,iin-1,jin)-y(kin,iin+1,jin))-y(kin,iin+2,jin))*c12i
          mm(3,1) = (z(kin,iin-2,jin)-8.*(z(kin,iin-1,jin)-z(kin,iin+1,jin))-z(kin,iin+2,jin))*c12i
       elseif(iin.eq.1) then
          mm(1,1) = -25./12.*x(kin,1,jin)+4.*x(kin,2,jin)-3.*x(kin,3,jin)+4./3.*x(kin,4,jin)-0.25*x(kin,5,jin)
          mm(2,1) = -25./12.*y(kin,1,jin)+4.*y(kin,2,jin)-3.*y(kin,3,jin)+4./3.*y(kin,4,jin)-0.25*y(kin,5,jin)
          mm(3,1) = -25./12.*z(kin,1,jin)+4.*z(kin,2,jin)-3.*z(kin,3,jin)+4./3.*z(kin,4,jin)-0.25*z(kin,5,jin)
       elseif(iin.eq.2) then
          mm(1,1) = -0.25*x(kin,1,jin)-5./6.*x(kin,2,jin)+1.5*x(kin,3,jin)-0.5*x(kin,4,jin)+1./12.*x(kin,5,jin)
          mm(2,1) = -0.25*y(kin,1,jin)-5./6.*y(kin,2,jin)+1.5*y(kin,3,jin)-0.5*y(kin,4,jin)+1./12.*y(kin,5,jin)
          mm(3,1) = -0.25*z(kin,1,jin)-5./6.*z(kin,2,jin)+1.5*z(kin,3,jin)-0.5*z(kin,4,jin)+1./12.*z(kin,5,jin)
       elseif(iin.eq.idim-1) then
          mm(1,1) = -1./12.*x(kin,idim-4,jin)+0.5*x(kin,idim-3,jin)-1.5*x(kin,idim-2,jin)&
                          + 5./6.*x(kin,idim-1,jin)+0.25*x(kin,idim,jin)
          mm(2,1) = -1./12.*y(kin,idim-4,jin)+0.5*y(kin,idim-3,jin)-1.5*y(kin,idim-2,jin)&
                          + 5./6.*y(kin,idim-1,jin)+0.25*y(kin,idim,jin)
          mm(3,1) = -1./12.*z(kin,idim-4,jin)+0.5*z(kin,idim-3,jin)-1.5*z(kin,idim-2,jin)&
                          + 5./6.*z(kin,idim-1,jin)+0.25*z(kin,idim,jin)
       elseif(iin.eq.idim) then
          mm(1,1) = 0.25*x(kin,idim-4,jin)-4./3.*x(kin,idim-3,jin)+3.*x(kin,idim-2,jin)&
                         - 4.*x(kin,idim-1,jin)+25./12.*x(kin,idim,jin)
          mm(2,1) = 0.25*y(kin,idim-4,jin)-4./3.*y(kin,idim-3,jin)+3.*y(kin,idim-2,jin)&
                         - 4.*y(kin,idim-1,jin)+25./12.*y(kin,idim,jin)
          mm(3,1) = 0.25*z(kin,idim-4,jin)-4./3.*z(kin,idim-3,jin)+3.*z(kin,idim-2,jin)&
                         - 4.*z(kin,idim-1,jin)+25./12.*z(kin,idim,jin)
       endif ! end if(iin)
    elseif(idim.eq.3.or.idim.eq.4) then
       if(iin.ge.2.and.iin.le.3) then
          mm(1,1) = -0.5*x(kin,iin-1,jin) + 0.5*x(kin,iin+1,jin)
          mm(2,1) = -0.5*y(kin,iin-1,jin) + 0.5*y(kin,iin+1,jin)
          mm(3,1) = -0.5*z(kin,iin-1,jin) + 0.5*z(kin,iin+1,jin)
       elseif(iin.eq.1) then
          mm(1,1) = -1.5*x(kin,1,jin) + 2.0*x(kin,2,jin) - 0.5*x(kin,3,jin)
          mm(2,1) = -1.5*y(kin,1,jin) + 2.0*y(kin,2,jin) - 0.5*y(kin,3,jin)
          mm(3,1) = -1.5*z(kin,1,jin) + 2.0*z(kin,2,jin) - 0.5*z(kin,3,jin)
       elseif(iin.eq.idim) then
          mm(1,1) = 0.5*x(kin,idim-2,jin) - 2.0*x(kin,idim-1,jin) + 1.5*x(kin,idim,jin)
          mm(2,1) = 0.5*y(kin,idim-2,jin) - 2.0*y(kin,idim-1,jin) + 1.5*y(kin,idim,jin)
          mm(3,1) = 0.5*z(kin,idim-2,jin) - 2.0*z(kin,idim-1,jin) + 1.5*z(kin,idim,jin)
       endif
    elseif(idim.eq.2) then
          mm(1,1) = x(kin,idim,jin) - x(kin,idim-1,jin)
          mm(2,1) = y(kin,idim,jin) - y(kin,idim-1,jin)
          mm(3,1) = z(kin,idim,jin) - z(kin,idim-1,jin)
    elseif(idim.eq.1) then
          mm(1,1) = 0.0
          mm(2,1) = 0.0
          mm(3,1) = 0.0
    endif

    if(jdim.ge.5) then
       if(jin.ge.3.and.jin.le.jdim-2) then
          mm(1,2) = (x(kin,iin,jin-2)-8.*(x(kin,iin,jin-1)-x(kin,iin,jin+1))-x(kin,iin,jin+2))*c12i
          mm(2,2) = (y(kin,iin,jin-2)-8.*(y(kin,iin,jin-1)-y(kin,iin,jin+1))-y(kin,iin,jin+2))*c12i
          mm(3,2) = (z(kin,iin,jin-2)-8.*(z(kin,iin,jin-1)-z(kin,iin,jin+1))-z(kin,iin,jin+2))*c12i
       elseif(jin.eq.1) then
          mm(1,2) = -25./12.*x(kin,iin,1)+4.*x(kin,iin,2)-3.*x(kin,iin,3)+4./3.*x(kin,iin,4)-0.25*x(kin,iin,5)
          mm(2,2) = -25./12.*y(kin,iin,1)+4.*y(kin,iin,2)-3.*y(kin,iin,3)+4./3.*y(kin,iin,4)-0.25*y(kin,iin,5)
          mm(3,2) = -25./12.*z(kin,iin,1)+4.*z(kin,iin,2)-3.*z(kin,iin,3)+4./3.*z(kin,iin,4)-0.25*z(kin,iin,5)
       elseif(jin.eq.2) then
          mm(1,2) = -0.25*x(kin,iin,1)-5./6.*x(kin,iin,2)+1.5*x(kin,iin,3)-0.5*x(kin,iin,4)+1./12.*x(kin,iin,5)
          mm(2,2) = -0.25*y(kin,iin,1)-5./6.*y(kin,iin,2)+1.5*y(kin,iin,3)-0.5*y(kin,iin,4)+1./12.*y(kin,iin,5)
          mm(3,2) = -0.25*z(kin,iin,1)-5./6.*z(kin,iin,2)+1.5*z(kin,iin,3)-0.5*z(kin,iin,4)+1./12.*z(kin,iin,5)
       elseif(jin.eq.jdim-1) then
          mm(1,2) = -1./12.*x(kin,iin,jdim-4)+0.5*x(kin,iin,jdim-3)-1.5*x(kin,iin,jdim-2)&
                          + 5./6.*x(kin,iin,jdim-1)+0.25*x(kin,iin,jdim)
          mm(2,2) = -1./12.*y(kin,iin,jdim-4)+0.5*y(kin,iin,jdim-3)-1.5*y(kin,iin,jdim-2)&
                          + 5./6.*y(kin,iin,jdim-1)+0.25*y(kin,iin,jdim)
          mm(3,2) = -1./12.*z(kin,iin,jdim-4)+0.5*z(kin,iin,jdim-3)-1.5*z(kin,iin,jdim-2)&
                          + 5./6.*z(kin,iin,jdim-1)+0.25*z(kin,iin,jdim)
       elseif(jin.eq.jdim) then
          mm(1,2) = 0.25*x(kin,iin,jdim-4)-4./3.*x(kin,iin,jdim-3)+3.*x(kin,iin,jdim-2)&
                        - 4.*x(kin,iin,jdim-1)+25./12.*x(kin,iin,jdim)
          mm(2,2) = 0.25*y(kin,iin,jdim-4)-4./3.*y(kin,iin,jdim-3)+3.*y(kin,iin,jdim-2)&
                        - 4.*y(kin,iin,jdim-1)+25./12.*y(kin,iin,jdim)
          mm(3,2) = 0.25*z(kin,iin,jdim-4)-4./3.*z(kin,iin,jdim-3)+3.*z(kin,iin,jdim-2)&
                        - 4.*z(kin,iin,jdim-1)+25./12.*z(kin,iin,jdim)
       endif
    elseif(jdim.eq.3.or.jdim.eq.4) then
       if(jin.ge.2.and.jin.le.3) then
          mm(1,2) = -0.5*x(kin,iin,jin-1) + 0.5*x(kin,iin,jin+1)
          mm(2,2) = -0.5*y(kin,iin,jin-1) + 0.5*y(kin,iin,jin+1)
          mm(3,2) = -0.5*z(kin,iin,jin-1) + 0.5*z(kin,iin,jin+1)
       elseif(jin.eq.1) then
          mm(1,2) = -1.5*x(kin,iin,1) + 2.0*x(kin,iin,2) - 0.5*x(kin,iin,3)
          mm(2,2) = -1.5*y(kin,iin,1) + 2.0*y(kin,iin,2) - 0.5*y(kin,iin,3)
          mm(3,2) = -1.5*z(kin,iin,1) + 2.0*z(kin,iin,2) - 0.5*z(kin,iin,3)
       elseif(jin.eq.jdim) then
          mm(1,2) = 0.5*x(kin,iin,jdim-2) - 2.0*x(kin,iin,jdim-1) + 1.5*x(kin,iin,jdim)
          mm(2,2) = 0.5*y(kin,iin,jdim-2) - 2.0*y(kin,iin,jdim-1) + 1.5*y(kin,iin,jdim)
          mm(3,2) = 0.5*z(kin,iin,jdim-2) - 2.0*z(kin,iin,jdim-1) + 1.5*z(kin,iin,jdim)
       endif
    elseif(jdim.eq.2) then
          mm(1,2) = x(kin,iin,jdim) - x(kin,iin,jdim-1)
          mm(2,2) = y(kin,iin,jdim) - y(kin,iin,jdim-1)
          mm(3,2) = z(kin,iin,jdim) - z(kin,iin,jdim-1)
    elseif(jdim.eq.1) then
          mm(1,2) = 0.0
          mm(2,2) = 0.0
          mm(3,2) = 0.0
    endif



    if(kdim.ge.5) then
       if(kin.ge.3.and.kin.le.kdim-2) then
          mm(1,3) = (x(kin-2,iin,jin)-8.*(x(kin-1,iin,jin)-x(kin+1,iin,jin))-x(kin+2,iin,jin))*c12i
          mm(2,3) = (y(kin-2,iin,jin)-8.*(y(kin-1,iin,jin)-y(kin+1,iin,jin))-y(kin+2,iin,jin))*c12i
          mm(3,3) = (z(kin-2,iin,jin)-8.*(z(kin-1,iin,jin)-z(kin+1,iin,jin))-z(kin+2,iin,jin))*c12i
       elseif(kin.eq.1) then
          mm(1,3) = -25./12.*x(1,iin,jin)+4.*x(2,iin,jin)-3.*x(3,iin,jin)+4./3.*x(4,iin,jin)-0.25*x(5,iin,jin)
          mm(2,3) = -25./12.*y(1,iin,jin)+4.*y(2,iin,jin)-3.*y(3,iin,jin)+4./3.*y(4,iin,jin)-0.25*y(5,iin,jin)
          mm(3,3) = -25./12.*z(1,iin,jin)+4.*z(2,iin,jin)-3.*z(3,iin,jin)+4./3.*z(4,iin,jin)-0.25*z(5,iin,jin)
       elseif(kin.eq.2) then
          mm(1,3) = -0.25*x(1,iin,jin)-5./6.*x(2,iin,jin)+1.5*x(3,iin,jin)-0.5*x(4,iin,jin)+1./12.*x(5,iin,jin)
          mm(2,3) = -0.25*y(1,iin,jin)-5./6.*y(2,iin,jin)+1.5*y(3,iin,jin)-0.5*y(4,iin,jin)+1./12.*y(5,iin,jin)
          mm(3,3) = -0.25*z(1,iin,jin)-5./6.*z(2,iin,jin)+1.5*z(3,iin,jin)-0.5*z(4,iin,jin)+1./12.*z(5,iin,jin)
       elseif(kin.eq.kdim-1) then
          mm(1,3) = -1./12.*x(kdim-4,iin,jin)+0.5*x(kdim-3,iin,jin)-1.5*x(kdim-2,iin,jin)&
                          + 5./6.*x(kdim-1,iin,jin)+0.25*x(kdim,iin,jin)
          mm(2,3) = -1./12.*y(kdim-4,iin,jin)+0.5*y(kdim-3,iin,jin)-1.5*y(kdim-2,iin,jin)&
                          + 5./6.*y(kdim-1,iin,jin)+0.25*y(kdim,iin,jin)
          mm(3,3) = -1./12.*z(kdim-4,iin,jin)+0.5*z(kdim-3,iin,jin)-1.5*z(kdim-2,iin,jin)&
                          + 5./6.*z(kdim-1,iin,jin)+0.25*z(kdim,iin,jin)
       elseif(kin.eq.kdim) then
          mm(1,3) = 0.25*x(kdim-4,iin,jin)-4./3.*x(kdim-3,iin,jin)+3.*x(kdim-2,iin,jin)&
                        - 4.*x(kdim-1,iin,jin)+25./12.*x(kdim,iin,jin)
          mm(2,3) = 0.25*y(kdim-4,iin,jin)-4./3.*y(kdim-3,iin,jin)+3.*y(kdim-2,iin,jin)&
                        - 4.*y(kdim-1,iin,jin)+25./12.*y(kdim,iin,jin)
          mm(3,3) = 0.25*z(kdim-4,iin,jin)-4./3.*z(kdim-3,iin,jin)+3.*z(kdim-2,iin,jin)&
                        - 4.*z(kdim-1,iin,jin)+25./12.*z(kdim,iin,jin)
       endif
    elseif(kdim.eq.3.or.kdim.eq.4) then
       if(kin.ge.2.and.kin.le.3) then
          mm(1,3) = -0.5*x(kin-1,iin,jin) + 0.5*x(kin+1,iin,jin)
          mm(2,3) = -0.5*y(kin-1,iin,jin) + 0.5*y(kin+1,iin,jin)
          mm(3,3) = -0.5*z(kin-1,iin,jin) + 0.5*z(kin+1,iin,jin)
       elseif(kin.eq.1) then
          mm(1,3) = -1.5*x(1,iin,jin) + 2.0*x(2,iin,jin) - 0.5*x(3,iin,jin)
          mm(2,3) = -1.5*y(1,iin,jin) + 2.0*y(2,iin,jin) - 0.5*y(3,iin,jin)
          mm(3,3) = -1.5*z(1,iin,jin) + 2.0*z(2,iin,jin) - 0.5*z(3,iin,jin)
       elseif(kin.eq.kdim) then
          mm(1,3) = 0.5*x(kdim-2,iin,jin) - 2.0*x(kdim-1,iin,jin) + 1.5*x(kdim,iin,jin)
          mm(2,3) = 0.5*y(kdim-2,iin,jin) - 2.0*y(kdim-1,iin,jin) + 1.5*y(kdim,iin,jin)
          mm(3,3) = 0.5*z(kdim-2,iin,jin) - 2.0*z(kdim-1,iin,jin) + 1.5*z(kdim,iin,jin)
       endif
    elseif(kdim.eq.2) then
          mm(1,3) = x(kdim,iin,jin) - x(kdim-1,iin,jin)
          mm(2,3) = y(kdim,iin,jin) - y(kdim-1,iin,jin)
          mm(3,3) = z(kdim,iin,jin) - z(kdim-1,iin,jin)
    elseif(kdim.eq.1) then
          mm(1,3) = 0.0
          mm(2,3) = 0.0
          mm(3,3) = 0.0
    endif

    CalmmShift = mm
    end function CalmmShift

    function Calmet_2D(mm)
    ! Calculate mesh metrics met
    ! met(2,2)
    ! met(1,1) = didx, met(1,2) = didy
    ! met(2,1) = djdx, met(2,2) = djdy
    !
    ! given
    !
    !  mm(2,2)
    !  mm(1,1) = dxdi, mm(1,2) = dxdj
    !  mm(2,1) = dydi, mm(2,2) = dydj

    real(8), intent(in) :: mm(2,2)
    real(8) :: met(2,2), Calmet_2D(2,2)
    real(8) :: vkapx, vkapy, vol

     ! Calculate didx, didy
     vkapx =   mm(2,2)
     vkapy = - mm(1,2)

     vol = + mm(1,1)*vkapx + mm(2,1)*vkapy

     met(1,1) = +vkapx/(vol+1.e-30)
     met(1,2) = +vkapy/(vol+1.e-30)


     ! Calculate djdx, djdy, djdz
     vkapx = - mm(2,1)
     vkapy =   mm(1,1)

     vol = mm(1,2)*vkapx + mm(2,2)*vkapy

     met(2,1) = +vkapx/(vol+1.e-30)
     met(2,2) = +vkapy/(vol+1.e-30)


     Calmet_2D = met

     end function Calmet_2D



    function Calmet(mm)
    ! Calculate mesh metrics met
    ! met(3,3)
    ! met(1,1) = didx, met(1,2) = didy, met(1,3) = didz
    ! met(2,1) = djdx, met(2,2) = djdy, met(2,3) = djdz
    ! met(3,1) = dkdx, met(3,2) = dkdy, met(3,3) = dkdz
    !
    ! given 
    !
    !  mm(3,3)
    !  mm(1,1) = dxdi, mm(1,2) = dxdj, mm(1,3) = dxdk
    !  mm(2,1) = dydi, mm(2,2) = dydj, mm(2,3) = dydk
    !  mm(3,1) = dzdi, mm(3,2) = dzdj, mm(3,3) = dzdk

    real(8), intent(in) :: mm(3,3)
    real(8) :: met(3,3), Calmet(3,3)
    real(8) :: vkapx, vkapy, vkapz, vol

     ! Calculate didx, didy, didz
     vkapx = mm(2,2)*mm(3,3) - mm(2,3)*mm(3,2)
     vkapy = mm(1,2)*mm(3,3) - mm(1,3)*mm(3,2)
     vkapz = mm(1,2)*mm(2,3) - mm(1,3)*mm(2,2)
     vol = + mm(1,1)*vkapx - mm(2,1)*vkapy + mm(3,1)*vkapz

     met(1,1) = +vkapx/(vol+1.e-30)
     met(1,2) = -vkapy/(vol+1.e-30)
     met(1,3) = +vkapz/(vol+1.e-30)
     
     ! Calculate djdx, djdy, djdz
     vkapx = mm(2,1)*mm(3,3) - mm(2,3)*mm(3,1)
     vkapy = mm(1,1)*mm(3,3) - mm(1,3)*mm(3,1)
     vkapz = mm(1,1)*mm(2,3) - mm(1,3)*mm(2,1)
     vol = - mm(1,2)*vkapx + mm(2,2)*vkapy - mm(3,2)*vkapz

     met(2,1) = -vkapx/(vol+1.e-30)
     met(2,2) = +vkapy/(vol+1.e-30)
     met(2,3) = -vkapz/(vol+1.e-30)

     ! Calculate dkdx, dkdy, dkdz
     vkapx = mm(2,1)*mm(3,2) - mm(2,2)*mm(3,1)
     vkapy = mm(1,1)*mm(3,2) - mm(1,2)*mm(3,1)
     vkapz = mm(1,1)*mm(2,2) - mm(1,2)*mm(2,1)
     vol = + mm(1,3)*vkapx - mm(2,3)*vkapy + mm(3,3)*vkapz

     met(3,1) = +vkapx/(vol+1.e-30)
     met(3,2) = -vkapy/(vol+1.e-30)
     met(3,3) = +vkapz/(vol+1.e-30)

     Calmet = met

!    For debug purposes
!     print *, 'Identity mm*met =', matmul(mm,met)

     end function Calmet







    function planeCalmm(idim,jdim,kdim, x, y, z)

    !  mm(3,3,:,:,:)
    !  mm(1,1,:,:,:) = dxdi, mm(1,2,:,:,:) = dxdj, mm(1,3,:,:,:) = dxdk
    !  mm(2,1,:,:,:) = dydi, mm(2,2,:,:,:) = dydj, mm(2,3,:,:,:) = dydk
    !  mm(3,1,:,:,:) = dzdi, mm(3,2,:,:,:) = dzdj, mm(3,3,:,:,:) = dzdk

    integer, intent(in) :: idim, jdim, kdim
    real(8), intent(in) :: x(:,:,:), y(:,:,:), z(:,:,:)
    integer :: i, j, k
    real(8),dimension(:,:,:,:,:), allocatable :: mm,planeCalmm

    allocate(mm(3,3,idim,jdim,kdim))
    allocate(planeCalmm(3,3,idim,jdim,kdim))

    if(idim.ge.5) then
      Do k=1,kdim
        Do j=1,jdim
          mm(1,1,1,j,k) = -25./12.*x(1,j,k)+4.*x(2,j,k)-3.*x(3,j,k)+4./3.*x(4,j,k)-0.25*x(5,j,k)
          mm(2,1,1,j,k) = -25./12.*y(1,j,k)+4.*y(2,j,k)-3.*y(3,j,k)+4./3.*y(4,j,k)-0.25*y(5,j,k)
          mm(3,1,1,j,k) = -25./12.*z(1,j,k)+4.*z(2,j,k)-3.*z(3,j,k)+4./3.*z(4,j,k)-0.25*z(5,j,k)

          mm(1,1,2,j,k) = -0.25*x(1,j,k)-5./6.*x(2,j,k)+1.5*x(3,j,k)-0.5*x(4,j,k)+1./12.*x(5,j,k)
          mm(2,1,2,j,k) = -0.25*y(1,j,k)-5./6.*y(2,j,k)+1.5*y(3,j,k)-0.5*y(4,j,k)+1./12.*y(5,j,k)
          mm(3,1,2,j,k) = -0.25*z(1,j,k)-5./6.*z(2,j,k)+1.5*z(3,j,k)-0.5*z(4,j,k)+1./12.*z(5,j,k)

          mm(1,1,idim-1,j,k) = -1./12.*x(idim-4,j,k)+0.5*x(idim-3,j,k)-1.5*x(idim-2,j,k)&
                          + 5./6.*x(idim-1,j,k)+0.25*x(idim,j,k)
          mm(2,1,idim-1,j,k) = -1./12.*y(idim-4,j,k)+0.5*y(idim-3,j,k)-1.5*y(idim-2,j,k)&
                          + 5./6.*y(idim-1,j,k)+0.25*y(idim,j,k)
          mm(3,1,idim-1,j,k) = -1./12.*z(idim-4,j,k)+0.5*z(idim-3,j,k)-1.5*z(idim-2,j,k)&
                          + 5./6.*z(idim-1,j,k)+0.25*z(idim,j,k)

          mm(1,1,idim,j,k) = 0.25*x(idim-4,j,k)-4./3.*x(idim-3,j,k)+3.*x(idim-2,j,k)&
                         - 4.*x(idim-1,j,k)+25./12.*x(idim,j,k)
          mm(2,1,idim,j,k) = 0.25*y(idim-4,j,k)-4./3.*y(idim-3,j,k)+3.*y(idim-2,j,k)&
                         - 4.*y(idim-1,j,k)+25./12.*y(idim,j,k)
          mm(3,1,idim,j,k) = 0.25*z(idim-4,j,k)-4./3.*z(idim-3,j,k)+3.*z(idim-2,j,k)&
                         - 4.*z(idim-1,j,k)+25./12.*z(idim,j,k)

          Do i=3,idim-2
            mm(1,1,i,j,k) = (x(i-2,j,k)-8.*(x(i-1,j,k)-x(i+1,j,k))-x(i+2,j,k))*c12i
            mm(2,1,i,j,k) = (y(i-2,j,k)-8.*(y(i-1,j,k)-y(i+1,j,k))-y(i+2,j,k))*c12i
            mm(3,1,i,j,k) = (z(i-2,j,k)-8.*(z(i-1,j,k)-z(i+1,j,k))-z(i+2,j,k))*c12i
          End Do
        End Do
       End Do

!    elseif(idim.eq.3.or.idim.eq.4) then
!      Do k=1,kdim
!        Do j=1,jdim
!          mm(1,1,1,j,k) = -1.5*x(1,j,k) + 2.0*x(2,j,k) - 0.5*x(3,j,k)
!          mm(2,1,1,j,k) = -1.5*y(1,j,k) + 2.0*y(2,j,k) - 0.5*y(3,j,k)
!          mm(3,1,1,j,k) = -1.5*z(1,j,k) + 2.0*z(2,j,k) - 0.5*z(3,j,k)

!          mm(1,1,idim,j,k) = 0.5*x(idim-2,j,k) - 2.0*x(idim-1,j,k) + 1.5*x(idim,j,k)
!          mm(2,1,idim,j,k) = 0.5*y(idim-2,j,k) - 2.0*y(idim-1,j,k) + 1.5*y(idim,j,k)
!          mm(3,1,idim,j,k) = 0.5*z(idim-2,j,k) - 2.0*z(idim-1,j,k) + 1.5*z(idim,j,k)

!          Do i=2,3
!          mm(1,1,i,j,k) = -0.5*x(i-1,j,k) + 0.5*x(i+1,j,k)
!          mm(2,1,i,j,k) = -0.5*y(i-1,j,k) + 0.5*y(i+1,j,k)
!          mm(3,1,i,j,k) = -0.5*z(i-1,j,k) + 0.5*z(i+1,j,k)
!          End Do


!    elseif(idim.eq.2) then




 !   elseif(idim.eq.1) then




     End If




!    if(idim.ge.5) then
!       if(iin.ge.3.and.iin.le.idim-2) then
!          mm(1,1) = (x(iin-2,jin,kin)-8.*(x(iin-1,jin,kin)-x(iin+1,jin,kin))-x(iin+2,jin,kin))*c12i
!          mm(2,1) = (y(iin-2,jin,kin)-8.*(y(iin-1,jin,kin)-y(iin+1,jin,kin))-y(iin+2,jin,kin))*c12i
!          mm(3,1) = (z(iin-2,jin,kin)-8.*(z(iin-1,jin,kin)-z(iin+1,jin,kin))-z(iin+2,jin,kin))*c12i
!       elseif(iin.eq.1) then
!          mm(1,1) = -25./12.*x(1,jin,kin)+4.*x(2,jin,kin)-3.*x(3,jin,kin)+4./3.*x(4,jin,kin)-0.25*x(5,jin,kin)
!          mm(2,1) = -25./12.*y(1,jin,kin)+4.*y(2,jin,kin)-3.*y(3,jin,kin)+4./3.*y(4,jin,kin)-0.25*y(5,jin,kin)
!          mm(3,1) = -25./12.*z(1,jin,kin)+4.*z(2,jin,kin)-3.*z(3,jin,kin)+4./3.*z(4,jin,kin)-0.25*z(5,jin,kin)
!       elseif(iin.eq.2) then
!          mm(1,1) = -0.25*x(1,jin,kin)-5./6.*x(2,jin,kin)+1.5*x(3,jin,kin)-0.5*x(4,jin,kin)+1./12.*x(5,jin,kin)
!          mm(2,1) = -0.25*y(1,jin,kin)-5./6.*y(2,jin,kin)+1.5*y(3,jin,kin)-0.5*y(4,jin,kin)+1./12.*y(5,jin,kin)
!          mm(3,1) = -0.25*z(1,jin,kin)-5./6.*z(2,jin,kin)+1.5*z(3,jin,kin)-0.5*z(4,jin,kin)+1./12.*z(5,jin,kin)
!       elseif(iin.eq.idim-1) then
!          mm(1,1) = -1./12.*x(idim-4,jin,kin)+0.5*x(idim-3,jin,kin)-1.5*x(idim-2,jin,kin)&
!                          + 5./6.*x(idim-1,jin,kin)+0.25*x(idim,jin,kin)
!          mm(2,1) = -1./12.*y(idim-4,jin,kin)+0.5*y(idim-3,jin,kin)-1.5*y(idim-2,jin,kin)&
!                          + 5./6.*y(idim-1,jin,kin)+0.25*y(idim,jin,kin)
!          mm(3,1) = -1./12.*z(idim-4,jin,kin)+0.5*z(idim-3,jin,kin)-1.5*z(idim-2,jin,kin)&
!                          + 5./6.*z(idim-1,jin,kin)+0.25*z(idim,jin,kin)
!       elseif(iin.eq.idim) then
!          mm(1,1) = 0.25*x(idim-4,jin,kin)-4./3.*x(idim-3,jin,kin)+3.*x(idim-2,jin,kin)&
!                         - 4.*x(idim-1,jin,kin)+25./12.*x(idim,jin,kin)
!          mm(2,1) = 0.25*y(idim-4,jin,kin)-4./3.*y(idim-3,jin,kin)+3.*y(idim-2,jin,kin)&
!                         - 4.*y(idim-1,jin,kin)+25./12.*y(idim,jin,kin)
!          mm(3,1) = 0.25*z(idim-4,jin,kin)-4./3.*z(idim-3,jin,kin)+3.*z(idim-2,jin,kin)&
!                         - 4.*z(idim-1,jin,kin)+25./12.*z(idim,jin,kin)
!       endif ! end if(iin)
!    elseif(idim.eq.3.or.idim.eq.4) then
!       if(iin.ge.2.and.iin.le.3) then
!          mm(1,1) = -0.5*x(iin-1,jin,kin) + 0.5*x(iin+1,jin,kin)
!          mm(2,1) = -0.5*y(iin-1,jin,kin) + 0.5*y(iin+1,jin,kin)
 !         mm(3,1) = -0.5*z(iin-1,jin,kin) + 0.5*z(iin+1,jin,kin)
!       elseif(iin.eq.1) then
!          mm(1,1) = -1.5*x(1,jin,kin) + 2.0*x(2,jin,kin) - 0.5*x(3,jin,kin)
!          mm(2,1) = -1.5*y(1,jin,kin) + 2.0*y(2,jin,kin) - 0.5*y(3,jin,kin)
!          mm(3,1) = -1.5*z(1,jin,kin) + 2.0*z(2,jin,kin) - 0.5*z(3,jin,kin)
!       elseif(iin.eq.idim) then
!          mm(1,1) = 0.5*x(idim-2,jin,kin) - 2.0*x(idim-1,jin,kin) + 1.5*x(idim,jin,kin)
!          mm(2,1) = 0.5*y(idim-2,jin,kin) - 2.0*y(idim-1,jin,kin) + 1.5*y(idim,jin,kin)
!          mm(3,1) = 0.5*z(idim-2,jin,kin) - 2.0*z(idim-1,jin,kin) + 1.5*z(idim,jin,kin)
!       endif
 !   elseif(idim.eq.2) then
!          mm(1,1) = x(idim,jin,kin) - x(idim-1,jin,kin)
!          mm(2,1) = y(idim,jin,kin) - y(idim-1,jin,kin)
!          mm(3,1) = z(idim,jin,kin) - z(idim-1,jin,kin)
!    elseif(idim.eq.1) then
!          mm(1,1) = 0.0
!          mm(2,1) = 0.0
!          mm(3,1) = 0.0
!    endif


    if(jdim.ge.5) then
      Do k=1,kdim
        Do i=1,idim
          mm(1,2,i,1,k) = -25./12.*x(i,1,k)+4.*x(i,2,k)-3.*x(i,3,k)+4./3.*x(i,4,k)-0.25*x(i,5,k)
          mm(2,2,i,1,k) = -25./12.*y(i,1,k)+4.*y(i,2,k)-3.*y(i,3,k)+4./3.*y(i,4,k)-0.25*y(i,5,k)
          mm(3,2,i,1,k) = -25./12.*z(i,1,k)+4.*z(i,2,k)-3.*z(i,3,k)+4./3.*z(i,4,k)-0.25*z(i,5,k)

          mm(1,2,i,2,k) = -0.25*x(i,1,k)-5./6.*x(i,2,k)+1.5*x(i,3,k)-0.5*x(i,4,k)+1./12.*x(i,5,k)
          mm(2,2,i,2,k) = -0.25*y(i,1,k)-5./6.*y(i,2,k)+1.5*y(i,3,k)-0.5*y(i,4,k)+1./12.*y(i,5,k)
          mm(3,2,i,2,k) = -0.25*z(i,1,k)-5./6.*z(i,2,k)+1.5*z(i,3,k)-0.5*z(i,4,k)+1./12.*z(i,5,k)

          mm(1,2,i,jdim-1,k) = -1./12.*x(i,jdim-4,k)+0.5*x(i,jdim-3,k)-1.5*x(i,jdim-2,k)&
                          + 5./6.*x(i,jdim-1,k)+0.25*x(i,jdim,k)
          mm(2,2,i,jdim-1,k) = -1./12.*y(i,jdim-4,k)+0.5*y(i,jdim-3,k)-1.5*y(i,jdim-2,k)&
                          + 5./6.*y(i,jdim-1,k)+0.25*y(i,jdim,k)
          mm(3,2,i,jdim-1,k) = -1./12.*z(i,jdim-4,k)+0.5*z(i,jdim-3,k)-1.5*z(i,jdim-2,k)&
                          + 5./6.*z(i,jdim-1,k)+0.25*z(i,jdim,k)

          mm(1,2,i,jdim,k) = 0.25*x(i,jdim-4,k)-4./3.*x(i,jdim-3,k)+3.*x(i,jdim-2,k)&
                        - 4.*x(i,jdim-1,k)+25./12.*x(i,jdim,k)
          mm(2,2,i,jdim,k) = 0.25*y(i,jdim-4,k)-4./3.*y(i,jdim-3,k)+3.*y(i,jdim-2,k)&
                        - 4.*y(i,jdim-1,k)+25./12.*y(i,jdim,k)
          mm(3,2,i,jdim,k) = 0.25*z(i,jdim-4,k)-4./3.*z(i,jdim-3,k)+3.*z(i,jdim-2,k)&
                        - 4.*z(i,jdim-1,k)+25./12.*z(i,jdim,k)

          Do j=3,jdim-2
          mm(1,2,i,j,k) = (x(i,j-2,k)-8.*(x(i,j-1,k)-x(i,j+1,k))-x(i,j+2,k))*c12i
          mm(2,2,i,j,k) = (y(i,j-2,k)-8.*(y(i,j-1,k)-y(i,j+1,k))-y(i,j+2,k))*c12i
          mm(3,2,i,j,k) = (z(i,j-2,k)-8.*(z(i,j-1,k)-z(i,j+1,k))-z(i,j+2,k))*c12i
          End Do
        End Do
       End Do
     End If


!    if(jdim.ge.5) then
!       if(jin.ge.3.and.jin.le.jdim-2) then
!          mm(1,2) = (x(iin,jin-2,kin)-8.*(x(iin,jin-1,kin)-x(iin,jin+1,kin))-x(iin,jin+2,kin))*c12i
!          mm(2,2) = (y(iin,jin-2,kin)-8.*(y(iin,jin-1,kin)-y(iin,jin+1,kin))-y(iin,jin+2,kin))*c12i
!          mm(3,2) = (z(iin,jin-2,kin)-8.*(z(iin,jin-1,kin)-z(iin,jin+1,kin))-z(iin,jin+2,kin))*c12i
!       elseif(jin.eq.1) then
!          mm(1,2) = -25./12.*x(iin,1,kin)+4.*x(iin,2,kin)-3.*x(iin,3,kin)+4./3.*x(iin,4,kin)-0.25*x(iin,5,kin)
!          mm(2,2) = -25./12.*y(iin,1,kin)+4.*y(iin,2,kin)-3.*y(iin,3,kin)+4./3.*y(iin,4,kin)-0.25*y(iin,5,kin)
!          mm(3,2) = -25./12.*z(iin,1,kin)+4.*z(iin,2,kin)-3.*z(iin,3,kin)+4./3.*z(iin,4,kin)-0.25*z(iin,5,kin)
!       elseif(jin.eq.2) then
!          mm(1,2) = -0.25*x(iin,1,kin)-5./6.*x(iin,2,kin)+1.5*x(iin,3,kin)-0.5*x(iin,4,kin)+1./12.*x(iin,5,kin)
!          mm(2,2) = -0.25*y(iin,1,kin)-5./6.*y(iin,2,kin)+1.5*y(iin,3,kin)-0.5*y(iin,4,kin)+1./12.*y(iin,5,kin)
!          mm(3,2) = -0.25*z(iin,1,kin)-5./6.*z(iin,2,kin)+1.5*z(iin,3,kin)-0.5*z(iin,4,kin)+1./12.*z(iin,5,kin)
!       elseif(jin.eq.jdim-1) then
!          mm(1,2) = -1./12.*x(iin,jdim-4,kin)+0.5*x(iin,jdim-3,kin)-1.5*x(iin,jdim-2,kin)&
!                          + 5./6.*x(iin,jdim-1,kin)+0.25*x(iin,jdim,kin)
!          mm(2,2) = -1./12.*y(iin,jdim-4,kin)+0.5*y(iin,jdim-3,kin)-1.5*y(iin,jdim-2,kin)&
!                          + 5./6.*y(iin,jdim-1,kin)+0.25*y(iin,jdim,kin)
!          mm(3,2) = -1./12.*z(iin,jdim-4,kin)+0.5*z(iin,jdim-3,kin)-1.5*z(iin,jdim-2,kin)&
!                          + 5./6.*z(iin,jdim-1,kin)+0.25*z(iin,jdim,kin)
!       elseif(jin.eq.jdim) then
!          mm(1,2) = 0.25*x(iin,jdim-4,kin)-4./3.*x(iin,jdim-3,kin)+3.*x(iin,jdim-2,kin)&
!                        - 4.*x(iin,jdim-1,kin)+25./12.*x(iin,jdim,kin)
!          mm(2,2) = 0.25*y(iin,jdim-4,kin)-4./3.*y(iin,jdim-3,kin)+3.*y(iin,jdim-2,kin)&
!                        - 4.*y(iin,jdim-1,kin)+25./12.*y(iin,jdim,kin)
!          mm(3,2) = 0.25*z(iin,jdim-4,kin)-4./3.*z(iin,jdim-3,kin)+3.*z(iin,jdim-2,kin)&
!                        - 4.*z(iin,jdim-1,kin)+25./12.*z(iin,jdim,kin)
!       endif
!    elseif(jdim.eq.3.or.jdim.eq.4) then
!       if(jin.ge.2.and.jin.le.3) then
!          mm(1,2) = -0.5*x(iin,jin-1,kin) + 0.5*x(iin,jin+1,kin)
!          mm(2,2) = -0.5*y(iin,jin-1,kin) + 0.5*y(iin,jin+1,kin)
!          mm(3,2) = -0.5*z(iin,jin-1,kin) + 0.5*z(iin,jin+1,kin)
!       elseif(jin.eq.1) then
!          mm(1,2) = -1.5*x(iin,1,kin) + 2.0*x(iin,2,kin) - 0.5*x(iin,3,kin)
!          mm(2,2) = -1.5*y(iin,1,kin) + 2.0*y(iin,2,kin) - 0.5*y(iin,3,kin)
!          mm(3,2) = -1.5*z(iin,1,kin) + 2.0*z(iin,2,kin) - 0.5*z(iin,3,kin)
!       elseif(jin.eq.jdim) then
!          mm(1,2) = 0.5*x(iin,jdim-2,kin) - 2.0*x(iin,jdim-1,kin) + 1.5*x(iin,jdim,kin)
!         mm(2,2) = 0.5*y(iin,jdim-2,kin) - 2.0*y(iin,jdim-1,kin) + 1.5*y(iin,jdim,kin)
!          mm(3,2) = 0.5*z(iin,jdim-2,kin) - 2.0*z(iin,jdim-1,kin) + 1.5*z(iin,jdim,kin)
!       endif
!    elseif(jdim.eq.2) then
!          mm(1,2) = x(iin,jdim,kin) - x(iin,jdim-1,kin)
!          mm(2,2) = y(iin,jdim,kin) - y(iin,jdim-1,kin)
!          mm(3,2) = z(iin,jdim,kin) - z(iin,jdim-1,kin)
!    elseif(jdim.eq.1) then
!          mm(1,2) = 0.0
!          mm(2,2) = 0.0
!          mm(3,2) = 0.0
!    endif


    if(kdim.ge.5) then
      Do j=1,jdim
        Do i=1,idim
          mm(1,3,i,j,1) = -25./12.*x(i,j,1)+4.*x(i,j,2)-3.*x(i,j,3)+4./3.*x(i,j,4)-0.25*x(i,j,5)
          mm(2,3,i,j,1) = -25./12.*y(i,j,1)+4.*y(i,j,2)-3.*y(i,j,3)+4./3.*y(i,j,4)-0.25*y(i,j,5)
          mm(3,3,i,j,1) = -25./12.*z(i,j,1)+4.*z(i,j,2)-3.*z(i,j,3)+4./3.*z(i,j,4)-0.25*z(i,j,5)

          mm(1,3,i,j,2) = -0.25*x(i,j,1)-5./6.*x(i,j,2)+1.5*x(i,j,3)-0.5*x(i,j,4)+1./12.*x(i,j,5)
          mm(2,3,i,j,2) = -0.25*y(i,j,1)-5./6.*y(i,j,2)+1.5*y(i,j,3)-0.5*y(i,j,4)+1./12.*y(i,j,5)
          mm(3,3,i,j,2) = -0.25*z(i,j,1)-5./6.*z(i,j,2)+1.5*z(i,j,3)-0.5*z(i,j,4)+1./12.*z(i,j,5)

          mm(1,3,i,j,kdim-1) = -1./12.*x(i,j,kdim-4)+0.5*x(i,j,kdim-3)-1.5*x(i,j,kdim-2)&
                          + 5./6.*x(i,j,kdim-1)+0.25*x(i,j,kdim)
          mm(2,3,i,j,kdim-1) = -1./12.*y(i,j,kdim-4)+0.5*y(i,j,kdim-3)-1.5*y(i,j,kdim-2)&
                          + 5./6.*y(i,j,kdim-1)+0.25*y(i,j,kdim)
          mm(3,3,i,j,kdim-1) = -1./12.*z(i,j,kdim-4)+0.5*z(i,j,kdim-3)-1.5*z(i,j,kdim-2)&
                          + 5./6.*z(i,j,kdim-1)+0.25*z(i,j,kdim)

          mm(1,3,i,j,kdim) = 0.25*x(i,j,kdim-4)-4./3.*x(i,j,kdim-3)+3.*x(i,j,kdim-2)&
                        - 4.*x(i,j,kdim-1)+25./12.*x(i,j,kdim)
          mm(2,3,i,j,kdim) = 0.25*y(i,j,kdim-4)-4./3.*y(i,j,kdim-3)+3.*y(i,j,kdim-2)&
                        - 4.*y(i,j,kdim-1)+25./12.*y(i,j,kdim)
          mm(3,3,i,j,kdim) = 0.25*z(i,j,kdim-4)-4./3.*z(i,j,kdim-3)+3.*z(i,j,kdim-2)&
                        - 4.*z(i,j,kdim-1)+25./12.*z(i,j,kdim)

          Do k=3,kdim-2
          mm(1,3,i,j,k) = (x(i,j,k-2)-8.*(x(i,j,k-1)-x(i,j,k+1))-x(i,j,k+2))*c12i
          mm(2,3,i,j,k) = (y(i,j,k-2)-8.*(y(i,j,k-1)-y(i,j,k+1))-y(i,j,k+2))*c12i
          mm(3,3,i,j,k) = (z(i,j,k-2)-8.*(z(i,j,k-1)-z(i,j,k+1))-z(i,j,k+2))*c12i
          End Do
        End Do
       End Do
     End If



!    if(kdim.ge.5) then
!       if(kin.ge.3.and.kin.le.kdim-2) then
!          mm(1,3) = (x(iin,jin,kin-2)-8.*(x(iin,jin,kin-1)-x(iin,jin,kin+1))-x(iin,jin,kin+2))*c12i
!          mm(2,3) = (y(iin,jin,kin-2)-8.*(y(iin,jin,kin-1)-y(iin,jin,kin+1))-y(iin,jin,kin+2))*c12i
!          mm(3,3) = (z(iin,jin,kin-2)-8.*(z(iin,jin,kin-1)-z(iin,jin,kin+1))-z(iin,jin,kin+2))*c12i
!       elseif(kin.eq.1) then
!          mm(1,3) = -25./12.*x(iin,jin,1)+4.*x(iin,jin,2)-3.*x(iin,jin,3)+4./3.*x(iin,jin,4)-0.25*x(iin,jin,5)
!          mm(2,3) = -25./12.*y(iin,jin,1)+4.*y(iin,jin,2)-3.*y(iin,jin,3)+4./3.*y(iin,jin,4)-0.25*y(iin,jin,5)
!          mm(3,3) = -25./12.*z(iin,jin,1)+4.*z(iin,jin,2)-3.*z(iin,jin,3)+4./3.*z(iin,jin,4)-0.25*z(iin,jin,5)
!       elseif(kin.eq.2) then
!          mm(1,3) = -0.25*x(iin,jin,1)-5./6.*x(iin,jin,2)+1.5*x(iin,jin,3)-0.5*x(iin,jin,4)+1./12.*x(iin,jin,5)
!          mm(2,3) = -0.25*y(iin,jin,1)-5./6.*y(iin,jin,2)+1.5*y(iin,jin,3)-0.5*y(iin,jin,4)+1./12.*y(iin,jin,5)
!          mm(3,3) = -0.25*z(iin,jin,1)-5./6.*z(iin,jin,2)+1.5*z(iin,jin,3)-0.5*z(iin,jin,4)+1./12.*z(iin,jin,5)
!       elseif(kin.eq.kdim-1) then
!          mm(1,3) = -1./12.*x(iin,jin,kdim-4)+0.5*x(iin,jin,kdim-3)-1.5*x(iin,jin,kdim-2)&
!                          + 5./6.*x(iin,jin,kdim-1)+0.25*x(iin,jin,kdim)
!          mm(2,3) = -1./12.*y(iin,jin,kdim-4)+0.5*y(iin,jin,kdim-3)-1.5*y(iin,jin,kdim-2)&
!                          + 5./6.*y(iin,jin,kdim-1)+0.25*y(iin,jin,kdim)
!          mm(3,3) = -1./12.*z(iin,jin,kdim-4)+0.5*z(iin,jin,kdim-3)-1.5*z(iin,jin,kdim-2)&
!                          + 5./6.*z(iin,jin,kdim-1)+0.25*z(iin,jin,kdim)
!       elseif(kin.eq.kdim) then
!          mm(1,3) = 0.25*x(iin,jin,kdim-4)-4./3.*x(iin,jin,kdim-3)+3.*x(iin,jin,kdim-2)&
!                        - 4.*x(iin,jin,kdim-1)+25./12.*x(iin,jin,kdim)
!          mm(2,3) = 0.25*y(iin,jin,kdim-4)-4./3.*y(iin,jin,kdim-3)+3.*y(iin,jin,kdim-2)&
!                        - 4.*y(iin,jin,kdim-1)+25./12.*y(iin,jin,kdim)
!          mm(3,3) = 0.25*z(iin,jin,kdim-4)-4./3.*z(iin,jin,kdim-3)+3.*z(iin,jin,kdim-2)&
!                        - 4.*z(iin,jin,kdim-1)+25./12.*z(iin,jin,kdim)
!       endif
!    elseif(kdim.eq.3.or.kdim.eq.4) then
 !      if(kin.ge.2.and.kin.le.3) then
!          mm(1,3) = -0.5*x(iin,jin,kin-1) + 0.5*x(iin,jin,kin+1)
!          mm(2,3) = -0.5*y(iin,jin,kin-1) + 0.5*y(iin,jin,kin+1)
!          mm(3,3) = -0.5*z(iin,jin,kin-1) + 0.5*z(iin,jin,kin+1)
!       elseif(kin.eq.1) then
!          mm(1,3) = -1.5*x(iin,jin,1) + 2.0*x(iin,jin,2) - 0.5*x(iin,jin,3)
!          mm(2,3) = -1.5*y(iin,jin,1) + 2.0*y(iin,jin,2) - 0.5*y(iin,jin,3)
!          mm(3,3) = -1.5*z(iin,jin,1) + 2.0*z(iin,jin,2) - 0.5*z(iin,jin,3)
!       elseif(kin.eq.kdim) then
!          mm(1,3) = 0.5*x(iin,jin,kdim-2) - 2.0*x(iin,jin,kdim-1) + 1.5*x(iin,jin,kdim)
!          mm(2,3) = 0.5*y(iin,jin,kdim-2) - 2.0*y(iin,jin,kdim-1) + 1.5*y(iin,jin,kdim)
!          mm(3,3) = 0.5*z(iin,jin,kdim-2) - 2.0*z(iin,jin,kdim-1) + 1.5*z(iin,jin,kdim)
!       endif
!    elseif(kdim.eq.2) then
!          mm(1,3) = x(iin,jin,kdim) - x(iin,jin,kdim-1)
!          mm(2,3) = y(iin,jin,kdim) - y(iin,jin,kdim-1)
!          mm(3,3) = z(iin,jin,kdim) - z(iin,jin,kdim-1)
!    elseif(kdim.eq.1) then
!          mm(1,3) = 0.0
!          mm(2,3) = 0.0
!          mm(3,3) = 0.0
!    endif

    planeCalmm = mm
    end function planeCalmm

    function planeCalmet(idim,jdim,kdim,mm)
    ! Calculate mesh metrics met
    ! met(3,3)
    ! met(1,1) = didx, met(1,2) = didy, met(1,3) = didz
    ! met(2,1) = djdx, met(2,2) = djdy, met(2,3) = djdz
    ! met(3,1) = dkdx, met(3,2) = dkdy, met(3,3) = dkdz
    !
    ! given 
    !
    !  mm(3,3)
    !  mm(1,1) = dxdi, mm(1,2) = dxdj, mm(1,3) = dxdk
    !  mm(2,1) = dydi, mm(2,2) = dydj, mm(2,3) = dydk
    !  mm(3,1) = dzdi, mm(3,2) = dzdj, mm(3,3) = dzdk
    Integer, intent(in) :: idim,jdim,kdim
    real(8), intent(in) :: mm(3,3,idim,jdim,kdim)
    real(8) :: met(3,3,idim,jdim,kdim), planeCalmet(3,3,idim,jdim,kdim)
    real(8) :: vkapx, vkapy, vkapz, vol
    Integer :: i,j,k


     Do k=1,kdim
       Do j=1,jdim
         Do i=1,idim
     ! Calculate didx, didy, didz
           vkapx = mm(2,2,i,j,k)*mm(3,3,i,j,k) - mm(2,3,i,j,k)*mm(3,2,i,j,k)
           vkapy = mm(1,2,i,j,k)*mm(3,3,i,j,k) - mm(1,3,i,j,k)*mm(3,2,i,j,k)
           vkapz = mm(1,2,i,j,k)*mm(2,3,i,j,k) - mm(1,3,i,j,k)*mm(2,2,i,j,k)
           vol = + mm(1,1,i,j,k)*vkapx - mm(2,1,i,j,k)*vkapy + mm(3,1,i,j,k)*vkapz

           met(1,1,i,j,k) = +vkapx/vol
           met(1,2,i,j,k) = -vkapy/vol
           met(1,3,i,j,k) = +vkapz/vol
     ! Calculate djdx, djdy, djdz
           vkapx = mm(2,1,i,j,k)*mm(3,3,i,j,k) - mm(2,3,i,j,k)*mm(3,1,i,j,k)
           vkapy = mm(1,1,i,j,k)*mm(3,3,i,j,k) - mm(1,3,i,j,k)*mm(3,1,i,j,k)
           vkapz = mm(1,1,i,j,k)*mm(2,3,i,j,k) - mm(1,3,i,j,k)*mm(2,1,i,j,k)
           vol = - mm(1,2,i,j,k)*vkapx + mm(2,2,i,j,k)*vkapy - mm(3,2,i,j,k)*vkapz

           met(2,1,i,j,k) = -vkapx/vol
           met(2,2,i,j,k) = +vkapy/vol
           met(2,3,i,j,k) = -vkapz/vol
      ! Calculate dkdx, dkdy, dkdz
           vkapx = mm(2,1,i,j,k)*mm(3,2,i,j,k) - mm(2,2,i,j,k)*mm(3,1,i,j,k)
           vkapy = mm(1,1,i,j,k)*mm(3,2,i,j,k) - mm(1,2,i,j,k)*mm(3,1,i,j,k)
           vkapz = mm(1,1,i,j,k)*mm(2,2,i,j,k) - mm(1,2,i,j,k)*mm(2,1,i,j,k)
           vol = + mm(1,3,i,j,k)*vkapx - mm(2,3,i,j,k)*vkapy + mm(3,3,i,j,k)*vkapz

           met(3,1,i,j,k) = +vkapx/vol
           met(3,2,i,j,k) = -vkapy/vol
           met(3,3,i,j,k) = +vkapz/vol

         End Do
       End Do
     End Do

     planeCalmet = met

!    For debug purposes
!     print *, 'Identity mm*met =', matmul(mm,met)

     end function planeCalmet












































end module modMetrics
