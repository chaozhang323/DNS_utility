module MFLOW
! contains subroutines for calculating molecular viscosity using either Sutherland law or Keyes law(both air and nitrogen).
implicit none
real(8) :: a0vis,a1vis,omegavis,a2vis

contains
 subroutine InitViscosity(ivis)
   integer,intent(in) :: ivis
   select case(ivis)
    case(1) ! sutherland(air)
     a0vis=1.458d-6; a1vis=110.4d0; omegavis=0.5d0; a2vis=0.d0
    case(2) ! keyes(air)
     a0vis=1.488d-6; a1vis=122.1d0; omegavis=0.5d0; a2vis=5.d0
    case(3) ! keyes(nitrogen)
     a0vis=1.418d-6; a1vis=116.4d0; omegavis=0.5d0; a2vis=5.d0
    case default
     print*,'Unknown iviscosity type.',ivis,' STOP!!!'
     stop
   end select
 end subroutine InitViscosity

 elemental real(8) function CalMu(tin)
   real(8),intent(in) :: tin
   real(8) :: tmp
   tmp=10.d0**(-a2vis/tin)
   calMu=a0vis*tin**(1.d0+omegavis)/(tin+a1vis*tmp)
 end function calMu

! ******************* calculate thermal conductivity assuming const Prandtl number **********************
  elemental real(8) function CalKappa(muin,cp,Prandtl)
      real(8), intent(in) :: muin,cp,Prandtl
      !if (iusePrandtl) then
        CalKappa = muin*cp/Prandtl
  end function CalKappa
!==============================================================================================
end module MFLOW
