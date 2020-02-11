!------How to compile?-------------------!!
!   1. mkdir build
!   2. cd build
!   3. cmake ..
!   4. make
!!---------------------------------------!!
! Definition of Cylindrical system:
!-----------------------------------------!
!             z ^     , r
!               |    /
!               |   /
!               |  /
!          theta|-/
!               |/
!               ---------------------> y
!            origin
!
!
!
!
!
!
!-----------------------------------------!
module coord_converter
implicit none

contains
!!-------------------------------------------------!!
subroutine coord_cart_to_cyl(y,z, theta, r, origin)
real(kind=8), dimension(:,:,:), intent(in) :: y, z
real(kind=8), dimension(:,:,:), intent(out) :: theta, r
real(kind=8), dimension(:), intent(in) :: origin

!theta = atan2(z-origin(2), y-origin(1))
theta = atan2(y-origin(1), z-origin(2))
r = sqrt((z-origin(2))**2+(y-origin(1))**2)

end subroutine
!!-------------------------------------------------!!
subroutine vel_cart_to_cyl(v,w, ut, ur, theta)
real(kind=8), dimension(:,:,:), intent(in) :: v, w, theta
real(kind=8), dimension(:,:,:), intent(out) :: ut, ur

!ut = -v*sin(theta) + w*cos(theta)
!ur =  v*cos(theta) + w*sin(theta)
ut = v*cos(theta) - w*sin(theta)
ur = v*sin(theta) + w*cos(theta)
end subroutine
!!-------------------------------------------------!!
subroutine coord_cyl_to_cart(theta, r, y, z, origin)
real(kind=8), dimension(:,:,:), intent(in) :: theta, r
real(kind=8), dimension(:,:,:), intent(out) :: y, z
real(kind=8), dimension(:), intent(in) :: origin

!z = origin(2) + r*sin(theta)
!y = origin(1) + r*cos(theta)
z = origin(2) + r*cos(theta)
y = origin(1) + r*sin(theta)


end subroutine
!!-------------------------------------------------!!
subroutine vel_cyl_to_cart(ut, ur, v, w, theta)
real(kind=8), dimension(:,:,:), intent(in) :: ut, ur, theta
real(kind=8), dimension(:,:,:), intent(out) :: v, w

!w = ur*sin(theta) + ut*cos(theta)
!v = ur*cos(theta) - ut*sin(theta)
w = ur*cos(theta) - ut*sin(theta)
v = ur*sin(theta) + ut*cos(theta)


end subroutine
!!-------------------------------------------------!!


end module
