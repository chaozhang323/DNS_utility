module mod_diff
implicit none
real(kind=8), parameter, private :: dp_1i12 = 1.d0 / 12.d0
real(kind=8), parameter, private :: dp_25i12 = 25.d0 / 12.d0
real(kind=8), parameter, private :: dp_4i3 = 4.d0 / 3.d0
real(kind=8), parameter, private :: dp_5i6 = 5.d0 / 6.d0
!public :: diff_m2p2, diff_m0p4, diff_m1p3, diff_m4p0, diff_m3p1
contains

subroutine diff_m2p2(f, dfdx)
real(kind=8), dimension(-2:2), intent(in) :: f
real(kind=8), intent(out) :: dfdx

dfdx = ( f(-2) - 8.d0*(f(-1)-f(1)) - f(2) ) * dp_1i12
return
end subroutine

subroutine diff_m0p4(f, dfdx)
real(kind=8), dimension(0:4), intent(in) :: f
real(kind=8), intent(out) :: dfdx

dfdx = -dp_25i12*f(0) + 4.d0*f(1) - 3.d0*f(2) + dp_4i3*f(3) - 0.25d0*f(4)
return
end subroutine

subroutine diff_m1p3(f, dfdx)
real(kind=8), dimension(-1:3), intent(in) :: f
real(kind=8), intent(out) :: dfdx

dfdx = -0.25d0*f(-1) - dp_5i6*f(0) + 1.5d0*f(1) - 0.5d0*f(2) + dp_1i12*f(3)
return
end subroutine

subroutine diff_m4p0(f, dfdx)
real(kind=8), dimension(-4:0), intent(in) :: f
real(kind=8), intent(out) :: dfdx

dfdx = dp_25i12*f(0) - 4.d0*f(-1) + 3.d0*f(-2) - dp_4i3*f(-3) + 0.25d0*f(-4)
return
end subroutine

subroutine diff_m3p1(f, dfdx)
real(kind=8), dimension(-3:1), intent(in) :: f
real(kind=8), intent(out) :: dfdx

dfdx = 0.25d0*f(1) + dp_5i6*f(0) - 1.5d0*f(-1) + 0.5d0*f(-2) - dp_1i12*f(-3)
return
end subroutine
!-------------------------------!
subroutine easydiff(f, df, nx)
integer, intent(in) :: nx
real(kind=8), dimension(nx), intent(in) :: f
real(kind=8), dimension(nx), intent(out) :: df
integer :: i

df = 0.d0
if(nx==1) then
    return
elseif(nx<5) then
    print*, 'In easydiff, nx must be greater than 5 or equal to 1.'
    stop
endif

call diff_m0p4(f(1:5), df(1))
call diff_m1p3(f(1:5), df(2))
do i = 3, nx-2
    call diff_m2p2(f(i-2:i+2), df(i))
enddo
call diff_m3p1(f(nx-4:nx), df(nx-1))
call diff_m4p0(f(nx-4:nx), df(nx))

return
end subroutine


end module
