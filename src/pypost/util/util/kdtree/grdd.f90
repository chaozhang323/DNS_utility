module grdd
implicit none

contains
!! ----------------- !!
subroutine grdd_update(x_in, dfdx, lambda, x_out, minimize)
real(kind=8), dimension(:), intent(in) :: x_in, dfdx
real(kind=8), intent(in) :: lambda
real(kind=8), dimension(:), intent(out) :: x_out 
integer, intent(in), optional :: minimize
real(kind=8) :: norm_dfdx

norm_dfdx =  sqrt(sum(dfdx*dfdx))
if(norm_dfdx/=0.d0) then
    x_out = dfdx / norm_dfdx
else
    x_out = 0.d0
endif

if(present(minimize) .and. minimize==1) then
    x_out = x_in - lambda*x_out
else
    x_out = x_in + lambda*x_out
endif

end subroutine
!! ----------------- !!
subroutine grdd_lower_bound(x_in, dfdx, lambda, lb, x_out)
real(kind=8), dimension(:), intent(in) :: x_in, dfdx
real(kind=8), dimension(:), intent(out) :: x_out
real(kind=8), intent(in) :: lb
real(kind=8), intent(inout):: lambda

x_out = lb - x_in
where(dfdx/=0.d0)
    x_out = x_out / dfdx
elsewhere
    x_out = 0.d0
endwhere

lambda = minval(x_out, mask=(x_out>0.d0))

end subroutine

end module

