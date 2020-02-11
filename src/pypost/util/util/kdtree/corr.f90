module corr
implicit none
contains
!! ------------------!
subroutine construct_R(x, R, theta, p)
real(kind=8), dimension(:,:), intent(in) :: x
real(kind=8), dimension(:,:), intent(inout) :: R
real(kind=8), dimension(:), intent(in) :: theta
real(kind=8), dimension(:), intent(in) :: p
integer :: i, j, nx

nx = size(x, 2)

R = 1.d0
do j = 1, nx
do i = j+1, nx
!    f_work = 0.d0
!    do n = 1, rank
!        f_work = f_work - theta(n)*(abs(x(i,n)-x(j,n))**p(n))
!    enddo
!    R(i,j) = exp(f_work)
    call corr_kriging(x(:,i), x(:,j), R(i,j), theta, p)
    R(j,i) = R(i,j)
enddo
enddo
end subroutine
!! ---------------------- 
subroutine corr_kriging(x, y, corr, theta, p)
real(kind=8), dimension(:), intent(in) :: x, y, theta, p
real(kind=8), intent(inout) :: corr
integer :: rank
integer :: n

rank = size(x, 1)

corr = 0.d0
do n = 1, rank
    corr = corr - theta(n)*(abs(x(n)-y(n))**p(n))
enddo
corr = exp(corr)
end subroutine


end module
