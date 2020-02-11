module cov
implicit none

contains
!! --------------
subroutine cal_cov(f1, f2, cov_f, n, rank)
integer, intent(in) :: n, rank
real(kind=8), dimension(n, rank), intent(in) :: f1, f2
real(kind=8), intent(out) :: cov_f
real(kind=8), dimension(rank) :: mean_f1, mean_f2
real(kind=8) :: real_n
real(kind=8), dimension(n, rank) :: f1_work, f2_work
integer :: i

real_n = dble(n)

! average
mean_f1 = sum(f1, 1) / real_n
mean_f2 = sum(f2, 1) / real_n
do i = 1, rank
    f1_work(:,i) = f1(:,i) - mean_f1(i)
    f2_work(:,i) = f2(:,i) - mean_f2(i)
enddo


cov_f = sum(f1_work * f2_work) / real_n

end subroutine
!!-----------------

end module
