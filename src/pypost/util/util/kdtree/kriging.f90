module kriging
use corr
use grdd
use det
implicit none
contains
!! ---------------------------- !!
subroutine init_params_naive(theta, p)
real(kind=8), dimension(:,:,:), intent(out) :: theta, p

theta = 1.d0
p = 2.d0
print*, "Distance parameters initiated with NAIVE assumptions."

end subroutine
!! ----------------------------- !!
subroutine weights_ok(x, indexes_nearest, x0, w, theta, p)
real(kind=8), dimension(:,:), intent(in) :: x
integer, dimension(:,:), intent(in) :: indexes_nearest
real(kind=8), dimension(:,:), intent(in) :: x0
real(kind=8), dimension(:,:,:), intent(out) :: w
real(kind=8), dimension(:,:,:), intent(in) :: theta, p
real(kind=8), dimension(:,:), allocatable :: a_work, b_work
real(kind=8), dimension(:), allocatable :: la_work
integer, dimension(:), allocatable :: IPIV_work
integer :: i, j, rank, nx0, n_nearest, nvar

rank = size(x, 1)
n_nearest = size(indexes_nearest, 1)
nx0 = size(x0, 2)
nvar = size(theta, 2)

allocate(a_work(n_nearest+1,n_nearest+1), b_work(n_nearest+1,1), la_work(n_nearest+1), IPIV_work(n_nearest+1))

do j = 1, nx0
do i = 1, nvar
   call kernel_weights_ok(x(:,indexes_nearest(:,j)), x0(:,j), w(:,i,j), theta(:,i,j), p(:,i,j), a_work, b_work, la_work, IPIV_work)
enddo
enddo

deallocate(a_work, b_work, la_work, IPIV_work)

end subroutine
!! ---------------------------- !!
subroutine kernel_weights_ok(x, x0, w, theta, p, a_work, b_work, la_work, IPIV_work)
real(kind=8), dimension(:,:), intent(in) :: x
real(kind=8), dimension(:), intent(in) :: x0
real(kind=8), dimension(:), intent(out) :: w
real(kind=8), dimension(:), intent(in) :: theta
real(kind=8), dimension(:), intent(in) :: p
real(kind=8), intent(inout) :: a_work(:,:), b_work(:,:), la_work(:)
integer, dimension(:), intent(inout) :: IPIV_work
integer :: n_nearest, rank

integer :: i, info

n_nearest = size(x, 2)
rank = size(x, 1)
!do i = 1, n_nearest; print*, x(:,i); enddo
call construct_R(x, a_work(1:n_nearest,1:n_nearest), theta, p)
a_work(1:n_nearest,n_nearest+1) = 1.d0
a_work(n_nearest+1,1:n_nearest) = 1.d0
a_work(n_nearest+1,n_nearest+1) = 0.d0

do i = 1, n_nearest
    call corr_kriging(x(:,i), x0, b_work(i,1), theta, p)
enddo
b_work(n_nearest+1, 1) = 1.d0

!do i=1,n_nearest+1; print*, a_work(i,:),'--->   ',b_work(i,1);enddo
!print*,'\n'

call DSYSV("U", n_nearest+1, 1, a_work, n_nearest+1, IPIV_work, b_work, n_nearest+1, la_work, size(la_work,1), info)
if(info/=0) then
    print*, 'ERROR in SOLVING LINEAR SYSTEM!!!'
    stop
endif
!print*, b_work
w = b_work(1:n_nearest, 1)


end subroutine
!! ----------------------------- !!
subroutine interp_ok(f, indexes_nearest, w, f0)
real(kind=8), dimension(:,:), intent(in) :: f
integer, dimension(:,:), intent(in) :: indexes_nearest
real(kind=8), dimension(:,:,:), intent(in) :: w
real(kind=8), dimension(:,:), intent(out) :: f0
integer :: i, j, nx0, nvar

nvar = size(w, 2)
nx0 = size(w, 3)

do j = 1, nx0
do i = 1, nvar
   call kernel_interp_ok(f(i, indexes_nearest(:,j)), w(:,i,j), f0(i,j)) 
enddo
enddo

end subroutine
!! ------------------------ !!
subroutine kernel_interp_ok(f, w, f0)
real(kind=8), dimension(:), intent(in) :: f
real(kind=8), dimension(:), intent(in) :: w
real(kind=8), intent(out) :: f0

f0 = dot_product(f, w)

end subroutine
!! ------------------- !!
subroutine get_Rinv(R, Rinv, IPIV_work, work)
real(kind=8), dimension(:,:), intent(in) :: R
real(kind=8), dimension(:,:), intent(out) :: Rinv
integer, dimension(:) :: IPIV_work
real(kind=8), dimension(:) :: work

integer :: n, info, i, j


n = size(R, 1)
Rinv = R

call DSYTRF("U", n, Rinv, n, IPIV_work, work, size(work,1), info)
print*, Rinv
print*, IPIV_work
call DSYTRI("U", n, Rinv, n, IPIV_work, work, info)
do j = 1, n
do i = j+1, n
    Rinv(i,j) = Rinv(j,i)
enddo
enddo
end subroutine

!! ------------------- !!
subroutine get_Rinv_and_detR(R, Rinv, detR, IPIV_work, work)
real(kind=8), dimension(:,:), intent(in) :: R
real(kind=8), dimension(:,:), intent(out) :: Rinv
real(kind=8), intent(out) :: detR
integer, dimension(:) :: IPIV_work
real(kind=8), dimension(:) :: work

integer :: n, info, i, j


n = size(R, 1)
Rinv = R

call DSYTRF("U", n, Rinv, n, IPIV_work, work, size(work,1), info)
if(info/=0) then
    print*, 'ERROR in LU factorization!!! info:', info
    stop
endif

!print*, 'Here some results'
!print*, R
!print*, Rinv
!print*, IPIV_work
!! --determinant
call dsydet(Rinv, IPIV_work, detR) ! Note: Rinv here is only a intermediate result from last step

!! --invert
call DSYTRI("U", n, Rinv, n, IPIV_work, work, info)
if(info/=0) then
    print*, 'ERROR in matrix inversion!!!'
    stop
endif

do j = 1, n
do i = j+1, n
    Rinv(i,j) = Rinv(j,i)
enddo
enddo
end subroutine
!! ------------------------- !!
subroutine get_mu(Rinv, y, mu, work)
real(kind=8), dimension(:,:), intent(in) :: Rinv
real(kind=8), dimension(:), intent(in) :: y
real(kind=8), intent(out) :: mu
real(kind=8), dimension(:) :: work
real(kind=8) :: Rinv11

work = matmul(Rinv, y)
!print*, sum(Rinv), '/.', sum(work)
!print*, Rinv, sum(Rinv)
Rinv11 = sum(Rinv)
if(Rinv11 /= 0.d0) then
    mu = sum(work) / Rinv11
else
    !mu = sum(y) / dble(size(y,1))
    mu = sum(y) / dble(size(y,1))
endif

end subroutine
!! -------------------------- !!
subroutine get_sigma2(Rinv, y, mu, sigma2, work1, work2)
real(kind=8), dimension(:,:), intent(in) :: Rinv
real(kind=8), dimension(:), intent(in) :: y
real(kind=8), intent(in) :: mu
real(kind=8), intent(out) :: sigma2
real(kind=8), dimension(:) :: work1, work2
integer :: n 

n = size(Rinv, 1)
work1 = y - mu
work2 = matmul(Rinv, work1)
sigma2 = dot_product(work1, work2) / dble(n)

end subroutine
!! --------------------- !!
subroutine get_mle(x, y,theta,p, R, Rinv, ml, IPIV_work, work1, work2)
real(kind=8), dimension(:,:), intent(in) :: x
real(kind=8), dimension(:), intent(in) :: y, theta, p
real(kind=8), dimension(:,:), intent(out) :: R, Rinv
real(kind=8), intent(out) :: ml
integer, dimension(:) :: IPIV_work
real(kind=8), dimension(:) :: work1, work2
real(kind=8) :: mu, sigma2, detR
integer :: n

call construct_R(x, R, theta, p)
call get_Rinv_and_detR(R, Rinv, detR, IPIV_work, work1)
call get_mu(Rinv, y, mu, work1)
call get_sigma2(Rinv, y, mu, sigma2, work1, work2)
n = size(Rinv, 1)
ml = -dble(n)*.5d0* log(abs(sigma2)) - .5d0*log(abs(detR))
#ifdef VERBOSE
!print*, Rinv
print*, 'While Calculating MLE, ml = ', ml, 'mu = ', mu,'sigma2 = ', sigma2, 'detR = ', detR
#endif
end subroutine
!! ------------------- !!
subroutine get_gradml(x,y, ml, h, theta,p, gradml, &
                      R_work, Rinv_work, IPIV_work, work1, work2, theta_work, p_work)
real(kind=8), dimension(:,:), intent(in) :: x
real(kind=8), dimension(:), intent(in) :: y
real(kind=8), intent(in) :: ml, h
real(kind=8), dimension(:), intent(in) :: theta, p
real(kind=8), dimension(:), intent(out) :: gradml
integer, dimension(:) :: IPIV_work
real(kind=8), dimension(:,:) :: R_work, Rinv_work
real(kind=8), dimension(:) :: work1, work2
real(kind=8), dimension(:) :: theta_work, p_work
real(kind=8) :: ml_work
integer :: n, i

n = size(theta, 1)

do i = 1, n
    theta_work = theta
    theta_work(i) = theta_work(i) + h
    call get_mle(x,y, theta_work, p, R_work, Rinv_work, ml_work, IPIV_work, work1, work2)
    gradml(i) = (ml_work - ml) / h
enddo

end subroutine 
!! -----------------------------------------------!!
subroutine grdd_kriging(x, y, indexes_nearest, theta, p, h, tol, maxcount)
real(kind=8), dimension(:,:), intent(in) :: x
real(kind=8), dimension(:,:), intent(in) :: y
integer, dimension(:,:), intent(in) :: indexes_nearest
real(kind=8), dimension(:,:,:), intent(inout) :: theta, p
real(kind=8), intent(in) :: h, tol
integer, intent(in) :: maxcount
integer :: n_nearest, nx0, rank, nvar, i, j
real(kind=8), dimension(:,:), allocatable :: R_work, Rinv_work
real(kind=8), dimension(:), allocatable :: gradml, work1, work2, theta_work, p_work
integer, dimension(:), allocatable :: IPIV_work

n_nearest = size(indexes_nearest, 1)
nx0 = size(indexes_nearest, 2)
rank = size(x, 1)
nvar = size(y, 1)

allocate(R_work(n_nearest,n_nearest), Rinv_work(n_nearest,n_nearest), &
         work1(n_nearest), work2(n_nearest), gradml(rank), &
         theta_work(rank), p_work(rank))
allocate(IPIV_work(n_nearest))

do j = 1, nx0
do i = 1, nvar
    call kernel_grdd_kriging(x(:,indexes_nearest(:,j)), y(i,indexes_nearest(:,j)), &
                             theta(:,i,j), p(:,i,j),  h, tol, maxcount, &
                             R_work, Rinv_work, IPIV_work, gradml, &
                             work1, work2, theta_work, p_work)
enddo
enddo

end subroutine
!! ------------------------------ !!
subroutine kernel_grdd_kriging(x,y, theta, p,  h, tol, maxcount, R_work, Rinv_work, IPIV_work, gradml, work1, work2, theta_work, p_work)
real(kind=8), dimension(:,:), intent(in) :: x
real(kind=8), dimension(:), intent(in) :: y
real(kind=8), dimension(:), intent(inout) :: theta, p
real(kind=8), intent(in) :: h, tol
integer, intent(in) :: maxcount
real(kind=8), dimension(:,:) :: R_work, Rinv_work
integer, dimension(:) :: IPIV_work
real(kind=8), dimension(:) :: gradml, work1, work2, theta_work, p_work
real(kind=8) :: lambda, ml, ml_star, diff_ml
integer :: grdd_info, count_,  rank, n_nearest

rank = size(theta, 1)
n_nearest = size(x, 1)
lambda = tol*2048.d0

call get_mle(x, y, theta, p, R_work, Rinv_work, ml, IPIV_work, work1, work2)

do count_ = 1, maxcount
    lambda = lambda * 2.d0
    call get_gradml(x,y, ml, h, theta,p, gradml, R_work, Rinv_work, IPIV_work, work1, work2, theta_work, p_work)

    grdd_info = 1
    do while(grdd_info /= 0)
        call grdd_update(theta, gradml, lambda, theta_work)
        call get_mle(x,y, theta_work, p, R_work, Rinv_work, ml_star, IPIV_work, work1, work2)
#ifdef VERBOSE
        print*, 'Got likelihood = ', ml_star, 'with theta = ', theta_work, 'and lambda = ', lambda
#endif
        if(ml_star < ml) then
#ifdef VERBOSE
            print*, 'Likelihood decreased with lambda =', lambda
#endif
            lambda = lambda * .5d0
            grdd_info = -1
        elseif(any(theta_work<0.d0) .eqv. .true.) then
#ifdef VERBOSE
            print*, 'Found theta below 0. likelihood = ', ml_star, 'theta = ', theta_work
#endif               
            call grdd_lower_bound(theta, gradml, lambda, 0.d0, theta_work)
            grdd_info = -2
!                where(theta_work<0.d0)
!                    theta = 0.d0
!                endwhere
!                theta_work = theta
!                call get_mle(x,y, theta_work, p, R_work, Rinv_work, ml_star, IPIV_work, work1, work2)
        else
            diff_ml = abs(ml - ml_star)
            ml = ml_star
            theta = theta_work
            grdd_info = 0
        endif
    enddo
#ifdef VERBOSE
    print*, 'Count = ', count_, ': likelihood = ', ml, 'theta = ', theta, 'lambda = ', lambda
    print*, '==================================================================='
#endif
    if( diff_ml < tol )  then
#ifdef VERBOSE
        print*, 'Convergence reached! theta = ', theta
#endif
        exit
    endif
enddo
end subroutine
!! -------------------- !!

end module
