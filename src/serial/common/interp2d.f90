module interp2d
use omp_lib
implicit none
contains
!------------------------------------------!
!   ge_solver: Gaussian_Elimination solver for linear system a*x = b
!   In: a_in: coefficient matrix (n by n)
!       b_in: vector on right hand side (n)
!       n:  dimension of the system
!   Out:    x
!-------------------------------------------!
subroutine ge_solver(a_in, b_in, x, n)
integer, intent(in) :: n
real(8), intent(in) :: a_in(n,n), b_in(n)
real(8) :: a(n,n), b(n)
real(8), intent(out) :: x(n)
integer :: k, i
real(8) :: w

a = a_in
b = b_in
do k = 2, n
    do i = k, n
        w = a(i,k-1) / a(k-1,k-1)
        a(i,k:n) = a(i,k:n) - w*a(k-1,k:n)
        b(i) = b(i) - w*b(k-1)
    enddo
enddo

!print*, a

x(n) = b(n) / a(n,n)
do k = n-1, 1, -1
    x(k) = (b(k) - sum( a(k,k+1:n)*x(k+1:n) )) / a(k,k)
enddo

end subroutine
!------------------------------------------!
!   get_intersect:  judge if two straight lines on a plane have a intersection
!   In: p00,p01:    two end points of line 1
!       p10,p11:    two end points of line 2
!   Out:    inter:  0 if there's no intersection, 1 if there is
!------------------------------------------!
subroutine get_intersect(p00, p01, p10, p11, inter)
real(8), dimension(2), intent(in) :: p00, p01, p10, p11
integer, intent(out) :: inter
real(8), dimension(2) :: vec0, vec1
real(8) :: one_over_det, l(2), v,t, det

vec0 = p01 - p00
vec1 = p11 - p10
det = -vec0(1)*vec1(2) + vec0(2)*vec1(1)
if( abs(det)<1.d-6 ) then
    inter = 0
    return
endif

one_over_det = 1.d0 / det!(-vec0(1)*vec1(2) + vec0(2)*vec1(1))
l = p10 - p00
t = (-vec1(2)*l(1) + vec1(1)*l(2) )*one_over_det
v = (-vec0(2)*l(1) + vec0(1)*l(2) )*one_over_det
!print*, t,v
if( t>=0.d0 .and. t<=1.d0 .and. v>=0.d0 .and. v<=1.d0) then
!   print*,'tv',t,v
!   print*, p00
!   print*, p01
!   print*, p10
!   print*, p11
    inter = 1
    return
else
    inter = 0
    return
endif

end subroutine
!--------------------------------------!
!   get_in_cell:    judge if a point is in current cell
!   In: x,y:    coordinate of current cell votices
!       p:  the point
!       p_far:  a point outside the current cell
!   Out:    nintersect: number of intersection the line joined by p & p_far make with the cell faces
!--------------------------------------!
subroutine get_in_cell(x, y, p, p_far, nintersect)
real(8), dimension(2,2), intent(in) :: x, y
real(8), dimension(2), intent(in) :: p, p_far
integer, intent(out) :: nintersect
integer :: inter

nintersect = 0
call get_intersect(p, p_far, (/x(1,1),y(1,1)/), (/x(2,1),y(2,1)/), inter)
nintersect = nintersect + inter
call get_intersect(p, p_far, (/x(2,1),y(2,1)/), (/x(2,2),y(2,2)/), inter)
nintersect = nintersect + inter
call get_intersect(p, p_far, (/x(2,2),y(2,2)/), (/x(1,2),y(1,2)/), inter)
nintersect = nintersect + inter
call get_intersect(p, p_far, (/x(1,2),y(1,2)/), (/x(1,1),y(1,1)/), inter)
nintersect = nintersect + inter
!print*, nintersect
return

end subroutine
!------------------------------------------!
!   get_current_cell:   find the current cell in which the point is
!   In: grid_x, grid_y: the structured grid
!       nx,ny:  the dimension of the grid
!       p: the point
!       p_far: the point outside of the whole set of grid
!       iguess_in,jguess_in:    the presumed index to start searching from
!   Out:    idx_i,idx_j:    the current cell index
!           ifnd:   1 if the cell if found, 0 if not
!------------------------------------------!
subroutine get_current_cell(grid_x, grid_y, nx,ny, p, p_far, iguess_in, jguess_in, idx_i, idx_j, ifnd)
integer, intent(in) :: nx, ny, iguess_in, jguess_in
real(8), dimension(nx,ny), intent(in) :: grid_x, grid_y
real(8), dimension(2), intent(in) :: p, p_far
integer, dimension(2), intent(out) :: idx_i, idx_j
integer, intent(out) :: ifnd
integer :: i, j, nintersect, iguess, jguess


ifnd = 0
!---guess---!
iguess = iguess_in
jguess = jguess_in
call get_in_cell(grid_x(iguess:iguess+1,jguess:jguess+1), &
                grid_y(iguess:iguess+1,jguess:jguess+1), p, p_far, nintersect)
if( mod(nintersect,2)==1 ) then
    idx_i = (/iguess, iguess+1/)
    idx_j = (/jguess, jguess+1/)
    ifnd = 1
    return
endif

iguess = iguess + 1
if( iguess>nx-1 .and. jguess<ny-1) then
    iguess = 1
    jguess = jguess + 1
endif
call get_in_cell(grid_x(iguess:iguess+1,jguess:jguess+1), &
                grid_y(iguess:iguess+1,jguess:jguess+1), p, p_far, nintersect)
if( mod(nintersect,2)==1 ) then
    idx_i = (/iguess, iguess+1/)
    idx_j = (/jguess, jguess+1/)
    ifnd = 1
    return
endif


do j = 1, ny-1
do i = 1, nx-1
    call get_in_cell(grid_x(i:i+1,j:j+1), grid_y(i:i+1,j:j+1), p, p_far, nintersect)
    !print*, i,j, nintersect
    if( mod(nintersect,2)==1 ) then
        idx_i = (/i, i+1/)
        idx_j = (/j, j+1/)
        ifnd = 1
        return
    endif
enddo
enddo

end subroutine
!_-----------------------------------------!
!   engine_linear:  linear interpolation engine
!   In: p1,p2,p3,p4:    4 points the define the cell
!       f1,f2,f3,f4:    4 values on the 4 points
!       p:  the point to be interpolated
!   out:    f:  the value on point p
!------------------------------------------!
subroutine engine_linear(p1,p2,p3,p4, f1,f2,f3,f4, p,f)
real(8), dimension(2), intent(in) :: p1,p2,p3,p4, p
real(8), intent(in) :: f1,f2,f3,f4
real(8), intent(out) :: f
real(8), dimension(4) :: a, poly
real(8), dimension(4,4) :: mat_lhs
real(8), dimension(4) :: vec_rhs
integer i
real(8) :: scalex, scaley

mat_lhs(:,1) = 1.
mat_lhs(:,2) = (/p1(1), p2(1), p3(1), p4(1)/) - p1(1)
scalex = maxval(mat_lhs(:,2))
mat_lhs(:,2) = mat_lhs(:,2) / scalex
mat_lhs(:,3) = (/p1(2), p2(2), p3(2), p4(2)/) - p1(2)
scaley = maxval(mat_lhs(:,3))
mat_lhs(:,3) = mat_lhs(:,3) / scaley
mat_lhs(:,4) = mat_lhs(:,2)*mat_lhs(:,3)
vec_rhs(1) = f1
vec_rhs(2) = f2
vec_rhs(3) = f3
vec_rhs(4) = f4

!print*, mat_lhs(1,:)
!print*, mat_lhs(2,:)
!print*, mat_lhs(3,:)
!print*, mat_lhs(4,:)
!print*, vec_rhs
call ge_solver(mat_lhs, vec_rhs, a, 4)
!print*, a

poly(1) = 1.d0
poly(2) = (p(1) - p1(1)) / scalex
poly(3) = (p(2) - p1(2)) / scaley
poly(4) = poly(2)*poly(3)
f = sum( a*poly )

!print*, poly
!print*, p
!print*,p1
!print*,p2
!print*,p3
!print*,p4


end subroutine
!------------------------------------------!
! Initialize interp2d_linear
!------------------------------------------!
subroutine init_interp2d_linear(grid_old_x, grid_old_y, grid_new_x, grid_new_y, nx_old,ny_old, nx_new,ny_new, p_far, indx)
integer, intent(in) :: nx_old, ny_old, nx_new, ny_new
real(8), dimension(nx_old, ny_old), intent(in) :: grid_old_x, grid_old_y
real(8), dimension(nx_new, ny_new), intent(in) :: grid_new_x, grid_new_y
integer, dimension(2, nx_new, ny_new), intent(out) :: indx
real(8), dimension(2), intent(in) :: p_far
integer :: i, j, idx_i(2), idx_j(2), ifnd
real(8), dimension(2) :: p
integer :: cr_t,cm_t,c1_t,c2_t
real(8) :: rate
print*,'initialize interp2d_linear ...'

idx_i = 1
idx_j = 1

! initialize the system clock
call system_clock(count_rate=cr_t)
call system_clock(count_max=cm_t)
rate = dble(cr_t)
call system_clock(c1_t)
!!$OMP PARALLEL &
!!$OMP private(i,j,p,ifnd)
!!$OMP DO
do j = 1, ny_new
do i = 1, nx_new
    p = (/ grid_new_x(i,j), grid_new_y(i,j) /)
    call get_current_cell(grid_old_x, grid_old_y, nx_old, ny_old, &
                        p, p_far, idx_i(1),idx_j(1), idx_i, idx_j, ifnd)
    if( ifnd==1 ) then
        !print*, idx_i, idx_j

        indx(1,i,j) = idx_i(1)
        indx(2,i,j) = idx_j(1)
    else
        indx(:,i,j) = -1
    endif
enddo
enddo
!!$OMP END DO
!!$OMP END PARALLEL
print*,' '
call system_clock(c2_t)
print*,'init_interp2d_linear time [sec]: ',dble(c2_t-c1_t)/rate
print*,' '
end subroutine
!------------------------------------------!
!   interp2d_linear:    linear interpolation
!   In:
!   Out:    data_new:new data interpolated on new grid
!------------------------------------------!
subroutine do_interp2d_linear(grid_old_x, grid_old_y, grid_new_x, grid_new_y, &
                              indx, data_old, data_new, nx_old,ny_old, nx_new, ny_new)
integer, intent(in) :: nx_old, ny_old, nx_new, ny_new
integer, dimension(2, nx_new, ny_new), intent(in) :: indx
real(8), dimension(nx_old, ny_old), intent(in) :: grid_old_x, grid_old_y, data_old
real(8), dimension(nx_new, ny_new), intent(in) :: grid_new_x, grid_new_y
real(8), dimension(nx_new, ny_new), intent(out) :: data_new
integer :: i, j, idx_i(2), idx_j(2), ifnd, inan, inan_tmp
integer, dimension(nx_new,ny_new) :: nan_map
real(8), dimension(2) :: p1,p2,p3,p4, p
integer :: cr_t,cm_t,c1_t,c2_t
real(8) :: rate

inan = 0
nan_map = 0
data_new = 0.d0

! initialize the system clock
call system_clock(count_rate=cr_t)
call system_clock(count_max=cm_t)
rate = dble(cr_t)
call system_clock(c1_t)
!!$OMP PARALLEL &
!!$OMP private(i,j,p,ifnd,idx_i,idx_j,p1,p2,p3,p4)
!!$OMP DO
do j = 1, ny_new
do i = 1, nx_new
    p = (/ grid_new_x(i,j), grid_new_y(i,j) /)
    if( any(indx(:,i,j).ne.-1) ) then
        !print*, idx_i, idx_j
        idx_i(1) = indx(1,i,j); idx_i(2) = indx(1,i,j)+1
        idx_j(1) = indx(2,i,j); idx_j(2) = indx(2,i,j)+1
        p1 = (/ grid_old_x(idx_i(1),idx_j(1)), grid_old_y(idx_i(1),idx_j(1)) /)
        p2 = (/ grid_old_x(idx_i(2),idx_j(1)), grid_old_y(idx_i(2),idx_j(1)) /)
        p3 = (/ grid_old_x(idx_i(2),idx_j(2)), grid_old_y(idx_i(2),idx_j(2)) /)
        p4 = (/ grid_old_x(idx_i(1),idx_j(2)), grid_old_y(idx_i(1),idx_j(2)) /)
        call engine_linear(p1,p2,p3,p4, data_old(idx_i(1),idx_j(1)), data_old(idx_i(2),idx_j(1)), &
                        data_old(idx_i(2),idx_j(2)), data_old(idx_i(1),idx_j(2)), p, data_new(i,j))
    else
        inan = inan + 1
        nan_map(i,j) = 1
    endif
enddo
enddo
!!$OMP END DO
!!$OMP END PARALLEL
print*,' '
call system_clock(c2_t)
print*,'do_interp2d_linear time [sec]: ',dble(c2_t-c1_t)/rate
print*,' '

do while( inan/=0 )
    inan = 0
    do j = 1, ny_new
    do i = 1, nx_new
        if( nan_map(i,j)==1 ) then
            call copy_nearest(i,j, nan_map, data_new, nx_new, ny_new, inan_tmp)
            inan = inan + inan_tmp
        endif
    enddo
    enddo
enddo

end subroutine

!------------------------------------------!
!   interp2d_linear:    linear interpolation
!   In: grid_old_x, grid_old_y: old grid coordinate
!       grid_new_x, grid_new_y: new grid
!       data_old:   old data on the old grid
!       nx_old,ny_old:  dimension of old grid
!       nx_new,ny_new:  dimension of new grid
!       p_far:  the point that is outside whole set of grid
!   Out:    data_new:new data interpolated on new grid
!------------------------------------------!
subroutine interp2d_linear(grid_old_x, grid_old_y, grid_new_x, grid_new_y, &
                            data_old, data_new, nx_old,ny_old, nx_new, ny_new, p_far)
integer, intent(in) :: nx_old, ny_old, nx_new, ny_new
real(8), dimension(nx_old, ny_old), intent(in) :: grid_old_x, grid_old_y, data_old
real(8), dimension(nx_new, ny_new), intent(in) :: grid_new_x, grid_new_y
real(8), dimension(nx_new, ny_new), intent(out) :: data_new
real(8), dimension(2), intent(in) :: p_far
integer :: i, j, idx_i(2), idx_j(2), ifnd, inan, inan_tmp
integer, dimension(nx_new,ny_new) :: nan_map
real(8), dimension(2) :: p1,p2,p3,p4, p
real(8) :: time_begin, time_end

idx_i = 1
idx_j = 1
inan = 0
nan_map = 0
data_new = 0.d0

time_begin = omp_get_wtime()
!!$OMP PARALLEL &
!!$OMP private(i,j,p,ifnd,p1,p2,p3,p4)
!!$OMP DO
do j = 1, ny_new
do i = 1, nx_new
    p = (/ grid_new_x(i,j), grid_new_y(i,j) /)
    call get_current_cell(grid_old_x, grid_old_y, nx_old, ny_old, &
                        p, p_far, idx_i(1),idx_j(1), idx_i, idx_j, ifnd)
    if( ifnd==1 ) then
        !print*, idx_i, idx_j

        p1 = (/ grid_old_x(idx_i(1),idx_j(1)), grid_old_y(idx_i(1),idx_j(1)) /)
        p2 = (/ grid_old_x(idx_i(2),idx_j(1)), grid_old_y(idx_i(2),idx_j(1)) /)
        p3 = (/ grid_old_x(idx_i(2),idx_j(2)), grid_old_y(idx_i(2),idx_j(2)) /)
        p4 = (/ grid_old_x(idx_i(1),idx_j(2)), grid_old_y(idx_i(1),idx_j(2)) /)
        call engine_linear(p1,p2,p3,p4, data_old(idx_i(1),idx_j(1)), data_old(idx_i(2),idx_j(1)), &
                        data_old(idx_i(2),idx_j(2)), data_old(idx_i(1),idx_j(2)), p, data_new(i,j))
    else
        inan = inan + 1
        nan_map(i,j) = 1
    endif
enddo
enddo
!!$OMP END DO
!!$OMP END PARALLEL
print*,' '
print*,'interp2d_linear OMP time [sec]: ',omp_get_wtime()-time_begin
print*,' '

do while( inan/=0 )
    inan = 0
    do j = 1, ny_new
    do i = 1, nx_new
        if( nan_map(i,j)==1 ) then
            call copy_nearest(i,j, nan_map, data_new, nx_new, ny_new, inan_tmp)
            inan = inan + inan_tmp
        endif

    enddo
    enddo
enddo


end subroutine
!-----------------------------------------------------------!
!   copy_nearest:   copy nearst data from the new dataset, if not interpolated
!   In: i,j:    index of point that interpolation failed
!   Inout:  nan_map:    map that marks failure of interpolation
!           dat:    dataset
!   Out: inan:  0 if no copy happens, otherwise positive
!
!-----------------------------------------------------------!
subroutine copy_nearest(i,j, nan_map, dat, nx, ny, inan)
integer, intent(in) :: nx,ny,i,j
integer, dimension(nx,ny), intent(inout) :: nan_map
real(8), dimension(nx,ny), intent(inout) :: dat
integer, intent(out) :: inan
integer i1,j1

inan = 0

j1 = j + 1
if(j1<=ny) then
    if( nan_map(i,j1)==0 ) then
        dat(i,j) = dat(i,j1)
        nan_map(i,j) = 0
        inan = inan + 1
        return
    endif
endif
j1 = j - 1
if(j1>=1) then
    if( nan_map(i,j1)==0 ) then
        dat(i,j) = dat(i,j1)
        nan_map(i,j) = 0
        inan = inan + 1
        return
    endif
endif
i1 = i + 1
if(i1<=nx) then
    if( nan_map(i1,j)==0 ) then
        dat(i,j) = dat(i1,j)
        nan_map(i,j) = 0
        inan = inan + 1
        return
    endif
endif
i1 = i - 1
if(i1>=1) then
    if( nan_map(i1,j)==0 ) then
        dat(i,j) = dat(i1,j)
        nan_map(i,j) = 0
        inan = inan + 1
        return
    endif
endif

end subroutine


end module
