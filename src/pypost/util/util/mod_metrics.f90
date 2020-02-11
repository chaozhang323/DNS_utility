module mod_metrics
use mod_diff
implicit none
real(8), parameter, private :: c12i = 1.d0/12.d0
!public :: extend_grid, compute_grid_derivative, compute_mesh_metrics
contains
!------------------------------------------------!
subroutine extend_1d(x, x_out, nx, stencilsize, bc)
integer, intent(in) :: nx, stencilsize, bc(2)
real(kind=8), intent(in) :: x(nx)
real(kind=8), intent(out) :: x_out(1-stencilsize:nx+stencilsize)
real(kind=8) :: lx
integer :: n

!! init
x_out = 0.d0
x_out(1:nx) = x
if(nx==1) then
    do n = 1-stencilsize, nx+stencilsize
        x_out(n) = x(1)
    enddo
    return
endif

!! extend
selectcase(bc(1))
case(0)! extrapolate
    lx = x_out(2) - x_out(1)
    do n = 1, stencilsize
        x_out(1-n) = x_out(2-n) - lx
    enddo
case(1)! mirror centered at n=1
    do n = 1, stencilsize
        x_out(1-n) = x_out(1+n)
    enddo
case(-1)! reflect centered at n=1
    if(x_out(1)/=0.d0) then
        print*, 'Error: reflect not at 0 axis!'
        stop
    endif
 do n = 1, stencilsize
        x_out(1-n) = -x_out(1+n)
    enddo
case(2)! mirror centered at n=0
    do n = 1, stencilsize
        x_out(1-n) = x_out(n)
    enddo
case(-2)! reflect centered at n=0
    if(x_out(1)==0.d0) then
        print*, 'Error: reflect at 0 axis!'
        stop
    endif
    do n = 1, stencilsize
        x_out(1-n) = -x_out(n)
    enddo
case default
    print*, 'Not implemented!'
    stop
endselect

selectcase(bc(2))
case(0)! extrapolate
    lx = x_out(nx) - x_out(nx-1)
    do n = 1, stencilsize
        x_out(nx+n) = x_out(nx-1+n) + lx
    enddo
case(1)! mirror centered at n=nx
    do n = 1, stencilsize
        x_out(nx+n) = x_out(nx-n)
    enddo
case(-1)! reflect centered at n=nx
    do n = 1, stencilsize
        x_out(nx+n) = -x_out(nx-n)
    enddo
case(2)! mirror centered at n=nx+1
    do n = 1, stencilsize
        x_out(nx+n) = x_out(nx-n+1)
    enddo
case(-2)! reflect centered at n=nx+1
    do n = 1, stencilsize
        x_out(nx+n) = -x_out(nx-n+1)
    enddo
case default
    print*, 'Not implemented!'
    stop
endselect

end subroutine


!------------------------------------------------!
subroutine extend_x(x, x_out, nx,ny,nz, stencilsize)
integer, intent(in) :: nx,ny,nz, stencilsize
real(kind=8), dimension(nx,ny,nz), intent(in) :: x
real(kind=8), dimension(1-stencilsize:nx+stencilsize, &
                        1-stencilsize:ny+stencilsize, &
                        1-stencilsize:nz+stencilsize), intent(out) :: x_out
real(kind=8), allocatable :: lx(:,:,:)
integer :: n

x_out = 0.d0
x_out(1:nx,1:ny,1:nz) = x

allocate(lx(1:ny,1:nz,2))
lx(:,:,1) = x_out(2, 1:ny,1:nz) - x_out(1,   1:ny,1:nz)
lx(:,:,2) = x_out(nx,1:ny,1:nz) - x_out(nx-1,1:ny,1:nz)
do n = 1, stencilsize
    x_out(1-n, 1:ny,1:nz) = x_out(2-n,   1:ny,1:nz) - lx(:,:,1)
    x_out(nx+n,1:ny,1:nz) = x_out(nx-1+n,1:ny,1:nz) + lx(:,:,2)
enddo
deallocate(lx)
allocate(lx(1-stencilsize:nx+stencilsize,1:nz,2))
lx(:,:,1) = x_out(1-stencilsize:nx+stencilsize,2, 1:nz) - x_out(1-stencilsize:nx+stencilsize,1,   1:nz) 
lx(:,:,2) = x_out(1-stencilsize:nx+stencilsize,ny,1:nz) - x_out(1-stencilsize:nx+stencilsize,ny-1,1:nz) 
do n = 1, stencilsize
    x_out(1-stencilsize:nx+stencilsize,1-n, 1:nz) = x_out(1-stencilsize:nx+stencilsize,2-n,   1:nz) - lx(:,:,1)
    x_out(1-stencilsize:nx+stencilsize,ny+n,1:nz) = x_out(1-stencilsize:nx+stencilsize,ny-1+n,1:nz) + lx(:,:,2)
enddo
deallocate(lx)
allocate(lx(1-stencilsize:nx+stencilsize,1-stencilsize:ny+stencilsize,2))
lx(:,:,1) = x_out(1-stencilsize:nx+stencilsize,1-stencilsize:ny+stencilsize, 2) &
          - x_out(1-stencilsize:nx+stencilsize,1-stencilsize:ny+stencilsize, 1) 
lx(:,:,2) = x_out(1-stencilsize:nx+stencilsize,1-stencilsize:ny+stencilsize,nz) &
          - x_out(1-stencilsize:nx+stencilsize,1-stencilsize:ny+stencilsize,nz-1) 
do n = 1, stencilsize
    x_out(1-stencilsize:nx+stencilsize,1-stencilsize:ny+stencilsize,   1-n) = &
    x_out(1-stencilsize:nx+stencilsize,1-stencilsize:ny+stencilsize,   2-n) - lx(:,:,1)
    x_out(1-stencilsize:nx+stencilsize,1-stencilsize:ny+stencilsize,  nz+n) = &
    x_out(1-stencilsize:nx+stencilsize,1-stencilsize:ny+stencilsize,nz-1+n) + lx(:,:,2)
enddo
deallocate(lx)

return
end subroutine

!------------------------------------------------!
subroutine extend_grid(x,y,z, x_out,y_out,z_out, nx,ny,nz,stencilsize, bc)
integer, intent(in) :: nx,ny,nz, stencilsize, bc(6)
real(kind=8), dimension(nx,ny,nz), intent(in) :: x,y,z
real(kind=8), dimension(1-stencilsize:nx+stencilsize, &
                        1-stencilsize:ny+stencilsize, &
                        1-stencilsize:nz+stencilsize), intent(out) :: x_out,y_out,z_out
real(kind=8) :: lx(ny,nz,3,2), ly(nx,nz,3,2), lz(nx,ny,3,2)
integer :: i, j, k

do k = 1, nz
do j = 1, ny
    call extend_1d(x(:,j,k), x_out(:,j,k), nx, stencilsize, -bc(1:2))
    call extend_1d(y(:,j,k), y_out(:,j,k), nx, stencilsize, bc(1:2))
    call extend_1d(z(:,j,k), z_out(:,j,k), nx, stencilsize, bc(1:2))
enddo
enddo


do k = 1, nz
do i = 1-stencilsize, nx+stencilsize
    call extend_1d(x_out(i,1:ny,k), x_out(i,:,k), ny, stencilsize, bc(3:4))
    call extend_1d(y_out(i,1:ny,k), y_out(i,:,k), ny, stencilsize, -bc(3:4))
    call extend_1d(z_out(i,1:ny,k), z_out(i,:,k), ny, stencilsize, bc(3:4))
enddo
enddo


do j = 1-stencilsize, ny+stencilsize
do i = 1-stencilsize, nx+stencilsize
    call extend_1d(x_out(i,j,1:nz), x_out(i,j,:), nz, stencilsize, bc(5:6))
    call extend_1d(y_out(i,j,1:nz), y_out(i,j,:), nz, stencilsize, bc(5:6))
    call extend_1d(z_out(i,j,1:nz), z_out(i,j,:), nz, stencilsize, -bc(5:6))
enddo
enddo
!call extend_x(x, x_out, nx,ny,nz,stencilsize)
!call extend_x(y, y_out, nx,ny,nz,stencilsize)
!call extend_x(z, z_out, nx,ny,nz,stencilsize)

return
end subroutine


!-------------------------------------------------------------!
subroutine compute_grid_derivative(x, y, z, grid_derivative, nx,ny,nz, stencilsize, bc)
integer, intent(in) :: nx,ny,nz, stencilsize, bc(6)
real(8), dimension(1-stencilsize:nx+stencilsize, &
                   1-stencilsize:ny+stencilsize, &
                   1-stencilsize:nz+stencilsize), intent(in) :: x,y,z
real(8), dimension(3,3, 1-stencilsize:nx+stencilsize, &
                        1-stencilsize:ny+stencilsize, &
                        1-stencilsize:nz+stencilsize), intent(out) :: grid_derivative
integer :: i, j, k, n
integer :: ist, ien, jst, jen, kst, ken

!! init 
grid_derivative = 0.d0

!! define interior points
ist = 1-stencilsize+2
ien = nx+stencilsize-2
jst = 1-stencilsize+2
jen = ny+stencilsize-2
kst = 1-stencilsize+2
ken = nz+stencilsize-2

!! interior points ( along with periodic boundary points)
do k = kst, ken
do j = jst, jen
    do i = ist, ien
        call diff_m2p2(x(i-2:i+2,j,k), grid_derivative(1,1,i,j,k)) 
        call diff_m2p2(y(i-2:i+2,j,k), grid_derivative(2,1,i,j,k))
        call diff_m2p2(z(i-2:i+2,j,k), grid_derivative(3,1,i,j,k))
    enddo
enddo
enddo
do k = kst, ken
do i = ist, ien
    do j = jst, jen
        call diff_m2p2(x(i,j-2:j+2,k), grid_derivative(1,2,i,j,k)) 
        call diff_m2p2(y(i,j-2:j+2,k), grid_derivative(2,2,i,j,k))
        call diff_m2p2(z(i,j-2:j+2,k), grid_derivative(3,2,i,j,k))
    enddo
enddo
enddo
do j = jst, jen
do i = ist, ien
    do k = kst, ken
        call diff_m2p2(x(i,j,k-2:k+2), grid_derivative(1,3,i,j,k)) 
        call diff_m2p2(y(i,j,k-2:k+2), grid_derivative(2,3,i,j,k))
        call diff_m2p2(z(i,j,k-2:k+2), grid_derivative(3,3,i,j,k))
    enddo
enddo
enddo

!! boundary points with one-sided finite difference
if(bc(1)/=0) then
    do k = kst, ken
    do j = jst, jen
        call diff_m0p4(x(1:5,j,k), grid_derivative(1,1,1,j,k))
        call diff_m0p4(y(1:5,j,k), grid_derivative(2,1,1,j,k))
        call diff_m0p4(z(1:5,j,k), grid_derivative(3,1,1,j,k))
        call diff_m1p3(x(1:5,j,k), grid_derivative(1,1,2,j,k))
        call diff_m1p3(y(1:5,j,k), grid_derivative(2,1,2,j,k))
        call diff_m1p3(z(1:5,j,k), grid_derivative(3,1,2,j,k))
    enddo
    enddo
endif
if(bc(2)/=0) then
    do k = kst, ken
    do j = jst, jen
        call diff_m3p1(x(nx-4:nx,j,k), grid_derivative(1,1,nx-1,j,k))
        call diff_m3p1(y(nx-4:nx,j,k), grid_derivative(2,1,nx-1,j,k))
        call diff_m3p1(z(nx-4:nx,j,k), grid_derivative(3,1,nx-1,j,k))
        call diff_m4p0(x(nx-4:nx,j,k), grid_derivative(1,1,nx,j,k))
        call diff_m4p0(y(nx-4:nx,j,k), grid_derivative(2,1,nx,j,k))
        call diff_m4p0(z(nx-4:nx,j,k), grid_derivative(3,1,nx,j,k))
    enddo
    enddo
endif
if(bc(3)/=0) then
    do k = kst, ken
    do i = ist, ien
        call diff_m0p4(x(i,1:5,k), grid_derivative(1,2,i,1,k))
        call diff_m0p4(y(i,1:5,k), grid_derivative(2,2,i,1,k))
        call diff_m0p4(z(i,1:5,k), grid_derivative(3,2,i,1,k))
        call diff_m1p3(x(i,1:5,k), grid_derivative(1,2,i,2,k))
        call diff_m1p3(y(i,1:5,k), grid_derivative(2,2,i,2,k))
        call diff_m1p3(z(i,1:5,k), grid_derivative(3,2,i,2,k))
    enddo
    enddo
endif
if(bc(4)/=0) then
    do k = kst, ken
    do i = ist, ien
        call diff_m3p1(x(i,ny-4:ny,k), grid_derivative(1,2,i,ny-1,k))
        call diff_m3p1(y(i,ny-4:ny,k), grid_derivative(2,2,i,ny-1,k))
        call diff_m3p1(z(i,ny-4:ny,k), grid_derivative(3,2,i,ny-1,k))
        call diff_m4p0(x(i,ny-4:ny,k), grid_derivative(1,2,i,ny,k))
        call diff_m4p0(y(i,ny-4:ny,k), grid_derivative(2,2,i,ny,k))
        call diff_m4p0(z(i,ny-4:ny,k), grid_derivative(3,2,i,ny,k))
    enddo
    enddo
endif
if(bc(5)/=0) then
    do j = jst, jen
    do i = ist, ien
        call diff_m0p4(x(i,j,1:5), grid_derivative(1,3,i,j,1))
        call diff_m0p4(y(i,j,1:5), grid_derivative(2,3,i,j,1))
        call diff_m0p4(z(i,j,1:5), grid_derivative(3,3,i,j,1))
        call diff_m1p3(x(i,j,1:5), grid_derivative(1,3,i,j,2))
        call diff_m1p3(y(i,j,1:5), grid_derivative(2,3,i,j,2))
        call diff_m1p3(z(i,j,1:5), grid_derivative(3,3,i,j,2))
    enddo
    enddo
endif
if(bc(6)/=0) then
    do j = jst, jen
    do i = ist, ien
        call diff_m3p1(x(i,j,nz-4:nz), grid_derivative(1,3,i,j,nz-1))
        call diff_m3p1(y(i,j,nz-4:nz), grid_derivative(2,3,i,j,nz-1))
        call diff_m3p1(z(i,j,nz-4:nz), grid_derivative(3,3,i,j,nz-1))
        call diff_m4p0(x(i,j,nz-4:nz), grid_derivative(1,3,i,j,nz))
        call diff_m4p0(y(i,j,nz-4:nz), grid_derivative(2,3,i,j,nz))
        call diff_m4p0(z(i,j,nz-4:nz), grid_derivative(3,3,i,j,nz))
    enddo
    enddo
endif


end subroutine
!! --------------------------------------------------------------- !!
subroutine diff_metrics(f, df, nx, stencilsize, bc)
integer, intent(in) :: nx, stencilsize, bc(2)
real(kind=8), intent(in) :: f(1-stencilsize:nx+stencilsize)
real(kind=8), intent(out) :: df(1-stencilsize:nx+stencilsize)
integer i, ist, ien

!! init
df = 0.d0

!! nx==1
if(nx==1) then
    df = 1.d0
    return
endif

!! define interior points
ist = 1-stencilsize+2
ien = nx+stencilsize-2

!! diff interior
do i = ist, ien
   call diff_m2p2(f(i-2:i+2), df(i))
enddo

!! diff boundary
if(bc(1)/=0) then
    call diff_m0p4(f(1:5), df(1))
    call diff_m1p3(f(1:5), df(2))
endif
if(bc(2)/=0) then
    call diff_m3p1(f(nx-4:nx), df(nx-1))
    call diff_m4p0(f(nx-4:nx), df(nx))
endif

end subroutine
!! --------------------------------------------------------------- !!
subroutine diff_metrics_3d(f, df, nx, ny, nz, stencilsize, bc)
integer, intent(in) :: nx,ny,nz, stencilsize, bc(6)
real(kind=8), intent(in) :: f(1-stencilsize:nx+stencilsize, 1-stencilsize:ny+stencilsize, 1-stencilsize:nz+stencilsize)
real(kind=8), intent(out) :: df(3, 1-stencilsize:nx+stencilsize, 1-stencilsize:ny+stencilsize, 1-stencilsize:nz+stencilsize)
integer :: i, j, k

do k = 1, nz
do j = 1, ny
    call diff_metrics(f(:,j,k), df(1,:,j,k), nx, stencilsize, bc(1:2))
enddo
enddo
do k = 1, nz
do i = 1, nx
    call diff_metrics(f(i,:,k), df(2,i,:,k), ny, stencilsize, bc(3:4))
enddo
enddo
do j = 1, ny
do i = 1, nx
    call diff_metrics(f(i,j,:), df(3,i,j,:), nz, stencilsize, bc(5:6))
enddo
enddo


end subroutine

!---------------------------------------------------------------!
subroutine compute_metric_identities(mm, vol, identity, nx,ny,nz, stencilsize, bc)
integer, intent(in) :: nx,ny,nz, stencilsize, bc(6)
real(kind=8), intent(in) :: mm(3,4, 1-stencilsize:nx+stencilsize, &
                                    1-stencilsize:ny+stencilsize, &
                                    1-stencilsize:nz+stencilsize)
real(kind=8), intent(in) :: vol(3, 1-stencilsize:nx+stencilsize, &
                                   1-stencilsize:ny+stencilsize, &
                                   1-stencilsize:nz+stencilsize)
real(kind=8), intent(out) :: identity(3, 1-stencilsize:nx+stencilsize, &
                                   1-stencilsize:ny+stencilsize, &
                                   1-stencilsize:nz+stencilsize)
real(kind=8) :: mm_work(3,3, 1-stencilsize:nx+stencilsize, &
                                         1-stencilsize:ny+stencilsize, &
                                         1-stencilsize:nz+stencilsize)
real(kind=8) :: a_work(3, 1-stencilsize:nx+stencilsize, &
                          1-stencilsize:ny+stencilsize, &
                          1-stencilsize:nz+stencilsize)
integer :: i,j,k
integer :: ist, ien, jst, jen, kst, ken

!! init 
identity = 0.d0
do j = 1, 3
do i = 1, 3
    mm_work(i,j,:,:,:) = mm(i,j,:,:,:)*vol(1,:,:,:)
enddo
enddo

!! define interior points
ist = 1-stencilsize+2
ien = nx+stencilsize-2
jst = 1-stencilsize+2
jen = ny+stencilsize-2
kst = 1-stencilsize+2
ken = nz+stencilsize-2

!! xi
do k = kst, ken
do j = jst, jen
    call diff_metrics(mm_work(1,1,:,j,k), a_work(1,:,j,k), nx, stencilsize, bc(1:2))
    call diff_metrics(mm_work(1,2,:,j,k), a_work(2,:,j,k), nx, stencilsize, bc(1:2))
    call diff_metrics(mm_work(1,3,:,j,k), a_work(3,:,j,k), nx, stencilsize, bc(1:2))
enddo
enddo
identity = identity + a_work

!! eta
do k = kst, ken
do i = ist, ien
    call diff_metrics(mm_work(2,1,i,:,k), a_work(1,i,:,k), ny, stencilsize, bc(3:4))
    call diff_metrics(mm_work(2,2,i,:,k), a_work(2,i,:,k), ny, stencilsize, bc(3:4))
    call diff_metrics(mm_work(2,3,i,:,k), a_work(3,i,:,k), ny, stencilsize, bc(3:4))
enddo
enddo
identity = identity + a_work

!! zeta
do j = jst, jen
do i = ist, ien
    call diff_metrics(mm_work(3,1,i,j,:), a_work(1,i,j,:), nz, stencilsize, bc(5:6))
    call diff_metrics(mm_work(3,2,i,j,:), a_work(2,i,j,:), nz, stencilsize, bc(5:6))
    call diff_metrics(mm_work(3,3,i,j,:), a_work(3,i,j,:), nz, stencilsize, bc(5:6))
enddo
enddo
identity = identity + a_work

identity(:,1-stencilsize:0,:,:) = 0.d0
identity(:,:,1-stencilsize:0,:) = 0.d0
identity(:,:,:,1-stencilsize:0) = 0.d0
identity(:,nx+1:nx+stencilsize,:,:) = 0.d0
identity(:,:,ny+1:ny+stencilsize,:) = 0.d0
identity(:,:,:,nz+1:nz+stencilsize) = 0.d0
end subroutine
!---------------------------------------------------------------!
subroutine compute_mmx(grid_derivative, mm, vol, nx,ny,nz)
integer, intent(in) :: nx,ny,nz
real(kind=8), dimension(3,3,nx,ny,nz), intent(in) :: grid_derivative
real(kind=8), dimension(4,nx,ny,nz), intent(out) :: mm
real(kind=8), dimension(nx,ny,nz), intent(out) :: vol
real(kind=8), dimension(nx,ny,nz) :: vkapx,vkapy,vkapz

mm = 0.d0

vkapx = grid_derivative(2,2,:,:,:)*grid_derivative(3,3,:,:,:) - grid_derivative(2,3,:,:,:)*grid_derivative(3,2,:,:,:)
vkapy = grid_derivative(1,2,:,:,:)*grid_derivative(3,3,:,:,:) - grid_derivative(1,3,:,:,:)*grid_derivative(3,2,:,:,:)
vkapz = grid_derivative(1,2,:,:,:)*grid_derivative(2,3,:,:,:) - grid_derivative(1,3,:,:,:)*grid_derivative(2,2,:,:,:)
vol = grid_derivative(1,1,:,:,:)*vkapx - grid_derivative(2,1,:,:,:)*vkapy + grid_derivative(3,1,:,:,:)*vkapz
where(vol(:,:,:)/=0.d0)
    mm(1,:,:,:) = vkapx/vol
    mm(2,:,:,:) = -vkapy/vol
    mm(3,:,:,:) = vkapz/vol
endwhere
mm(4,:,:,:) = sqrt(sum(mm(1:3,:,:,:)**2,1))

return
end subroutine
!----------------------------------------------------------------!
subroutine compute_mmy(grid_derivative, mm, vol, nx,ny,nz)
integer, intent(in) :: nx,ny,nz
real(kind=8), dimension(3,3,nx,ny,nz), intent(in) :: grid_derivative
real(kind=8), dimension(4,nx,ny,nz), intent(out) :: mm
real(kind=8), dimension(nx,ny,nz), intent(out) :: vol
real(kind=8), dimension(nx,ny,nz) :: vkapx,vkapy,vkapz

mm = 0.d0

vkapx = grid_derivative(2,1,:,:,:)*grid_derivative(3,3,:,:,:) - grid_derivative(2,3,:,:,:)*grid_derivative(3,1,:,:,:)
vkapy = grid_derivative(1,1,:,:,:)*grid_derivative(3,3,:,:,:) - grid_derivative(1,3,:,:,:)*grid_derivative(3,1,:,:,:)
vkapz = grid_derivative(1,1,:,:,:)*grid_derivative(2,3,:,:,:) - grid_derivative(1,3,:,:,:)*grid_derivative(2,1,:,:,:)
vol(:,:,:) = - grid_derivative(1,2,:,:,:)*vkapx + grid_derivative(2,2,:,:,:)*vkapy - grid_derivative(3,2,:,:,:)*vkapz
where(vol(:,:,:)/=0.d0)
    mm(1,:,:,:) = -vkapx/vol(:,:,:)
    mm(2,:,:,:) = vkapy/vol(:,:,:)
    mm(3,:,:,:) = -vkapz/vol(:,:,:)
endwhere
mm(4,:,:,:) = sqrt(sum(mm(1:3,:,:,:)**2,1))

return
end subroutine
!----------------------------------------------------------------!
subroutine compute_mmz(grid_derivative, mm, vol, nx,ny,nz)
integer, intent(in) :: nx,ny,nz
real(kind=8), dimension(3,3,nx,ny,nz), intent(in) :: grid_derivative
real(kind=8), dimension(4,nx,ny,nz), intent(out) :: mm
real(kind=8), dimension(nx,ny,nz), intent(out) :: vol
real(kind=8), dimension(nx,ny,nz) :: vkapx,vkapy,vkapz

mm = 0.d0

vkapx = grid_derivative(2,1,:,:,:)*grid_derivative(3,2,:,:,:) - grid_derivative(2,2,:,:,:)*grid_derivative(3,1,:,:,:)
vkapy = grid_derivative(1,1,:,:,:)*grid_derivative(3,2,:,:,:) - grid_derivative(1,2,:,:,:)*grid_derivative(3,1,:,:,:)
vkapz = grid_derivative(1,1,:,:,:)*grid_derivative(2,2,:,:,:) - grid_derivative(1,2,:,:,:)*grid_derivative(2,1,:,:,:)
vol(:,:,:) = grid_derivative(1,3,:,:,:)*vkapx - grid_derivative(2,3,:,:,:)*vkapy + grid_derivative(3,3,:,:,:)*vkapz
where(vol(:,:,:)/=0.d0)
    mm(1,:,:,:) = vkapx/vol(:,:,:)
    mm(2,:,:,:) = -vkapy/vol(:,:,:)
    mm(3,:,:,:) = vkapz/vol(:,:,:)
endwhere
mm(4,:,:,:) = sqrt(sum(mm(1:3,:,:,:)**2,1))

return
end subroutine
!----------------------------------------------------------------!
subroutine compute_mesh_metrics(grid_derivative, mm,vol, nx,ny,nz)
integer, intent(in) :: nx,ny,nz
real(kind=8), dimension(3,3,nx,ny,nz), intent(in) :: grid_derivative
real(kind=8), dimension(3,4,nx,ny,nz), intent(out) :: mm
real(kind=8), dimension(3,nx,ny,nz), intent(out) :: vol

call compute_mmx(grid_derivative, mm(1,:,:,:,:), vol(1,:,:,:), nx,ny,nz)
call compute_mmy(grid_derivative, mm(2,:,:,:,:), vol(2,:,:,:), nx,ny,nz)
call compute_mmz(grid_derivative, mm(3,:,:,:,:), vol(3,:,:,:), nx,ny,nz)

return
end subroutine
!---------------------------------------------------------------!
!----------------------------------------------------------------!
subroutine compute_mesh_metrics_GCL(x,y,z, grid_derivative, mm,vol, nx,ny,nz, stencilsize, bc)
integer, intent(in) :: nx,ny,nz, stencilsize, bc(6)
real(kind=8), dimension(1-stencilsize:nx+stencilsize, &
                        1-stencilsize:ny+stencilsize, &
                        1-stencilsize:nz+stencilsize), intent(in) :: x,y,z
real(kind=8), dimension(3,3,1-stencilsize:nx+stencilsize, &
                        1-stencilsize:ny+stencilsize, &
                        1-stencilsize:nz+stencilsize), intent(in) :: grid_derivative
real(kind=8), dimension(3,4, 1-stencilsize:nx+stencilsize, &
                        1-stencilsize:ny+stencilsize, &
                        1-stencilsize:nz+stencilsize), intent(out) :: mm
real(kind=8), dimension(3, 1-stencilsize:nx+stencilsize, &
                        1-stencilsize:ny+stencilsize, &
                        1-stencilsize:nz+stencilsize), intent(out) :: vol
integer :: ist, ien, jst, jen, kst, ken, i, j, k
real(kind=8), dimension(1-stencilsize:nx+stencilsize, &
                        1-stencilsize:ny+stencilsize, &
                        1-stencilsize:nz+stencilsize) :: a_work, b_work

call compute_mesh_metrics(grid_derivative, mm,vol, nx+2*stencilsize,ny+2*stencilsize,nz+2*stencilsize)

!! define interior points
ist = 1-stencilsize+2
ien = nx+stencilsize-2
jst = 1-stencilsize+2
jen = ny+stencilsize-2
kst = 1-stencilsize+2
ken = nz+stencilsize-2

do j = jst, jen
do i = ist, ien
    call diff_metrics(grid_derivative(2,2,i,j,:)*z(i,j,:), a_work(i,j,:), nz, stencilsize, bc(5:6))
enddo
enddo
do k = kst, ken
do i = ist, ien
    call diff_metrics(grid_derivative(2,3,i,:,k)*z(i,:,k), b_work(i,:,k), ny, stencilsize, bc(3:4))
enddo
enddo
mm(1,1,:,:,:) = a_work - b_work

do j = jst, jen
do i = ist, ien
    call diff_metrics(grid_derivative(3,2,i,j,:)*x(i,j,:), a_work(i,j,:), nz, stencilsize, bc(5:6))
enddo
enddo
do k = kst, ken
do i = ist, ien
    call diff_metrics(grid_derivative(3,3,i,:,k)*x(i,:,k), b_work(i,:,k), ny, stencilsize, bc(3:4))
enddo
enddo
mm(1,2,:,:,:) = a_work - b_work

do j = jst, jen
do i = ist, ien
    call diff_metrics(grid_derivative(1,2,i,j,:)*y(i,j,:), a_work(i,j,:), nz, stencilsize, bc(5:6))
enddo
enddo
do k = kst, ken
do i = ist, ien
    call diff_metrics(grid_derivative(1,3,i,:,k)*y(i,:,k), b_work(i,:,k), ny, stencilsize, bc(3:4))
enddo
enddo
mm(1,3,:,:,:) = a_work - b_work

do k = kst, ken
do j = jst, jen
    call diff_metrics(grid_derivative(2,3,:,j,k)*z(:,j,k), a_work(:,j,k), nx, stencilsize, bc(1:2))
enddo
enddo
do j = jst, jen
do i = ist, ien
    call diff_metrics(grid_derivative(2,1,i,j,:)*z(i,j,:), b_work(i,j,:), nz, stencilsize, bc(5:6))
enddo
enddo
mm(2,1,:,:,:) = a_work - b_work

do k = kst, ken
do j = jst, jen
    call diff_metrics(grid_derivative(3,3,:,j,k)*x(:,j,k), a_work(:,j,k), nx, stencilsize, bc(1:2))
enddo
enddo
do j = jst, jen
do i = ist, ien
    call diff_metrics(grid_derivative(3,1,i,j,:)*x(i,j,:), b_work(i,j,:), nz, stencilsize, bc(5:6))
enddo
enddo
mm(2,2,:,:,:) = a_work - b_work

do k = kst, ken
do j = jst, jen
    call diff_metrics(grid_derivative(1,3,:,j,k)*y(:,j,k), a_work(:,j,k), nx, stencilsize, bc(1:2))
enddo
enddo
do j = jst, jen
do i = ist, ien
    call diff_metrics(grid_derivative(1,1,i,j,:)*y(i,j,:), b_work(i,j,:), nz, stencilsize, bc(5:6))
enddo
enddo
mm(2,3,:,:,:) = a_work - b_work

do k = kst, ken
do i = ist, ien
    call diff_metrics(grid_derivative(2,1,i,:,k)*z(i,:,k), a_work(i,:,k), ny, stencilsize, bc(3:4))
enddo
enddo
do k = kst, ken
do j = jst, jen
    call diff_metrics(grid_derivative(2,2,:,j,k)*z(:,j,k), b_work(:,j,k), nx, stencilsize, bc(1:2))
enddo
enddo
mm(3,1,:,:,:) = a_work - b_work

do k = kst, ken
do i = ist, ien
    call diff_metrics(grid_derivative(3,1,i,:,k)*x(i,:,k), a_work(i,:,k), ny, stencilsize, bc(3:4))
enddo
enddo
do k = kst, ken
do j = jst, jen
    call diff_metrics(grid_derivative(3,2,:,j,k)*x(:,j,k), b_work(:,j,k), nx, stencilsize, bc(1:2))
enddo
enddo
mm(3,2,:,:,:) = a_work - b_work

do k = kst, ken
do i = ist, ien
    call diff_metrics(grid_derivative(1,1,i,:,k)*y(i,:,k), a_work(i,:,k), ny, stencilsize, bc(3:4))
enddo
enddo
do k = kst, ken
do j = jst, jen
    call diff_metrics(grid_derivative(1,2,:,j,k)*y(:,j,k), b_work(:,j,k), nx, stencilsize, bc(1:2))
enddo
enddo
mm(3,3,:,:,:) = a_work - b_work

do j = 1, 3
do i = 1, 3
    where(vol(1,:,:,:)/=0.d0)
        mm(i,j,:,:,:) = mm(i,j,:,:,:) / vol(1,:,:,:)
    endwhere
enddo
enddo




return
end subroutine

!! ----------------------------------------------------- !!
subroutine cylindrical_mesh_metrics(mm, z, mm_cyl, nx,ny,nz)
integer, intent(in) :: nx,ny,nz
real(kind=8), dimension(3,4,nx,ny,nz), intent(in) :: mm
real(kind=8), dimension(nx,ny,nz), intent(in) :: z
real(kind=8), dimension(3,4,nx,ny,nz), intent(out) :: mm_cyl
integer i

mm_cyl = mm

do i = 1, 3
where(z/=0.d0)
    mm_cyl(i,2,:,:,:) = mm(i,2,:,:,:) / z
endwhere
enddo

end subroutine

!--------------------------------------------------------!
subroutine compute_mesh_metrics_halfpoints(grid_derivative, mmx,mmy,mmz, volx,voly,volz, nx,ny,nz)
integer, intent(in) :: nx,ny,nz
real(kind=8), dimension(3,3,nx,ny,nz), intent(in) :: grid_derivative
real(kind=8), dimension(4,nx-1,ny,nz), intent(out) :: mmx
real(kind=8), dimension(4,nx,ny-1,nz), intent(out) :: mmy
real(kind=8), dimension(4,nx,ny,nz-1), intent(out) :: mmz
real(kind=8), dimension(nx-1,ny,nz), intent(out) :: volx
real(kind=8), dimension(nx,ny-1,nz), intent(out) :: voly
real(kind=8), dimension(nx,ny,nz-1), intent(out) :: volz
real(kind=8), dimension(3,3,nx-1,ny,nz) :: gdx
real(kind=8), dimension(3,3,nx,ny-1,nz) :: gdy
real(kind=8), dimension(3,3,nx,ny,nz-1) :: gdz

gdx = 5.d-1 * (grid_derivative(:,:,1:nx-1,:,:)+grid_derivative(:,:,2:nx,:,:))
gdy = 5.d-1 * (grid_derivative(:,:,:,1:ny-1,:)+grid_derivative(:,:,:,2:ny,:))
gdz = 5.d-1 * (grid_derivative(:,:,:,:,1:nz-1)+grid_derivative(:,:,:,:,2:nz))

call compute_mmx(gdx, mmx, volx, nx-1,ny,nz)
call compute_mmy(gdy, mmy, voly, nx,ny-1,nz)
call compute_mmz(gdz, mmz, volz, nx,ny,nz-1)

return
end subroutine
!------------------------------------------------!
subroutine mm_calculator_2d(x, y, mm, nx, ny)
real(8), dimension(nx,ny), intent(in) :: x, y
integer nx, ny
real(8), dimension(2,2,nx,ny), intent(out) :: mm
integer i, j

do j = 1, ny
do i = 1, nx
    mm(:,:,i,j) = calmm_2d(i,j,x,y)
enddo
enddo

end subroutine
!----------------------------------------------!





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












































end module 
