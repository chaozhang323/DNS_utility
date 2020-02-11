module flow
implicit none

contains
subroutine rortex_2d(dudx, dudy, dvdx, dvdy, rortex, nx,ny)
integer, intent(in) :: nx,ny
real(8), intent(in), dimension(nx,ny) :: dudx,dudy,dvdx,dvdy
real(8), intent(out), dimension(nx,ny) :: rortex
real(8), dimension(nx,ny) :: alpha, beta, alpha_sq, beta_sq

alpha_sq = 0.25d0 * ((dvdy-dudx)**2+(dvdx+dudy)**2)
beta = 0.5d0 * (dvdx-dudy)
beta_sq = beta**2

where(beta_sq>alpha_sq)
    alpha = sqrt(alpha_sq)
    where(beta>0.d0)
        rortex = 2.d0*(beta-alpha)
    elsewhere
        rortex = -2.d0*(beta+alpha)
    endwhere
elsewhere
    rortex = 0.d0
endwhere
end subroutine
!-----------------------------------------------------------!
subroutine dudx_3d(mm, u, dudx, nx,ny,nz)
use mod_diff, only: easydiff
integer, intent(in) :: nx, ny, nz
real(kind=8), dimension(3,nx,ny,nz), intent(in) :: mm
real(kind=8), dimension(nx,ny,nz), intent(in) :: u
real(kind=8), dimension(nx,ny,nz), intent(out) ::dudx 
real(kind=8), dimension(nx,ny,nz) :: a_work, b_work, c_work
integer :: i, j, k

!! dudx
do k = 1, nz; do j = 1, ny
    call easydiff(u(:,j,k), a_work(:,j,k), nx)
enddo; enddo
do k = 1, nz; do i = 1, nx
    call easydiff(u(i,:,k), b_work(i,:,k), ny)
enddo; enddo
do j = 1, ny; do i = 1, nx
    call easydiff(u(i,j,:), c_work(i,j,:), nz)
enddo; enddo
dudx = mm(1,:,:,:)*a_work + mm(2,:,:,:)*b_work + mm(3,:,:,:)*c_work
end subroutine

!-----------------------------------------------------------!
subroutine div_3d(mm, u, v, w, div, nx, ny, nz)
use mod_diff, only: easydiff
integer, intent(in) :: nx, ny, nz
real(kind=8), dimension(3,4,nx,ny,nz), intent(in) :: mm
real(kind=8), dimension(nx,ny,nz), intent(in) :: u, v, w
real(kind=8), dimension(nx,ny,nz), intent(out) :: div
real(kind=8), dimension(nx,ny,nz) :: dudx, dvdy, dwdz

!!
call dudx_3d(mm(:,1,:,:,:), u, dudx, nx,ny,nz)
call dudx_3d(mm(:,2,:,:,:), v, dvdy, nx,ny,nz)
call dudx_3d(mm(:,3,:,:,:), w, dwdz, nx,ny,nz)
div = dudx+dvdy+dwdz

end subroutine

!-----------------------------------------------------------!
subroutine curl_3d(mm, u, v, w, curl, nx, ny, nz)
integer, intent(in) :: nx, ny, nz
real(kind=8), dimension(3,4,nx,ny,nz), intent(in) :: mm
real(kind=8), dimension(nx,ny,nz), intent(in) :: u, v, w
real(kind=8), dimension(3, nx,ny,nz), intent(out) :: curl
real(kind=8), dimension(nx,ny,nz) :: a_work, b_work

!!
call dudx_3d(mm(:,2,:,:,:), w, a_work, nx,ny,nz)! dwdy
call dudx_3d(mm(:,3,:,:,:), v, b_work, nx,ny,nz)! dvdz
curl(1,:,:,:) = a_work - b_work
call dudx_3d(mm(:,1,:,:,:), w, a_work, nx,ny,nz)! dwdx
call dudx_3d(mm(:,3,:,:,:), u, b_work, nx,ny,nz)! dudz
curl(2,:,:,:) = b_work - a_work
call dudx_3d(mm(:,1,:,:,:), w, a_work, nx,ny,nz)! dwdx
call dudx_3d(mm(:,2,:,:,:), u, b_work, nx,ny,nz)! dudy
curl(3,:,:,:) = a_work - b_work

end subroutine




end module


