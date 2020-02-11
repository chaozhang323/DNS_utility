module mod_time
implicit none
contains
!! ----------------------------------------- !!
!!  Will only work with Colonius
!! ----------------------------------------- !!
subroutine CalCFL(kapxi, kapeta, kapzeta, &
                  u, v, w, a, mu, rho, gamma_, pr, &
                  cfl, cflv, &
                  dtn )
real(kind=8), dimension(:,:,:,:), intent(in) :: kapxi, kapeta, kapzeta
real(kind=8), dimension(:,:,:), intent(in) :: u, v, w, a, mu, rho
real(kind=8), intent(in) :: gamma_, pr, cfl, cflv
real(kind=8), intent(out) :: dtn
integer i,j,k, ilen,jlen,klen
real(8) :: dtconv,dtvis
real(8) :: fctr,ccc
real(kind=8) :: dti,dtj,dtk,dtvijk,dtcijk
real(8) :: ub,vb,wb
real(8), dimension(4) :: kxi, keta, kzeta
integer :: ktmp
integer, dimension(:), allocatable :: cr

ilen = size(u, 2)
jlen = size(u, 3)
klen = size(u, 1)
print*, "Your flow field is sized )", klen, ilen, jlen, "(jik)"
allocate(cr(klen))
cr = 1

dtconv=1000.
dtvis=1000.
fctr =1.0
ccc=1.0d8

do k=1,klen
    ktmp=k
    do j=1,jlen
    do i=1,ilen
        kxi = kapxi(:,i,j,k)
        kzeta = kapzeta(:,k,i,j)
        keta = kapeta(:,j,i,ktmp)/dble(cr(k))  ! ktmp = 2 for k = 1 at Polar axis; ktmp = k otherwise
        keta = kapeta(:,j,i,k)
        keta(2) = kapeta(2,j,i,k)/cr(k)
        keta(4) = sqrt(keta(1)**2 + keta(2)**2 + keta(3)**2)

        !viscous constrain
        ccc = 1./(kxi(4)*kxi(4)+keta(4)*keta(4)+kzeta(4)*kzeta(4))

        ub=u(k,i,j)*kxi(1)+v(k,i,j)*kxi(2)+w(k,i,j)*kxi(3)
        vb=u(k,i,j)*keta(1)+v(k,i,j)*keta(2)+w(k,i,j)*keta(3)
        wb=u(k,i,j)*kzeta(1)+v(k,i,j)*kzeta(2)+w(k,i,j)*kzeta(3)

        dti = abs(ub)+a(k,i,j)*kxi(4)
        dtj = abs(vb)+a(k,i,j)*keta(4)
        dtk = abs(wb)+a(k,i,j)*kzeta(4)

        dtcijk = 1.d0/(dti+dtj+dtk)
        dtconv = min(dtconv,dtcijk) ! due to convection constraint
        dtvijk = ccc*rho(k,i,j)/(mu(k,i,j))  ! due to viscous constraint
  
        dtvis  = min(dtvis,dtvijk)
        !pr     = (rbar+cv)*mu(k,i,j)/kappa(k,i,j)
        fctr   = max(gamma_/pr,fctr)
    end do
    end do
end do

dtconv=cfl*dtconv
dtvis=cflv*dtvis/fctr
dtn = min(dtconv,dtvis)
dtn = max(dtn,1.d-30) ! Avoid floating point error

return
end subroutine

end module
