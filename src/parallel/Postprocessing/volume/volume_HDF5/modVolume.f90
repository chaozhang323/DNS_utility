module modVolume
use MDerivative
use decomp_2d
use decomp_2d_fft
implicit none
type index_output
    integer :: ist, iend, isp
    integer :: jst, jend, jsp
    integer :: kst, kend, ksp
end type
contains 

    subroutine Callambda2(idx,xin,yin,zin,uin,vin,win,Q)
       type(index_output), intent(in) :: idx
       real(8), intent(in) :: xin(:,:,:), yin(:,:,:), zin(:,:,:)
       real(8), intent(in) :: uin(:,:,:), vin(:,:,:), win(:,:,:)
       real(8), intent(out), dimension(:,:,:), allocatable :: Q
       integer :: i,j,k,ilen,jlen,klen,iin,jin,kin
       real(8) :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz

       ilen = (idx%iend-idx%ist)/idx%isp+1
       jlen = (idx%jend-idx%jst)/idx%jsp+1
       klen = (idx%kend-idx%kst)/idx%ksp+1
       if(allocated(Q)) deallocate(Q)
       allocate(Q(ilen,jlen,klen))
       Q = 0.d0
       do k=1,klen
         do j=1,jlen
            do i=1,ilen
               dudx= ddgrid(iin,jin,kin,xin,yin,zin,uin,1)
               dudy= ddgrid(iin,jin,kin,xin,yin,zin,uin,2)
               dudz= ddgrid(iin,jin,kin,xin,yin,zin,uin,3)
               dvdx= ddgrid(iin,jin,kin,xin,yin,zin,vin,1)
               dvdy= ddgrid(iin,jin,kin,xin,yin,zin,vin,2)
               dvdz= ddgrid(iin,jin,kin,xin,yin,zin,vin,3)
               dwdx= ddgrid(iin,jin,kin,xin,yin,zin,win,1)
               dwdy= ddgrid(iin,jin,kin,xin,yin,zin,win,2)
               dwdz= ddgrid(iin,jin,kin,xin,yin,zin,win,3)
               Q(i,j,k) = -0.5*(dudx**2 + dvdy**2 + dwdz**2)-(dvdx*dudy+dwdy*dvdz+dudz*dwdx)
            enddo
          enddo
        enddo

    end subroutine Callambda2



    subroutine CalVorticity(idx,idim,jdim,kdim,xin,yin,zin,uin,vin,win,omx,omy,omz)
       type(index_output), intent(in) :: idx
       integer, intent(in) :: idim, jdim, kdim
       real(8), intent(in) :: xin(1:kdim,-1:idim-2,-1:jdim-2), yin(1:kdim,-1:idim-2,-1:jdim-2), zin(1:kdim,-1:idim-2,-1:jdim-2)
       real(8), intent(in) :: uin(1:kdim,-1:idim-2,-1:jdim-2), vin(1:kdim,-1:idim-2,-1:jdim-2), win(1:kdim,-1:idim-2,-1:jdim-2)
       real(8), intent(out), dimension(:,:,:), allocatable :: omx, omy, omz
       integer :: i,j,k,ilen,jlen,klen,iin,jin,kin
       real(8) :: dudy,dudz,dvdx,dvdz,dwdx,dwdy

!       ilen = (idx%iend-idx%ist)/idx%isp+1
!       jlen = (idx%jend-idx%jst)/idx%jsp+1
!       klen = (idx%kend-idx%kst)/idx%ksp+1
       ilen = idim
       jlen = jdim
       klen = kdim

       if(allocated(omx)) deallocate(omx)
       if(allocated(omy)) deallocate(omy)
       if(allocated(omz)) deallocate(omz)
     !  allocate(omx(ilen,jlen,klen), omy(ilen,jlen,klen),omz(ilen,jlen,klen))
       allocate(omx(klen,ilen-4,jlen-4), omy(klen,ilen-4,jlen-4),omz(klen,ilen-4,jlen-4))

       omx = 0.d0; omy = 0.d0; omz = 0.d0
       do k=1,klen
         do j=1,jlen-4
            do i=1,ilen-4
               iin = i
               jin = j
               kin = k
               dudy= ddgrid_P(iin,jin,kin,idim,jdim,kdim,xin(1:kdim,-1:idim-2,-1:jdim-2),yin(1:kdim,-1:idim-2,-1:jdim-2), &
                                   zin(1:kdim,-1:idim-2,-1:jdim-2),uin(1:kdim,-1:idim-2,-1:jdim-2),2)
               dudz= ddgrid_P(iin,jin,kin,idim,jdim,kdim,xin(1:kdim,-1:idim-2,-1:jdim-2),yin(1:kdim,-1:idim-2,-1:jdim-2), &
                                   zin(1:kdim,-1:idim-2,-1:jdim-2),uin(1:kdim,-1:idim-2,-1:jdim-2),3)
               dvdx= ddgrid_P(iin,jin,kin,idim,jdim,kdim,xin(1:kdim,-1:idim-2,-1:jdim-2),yin(1:kdim,-1:idim-2,-1:jdim-2), &
                                   zin(1:kdim,-1:idim-2,-1:jdim-2),vin(1:kdim,-1:idim-2,-1:jdim-2),1)
               dvdz= ddgrid_P(iin,jin,kin,idim,jdim,kdim,xin(1:kdim,-1:idim-2,-1:jdim-2),yin(1:kdim,-1:idim-2,-1:jdim-2), &
                                   zin(1:kdim,-1:idim-2,-1:jdim-2),vin(1:kdim,-1:idim-2,-1:jdim-2),3)
               dwdx= ddgrid_P(iin,jin,kin,idim,jdim,kdim,xin(1:kdim,-1:idim-2,-1:jdim-2),yin(1:kdim,-1:idim-2,-1:jdim-2), &
                                   zin(1:kdim,-1:idim-2,-1:jdim-2),win(1:kdim,-1:idim-2,-1:jdim-2),1)
               dwdy= ddgrid_P(iin,jin,kin,idim,jdim,kdim,xin(1:kdim,-1:idim-2,-1:jdim-2),yin(1:kdim,-1:idim-2,-1:jdim-2), &
                                   zin(1:kdim,-1:idim-2,-1:jdim-2),win(1:kdim,-1:idim-2,-1:jdim-2),2)
               omx(k,i,j) = dwdy-dvdz
               omy(k,i,j) = dudz-dwdx
               omz(k,i,j) = dvdx-dudy
            enddo
          enddo
        enddo
     end subroutine CalVorticity

    subroutine CalDivergence(idx,idim,jdim,kdim,xin,yin,zin,uin,vin,win,div)
       type(index_output), intent(in) :: idx
       integer, intent(in) :: idim, jdim, kdim
       real(8), intent(in) :: xin(1:kdim,-1:idim-2,-1:jdim-2), yin(1:kdim,-1:idim-2,-1:jdim-2), zin(1:kdim,-1:idim-2,-1:jdim-2)
       real(8), intent(in) :: uin(1:kdim,-1:idim-2,-1:jdim-2), vin(1:kdim,-1:idim-2,-1:jdim-2), win(1:kdim,-1:idim-2,-1:jdim-2)
!       real(8), intent(in) :: xin(:,:,:), yin(:,:,:), zin(:,:,:)
!       real(8), intent(in) :: uin(:,:,:), vin(:,:,:), win(:,:,:)
       real(8), intent(out), dimension(:,:,:), allocatable :: div
       integer :: i, j, k, ilen, jlen, klen, iin, jin, kin
       real(8) :: dudx, dvdy, dwdz

!       ilen = (idx%iend-idx%ist)/idx%isp+1
!       jlen = (idx%jend-idx%jst)/idx%jsp+1
!       klen = (idx%kend-idx%kst)/idx%ksp+1
       ilen = idim
       jlen = jdim
       klen = kdim

       if(allocated(div)) deallocate(div)
!       allocate(div(ilen,jlen,klen))
       allocate(div(klen,ilen-4,jlen-4))
       div = 0.d0
       do k=1,klen
         do j=1,jlen-4
            do i=1,ilen-4
               iin = i
               jin = j
               kin = k
               dudx= ddgrid_P(iin,jin,kin,idim,jdim,kdim,xin(1:kdim,-1:idim-2,-1:jdim-2),yin(1:kdim,-1:idim-2,-1:jdim-2), &
                                   zin(1:kdim,-1:idim-2,-1:jdim-2),uin(1:kdim,-1:idim-2,-1:jdim-2),1)
               dvdy= ddgrid_P(iin,jin,kin,idim,jdim,kdim,xin(1:kdim,-1:idim-2,-1:jdim-2),yin(1:kdim,-1:idim-2,-1:jdim-2), &
                                   zin(1:kdim,-1:idim-2,-1:jdim-2),vin(1:kdim,-1:idim-2,-1:jdim-2),2)
               dwdz= ddgrid_P(iin,jin,kin,idim,jdim,kdim,xin(1:kdim,-1:idim-2,-1:jdim-2),yin(1:kdim,-1:idim-2,-1:jdim-2), &
                                   zin(1:kdim,-1:idim-2,-1:jdim-2),win(1:kdim,-1:idim-2,-1:jdim-2),3)
               div(k,i,j) = dudx + dvdy + dwdz
            enddo
          enddo
        enddo
    end subroutine CalDivergence

    ! parallel version
    subroutine CalGradient(idx,idim,jdim,kdim,xin,yin,zin,vin,grad_v)
       type(index_output), intent(in) :: idx
       integer, intent(in) :: idim, jdim, kdim
       real(8), intent(in) :: xin(1:kdim,-1:idim-2,-1:jdim-2), yin(1:kdim,-1:idim-2,-1:jdim-2), zin(1:kdim,-1:idim-2,-1:jdim-2)
       real(8), intent(in) :: vin(1:kdim,-1:idim-2,-1:jdim-2)
       real(8), intent(out), dimension(:,:,:), allocatable :: grad_v
       integer :: i, j, k, ilen, jlen, klen, iin, jin, kin
       real(8) :: dvdx, dvdy, dvdz

!       ilen = (idx%iend-idx%ist)/idx%isp+1
!       jlen = (idx%jend-idx%jst)/idx%jsp+1
!       klen = (idx%kend-idx%kst)/idx%ksp+1
       ilen = idim
       jlen = jdim
       klen = kdim

       if(allocated(grad_v)) deallocate(grad_v)
       allocate( grad_v(klen,ilen-4,jlen-4) )

       grad_v = 0.d0
       do k=1,klen
         do j=1,jlen-4
            do i=1,ilen-4
               iin = i
               jin = j
               kin = k
               dvdx= ddgrid_P(iin,jin,kin,idim,jdim,kdim,xin(1:kdim,-1:idim-2,-1:jdim-2),yin(1:kdim,-1:idim-2,-1:jdim-2), &
                                        zin(1:kdim,-1:idim-2,-1:jdim-2),vin(1:kdim,-1:idim-2,-1:jdim-2),1)
               dvdy= ddgrid_P(iin,jin,kin,idim,jdim,kdim,xin(1:kdim,-1:idim-2,-1:jdim-2),yin(1:kdim,-1:idim-2,-1:jdim-2), &
                                        zin(1:kdim,-1:idim-2,-1:jdim-2),vin(1:kdim,-1:idim-2,-1:jdim-2),2)
               dvdz= ddgrid_P(iin,jin,kin,idim,jdim,kdim,xin(1:kdim,-1:idim-2,-1:jdim-2),yin(1:kdim,-1:idim-2,-1:jdim-2), &
                                        zin(1:kdim,-1:idim-2,-1:jdim-2),vin(1:kdim,-1:idim-2,-1:jdim-2),3)
               grad_v(k,i,j) = sqrt(dvdx**2+dvdy**2+dvdz**2)
            enddo
          enddo
        enddo

    end subroutine CalGradient

    subroutine CalSwirl(idx,idim,jdim,kdim,xin,yin,zin,uin,vin,win,swirl)
       type(index_output), intent(in) :: idx
       integer, intent(in) :: idim, jdim, kdim
       real(8), intent(in) :: xin(1:kdim,-1:idim-2,-1:jdim-2), yin(1:kdim,-1:idim-2,-1:jdim-2), zin(1:kdim,-1:idim-2,-1:jdim-2)
       real(8), intent(in) :: uin(1:kdim,-1:idim-2,-1:jdim-2), vin(1:kdim,-1:idim-2,-1:jdim-2), win(1:kdim,-1:idim-2,-1:jdim-2)
       real(8), intent(out), dimension(:,:,:), allocatable :: swirl
       integer :: i, j, k, ilen, jlen, klen, iin, jin, kin
       real(8) :: dux,dvx,dwx,duz,dvz,dwz,duy,dvy,dwy
       real :: s11,s12,s13,s22,s23,s33,tripleS
       real :: g11,g12,g13,g21,g22,g23,g31,g32,g33,sigmaS
       real :: PP,Q1,QQ,R1,RR,dis,verde,azul
       real :: swirlmax

!       ilen = (idx%iend-idx%ist)/idx%isp+1
!       jlen = (idx%jend-idx%jst)/idx%jsp+1
!       klen = (idx%kend-idx%kst)/idx%ksp+1
       ilen = idim
       jlen = jdim
       klen = kdim

       if(allocated(swirl)) deallocate(swirl)
!       allocate(swirl(ilen,jlen,klen))
       allocate(swirl(klen,ilen-4,jlen-4))
       swirl = 0.d0
       swirlmax = 0.d0
       do k=1,klen
         do j=1,jlen-4
            do i=1,ilen-4
               iin = i
               jin = j
               kin = k
               dux= ddgrid_P(iin,jin,kin,idim,jdim,kdim,xin(1:kdim,-1:idim-2,-1:jdim-2),yin(1:kdim,-1:idim-2,-1:jdim-2), &
                                   zin(1:kdim,-1:idim-2,-1:jdim-2),uin(1:kdim,-1:idim-2,-1:jdim-2),1)
               dvx= ddgrid_P(iin,jin,kin,idim,jdim,kdim,xin(1:kdim,-1:idim-2,-1:jdim-2),yin(1:kdim,-1:idim-2,-1:jdim-2), &
                                   zin(1:kdim,-1:idim-2,-1:jdim-2),vin(1:kdim,-1:idim-2,-1:jdim-2),1)
               dwx= ddgrid_P(iin,jin,kin,idim,jdim,kdim,xin(1:kdim,-1:idim-2,-1:jdim-2),yin(1:kdim,-1:idim-2,-1:jdim-2), &
                                   zin(1:kdim,-1:idim-2,-1:jdim-2),win(1:kdim,-1:idim-2,-1:jdim-2),1)
               duz= ddgrid_P(iin,jin,kin,idim,jdim,kdim,xin(1:kdim,-1:idim-2,-1:jdim-2),yin(1:kdim,-1:idim-2,-1:jdim-2), &
                                   zin(1:kdim,-1:idim-2,-1:jdim-2),uin(1:kdim,-1:idim-2,-1:jdim-2),3)
               dvz= ddgrid_P(iin,jin,kin,idim,jdim,kdim,xin(1:kdim,-1:idim-2,-1:jdim-2),yin(1:kdim,-1:idim-2,-1:jdim-2), &
                                   zin(1:kdim,-1:idim-2,-1:jdim-2),vin(1:kdim,-1:idim-2,-1:jdim-2),3)
               dwz= ddgrid_P(iin,jin,kin,idim,jdim,kdim,xin(1:kdim,-1:idim-2,-1:jdim-2),yin(1:kdim,-1:idim-2,-1:jdim-2), &
                                   zin(1:kdim,-1:idim-2,-1:jdim-2),win(1:kdim,-1:idim-2,-1:jdim-2),3)
               duy= ddgrid_P(iin,jin,kin,idim,jdim,kdim,xin(1:kdim,-1:idim-2,-1:jdim-2),yin(1:kdim,-1:idim-2,-1:jdim-2), &
                                   zin(1:kdim,-1:idim-2,-1:jdim-2),uin(1:kdim,-1:idim-2,-1:jdim-2),2)
               dvy= ddgrid_P(iin,jin,kin,idim,jdim,kdim,xin(1:kdim,-1:idim-2,-1:jdim-2),yin(1:kdim,-1:idim-2,-1:jdim-2), &
                                   zin(1:kdim,-1:idim-2,-1:jdim-2),vin(1:kdim,-1:idim-2,-1:jdim-2),2)
               dwy= ddgrid_P(iin,jin,kin,idim,jdim,kdim,xin(1:kdim,-1:idim-2,-1:jdim-2),yin(1:kdim,-1:idim-2,-1:jdim-2), &
                                   zin(1:kdim,-1:idim-2,-1:jdim-2),win(1:kdim,-1:idim-2,-1:jdim-2),2)
               PP=-(dux+dvy+dwz)
               s11=dux
               s12=0.5*(duy+dvx)
               s13=0.5*(duz+dwx)
               s22=dvy
               s23=0.5*(dvz+dwy)
               s33=dwz
               tripleS=   s11*(s11*s11+s12*s12+s13*s13) &
                      +   s22*(s12*s12+s22*s22+s23*s23) &
                      +   s33*(s13*s13+s23*s23+s33*s33) &
                      +2.*s12*(s11*s12+s22*s12+s23*s13) &
                      +2.*s13*(s11*s13+s12*s23+s13*s33) &
                      +2.*s23*(s12*s13+s22*s23+s23*s33)
               g11= 0.0
               g12=-0.5*(duy-dvx)
               g13=-0.5*(duz-dwx)
               g21=-0.5*(dvx-duy)
               g22= 0.0
               g23=-0.5*(dvz-dwy)
               g31=-0.5*(dwx-duz)
               g32=-0.5*(dwy-dvz)
               g33= 0.0
               sigmaS=g11*g11*s11+g11*g12*s12+g11*g13*s13 &
                     +g12*g21*s11+g12*g22*s12+g12*g23*s13 &
                     +g13*g31*s11+g13*g32*s12+g13*g33*s13 &
                     +g21*g11*s12+g21*g12*s22+g21*g13*s23 &
                     +g22*g21*s12+g22*g22*s22+g22*g23*s23 &
                     +g23*g31*s12+g23*g32*s22+g23*g33*s23 &
                     +g31*g11*s13+g31*g12*s23+g31*g13*s33 &
                     +g32*g21*s13+g32*g22*s23+g32*g23*s33 &
                     +g33*g31*s13+g33*g32*s23+g33*g33*s33
               Q1=0.5*(PP**2 &
                 -(s11*s11+s22*s22+s33*s33+2.*(s12*s12+s13*s13+s23*s23)) &
                 -(g11*g11+g22*g22+g33*g33+2.*(g12*g21+g13*g31+g23*g32)))
               R1=(-PP**3+3.*PP*Q1-tripleS-3.*sigmaS)/3.
               RR=R1+2.*(PP**3)/27.-PP*Q1/3.
               QQ=Q1-(PP**2)/3.
               dis=0.25*RR**2+(QQ**3)/27.
               if(dis.gt.0.0)then
!                 verde=-0.5*RR+sqrt(dis)
!                 azul =-0.5*RR-sqrt(dis)
                  swirl(k,i,j)=dis
               else
                  swirl(k,i,j)=1.d-30
               endif
               swirlmax=max(swirlmax,swirl(k,i,j))
            enddo
          enddo
        enddo
        swirl = swirl/swirlmax
    end subroutine CalSwirl

    subroutine CalWallStat(idx,xin,yin,zin,uin,tin,dudn,dtdn,ds)
       type(index_output), intent(in) :: idx
       real(8), intent(in) :: xin(:,:,:), yin(:,:,:), zin(:,:,:)
       real(8), intent(in) :: uin(:,:,:), tin(:,:,:)
       real(8), intent(out), dimension(:,:,:), allocatable :: dudn, dtdn, ds
       integer :: i, j, ilen, jlen, iin, jin
       real(8) :: dux, duy, duz, dtx, dty, dtz, len, met(3,3)

       ilen = (idx%iend-idx%ist)/idx%isp+1
       jlen = (idx%jend-idx%jst)/idx%jsp+1
       if(allocated(dudn)) deallocate(dudn)
       if(allocated(dtdn)) deallocate(dtdn)
       if(allocated(ds)) deallocate(ds)
       allocate(dudn(ilen,jlen,1), dtdn(ilen,jlen,1), ds(ilen,jlen,1))
       dudn = 0.d0; dtdn = 0.d0; ds = 0.d0
       do j=1,jlen
          do i=1,ilen
             iin = idx%ist+(i-1)*idx%isp
             jin = idx%jst+(j-1)*idx%jsp
             dux= ddgrid(iin,jin,1,xin,yin,zin,uin,1)
             duy= ddgrid(iin,jin,1,xin,yin,zin,uin,2)
             duz= ddgrid(iin,jin,1,xin,yin,zin,uin,3)
             dtx= ddgrid(iin,jin,1,xin,yin,zin,tin,1)
             dty= ddgrid(iin,jin,1,xin,yin,zin,tin,2)
             dtz= ddgrid(iin,jin,1,xin,yin,zin,tin,3)
!             met = Calmet(Calmm(iin,jin,1,xin,yin,zin))
             len = sqrt(met(3,1)**2+met(3,2)**2+met(3,3)**2)
             dudn(i,j,1) = ( dux*met(3,1)+duy*met(3,2)+duz*met(3,3) )/len
             dtdn(i,j,1) = ( dtx*met(3,1)+dty*met(3,2)+dtz*met(3,3) )/len
             if(met(3,3).ne.0.d0) ds(i,j,1) = len/abs(met(3,3))
          enddo
       enddo

    end subroutine CalWallStat

end module modVolume
