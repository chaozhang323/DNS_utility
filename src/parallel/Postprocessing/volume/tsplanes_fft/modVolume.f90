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



    subroutine CalVorticity_tsjplane(idx,imax,kmax,buffer_grid,uin,vin,win,dudj,dvdj,dwdj,omx,omy,omz,ist_tsjplane,kst_tsjplane)
      type(index_output), intent(in) :: idx
      integer, intent(in) :: imax, kmax, ist_tsjplane,kst_tsjplane
      real(8), intent(in) :: buffer_grid(:,:,:,:)
      real(8), intent(in) :: uin(:,:), vin(:,:), win(:,:)
      real(8), intent(in) :: dudj(:,:), dvdj(:,:), dwdj(:,:)
      real(8), intent(out), dimension(:,:,:), allocatable :: omx, omy, omz
      integer :: i,j,k,ilen,jlen,klen,iin,jin,kin
      real(8) :: dudy,dudz,dvdx,dvdz,dwdx,dwdy

       ilen = (idx%iend-idx%ist)/idx%isp+1
       jlen = (idx%jend-idx%jst)/idx%jsp+1
       klen = (idx%kend-idx%kst)/idx%ksp+1
       if(allocated(omx)) deallocate(omx)
       if(allocated(omy)) deallocate(omy)
       if(allocated(omz)) deallocate(omz)
       allocate(omx(ilen,jlen,klen), omy(ilen,jlen,klen),omz(ilen,jlen,klen))
       omx = 0.d0; omy = 0.d0; omz = 0.d0
       do k=1, klen
         do j=1, jlen
           do i=1, ilen
             iin = idx%ist+(i-1)*idx%isp-ist_tsjplane+1
             jin = 1
             kin = idx%kst+(k-1)*idx%ksp-kst_tsjplane+1
             dudy= ddgrid_tsjplane(iin,jin,kin,imax,kmax,buffer_grid(iin,kin,1,1:12),uin,dudj(iin,kin),2)
             dudz= ddgrid_tsjplane(iin,jin,kin,imax,kmax,buffer_grid(iin,kin,1,1:12),uin,dudj(iin,kin),3)
             dvdx= ddgrid_tsjplane(iin,jin,kin,imax,kmax,buffer_grid(iin,kin,1,1:12),vin,dvdj(iin,kin),1)
             dvdz= ddgrid_tsjplane(iin,jin,kin,imax,kmax,buffer_grid(iin,kin,1,1:12),vin,dvdj(iin,kin),3)
             dwdx= ddgrid_tsjplane(iin,jin,kin,imax,kmax,buffer_grid(iin,kin,1,1:12),win,dwdj(iin,kin),1)
             dwdx= ddgrid_tsjplane(iin,jin,kin,imax,kmax,buffer_grid(iin,kin,1,1:12),win,dwdj(iin,kin),2)
             omx(i,j,k) = dwdy-dvdz
             omy(i,j,k) = dudz-dwdx
             omz(i,j,k) = dvdx-dudy
           enddo
         enddo
       enddo

    end subroutine CalVorticity_tsjplane






    subroutine CalVorticity(idx,xin,yin,zin,uin,vin,win,omx,omy,omz)
       type(index_output), intent(in) :: idx
       real(8), intent(in) :: xin(:,:,:), yin(:,:,:), zin(:,:,:)
       real(8), intent(in) :: uin(:,:,:), vin(:,:,:), win(:,:,:)
       real(8), intent(out), dimension(:,:,:), allocatable :: omx, omy, omz
       integer :: i,j,k,ilen,jlen,klen,iin,jin,kin
       real(8) :: dudy,dudz,dvdx,dvdz,dwdx,dwdy

       ilen = (idx%iend-idx%ist)/idx%isp+1
       jlen = (idx%jend-idx%jst)/idx%jsp+1
       klen = (idx%kend-idx%kst)/idx%ksp+1
       if(allocated(omx)) deallocate(omx)
       if(allocated(omy)) deallocate(omy)
       if(allocated(omz)) deallocate(omz)
       allocate(omx(ilen,jlen,klen), omy(ilen,jlen,klen),omz(ilen,jlen,klen))
       omx = 0.d0; omy = 0.d0; omz = 0.d0
       do k=1,klen
         do j=1,jlen
            do i=1,ilen
               iin = idx%ist+(i-1)*idx%isp
               jin = idx%jst+(j-1)*idx%jsp
               kin = idx%kst+(k-1)*idx%ksp
               dudy= ddgrid(iin,jin,kin,xin,yin,zin,uin,2)
               dudz= ddgrid(iin,jin,kin,xin,yin,zin,uin,3)
               dvdx= ddgrid(iin,jin,kin,xin,yin,zin,vin,1)
               dvdz= ddgrid(iin,jin,kin,xin,yin,zin,vin,3)
               dwdx= ddgrid(iin,jin,kin,xin,yin,zin,win,1)
               dwdy= ddgrid(iin,jin,kin,xin,yin,zin,win,2)
               omx(i,j,k) = dwdy-dvdz
               omy(i,j,k) = dudz-dwdx
               omz(i,j,k) = dvdx-dudy
            enddo
          enddo
        enddo
     end subroutine CalVorticity

     subroutine CalDivergence_tsjplane(idx,imax,kmax,buffer_grid,uin,vin,win,dudj,dvdj,dwdj,div,ist_tsjplane,kst_tsjplane)
       type(index_output), intent(in) :: idx
       integer, intent(in) :: imax, kmax, ist_tsjplane,kst_tsjplane
       real(8), intent(in) :: buffer_grid(:,:,:,:)
       real(8), intent(in) :: uin(:,:), vin(:,:), win(:,:)
       real(8), intent(in) :: dudj(:,:), dvdj(:,:), dwdj(:,:)
       real(8), intent(out), dimension(:,:,:), allocatable :: div
       integer :: i,j,k,ilen,jlen,klen,iin,jin,kin
       real(8) :: dudx, dvdy, dwdz

       ilen = (idx%iend-idx%ist)/idx%isp+1
       jlen = (idx%jend-idx%jst)/idx%jsp+1
       klen = (idx%kend-idx%kst)/idx%ksp+1
       if(allocated(div)) deallocate(div)
       allocate(div(ilen,jlen,klen))
       div = 0.d0
       do k=1, klen
         do j=1, jlen
           do i=1, ilen
             iin = idx%ist+(i-1)*idx%isp-ist_tsjplane+1
             jin = 1
             kin = idx%kst+(k-1)*idx%ksp-kst_tsjplane+1
             dudx= ddgrid_tsjplane(iin,jin,kin,imax,kmax,buffer_grid(iin,kin,1,1:12),uin,dudj(iin,kin),1)
             dvdy= ddgrid_tsjplane(iin,jin,kin,imax,kmax,buffer_grid(iin,kin,1,1:12),vin,dvdj(iin,kin),2)
             dwdz= ddgrid_tsjplane(iin,jin,kin,imax,kmax,buffer_grid(iin,kin,1,1:12),win,dwdj(iin,kin),3)
             div(i,j,k) = dudx + dvdy + dwdz
           enddo
         enddo
       enddo

     end subroutine CalDivergence_tsjplane



    subroutine CalDivergence(idx,xin,yin,zin,uin,vin,win,div)
       type(index_output), intent(in) :: idx
       real(8), intent(in) :: xin(:,:,:), yin(:,:,:), zin(:,:,:)
       real(8), intent(in) :: uin(:,:,:), vin(:,:,:), win(:,:,:)
       real(8), intent(out), dimension(:,:,:), allocatable :: div
       integer :: i, j, k, ilen, jlen, klen, iin, jin, kin
       real(8) :: dudx, dvdy, dwdz

       ilen = (idx%iend-idx%ist)/idx%isp+1
       jlen = (idx%jend-idx%jst)/idx%jsp+1
       klen = (idx%kend-idx%kst)/idx%ksp+1
       if(allocated(div)) deallocate(div)
       allocate(div(ilen,jlen,klen))
       div = 0.d0
       do k=1,klen
         do j=1,jlen
            do i=1,ilen
               iin = idx%ist+(i-1)*idx%isp
               jin = idx%jst+(j-1)*idx%jsp
               kin = idx%kst+(k-1)*idx%ksp
               dudx= ddgrid(iin,jin,kin,xin,yin,zin,uin,1)
               dvdy= ddgrid(iin,jin,kin,xin,yin,zin,vin,2)
               dwdz= ddgrid(iin,jin,kin,xin,yin,zin,win,3)
               div(i,j,k) = dudx + dvdy + dwdz
            enddo
          enddo
        enddo
    end subroutine CalDivergence

    subroutine CalGradient(idx,xin,yin,zin,vin,grad_v)
       type(index_output), intent(in) :: idx
       real(8), intent(in) :: xin(:,:,:), yin(:,:,:), zin(:,:,:)
       real(8), intent(in) :: vin(:,:,:)
       real(8), intent(out), dimension(:,:,:), allocatable :: grad_v
       integer :: i, j, k, ilen, jlen, klen, iin, jin, kin
       real(8) :: dvdx, dvdy, dvdz

       ilen = (idx%iend-idx%ist)/idx%isp+1
       jlen = (idx%jend-idx%jst)/idx%jsp+1
       klen = (idx%kend-idx%kst)/idx%ksp+1
       if(allocated(grad_v)) deallocate(grad_v)
       allocate(grad_v(ilen,jlen,klen))
       grad_v = 0.d0
       do k=1,klen
         do j=1,jlen
            do i=1,ilen
               iin = idx%ist+(i-1)*idx%isp
               jin = idx%jst+(j-1)*idx%jsp
               kin = idx%kst+(k-1)*idx%ksp
               dvdx= ddgrid(iin,jin,kin,xin,yin,zin,vin,1)
               dvdy= ddgrid(iin,jin,kin,xin,yin,zin,vin,2)
               dvdz= ddgrid(iin,jin,kin,xin,yin,zin,vin,3)
               grad_v(i,j,k) = sqrt(dvdx**2+dvdy**2+dvdz**2)
            enddo
          enddo
        enddo
    end subroutine CalGradient

!    subroutine CalGradient_tsjplane(idx,x,y,z,didx,djdx,dkdx,didy,djdy,dkdy,didz,djdz,dkdz,vin,dvindj,grad_v,ist_tsjplane,kst_tsjplane)
    subroutine CalGradient_tsjplane(idx,imax,kmax,buffer_grid,vin,dvindj,grad_v,ist_tsjplane,kst_tsjplane)
       type(index_output), intent(in) :: idx
       real(8), intent(in) :: buffer_grid(:,:,:,:)
       real(8), intent(in) :: vin(:,:)                                 !!!!!!!
       real(8), intent(in) :: dvindj(:,:)
       real(8), intent(out), dimension(:,:,:), allocatable :: grad_v
       integer, intent(in) :: imax, kmax
       integer, intent(in) :: ist_tsjplane, kst_tsjplane
       integer :: i, j, k, ilen, jlen, klen, iin, jin, kin
       real(8) :: dvdx, dvdy, dvdz

       ilen = (idx%iend-idx%ist)/idx%isp+1
       jlen = (idx%jend-idx%jst)/idx%jsp+1
       klen = (idx%kend-idx%kst)/idx%ksp+1
       if(allocated(grad_v)) deallocate(grad_v)
       allocate(grad_v(ilen,jlen,klen))
       grad_v = 0.d0

 !   print *, 'dvindj = ', dvindj(1,1)

       do k=1, klen
         do j=1, jlen
           do i=1, ilen
               iin = idx%ist+(i-1)*idx%isp-ist_tsjplane+1
               jin = 1
               kin = idx%kst+(k-1)*idx%ksp-kst_tsjplane+1
               dvdx= ddgrid_tsjplane(iin,jin,kin,imax,kmax,buffer_grid(iin,kin,1,1:12),vin,dvindj(iin,kin),1)
               dvdy= ddgrid_tsjplane(iin,jin,kin,imax,kmax,buffer_grid(iin,kin,1,1:12),vin,dvindj(iin,kin),2)
               dvdz= ddgrid_tsjplane(iin,jin,kin,imax,kmax,buffer_grid(iin,kin,1,1:12),vin,dvindj(iin,kin),3)
               grad_v(i,j,k) = sqrt(dvdx**2+dvdy**2+dvdz**2)
           enddo
         enddo
       enddo
    end subroutine CalGradient_tsjplane


    subroutine CalGradient_tsiplane_P(idx,jmax,kmax,buffer_grid,vin,dvindi,grad_v,jst_tsiplane,kst_tsiplane)
       type(index_output), intent(in) :: idx
       real(8), intent(in) :: buffer_grid(1:jmax,-1:kmax-2,1:1,1:12)
       real(8), intent(in) :: vin(jmax,-1:kmax-2)
       real(8), intent(in) :: dvindi(jmax,-1:kmax-2)
       real(8), intent(out), dimension(:,:,:), allocatable :: grad_v
       integer, intent(in) :: jmax, kmax
       integer, intent(in) :: jst_tsiplane, kst_tsiplane
       integer :: i, j, k, ilen, jlen, klen, iin, jin, kin
       real(8) :: dvdx, dvdz
       real(8), dimension(:,:,:), allocatable :: dvdy
       real(8) :: dvdy1
       real(8) :: Lx



       Lx = buffer_grid(jmax,1,1,2) - buffer_grid(1,1,1,2)

!       ilen = (idx%iend-idx%ist)/idx%isp+1
!       jlen = (idx%jend-idx%jst)/idx%jsp+1
!       klen = (idx%kend-idx%kst)/idx%ksp+1
       ilen = 1
       jlen = jmax
       klen = kmax

       if(allocated(grad_v)) deallocate(grad_v)
       allocate(grad_v(ilen,jlen,klen))
       grad_v = 0.d0
       allocate(dvdy(jmax,kmax,1))

       call ddy_fft1st(1,jmax,kmax,vin,Lx,dvdy)

       do k=1, klen
         do j=1, jlen
           do i=1, ilen
               iin = i
               jin = j
               kin = k
               dvdx  = ddgrid_tsiplane_P(iin,jin,kin,jmax,kmax,buffer_grid(jin,kin,1,1:12),vin,dvindi(jin,kin),1)
               dvdy1 = ddgrid_tsiplane_P(iin,jin,kin,jmax,kmax,buffer_grid(jin,kin,1,1:12),vin,dvindi(jin,kin),2)
               dvdz  = ddgrid_tsiplane_P(iin,jin,kin,jmax,kmax,buffer_grid(jin,kin,1,1:12),vin,dvindi(jin,kin),3)
!               dvdx = ddgrid_tsiplane(iin,jin,kin,jmax,kmax,buffer_grid(jin,kin,1,1:12),vin,dvindi(jin,kin),1)
!               dvdz = ddgrid_tsiplane(iin,jin,kin,jmax,kmax,buffer_grid(jin,kin,1,1:12),vin,dvindi(jin,kin),3)

               grad_v(i,j,k) = sqrt(dvdx**2+dvdy1**2+dvdz**2)
!               grad_v(i,j,k) = sqrt(dvdx**2+dvdy(j,k,1)**2+dvdz**2)
           enddo
         enddo
       enddo


       ! check the grid info
     !  print *, 'nrank = ', nrank
     !  print *, 'size(buffer_grid,1) = ', size(buffer_grid,1)
     !  print *, 'size(buffer_grid,2) = ', size(buffer_grid,2)
     !  print *, 'size(buffer_grid,3) = ', size(buffer_grid,3)
     !  print *, 'buffer_grid(1,1,1,3) = ', buffer_grid(1,1,1,3)
     !  print *, 'buffer_grid(1,2,1,3) = ', buffer_grid(1,2,1,3)
     !  print *, 'buffer_grid(1,0,1,3) = ', buffer_grid(1,0,1,3)
     !  print *, 'buffer_grid(1,-1,1,3) = ', buffer_grid(1,-1,1,3)
     !  print *, 'buffer_grid(1,231,1,3) = ', buffer_grid(1,231,1,3)
     !  print *, 'buffer_grid(1,232,1,3) = ', buffer_grid(1,232,1,3)


    end subroutine CalGradient_tsiplane_P

    subroutine ddy_fft1st(nx,ny,nz,in,Lx,out)
      implicit none
      integer, intent(in) :: nx, ny, nz
      real(8), intent(in) :: Lx
      real(8), dimension(:,:), intent(in) :: in
      real(8), dimension(:,:,:), intent(out) :: out
      real(8), dimension(:,:,:), allocatable :: in_2        !!!!
      complex(8), dimension(:,:,:), allocatable :: out_fft
      real(8), parameter :: pi=4.d0*dble(atan(1.0))
      complex(8), parameter :: ii = (0.0, 1.0)
      integer, dimension(3) :: fft_start, fft_end, fft_size
      integer :: i, j, k
      real(8) :: kk


      allocate(  in_2(size(in,1),size(in,2),1) )
      in_2(:,:,1) = in

      call decomp_2d_fft_init
      call decomp_2d_fft_get_size(fft_start,fft_end,fft_size)
      allocate (out_fft(fft_start(1):fft_end(1), &
                        fft_start(2):fft_end(2), &
                        fft_start(3):fft_end(3)))
      call decomp_2d_fft_3d(in_2,out_fft)

      do i=fft_start(1), fft_end(1)
        kk = 2.d0*pi*(i-1)/Lx
        out_fft(i,:,:) = out_fft(i,:,:)*ii*kk
      enddo

      call decomp_2d_fft_3d(out_fft,out)

      ! normalization
      out = out/real(nx)/real(ny)/real(nz)

      call decomp_2d_fft_finalize

    end subroutine ddy_fft1st













    subroutine CalSwirl(idx,xin,yin,zin,uin,vin,win,swirl)
       type(index_output), intent(in) :: idx
       real(8), intent(in) :: xin(:,:,:), yin(:,:,:), zin(:,:,:)
       real(8), intent(in) :: uin(:,:,:), vin(:,:,:), win(:,:,:)
       real(8), intent(out), dimension(:,:,:), allocatable :: swirl
       integer :: i, j, k, ilen, jlen, klen, iin, jin, kin
       real(8) :: dux,dvx,dwx,duz,dvz,dwz,duy,dvy,dwy
       real :: s11,s12,s13,s22,s23,s33,tripleS
       real :: g11,g12,g13,g21,g22,g23,g31,g32,g33,sigmaS
       real :: PP,Q1,QQ,R1,RR,dis,verde,azul
       real :: swirlmax

       ilen = (idx%iend-idx%ist)/idx%isp+1
       jlen = (idx%jend-idx%jst)/idx%jsp+1
       klen = (idx%kend-idx%kst)/idx%ksp+1
       if(allocated(swirl)) deallocate(swirl)
       allocate(swirl(ilen,jlen,klen))
       swirl = 0.d0
       swirlmax = 0.d0
       do k=1,klen
         do j=1,jlen
            do i=1,ilen
               iin = idx%ist+(i-1)*idx%isp
               jin = idx%jst+(j-1)*idx%jsp
               kin = idx%kst+(k-1)*idx%ksp
               dux= ddgrid(iin,jin,kin,xin,yin,zin,uin,1)
               dvx= ddgrid(iin,jin,kin,xin,yin,zin,vin,1)
               dwx= ddgrid(iin,jin,kin,xin,yin,zin,win,1)
               duz= ddgrid(iin,jin,kin,xin,yin,zin,uin,3)
               dvz= ddgrid(iin,jin,kin,xin,yin,zin,vin,3)
               dwz= ddgrid(iin,jin,kin,xin,yin,zin,win,3)
               duy= ddgrid(iin,jin,kin,xin,yin,zin,uin,2)
               dvy= ddgrid(iin,jin,kin,xin,yin,zin,vin,2)
               dwy= ddgrid(iin,jin,kin,xin,yin,zin,win,2)
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
                  swirl(i,j,k)=dis
               else
                  swirl(i,j,k)=1.d-30
               endif
               swirlmax=max(swirlmax,swirl(i,j,k))
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
             met = Calmet(Calmm(iin,jin,1,xin,yin,zin))
             len = sqrt(met(3,1)**2+met(3,2)**2+met(3,3)**2)
             dudn(i,j,1) = ( dux*met(3,1)+duy*met(3,2)+duz*met(3,3) )/len
             dtdn(i,j,1) = ( dtx*met(3,1)+dty*met(3,2)+dtz*met(3,3) )/len
             if(met(3,3).ne.0.d0) ds(i,j,1) = len/abs(met(3,3))
          enddo
       enddo

    end subroutine CalWallStat

end module modVolume
