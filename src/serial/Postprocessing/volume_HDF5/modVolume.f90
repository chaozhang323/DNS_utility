module modVolume
use MDerivative
implicit none
type index_output
    integer :: ist, iend, isp
    integer :: jst, jend, jsp
    integer :: kst, kend, ksp
end type
contains 

    subroutine Calp0(idx,uin,vin,win,pin,tin,p0)
       type(index_output), intent(in) :: idx
       real(8), intent(in) :: uin(:,:,:), vin(:,:,:), win(:,:,:)
       real(8), intent(in) :: pin(:,:,:), tin(:,:,:)
       real(8), intent(out), dimension(:,:,:), allocatable :: p0
       integer :: i,j,k,ilen,jlen,klen,iin,jin,kin
       real(8) :: gamma =1.4, rbar = 287.d0
       real(8) :: tmp

       ilen = (idx%iend-idx%ist)/idx%isp+1
       jlen = (idx%jend-idx%jst)/idx%jsp+1
       klen = (idx%kend-idx%kst)/idx%ksp+1
       if(allocated(p0)) deallocate(p0)
       allocate(p0(ilen,jlen,klen))
       p0 = 0.d0

       do k=1, klen
         do j=1, jlen
           do i=1, ilen
             iin = idx%ist+(i-1)*idx%isp
             jin = idx%jst+(j-1)*idx%jsp
             kin = idx%kst+(k-1)*idx%ksp
             tmp = uin(i,j,k)**2 + vin(i,j,k)**2 + win(i,j,k)**2
             p0(i,j,k) = pin(i,j,k)*( 1.d0 + (gamma-1.d0)/2.d0*(tmp/(gamma*rbar*tin(i,j,k))))**(gamma/(gamma-1.d0))
           enddo
         enddo
       enddo


    end subroutine Calp0


    ! ref: JFM: On the relationships between local vortex identification schemes
    ! By Pinaki Chakraborty, S. Balachandar and Ronald J. Adrian
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
               iin = idx%ist+(i-1)*idx%isp
               jin = idx%jst+(j-1)*idx%jsp
               kin = idx%kst+(k-1)*idx%ksp
               dudx= ddgrid(iin,jin,kin,xin,yin,zin,uin,1)
               dudy= ddgrid(iin,jin,kin,xin,yin,zin,uin,2)
               dudz= ddgrid(iin,jin,kin,xin,yin,zin,uin,3)
               dvdx= ddgrid(iin,jin,kin,xin,yin,zin,vin,1)
               dvdy= ddgrid(iin,jin,kin,xin,yin,zin,vin,2)
               dvdz= ddgrid(iin,jin,kin,xin,yin,zin,vin,3)
               dwdx= ddgrid(iin,jin,kin,xin,yin,zin,win,1)
               dwdy= ddgrid(iin,jin,kin,xin,yin,zin,win,2)
               dwdz= ddgrid(iin,jin,kin,xin,yin,zin,win,3)
               !Q(i,j,k) = -0.5*(dudx**2 + dvdy**2 + dwdz**2)-(dvdx*dudy+dwdy*dvdz+dudz*dwdx)
               Q(i,j,k) = 0.5*( (dudx+dvdy+dwdz)**2 - (dudx**2+dvdy**2+dwdz**2+2.d0*dudy*dvdx+2.d0*dudz*dwdx+2.d0*dvdz*dwdy) )
            enddo
          enddo
        enddo

    end subroutine Callambda2



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
       if(size(xin,dim=2).ne.1) then
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
        else
          do k=1,klen
            do j=1,jlen
              do i=1,ilen
                iin = idx%ist+(i-1)*idx%isp
                kin = idx%kst+(k-1)*idx%ksp
                dudz = ddgrid_2D( iin,kin,xin(:,1,:),zin(:,1,:),uin(:,1,:),2 )
                dwdx = ddgrid_2D( iin,kin,xin(:,1,:),zin(:,1,:),win(:,1,:),1 )
                omy(i,j,k) = dudz - dwdx
              enddo
            enddo
          enddo
        endif

     end subroutine CalVorticity


    subroutine CalVorticity_cylin(idx,xin,yin,zin,uin,vin,win,omx,omy,omz)
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
               omx(i,j,k) = dwdy/(zin(iin,jin,kin)+1.e-30) - dvdz - vin(iin,jin,kin)/(zin(iin,jin,kin)+1.e-30)
               omy(i,j,k) = dudz-dwdx
               omz(i,j,k) = dvdx-dudy/(zin(iin,jin,kin)+1.e-30)
            enddo
          enddo
        enddo
     end subroutine CalVorticity_cylin

    subroutine Calsrc_Philip(idx,xin,yin,zin,uin,vin,win,src)
       type(index_output), intent(in) :: idx
       real(8), intent(in) :: xin(:,:,:), yin(:,:,:), zin(:,:,:)
       real(8), intent(in) :: uin(:,:,:), vin(:,:,:), win(:,:,:)
       real(8), intent(out), dimension(:,:,:), allocatable :: src
       integer :: i, j, k, ilen, jlen, klen, iin, jin, kin
       real(8) :: dudx, dvdy, dwdz, dudy, dvdx, dudz, dwdx, dvdz, dwdy

       ilen = (idx%iend-idx%ist)/idx%isp+1
       jlen = (idx%jend-idx%jst)/idx%jsp+1
       klen = (idx%kend-idx%kst)/idx%ksp+1
       if(allocated(src)) deallocate(src)
       allocate(src(ilen,jlen,klen))
       src = 0.d0
       do k=1,klen
         do j=1,jlen
            do i=1,ilen
               iin = idx%ist+(i-1)*idx%isp
               jin = idx%jst+(j-1)*idx%jsp
               kin = idx%kst+(k-1)*idx%ksp
               dudx= ddgrid(iin,jin,kin,xin,yin,zin,uin,1)
               dvdy= ddgrid(iin,jin,kin,xin,yin,zin,vin,2)
               dwdz= ddgrid(iin,jin,kin,xin,yin,zin,win,3)
               dudy= ddgrid(iin,jin,kin,xin,yin,zin,uin,2)
               dvdx= ddgrid(iin,jin,kin,xin,yin,zin,vin,1)
               dudz= ddgrid(iin,jin,kin,xin,yin,zin,uin,3)
               dwdx= ddgrid(iin,jin,kin,xin,yin,zin,win,1)
               dvdz= ddgrid(iin,jin,kin,xin,yin,zin,vin,3)
               dwdy= ddgrid(iin,jin,kin,xin,yin,zin,win,2)
               src(i,j,k) = dudx**2 + dvdy**2 + dwdz**2 + 2.d0*(dudy*dvdx + dudz*dwdx + dvdz*dwdy)
            enddo
          enddo
        enddo
    end subroutine Calsrc_Philip


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

       if(size(xin,dim=2).eq.1) then
          do k=1,klen
            do j=1,jlen
              do i=1,ilen
                iin = idx%ist+(i-1)*idx%isp
                kin = idx%kst+(k-1)*idx%ksp
                !dudx= ddgrid(iin,jin,kin,xin,yin,zin,uin,1)
                !dwdz= ddgrid(iin,jin,kin,xin,yin,zin,win,3)
                dudx= ddgrid_2D(iin,kin,xin(:,1,:),zin(:,1,:),uin(:,1,:),1)
                dwdz= ddgrid_2D(iin,kin,xin(:,1,:),zin(:,1,:),win(:,1,:),2)
                div(i,j,k) = dudx + dwdz
               enddo
            enddo
          enddo
       elseif(size(xin,dim=3).eq.1) then
          do k=1,klen
            do j=1,jlen
              do i=1,ilen
                iin = idx%ist+(i-1)*idx%isp
                jin = idx%jst+(j-1)*idx%jsp
                !dudx= ddgrid(iin,jin,kin,xin,yin,zin,uin,1)
                !dwdz= ddgrid(iin,jin,kin,xin,yin,zin,win,3)
                dudx= ddgrid_2D(iin,jin,xin(:,:,1),yin(:,:,1),uin(:,:,1),1)
                dvdy= ddgrid_2D(iin,jin,xin(:,:,1),yin(:,:,1),vin(:,:,1),2)
                div(i,j,k) = dudx + dvdy
               enddo
            enddo
          enddo
       else
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
        endif

    end subroutine CalDivergence

    subroutine CalDivergence_cylin(idx,xin,yin,zin,uin,vin,win,div)
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
               div(i,j,k) = dudx + dvdy/(zin(iin,jin,kin)+1.e-30) + dwdz + win(iin,jin,kin)/(zin(iin,jin,kin)+1.e-30)
            enddo
          enddo
        enddo
    end subroutine CalDivergence_cylin

    !compute laplacian of density field for shocklet detection
    subroutine CalLaplacianGCC(idx,xin,yin,zin,rhoin,uin,vin,win,lap,shockloc)
        type(index_output), intent(in) :: idx
        real(8), intent(in) :: xin(:,:,:), yin(:,:,:), zin(:,:,:)
        real(8), intent(in) :: rhoin(:,:,:), uin(:,:,:), vin(:,:,:), win(:,:,:)
        real(8), intent(out), dimension(:,:,:), allocatable :: lap, shockloc
        integer, dimension(:,:,:), allocatable :: shock1, shock2
        integer :: i,j,k,ilen,jlen,klen,iin,jin,kin
        real(8), dimension(:,:,:), allocatable :: drhodx, drhody, drhodz, div
        real(8), dimension(:,:), allocatable :: divave, div2ave, divrms
        real(8) :: d2pdx, d2pdy, d2pdz, dudx, dvdy, dwdz

        ilen = (idx%iend-idx%ist)/idx%isp+1
        jlen = (idx%jend-idx%jst)/idx%jsp+1
        klen = (idx%kend-idx%kst)/idx%ksp+1
        if(allocated(lap)) deallocate(lap)
        allocate(lap(ilen,jlen,klen))
        if(allocated(drhodx)) deallocate(drhodx)
        allocate(drhodx(ilen,jlen,klen))
        if(allocated(drhody)) deallocate(drhody)
        allocate(drhody(ilen,jlen,klen))
        if(allocated(drhodz)) deallocate(drhodz)
        allocate(drhodz(ilen,jlen,klen))
        
        allocate(div(ilen,jlen,klen),shock1(ilen,jlen,klen),shock2(ilen,jlen,klen),shockloc(ilen,jlen,klen))
        allocate(divave(ilen,klen),div2ave(ilen,klen),divrms(ilen,klen))

        do k=1,klen
          do j=1,jlen
            do i=1,ilen
              iin = idx%ist+(i-1)*idx%isp
              jin = idx%jst+(j-1)*idx%jsp
              kin = idx%kst+(k-1)*idx%ksp
              drhodx(i,j,k) = ddgrid(iin,jin,kin,xin,yin,zin,rhoin,1)
              drhody(i,j,k) = ddgrid(iin,jin,kin,xin,yin,zin,rhoin,2)
              drhodz(i,j,k) = ddgrid(iin,jin,kin,xin,yin,zin,rhoin,3)
              dudx= ddgrid(iin,jin,kin,xin,yin,zin,uin,1)
              dvdy= ddgrid(iin,jin,kin,xin,yin,zin,vin,2)
              dwdz= ddgrid(iin,jin,kin,xin,yin,zin,win,3)
              div(i,j,k) = dudx + dvdy + dwdz
            enddo
          enddo
        enddo        

        do k=1,klen
          do j=1,jlen
            do i=1,ilen
              iin=idx%ist+(i-1)*idx%isp
              jin=idx%jst+(j-1)*idx%jsp
              kin=idx%kst+(k-1)*idx%ksp
              d2pdx = ddgrid(iin,jin,kin,xin,yin,zin,drhodx,1)
              d2pdy = ddgrid(iin,jin,kin,xin,yin,zin,drhody,2)
              d2pdz = ddgrid(iin,jin,kin,xin,yin,zin,drhodz,3)
              lap(i,j,k) = d2pdx + d2pdy + d2pdz
            enddo
           enddo
          enddo
          do i=1, ilen
            do k=1, klen
              divave(1:ilen,1:klen)  = sum(div(i,1:jlen,k))
              div2ave(1:ilen,1:klen) = sum(div(i,1:jlen,k)**2)
            enddo
          enddo
          divave  = divave/dble(jlen)
          div2ave = div2ave/dble(jlen)
          divrms = sqrt(abs(div2ave-divave**2)+1.e-30)

          shock1 = 0; shock2 = 0; shockloc = 0
          do j=1, jlen
            do k=1, klen
              do i=1, ilen
                if(-div(i,j,k).gt.3.*divrms(i,k)) then
                  shock2(i,j,k) = 1
                endif
              enddo
              do i=1, ilen-1
                if((lap(i,j,k).gt.0.and.lap(i+1,j,k).le.0).or.(lap(i,j,k).le.0.and.lap(i+1,j,k).gt.0) ) then
                  shock1(i,j,k) = 1
                endif
              enddo
            enddo
          enddo

          do i=1, ilen
            do j=1, jlen
              do k=1, klen
                if(shock1(i,j,k).eq.1.and.shock2(i,j,k).eq.1) then
                  shockloc(i,j,k) = 1
                endif
              enddo
            enddo
          enddo

        print *, "finished modVolume--CalLaplacian"
    end subroutine CalLaplacianGCC

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

       if(size(xin,dim=2).eq.1) then
         do k=1,klen
           do j=1,jlen
             do i=1,ilen
               iin = idx%ist+(i-1)*idx%isp
               kin = idx%kst+(k-1)*idx%ksp
               dvdx= ddgrid_2D(iin,kin,xin(:,1,:),zin(:,1,:),vin(:,1,:),1)
               dvdz= ddgrid_2D(iin,kin,xin(:,1,:),zin(:,1,:),vin(:,1,:),2)
               grad_v(i,j,k) = sqrt(dvdx**2+ dvdz**2)
             enddo
           enddo
         enddo
       elseif(size(xin,dim=3).eq.1) then
         do k=1,klen
           do j=1,jlen
             do i=1,ilen
               iin = idx%ist+(i-1)*idx%isp
               jin = idx%jst+(j-1)*idx%jsp
               dvdx= ddgrid_2D(iin,jin,xin(:,:,1),yin(:,:,1),vin(:,:,1),1)
               dvdy= ddgrid_2D(iin,jin,xin(:,:,1),yin(:,:,1),vin(:,:,1),2)
               grad_v(i,j,k) = sqrt(dvdx**2+ dvdy**2)
             enddo
           enddo
         enddo
       else
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

       endif

    end subroutine CalGradient


    subroutine CalGradient_cylin(idx,xin,yin,zin,vin,grad_v)
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
               dvdy= ddgrid(iin,jin,kin,xin,yin,zin,vin,2)/(zin(iin,jin,kin)+1.e-30)
               dvdz= ddgrid(iin,jin,kin,xin,yin,zin,vin,3)
               grad_v(i,j,k) = sqrt(dvdx**2+dvdy**2+dvdz**2)
            enddo
          enddo
        enddo
    end subroutine CalGradient_cylin


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
