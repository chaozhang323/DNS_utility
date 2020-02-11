    module MTSStat_iplane
      use MTSHDF5
      implicit none
      integer, parameter :: nv2d_ipInt = 13, nv2d_ipStat = 85+5+5+4
      real(8), parameter, private :: gamma = 1.4, Prandtl = 0.71  ! (cp = 3.5*rbar)
      real(8) :: rbar, cp
!      real(8) :: uinf, rhoinf

      type tp_stat_iplane
        integer :: nsample
        real(8), dimension(:,:), allocatable :: u, v, w, p, t
        real(8), dimension(:,:), allocatable :: dudx, dvdx, dwdx, dpdx, dtdx
        real(8), dimension(:,:), allocatable :: dudy, dvdy, dwdy, dpdy, dtdy
        real(8), dimension(:,:), allocatable :: dudz, dvdz, dwdz, dpdz, dtdz
        real(8), dimension(:,:), allocatable :: dudt, dvdt, dwdt, dpdt, dtdt

        real(8), dimension(:,:), allocatable :: dudi, dvdi, dwdi, dpdi, dtdi
        real(8), dimension(:,:), allocatable :: dudj, dvdj, dwdj, dpdj, dtdj
        real(8), dimension(:,:), allocatable :: dudk, dvdk, dwdk, dpdk, dtdk

      end type tp_stat_iplane

      type(tp_stat_iplane), private :: TSStat
      character(30), private :: dname2d_ipInt(nv2d_ipInt+1), dname2d_ipStat(nv2d_ipStat+1)
      parameter(dname2d_ipInt = (/'uave','vave','wave','pave','tave','rhoave','u2','v2','w2','p2','t2','rho2','tauw','z'/))
      parameter(dname2d_ipStat = (/'uave','vave','wave','pave','tave', 'rhoave', 'dpdx','dpdt','dpdxdpdt','divave','omxave','omyave','omzave','magomave','save',&
                           't0ave','p0ave','u2','v2','w2','p2','t2','rho2','dpdx2', 'dpdt2', 'div2', 'omx2', 'omy2', 'omz2', 'magom2', 's2', 't02', 'p02',&
                           'u3', 'u4', 'p3', 'p4', 'uv', 'uw', 'vw', 'ut', 'vt', 'wt', 'up', 'vp', 'wp', 'pt', 'pr', 'ps', 'som', 'pom', 'ruT0', 'ru', &
                           'rv', 'rw', 'ru2', 'rv2', 'rw2', 'ruv', 'ruw', 'rvw', 'r2u2', 'rut', 'rvt', 'rwt', 'uvd', 'dudi', 'dvdi', 'dwdi', 'dpdi', 'dtdi', &
                           'dudj', 'dvdj', 'dwdj', 'dpdj', 'dtdj', 'dudk', 'dvdk', 'dwdk', 'dpdk', 'dtdk', 'dudt', 'dvdt', 'dwdt', 'dtdt', 'dudxdudt', 'dudt2', &
                           'dudx2', 'dvdy2', 'dwdz2', 'dudxdvdy', 'dudxdwdz', 'dvdydwdz', 'M', 'M2', 'rho3', 'rho4', 't3', 't4', 'z'/) )

      real(8), dimension(:), allocatable, private :: dudn, dtdn
      real(8), dimension(:,:), allocatable, private :: mu
      type(tp_hyperslab) :: TSiout

   contains

      subroutine InitTSStat_iplane(ns,kdim,rbar_in)
        integer, intent(in) :: ns, kdim
        real(8), intent(in) :: rbar_in
         
        rbar = rbar_in
        cp = rbar*gamma/(gamma-1.0)

        TSiout%fname = 'none'

        TSStat%nsample = ns
        allocate(TSStat%u(ns,kdim))
        allocate(TSStat%v(ns,kdim))
        allocate(TSStat%w(ns,kdim))
        allocate(TSStat%p(ns,kdim))
        allocate(TSStat%t(ns,kdim))
        allocate(TSStat%dudx(ns,kdim), TSStat%dudy(ns,kdim), TSStat%dudz(ns,kdim))
        allocate(TSStat%dvdx(ns,kdim), TSStat%dvdy(ns,kdim), TSStat%dvdz(ns,kdim))
        allocate(TSStat%dwdx(ns,kdim), TSStat%dwdy(ns,kdim), TSStat%dwdz(ns,kdim))
        allocate(TSStat%dpdx(ns,kdim), TSStat%dpdt(ns,kdim) )
        allocate(TSStat%dudt(ns,kdim), TSStat%dvdt(ns,kdim), TSStat%dwdt(ns,kdim))
        allocate(TSStat%dtdt(ns,kdim))

        allocate(TSStat%dudi(ns,kdim), TSStat%dvdi(ns,kdim), TSStat%dwdi(ns,kdim), TSStat%dpdi(ns,kdim), TSStat%dtdi(ns,kdim))
        allocate(TSStat%dudj(ns,kdim), TSStat%dvdj(ns,kdim), TSStat%dwdj(ns,kdim), TSStat%dpdj(ns,kdim), TSStat%dtdj(ns,kdim))
        allocate(TSStat%dudk(ns,kdim), TSStat%dvdk(ns,kdim), TSStat%dwdk(ns,kdim), TSStat%dpdk(ns,kdim), TSStat%dtdk(ns,kdim))

        allocate(dudn(TSStat%nsample), dtdn(TSStat%nsample))
        allocate(mu(TSStat%nsample,kdim))

      end subroutine


      subroutine CalTSDeriv_iplane(nt,ny,nz,nv,buffer,buffer_grid,dt)
        integer, intent(in) :: nt, ny, nz, nv
!        real(8), intent(in) :: buffer(nt,ny,nz,nv), buffer_grid(ny,9),dt
        real(8), intent(in) :: buffer(nt,ny,nz,nv), buffer_grid(ny,nz,9),dt
        integer :: i, j, k, n, nn
        real(8) :: didx, djdx, dkdx, didy, djdy, dkdy, didz, djdz, dkdz
!        real(8) :: dudi, dvdi, dwdi, dudj, dvdj, dwdj, dudk, dvdk, dwdk
!        real(8) :: dpdi, dpdj, dpdk
        real(8), parameter :: c12i = 1.d0/12.d0
        real(8) :: dti

        dti = 1.d0/dt

        do k=1, nz
          TSStat%u(:,k) = reshape(buffer(3:nt-2,3:ny-2,k,1),(/TSStat%nsample/))
          TSStat%v(:,k) = reshape(buffer(3:nt-2,3:ny-2,k,2),(/TSStat%nsample/))
          TSStat%w(:,k) = reshape(buffer(3:nt-2,3:ny-2,k,3),(/TSStat%nsample/))
          TSStat%p(:,k) = reshape(buffer(3:nt-2,3:ny-2,k,4),(/TSStat%nsample/))
          TSStat%t(:,k) = reshape(buffer(3:nt-2,3:ny-2,k,5),(/TSStat%nsample/))
        enddo

        do i=1, TSStat%nsample
          do k=1, nz
            mu(i,k) = CalMu(TSStat%t(i,k))
          enddo
        enddo

        do k=1, nz
          nn = 0
          do j=3, ny-2
            ! didx, djdx, dkdx, didy, djdy, dkdy, didz, djdz, dkdz
            didx = buffer_grid(j,k,1)
            djdx = buffer_grid(j,k,2)
            dkdx = buffer_grid(j,k,3)
            didy = buffer_grid(j,k,4)
            djdy = buffer_grid(j,k,5)
            dkdy = buffer_grid(j,k,6)
            didz = buffer_grid(j,k,7)
            djdz = buffer_grid(j,k,8)
            dkdz = buffer_grid(j,k,9)

            do n=3, nt-2
              nn = nn + 1

              TSStat%dudi(nn,k) = buffer(n,j,k,6)
              TSStat%dvdi(nn,k) = buffer(n,j,k,7)
              TSStat%dwdi(nn,k) = buffer(n,j,k,8)
              TSStat%dpdi(nn,k) = buffer(n,j,k,9)
              TSStat%dtdi(nn,k) = buffer(n,j,k,10)

              TSStat%dudj(nn,k) = ( buffer(n,j-2,k,1) - 8.0*(buffer(n,j-1,k,1) - buffer(n,j+1,k,1)) - buffer(n,j+2,k,1) )*c12i
              TSStat%dvdj(nn,k) = ( buffer(n,j-2,k,2) - 8.0*(buffer(n,j-1,k,2) - buffer(n,j+1,k,2)) - buffer(n,j+2,k,2) )*c12i
              TSStat%dwdj(nn,k) = ( buffer(n,j-2,k,3) - 8.0*(buffer(n,j-1,k,3) - buffer(n,j+1,k,3)) - buffer(n,j+2,k,3) )*c12i
              TSStat%dpdj(nn,k) = ( buffer(n,j-2,k,4) - 8.0*(buffer(n,j-1,k,4) - buffer(n,j+1,k,4)) - buffer(n,j+2,k,4) )*c12i
              TSStat%dtdj(nn,k) = ( buffer(n,j-2,k,5) - 8.0*(buffer(n,j-1,k,5) - buffer(n,j+1,k,5)) - buffer(n,j+2,k,5) )*c12i

              if(k.ge.3.and.k.le.nz-2) then
                TSStat%dudk(nn,k) = (buffer(n,j,k-2,1) - 8.0*(buffer(n,j,k-1,1) - buffer(n,j,k+1,1)) - buffer(n,j,k+2,1))*c12i
                TSStat%dvdk(nn,k) = (buffer(n,j,k-2,2) - 8.0*(buffer(n,j,k-1,2) - buffer(n,j,k+1,2)) - buffer(n,j,k+2,2))*c12i
                TSStat%dwdk(nn,k) = (buffer(n,j,k-2,3) - 8.0*(buffer(n,j,k-1,3) - buffer(n,j,k+1,3)) - buffer(n,j,k+2,3))*c12i
                TSStat%dpdk(nn,k) = (buffer(n,j,k-2,4) - 8.0*(buffer(n,j,k-1,4) - buffer(n,j,k+1,4)) - buffer(n,j,k+2,4))*c12i
                TSStat%dtdk(nn,k) = (buffer(n,j,k-2,5) - 8.0*(buffer(n,j,k-1,5) - buffer(n,j,k+1,5)) - buffer(n,j,k+2,5))*c12i
              elseif(k.eq.1) then
                TSStat%dudk(nn,k) = -25./12.*buffer(n,j,1,1) + 4.*buffer(n,j,2,1) - 3.*buffer(n,j,3,1) + 4./3.*buffer(n,j,4,1) - 0.25*buffer(n,j,5,1)
                TSStat%dvdk(nn,k) = -25./12.*buffer(n,j,1,2) + 4.*buffer(n,j,2,2) - 3.*buffer(n,j,3,2) + 4./3.*buffer(n,j,4,2) - 0.25*buffer(n,j,5,2)
                TSStat%dwdk(nn,k) = -25./12.*buffer(n,j,1,3) + 4.*buffer(n,j,2,3) - 3.*buffer(n,j,3,3) + 4./3.*buffer(n,j,4,3) - 0.25*buffer(n,j,5,3)
                TSStat%dpdk(nn,k) = -25./12.*buffer(n,j,1,4) + 4.*buffer(n,j,2,4) - 3.*buffer(n,j,3,4) + 4./3.*buffer(n,j,4,4) - 0.25*buffer(n,j,5,4)
                TSStat%dtdk(nn,k) = -25./12.*buffer(n,j,1,5) + 4.*buffer(n,j,2,5) - 3.*buffer(n,j,3,5) + 4./3.*buffer(n,j,4,5) - 0.25*buffer(n,j,5,5)
              elseif(k.eq.2) then
                TSStat%dudk(nn,k) = -0.25*buffer(n,j,1,1) - 5./6.*buffer(n,j,2,1) + 1.5*buffer(n,j,3,1) - 0.5*buffer(n,j,4,1) + 1./12.*buffer(n,j,5,1)
                TSStat%dvdk(nn,k) = -0.25*buffer(n,j,1,2) - 5./6.*buffer(n,j,2,2) + 1.5*buffer(n,j,3,2) - 0.5*buffer(n,j,4,2) + 1./12.*buffer(n,j,5,2)
                TSStat%dwdk(nn,k) = -0.25*buffer(n,j,1,3) - 5./6.*buffer(n,j,2,3) + 1.5*buffer(n,j,3,3) - 0.5*buffer(n,j,4,3) + 1./12.*buffer(n,j,5,3)
                TSStat%dpdk(nn,k) = -0.25*buffer(n,j,1,4) - 5./6.*buffer(n,j,2,4) + 1.5*buffer(n,j,3,4) - 0.5*buffer(n,j,4,4) + 1./12.*buffer(n,j,5,4)
                TSStat%dtdk(nn,k) = -0.25*buffer(n,j,1,5) - 5./6.*buffer(n,j,2,5) + 1.5*buffer(n,j,3,5) - 0.5*buffer(n,j,4,5) + 1./12.*buffer(n,j,5,5)
              elseif(k.eq.nz-1) then
                TSStat%dudk(nn,k) = -1./12.*buffer(n,j,k-4,1) + 0.5*buffer(n,j,k-3,1) - 1.5*buffer(n,j,k-2,1) + 5./6.*buffer(n,j,k-1,1) + 0.25*buffer(n,j,k,1)
                TSStat%dvdk(nn,k) = -1./12.*buffer(n,j,k-4,2) + 0.5*buffer(n,j,k-3,2) - 1.5*buffer(n,j,k-2,2) + 5./6.*buffer(n,j,k-1,2) + 0.25*buffer(n,j,k,2)
                TSStat%dwdk(nn,k) = -1./12.*buffer(n,j,k-4,3) + 0.5*buffer(n,j,k-3,3) - 1.5*buffer(n,j,k-2,3) + 5./6.*buffer(n,j,k-1,3) + 0.25*buffer(n,j,k,3)
                TSStat%dpdk(nn,k) = -1./12.*buffer(n,j,k-4,4) + 0.5*buffer(n,j,k-3,4) - 1.5*buffer(n,j,k-2,4) + 5./6.*buffer(n,j,k-1,4) + 0.25*buffer(n,j,k,4)
                TSStat%dtdk(nn,k) = -1./12.*buffer(n,j,k-4,5) + 0.5*buffer(n,j,k-3,5) - 1.5*buffer(n,j,k-2,5) + 5./6.*buffer(n,j,k-1,5) + 0.25*buffer(n,j,k,5)
              elseif(k.eq.nz) then
                TSStat%dudk(nn,k) = 0.25*buffer(n,j,k-4,1) - 4./3.*buffer(n,j,k-3,1) + 3.*buffer(n,j,k-2,1) - 4.*buffer(n,j,k-1,1) + 25./12.*buffer(n,j,k,1)
                TSStat%dvdk(nn,k) = 0.25*buffer(n,j,k-4,2) - 4./3.*buffer(n,j,k-3,2) + 3.*buffer(n,j,k-2,2) - 4.*buffer(n,j,k-1,2) + 25./12.*buffer(n,j,k,2)
                TSStat%dwdk(nn,k) = 0.25*buffer(n,j,k-4,3) - 4./3.*buffer(n,j,k-3,3) + 3.*buffer(n,j,k-2,3) - 4.*buffer(n,j,k-1,3) + 25./12.*buffer(n,j,k,3)
                TSStat%dpdk(nn,k) = 0.25*buffer(n,j,k-4,4) - 4./3.*buffer(n,j,k-3,4) + 3.*buffer(n,j,k-2,4) - 4.*buffer(n,j,k-1,4) + 25./12.*buffer(n,j,k,4)
                TSStat%dtdk(nn,k) = 0.25*buffer(n,j,k-4,5) - 4./3.*buffer(n,j,k-3,5) + 3.*buffer(n,j,k-2,5) - 4.*buffer(n,j,k-1,5) + 25./12.*buffer(n,j,k,5)
              endif
              TSStat%dudx(nn,k) = TSStat%dudi(nn,k)*didx + TSStat%dudj(nn,k)*djdx + TSStat%dudk(nn,k)*dkdx
              TSStat%dvdx(nn,k) = TSStat%dvdi(nn,k)*didx + TSStat%dvdj(nn,k)*djdx + TSStat%dvdk(nn,k)*dkdx
              TSStat%dwdx(nn,k) = TSStat%dwdi(nn,k)*didx + TSStat%dwdj(nn,k)*djdx + TSStat%dwdk(nn,k)*dkdx
              TSStat%dpdx(nn,k) = TSStat%dpdi(nn,k)*didx + TSStat%dpdj(nn,k)*djdx + TSStat%dpdk(nn,k)*dkdx

              TSStat%dudy(nn,k) = TSStat%dudi(nn,k)*didy + TSStat%dudj(nn,k)*djdy + TSStat%dudk(nn,k)*dkdy
              TSStat%dvdy(nn,k) = TSStat%dvdi(nn,k)*didy + TSStat%dvdj(nn,k)*djdy + TSStat%dvdk(nn,k)*dkdy
              TSStat%dwdy(nn,k) = TSStat%dwdi(nn,k)*didy + TSStat%dwdj(nn,k)*djdy + TSStat%dwdk(nn,k)*dkdy

              TSStat%dudz(nn,k) = TSStat%dudi(nn,k)*didz + TSStat%dudj(nn,k)*djdz + TSStat%dudk(nn,k)*dkdz
              TSStat%dvdz(nn,k) = TSStat%dvdi(nn,k)*didz + TSStat%dvdj(nn,k)*djdz + TSStat%dvdk(nn,k)*dkdz
              TSStat%dwdz(nn,k) = TSStat%dwdi(nn,k)*didz + TSStat%dwdj(nn,k)*djdz + TSStat%dwdk(nn,k)*dkdz

              TSStat%dpdt(nn,k) = ( buffer(n-2,j,k,4) - 8.0*(buffer(n-1,j,k,4) - buffer(n+1,j,k,4)) - buffer(n+2,j,k,4) )*c12i*dti
              TSStat%dudt(nn,k) = ( buffer(n-2,j,k,1) - 8.0*(buffer(n-1,j,k,1) - buffer(n+1,j,k,1)) - buffer(n+2,j,k,1) )*c12i*dti
              TSStat%dvdt(nn,k) = ( buffer(n-2,j,k,2) - 8.0*(buffer(n-1,j,k,2) - buffer(n+1,j,k,2)) - buffer(n+2,j,k,2) )*c12i*dti
              TSStat%dwdt(nn,k) = ( buffer(n-2,j,k,3) - 8.0*(buffer(n-1,j,k,3) - buffer(n+1,j,k,3)) - buffer(n+2,j,k,3) )*c12i*dti
              TSStat%dtdt(nn,k) = ( buffer(n-2,j,k,5) - 8.0*(buffer(n-1,j,k,5) - buffer(n+1,j,k,5)) - buffer(n+2,j,k,5) )*c12i*dti
            enddo ! end n loop
          enddo ! end j loop
        enddo ! end k loop

      end subroutine CalTSDeriv_iplane

      subroutine CalTSStat_ipStat(nz,buffer)
      integer, intent(in) :: nz
      real(8), intent(out) :: buffer(nz,nv2d_ipStat)
      integer :: i, k
      real(8) :: num_dble, rho, div, omx, omy, omz, magom, s, p0, t0
      real(8) :: div_tmp

      num_dble   = dble(TSStat%nsample)
      do k=1, nz
        buffer(k,1)  = sum(TSStat%u(:,k))/num_dble  ! uave
        buffer(k,2)  = sum(TSStat%v(:,k))/num_dble  ! vave
        buffer(k,3)  = sum(TSStat%w(:,k))/num_dble  ! wave
        buffer(k,4)  = sum(TSStat%p(:,k))/num_dble  ! pave
        buffer(k,5)  = sum(TSStat%t(:,k))/num_dble  ! tave

        buffer(k,7)  = sum(TSStat%dpdx(:,k))/num_dble ! dpdxave
        buffer(k,8)  = sum(TSStat%dpdt(:,k))/num_dble ! dpdtave
        buffer(k,9)  = sum(TSStat%dpdx(:,k)*TSStat%dpdt(:,k))/num_dble ! dpdxdpdtave

        buffer(k,18) = sum(TSStat%u(:,k)**2)/num_dble   ! u2ave
        buffer(k,19) = sum(TSStat%v(:,k)**2)/num_dble   ! v2ave
        buffer(k,20) = sum(TSStat%w(:,k)**2)/num_dble   ! w2ave
        buffer(k,21) = sum(TSStat%p(:,k)**2)/num_dble   ! p2ave
        buffer(k,22) = sum(TSStat%t(:,k)**2)/num_dble   ! t2ave

        buffer(k,24) = sum(TSStat%dpdx(:,k)**2)/num_dble ! dpdx2ave
        buffer(k,25) = sum(TSStat%dpdt(:,k)**2)/num_dble ! dpdt2ave

        buffer(k,34) = sum(TSStat%u(:,k)**3)/num_dble    ! u3ave
        buffer(k,35) = sum(TSStat%u(:,k)**4)/num_dble    ! u4ave
        buffer(k,36) = sum(TSStat%p(:,k)**3)/num_dble    ! p3ave
        buffer(k,37) = sum(TSStat%p(:,k)**4)/num_dble    ! p4ave

        buffer(k,38) = sum(TSStat%u(:,k)*TSStat%v(:,k))/num_dble ! uvave
        buffer(k,39) = sum(TSStat%u(:,k)*TSStat%w(:,k))/num_dble ! uwave
        buffer(k,40) = sum(TSStat%v(:,k)*TSStat%w(:,k))/num_dble ! vwave
        buffer(k,41) = sum(TSStat%u(:,k)*TSStat%t(:,k))/num_dble ! utave
        buffer(k,42) = sum(TSStat%v(:,k)*TSStat%t(:,k))/num_dble ! vtave
        buffer(k,43) = sum(TSStat%w(:,k)*TSStat%t(:,k))/num_dble ! wtave
        buffer(k,44) = sum(TSStat%u(:,k)*TSStat%p(:,k))/num_dble ! upave
        buffer(k,45) = sum(TSStat%v(:,k)*TSStat%p(:,k))/num_dble ! vpave
        buffer(k,46) = sum(TSStat%w(:,k)*TSStat%p(:,k))/num_dble ! wpave
        buffer(k,47) = sum(TSStat%t(:,k)*TSStat%p(:,k))/num_dble ! tpave

        buffer(k,67) = sum(TSStat%dudi(:,k))/num_dble            ! dudiave
        buffer(k,68) = sum(TSStat%dvdi(:,k))/num_dble            ! dvdiave
        buffer(k,69) = sum(TSStat%dwdi(:,k))/num_dble            ! dwdiave
        buffer(k,70) = sum(TSStat%dpdi(:,k))/num_dble            ! dpdiave
        buffer(k,71) = sum(TSStat%dtdi(:,k))/num_dble            ! dtdiave

        buffer(k,72) = sum(TSStat%dudj(:,k))/num_dble            ! dudjave
        buffer(k,73) = sum(TSStat%dvdj(:,k))/num_dble            ! dvdjave
        buffer(k,74) = sum(TSStat%dwdj(:,k))/num_dble            ! dwdjave
        buffer(k,75) = sum(TSStat%dpdj(:,k))/num_dble            ! dpdjave
        buffer(k,76) = sum(TSStat%dtdj(:,k))/num_dble            ! dtdjave

        buffer(k,77) = sum(TSStat%dudk(:,k))/num_dble            ! dudkave
        buffer(k,78) = sum(TSStat%dvdk(:,k))/num_dble            ! dvdkave
        buffer(k,79) = sum(TSStat%dwdk(:,k))/num_dble            ! dwdkave
        buffer(k,80) = sum(TSStat%dpdk(:,k))/num_dble            ! dpdkave
        buffer(k,81) = sum(TSStat%dtdk(:,k))/num_dble            ! dtdkave

        buffer(k,82) = sum(TSStat%dudt(:,k))/num_dble            ! dudtave
        buffer(k,83) = sum(TSStat%dvdt(:,k))/num_dble            ! dvdtave
        buffer(k,84) = sum(TSStat%dwdt(:,k))/num_dble            ! dwdtave
        buffer(k,85) = sum(TSStat%dtdt(:,k))/num_dble            ! dtdtave

        buffer(k,86) = sum(TSStat%dudx(:,k)*TSStat%dudt(:,k))/num_dble ! dudxdudtave
        buffer(k,87) = sum(TSStat%dudt(:,k)**2)/num_dble         ! dudt2
        buffer(k,88) = sum(TSStat%dudx(:,k)**2)/num_dble         ! dudx2
        buffer(k,89) = sum(TSStat%dvdy(:,k)**2)/num_dble         ! dvdy2
        buffer(k,90) = sum(TSStat%dwdz(:,k)**2)/num_dble         ! dwdz2
        buffer(k,91) = sum(TSStat%dudx(:,k)*TSStat%dvdy(:,k))/num_dble  ! dudxdvdy
        buffer(k,92) = sum(TSStat%dudx(:,k)*TSStat%dwdz(:,k))/num_dble  ! dudxdwdz
        buffer(k,93) = sum(TSStat%dvdy(:,k)*TSStat%dwdz(:,k))/num_dble  ! dvdydwdz


        buffer(k,94) = sum(TSStat%u(:,k)/sqrt(gamma*rbar*TSStat%t(:,k)))/num_dble      ! Mave
        buffer(k,95) = sum((TSStat%u(:,k)/sqrt(gamma*rbar*TSStat%t(:,k)))**2)/num_dble ! M2ave

        buffer(k,98) = sum(TSStat%t(:,k)**3)/num_dble   ! t3ave
        buffer(k,99) = sum(TSStat%t(:,k)**4)/num_dble   ! t4ave

        buffer(k,6) = 0.d0; buffer(k,23) = 0.d0
        buffer(k,10:17) = 0.d0
        buffer(k,26:33) = 0.d0
        buffer(k,48:65) = 0.d0

        buffer(k,96:97) = 0.d0

        do i=1, TSStat%nsample
          rho = TSStat%p(i,k)/rbar/TSStat%t(i,k)
          buffer(k,6)  = buffer(k,6)  + rho      ! rhoave
          buffer(k,23) = buffer(k,23) + rho**2   ! rho2ave
          buffer(k,96) = buffer(k,96) + rho**3   ! rho3ave
          buffer(k,97) = buffer(k,97) + rho**4   ! rho4ave

          div = TSStat%dudx(i,k) + TSStat%dvdy(i,k) + TSStat%dwdz(i,k)
          omx = TSStat%dwdy(i,k) - TSStat%dvdz(i,k)
          omy = TSStat%dudz(i,k) - TSStat%dwdx(i,k)
          omz = TSStat%dvdx(i,k) - TSStat%dudy(i,k)
          magom = sqrt(omx**2 + omy**2 + omz**2)

!          div_tmp =
!          buffer(k.91) = buffer(k,91) +

          buffer(k,10) = buffer(k,10) + div      ! divave
          buffer(k,11) = buffer(k,11) + omx      ! omxave
          buffer(k,12) = buffer(k,12) + omy      ! omyave
          buffer(k,13) = buffer(k,13) + omz      ! omzave
          buffer(k,14) = buffer(k,14) + magom    ! magomave

          buffer(k,26) = buffer(k,26) + div**2   ! div2ave
          buffer(k,27) = buffer(k,27) + omx**2   ! omxave
          buffer(k,28) = buffer(k,28) + omy**2   ! omyave
          buffer(k,29) = buffer(k,29) + omz**2   ! omzave
          buffer(k,30) = buffer(k,30) + magom**2 ! magomave

          s = cp*log( TSStat%t(i,k)/buffer(k,5) ) - rbar*log( TSStat%p(i,k)/buffer(k,4) )
          t0 = TSStat%t(i,k) + 0.5/cp*( TSStat%u(i,k)**2 + TSStat%v(i,k)**2 + TSStat%w(i,k)**2 )
          p0 = TSStat%p(i,k)*(t0/TSStat%t(i,k))**(gamma/(gamma-1.))
          buffer(k,15) = buffer(k,15) + s        ! save
          buffer(k,16) = buffer(k,16) + t0       ! t0ave
          buffer(k,17) = buffer(k,17) + p0       ! p0ave
          buffer(k,31) = buffer(k,31) + s**2     ! s2ave
          buffer(k,32) = buffer(k,32) + t0**2    ! t02ave
          buffer(k,33) = buffer(k,33) + p0**2    ! p02ave

          buffer(k,48) = buffer(k,48) + rho*TSStat%p(i,k)         ! prave
          buffer(k,49) = buffer(k,49) + s*TSStat%p(i,k)           ! psave
          buffer(k,50) = buffer(k,50) + s*magom                   ! somave
          buffer(k,51) = buffer(k,51) + TSStat%p(i,k)*magom       ! pomave
          buffer(k,52) = buffer(k,52) + rho*TSStat%u(i,k)*t0      ! ruT0ave
          buffer(k,53) = buffer(k,53) + rho*TSStat%u(i,k)         ! ruave
          buffer(k,54) = buffer(k,54) + rho*TSStat%v(i,k)         ! rvave
          buffer(k,55) = buffer(k,55) + rho*TSStat%w(i,k)         ! rwave
          buffer(k,56) = buffer(k,56) + rho*TSStat%u(i,k)**2      ! ru2ave
          buffer(k,57) = buffer(k,57) + rho*TSStat%v(i,k)**2      ! rv2ave
          buffer(k,58) = buffer(k,58) + rho*TSStat%w(i,k)**2      ! rw2ave
          buffer(k,59) = buffer(k,59) + rho*TSStat%u(i,k)*TSStat%v(i,k)   ! ruvave
          buffer(k,60) = buffer(k,60) + rho*TSStat%u(i,k)*TSStat%w(i,k)   ! ruwave
          buffer(k,61) = buffer(k,61) + rho*TSStat%v(i,k)*TSStat%w(i,k)   ! rvwave
          buffer(k,62) = buffer(k,62) + (rho*TSStat%u(i,k))**2            ! r2u2ave
          buffer(k,63) = buffer(k,63) + rho*TSStat%u(i,k)*TSStat%t(i,k)   ! rutave
          buffer(k,64) = buffer(k,64) + rho*TSStat%v(i,k)*TSStat%t(i,k)   ! rvtave
          buffer(k,65) = buffer(k,65) + rho*TSStat%w(i,k)*TSStat%t(i,k)   ! rwtave

        enddo
        buffer(k,6)  = buffer(k,6)/num_dble
        buffer(k,23) = buffer(k,23)/num_dble
        buffer(k,96) = buffer(k,96)/num_dble  ! rho3ave
        buffer(k,97) = buffer(k,97)/num_dble  ! rho4ave
        buffer(k,10:17) = buffer(k,10:17)/num_dble
        buffer(k,26:33) = buffer(k,26:33)/num_dble
        buffer(k,48:65) = buffer(k,48:65)/num_dble
        buffer(k,66) = 1.e30 !!!!
      enddo ! end k loop


      end subroutine CalTSStat_ipStat

      subroutine CalTSStat_ipInt(nz,zloc,buffer,buffer_grid)
        integer, intent(in) :: nz
        real(8), intent(in) :: zloc(nz), buffer_grid(nz,9)
        real(8), intent(out) :: buffer(nz,nv2d_ipInt)
        integer :: i, k
        real(8) :: num_dble, rho, tauwtmp
        integer :: temp(1)

        !!!!!!
        dudn = -25./12.*TSStat%u(:,1) + 4.*TSStat%u(:,2) - 3.*TSStat%u(:,3) + 4./3.*TSStat%u(:,4) &
               - 0.25*TSStat%u(:,5)
        dtdn = -25./12.*TSStat%t(:,1) + 4.*TSStat%t(:,2) - 3.*TSStat%t(:,3) + 4./3.*TSStat%t(:,4) &
               - 0.25*TSStat%t(:,5)
!        do i=1, TSStat%nsample
!          do k=1, nz
!            mu(i,k) = CalMu(TSStat%t(i,k))
!          enddo
!        enddo

        num_dble = dble(TSStat%nsample)
        do k=1, nz
          buffer(k,1)  = sum(TSStat%u(:,k),dim=1)/num_dble  ! uave
          buffer(k,2)  = sum(TSStat%v(:,k),dim=1)/num_dble  ! vave
          buffer(k,3)  = sum(TSStat%w(:,k),dim=1)/num_dble  ! wave
          buffer(k,4)  = sum(TSStat%p(:,k),dim=1)/num_dble  ! pave
          buffer(k,5)  = sum(TSStat%t(:,k),dim=1)/num_dble  ! tave

          buffer(k,7)  = sum(TSStat%u(:,k)**2,dim=1)/num_dble ! u2
          buffer(k,8)  = sum(TSStat%v(:,k)**2,dim=1)/num_dble ! v2
          buffer(k,9)  = sum(TSStat%w(:,k)**2,dim=1)/num_dble ! w2
          buffer(k,10) = sum(TSStat%p(:,k)**2,dim=1)/num_dble ! p2
          buffer(k,11) = sum(TSStat%t(:,k)**2,dim=1)/num_dble ! t2

          buffer(k,6)  = 0.d0; buffer(k,12) = 0.d0
          buffer(k,13) = 0.d0
          do i=1, TSStat%nsample
            rho = TSStat%p(i,k)/rbar/TSStat%t(i,k)
            buffer(k,6)  = buffer(k,6) + rho
            buffer(k,12) = buffer(k,12) + rho**2
            buffer(k,13) = buffer(k,13) + mu(i,k)*dudn(i)*buffer_grid(k,9)

          enddo ! end i loop
        enddo ! end k loop
        buffer(:,6)  = buffer(:,6)/num_dble   ! rhoave
        buffer(:,12) = buffer(:,12)/num_dble  ! rho2
        buffer(:,13) = buffer(:,13)/num_dble
        buffer(:,13) = buffer(1,13)           ! tauw

      end subroutine CalTSStat_ipInt


    subroutine WriteTSStat_ipStat(fn,iloc,nz,buffer)
      character(*), intent(in) :: fn
      integer, intent(in) :: nz, iloc
      real(8), intent(in) :: buffer(1,nz,nv2d_ipStat+1)
      integer :: i,n
      character(4) :: fnum
      character(500) :: filename

      write(unit=fnum,fmt='(I04.4)') iloc

      TSiout%fname = trim(fn)
      TSiout%gname = '/'//fnum//'Stat'
      TSiout%rank = 2
      TSiout%dnum = nv2d_ipStat+1

      if(allocated(TSiout%dname)) then
        deallocate(TSiout%dname,TSiout%dimsf,TSiout%dimsm)
        deallocate(TSiout%block,TSiout%count,TSiout%stride)
        deallocate(TSiout%offset)
      endif

      allocate(TSiout%dname(TSiout%dnum))
      allocate(TSiout%dimsf(TSiout%rank))
      allocate(TSiout%dimsm(TSiout%rank))
      allocate(TSiout%block(TSiout%rank))
      allocate(TSiout%count(TSiout%rank))
      allocate(TSiout%stride(TSiout%rank))
      allocate(TSiout%offset(TSiout%rank))
      TSiout%dname = dname2d_ipStat
      TSiout%dimsf(1) = 1
      TSiout%dimsf(2) = nz
      TSiout%dimsm  = TSiout%dimsf
      TSiout%block  = TSiout%dimsm
      TSiout%count  = 1
      TSiout%stride = 1
      TSiout%offset = 0
      TSiout%IsHSInitialized = .true.
      TSiout%IsMultiGroup = .true.
      call WriteTSHDF5_2D(TSiout, buffer)

      ! write tec360 ascii data file
      filename = 'AveAcoustic_iplane_'//fnum//'-timeave.dat'
      print*,'Writing file: ',trim(filename)
      open(33,file=filename,status='unknown')
      write(33,*) 'Variables = ',(dname2d_ipStat(n),n=1,TSiout%dnum)
      write(33,*) 'Zone I = ',size(buffer,dim=2)
      do i=1,size(buffer,dim=2)
          write(33,*) (buffer(1,i,n),n=1,TSiout%dnum)
      enddo
      close(33)

    end subroutine WriteTSStat_ipStat

    subroutine WriteTSStat_ipInt(fn,iloc,nz,buffer)
      character(*), intent(in) :: fn
      integer, intent(in) :: nz, iloc
      real(8), intent(in) :: buffer(1,nz,nv2d_ipInt+1)
      character(4) :: fnum

      write(unit=fnum,fmt='(I04.4)') iloc

      if(trim(TSiout%fname).eq.'none') then
        TSiout%IsMultiGroup = .false.
      elseif(trim(TSiout%fname).eq.trim(fn)) then
        TSiout%IsMultiGroup = .true.
      else
        TSiout%IsMultiGroup = .false.
      endif

      TSiout%fname = trim(fn)
      TSiout%gname = '/'//fnum//'Int'
      TSiout%rank = 2
      TSiout%dnum = nv2d_ipInt+1

      if(allocated(TSiout%dname)) then
        deallocate(TSiout%dname,TSiout%dimsf,TSiout%dimsm)
        deallocate(TSiout%block,TSiout%count,TSiout%stride)
        deallocate(TSiout%offset)
      endif

      allocate(TSiout%dname(TSiout%dnum))
      allocate(TSiout%dimsf(TSiout%rank))
      allocate(TSiout%dimsm(TSiout%rank))
      allocate(TSiout%block(TSiout%rank))
      allocate(TSiout%count(TSiout%rank))
      allocate(TSiout%stride(TSiout%rank))
      allocate(TSiout%offset(TSiout%rank))
      TSiout%dname = dname2d_ipInt
      TSiout%dimsf(1) = 1
      TSiout%dimsf(2) = nz
      TSiout%dimsm  = TSiout%dimsf
      TSiout%block  = TSiout%dimsm
      TSiout%count  = 1
      TSiout%stride = 1
      TSiout%offset = 0
      TSiout%IsHSInitialized = .true.
!      if( .not.TSiout%IsMultiGroup) then
        call WriteTSHDF5_2D(TSiout, buffer)
!        TSiout%IsMultiGroup = .true.
!      else
!        call WriteTSHDF5_2D(TSiout, buffer)
!      endif

    end subroutine WriteTSStat_ipInt

    elemental real(8) function CalMu(tin)
      real(8), intent(in) :: tin

      CalMu = (1.458E-6)*((tin)**1.5)/(tin+110.4)

    end function CalMu




    end module MTSStat_iplane
