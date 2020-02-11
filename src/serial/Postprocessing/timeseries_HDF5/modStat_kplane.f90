 module MTSStat_kplane
    use MTSHDF5
    implicit none
    integer, parameter :: nv2d_kpStat = 85+5+5, nv2d_kpInt = 5+2+4
    real(8), parameter :: rbar = 287.d0, gamma = 1.4, cp = 1004.5, Prandtl = 0.71 ! (cp = 3.5*rbar)
    type tp_stat_kplane
        integer :: nsample
        real(8), dimension(:), allocatable :: u, v, w, p, t
        real(8), dimension(:), allocatable :: dudx, dvdx, dwdx, dpdx, dtdx
        real(8), dimension(:), allocatable :: dudy, dvdy, dwdy, dpdy, dtdy
        real(8), dimension(:), allocatable :: dudz, dvdz, dwdz, dpdz, dtdz
        real(8), dimension(:), allocatable :: dudt, dvdt, dwdt, dpdt, dtdt

        real(8), dimension(:), allocatable :: dudi, dvdi, dwdi, dpdi, dtdi
        real(8), dimension(:), allocatable :: dudj, dvdj, dwdj, dpdj, dtdj
        real(8), dimension(:), allocatable :: dudk, dvdk, dwdk, dpdk, dtdk

    end type tp_stat_kplane
    type(tp_stat_kplane), private :: TSStat
    character(30), private :: dname2d_kpStat(nv2d_kpStat+2), dname2d_kpInt(nv2d_kpInt+1)
    parameter( dname2d_kpStat = (/'uave','vave','wave','pave','tave', 'rhoave', 'dpdx','dpdt','dpdxdpdt','divave','omxave','omyave','omzave','magomave','save',&
                           't0ave','p0ave','u2','v2','w2','p2','t2','rho2','dpdx2', 'dpdt2', 'div2', 'omx2', 'omy2', 'omz2', 'magom2', 's2', 't02', 'p02',&
                           'u3', 'u4', 'p3', 'p4', 'uv', 'uw', 'vw', 'ut', 'vt', 'wt', 'up', 'vp', 'wp', 'pt', 'pr', 'ps', 'som', 'pom', 'ruT0', 'ru', &
                           'rv', 'rw', 'ru2', 'rv2', 'rw2', 'ruv', 'ruw', 'rvw', 'r2u2', 'rut', 'rvt', 'rwt',  'uvd', 'dudi', 'dvdi', 'dwdi', 'dpdi', 'dtdi', &
                           'dudj', 'dvdj', 'dwdj', 'dpdj', 'dtdj', 'dudk', 'dvdk', 'dwdk', 'dpdk', 'dtdk', 'dudt', 'dvdt', 'dwdt', 'dtdt', 'dudxdudt', 'dudt2', &
                           'dudx2', 'dvdy2', 'dwdz2', 'dudxdvdy', 'dudxdwdz', 'dvdydwdz', 'M', 'M2', 'x', 'z'/) )

    parameter(dname2d_kpInt = (/'tauw', 'tauw2', 'Cf','utau','ztau', 'qw', 'qw2', 'tauw3', 'tauw4', 'qw3', 'qw4', 'x'/))

    real(8), dimension(:), allocatable, private :: dudn, dtdn
    real(8), dimension(:), allocatable, private :: mu

    type(tp_hyperslab), private :: TSkout

 contains
   subroutine InitTSStat(ns)
     integer, intent(in) :: ns

     TSkout%fname = 'none'

     TSStat%nsample = ns
     allocate(TSStat%u(ns))
     allocate(TSStat%v(ns))
     allocate(TSStat%w(ns))
     allocate(TSStat%p(ns))
     allocate(TSStat%t(ns))
     allocate(TSStat%dudx(ns), TSStat%dudy(ns), TSStat%dudz(ns))
     allocate(TSStat%dvdx(ns), TSStat%dvdy(ns), TSStat%dvdz(ns))
     allocate(TSStat%dwdx(ns), TSStat%dwdy(ns), TSStat%dwdz(ns))
     allocate(TSStat%dpdx(ns), TSStat%dpdt(ns) )
     allocate(TSStat%dudt(ns), TSStat%dvdt(ns), TSStat%dwdt(ns))
     allocate(TSStat%dtdt(ns))

     allocate(TSStat%dudi(ns), TSStat%dvdi(ns), TSStat%dwdi(ns), TSStat%dpdi(ns), TSStat%dtdi(ns))
     allocate(TSStat%dudj(ns), TSStat%dvdj(ns), TSStat%dwdj(ns), TSStat%dpdj(ns), TSStat%dtdj(ns))
     allocate(TSStat%dudk(ns), TSStat%dvdk(ns), TSStat%dwdk(ns), TSStat%dpdk(ns), TSStat%dtdk(ns))

     allocate(dudn(TSStat%nsample), dtdn(TSStat%nsample))
     allocate(mu(TSStat%nsample))

   end subroutine InitTSStat

   subroutine CalTSDeriv_kplane(nt, ny, nv, buffer, buffer_grid, dt)
     integer, intent(in) :: nt, ny, nv
     real(8), intent(in) :: buffer(nt,ny,-2:2,nv), buffer_grid(ny, 9), dt
!     type(tp_stat), intent(out) :: TSStat

     integer :: j, n, nn
     real(8) :: didx, djdx, dkdx, didy, djdy, dkdy, didz, djdz, dkdz
!     real(8) :: dudi, dvdi, dwdi, dudj, dvdj, dwdj, dudk, dvdk, dwdk
!     real(8) :: dpdi, dpdj, dpdk
     real(8), parameter :: c12i = 1.d0/12.d0
     real(8) :: dti

     dti = 1.d0/dt

     TSStat%u = reshape(buffer(3:nt-2,3:ny-2,0,1),(/TSStat%nsample/))
     TSStat%v = reshape(buffer(3:nt-2,3:ny-2,0,2),(/TSStat%nsample/))
     TSStat%w = reshape(buffer(3:nt-2,3:ny-2,0,3),(/TSStat%nsample/))
     TSStat%p = reshape(buffer(3:nt-2,3:ny-2,0,4),(/TSStat%nsample/))
     TSStat%t = reshape(buffer(3:nt-2,3:ny-2,0,5),(/TSStat%nsample/))
     dudn = reshape(buffer(3:nt-2,3:ny-2,1,6),(/TSStat%nsample/))
     dtdn = reshape(buffer(3:nt-2,3:ny-2,1,10),(/TSStat%nsample/))
     nn = 0
     do j = 3, ny-2
           ! didx, djdx, dkdx, didy, djdy, dkdy, didz, djdz, dkdz
           didx = buffer_grid(j,1)
           djdx = buffer_grid(j,2)
           dkdx = buffer_grid(j,3)
           didy = buffer_grid(j,4)
           djdy = buffer_grid(j,5)
           dkdy = buffer_grid(j,6)
           didz = buffer_grid(j,7)
           djdz = buffer_grid(j,8)
           dkdz = buffer_grid(j,9)
           do n = 3, nt-2
               nn = nn + 1
               TSStat%dudi(nn) = ( buffer(n,j,-2,1) - 8.0*(buffer(n,j,-1,1) - buffer(n,j,1,1)) - buffer(n,j,2,1) )*c12i
               TSStat%dvdi(nn) = ( buffer(n,j,-2,2) - 8.0*(buffer(n,j,-1,2) - buffer(n,j,1,2)) - buffer(n,j,2,2) )*c12i
               TSStat%dwdi(nn) = ( buffer(n,j,-2,3) - 8.0*(buffer(n,j,-1,3) - buffer(n,j,1,3)) - buffer(n,j,2,3) )*c12i
               TSStat%dpdi(nn) = ( buffer(n,j,-2,4) - 8.0*(buffer(n,j,-1,4) - buffer(n,j,1,4)) - buffer(n,j,2,4) )*c12i
               TSStat%dtdi(nn) = ( buffer(n,j,-2,5) - 8.0*(buffer(n,j,-1,5) - buffer(n,j,1,5)) - buffer(n,j,2,5) )*c12i

               TSStat%dudj(nn) = ( buffer(n,j-2,0,1) - 8.0*(buffer(n,j-1,0,1) - buffer(n,j+1,0,1)) - buffer(n,j+2,0,1) )*c12i
               TSStat%dvdj(nn) = ( buffer(n,j-2,0,2) - 8.0*(buffer(n,j-1,0,2) - buffer(n,j+1,0,2)) - buffer(n,j+2,0,2) )*c12i
               TSStat%dwdj(nn) = ( buffer(n,j-2,0,3) - 8.0*(buffer(n,j-1,0,3) - buffer(n,j+1,0,3)) - buffer(n,j+2,0,3) )*c12i
               TSStat%dpdj(nn) = ( buffer(n,j-2,0,4) - 8.0*(buffer(n,j-1,0,4) - buffer(n,j+1,0,4)) - buffer(n,j+2,0,4) )*c12i
               TSStat%dtdj(nn) = ( buffer(n,j-2,0,5) - 8.0*(buffer(n,j-1,0,5) - buffer(n,j+1,0,5)) - buffer(n,j+2,0,5) )*c12i

               TSStat%dudk(nn) = buffer(n,j,0,6)
               TSStat%dvdk(nn) = buffer(n,j,0,7)
               TSStat%dwdk(nn) = buffer(n,j,0,8)
               TSStat%dpdk(nn) = buffer(n,j,0,9)
               TSStat%dtdk(nn) = buffer(n,j,0,10)

               TSStat%dudx(nn) = TSStat%dudi(nn)*didx + TSStat%dudj(nn)*djdx + TSStat%dudk(nn)*dkdx
               TSStat%dvdx(nn) = TSStat%dvdi(nn)*didx + TSStat%dvdj(nn)*djdx + TSStat%dvdk(nn)*dkdx
               TSStat%dwdx(nn) = TSStat%dwdi(nn)*didx + TSStat%dwdj(nn)*djdx + TSStat%dwdk(nn)*dkdx
               TSStat%dpdx(nn) = TSStat%dpdi(nn)*didx + TSStat%dpdj(nn)*djdx + TSStat%dpdk(nn)*dkdx

               TSStat%dudy(nn) = TSStat%dudi(nn)*didy + TSStat%dudj(nn)*djdy + TSStat%dudk(nn)*dkdy
               TSStat%dvdy(nn) = TSStat%dvdi(nn)*didy + TSStat%dvdj(nn)*djdy + TSStat%dvdk(nn)*dkdy
               TSStat%dwdy(nn) = TSStat%dwdi(nn)*didy + TSStat%dwdj(nn)*djdy + TSStat%dwdk(nn)*dkdy
 !              TSStat%dpdy(nn) = dpdi*didy + dpdj*djdy + dpdk*dkdy

               TSStat%dudz(nn) = TSStat%dudi(nn)*didz + TSStat%dudj(nn)*djdz + TSStat%dudk(nn)*dkdz
               TSStat%dvdz(nn) = TSStat%dvdi(nn)*didz + TSStat%dvdj(nn)*djdz + TSStat%dvdk(nn)*dkdz
               TSStat%dwdz(nn) = TSStat%dwdi(nn)*didz + TSStat%dwdj(nn)*djdz + TSStat%dwdk(nn)*dkdz
!               TSStat%dpdz(nn) = dpdi*didz + dpdj*djdz + dpdk*dkdz

               TSStat%dpdt(nn) = ( buffer(n-2,j,0,4) - 8.0*(buffer(n-1,j,0,4) - buffer(n+1,j,0,4)) - buffer(n+2,j,0,4) )*c12i*dti
               TSStat%dudt(nn) = ( buffer(n-2,j,0,1) - 8.0*(buffer(n-1,j,0,1) - buffer(n+1,j,0,1)) - buffer(n+2,j,0,1) )*c12i*dti
               TSStat%dvdt(nn) = ( buffer(n-2,j,0,2) - 8.0*(buffer(n-1,j,0,2) - buffer(n+1,j,0,2)) - buffer(n+2,j,0,2) )*c12i*dti
               TSStat%dwdt(nn) = ( buffer(n-2,j,0,3) - 8.0*(buffer(n-1,j,0,3) - buffer(n+1,j,0,3)) - buffer(n+2,j,0,3) )*c12i*dti
               TSStat%dtdt(nn) = ( buffer(n-2,j,0,5) - 8.0*(buffer(n-1,j,0,5) - buffer(n+1,j,0,5)) - buffer(n+2,j,0,5) )*c12i*dti
          enddo
      enddo

   end subroutine CalTSDeriv_kplane

   subroutine CalTSStat_kpStat(buffer)
!      type(tp_stat), intent(in) :: TSStat
      real(8), intent(out) :: buffer(nv2d_kpStat)
      integer :: i
      real(8) :: num_dble, rho, div, omx, omy, omz, magom, s, p0, t0

     num_dble   = dble(TSStat%nsample)
     buffer(1)  = sum(TSStat%u)/num_dble  ! uave
     buffer(2)  = sum(TSStat%v)/num_dble  ! vave
     buffer(3)  = sum(TSStat%w)/num_dble  ! wave
     buffer(4)  = sum(TSStat%p)/num_dble  ! pave
     buffer(5)  = sum(TSStat%t)/num_dble  ! tave
!     buffer(6)  = sum(TSStat%p/rbar/TSStat%t) ! rhoave

     buffer(7)  = sum(TSStat%dpdx)/num_dble  ! dpdxave
     buffer(8)  = sum(TSStat%dpdt)/num_dble  ! dpdtave
     buffer(9)  = sum(TSStat%dpdx*TSStat%dpdt)/num_dble ! dpdxdpdtave
!     buffer(10) = sum(TSStat%dudx + TSStat%dvdy + TSStat%dwdz)/num_dble ! div

     buffer(18) = sum(TSStat%u**2)/num_dble ! u2ave
     buffer(19) = sum(TSStat%v**2)/num_dble ! v2ave
     buffer(20) = sum(TSStat%w**2)/num_dble ! w2ave
     buffer(21) = sum(TSStat%p**2)/num_dble ! p2ave
     buffer(22) = sum(TSStat%t**2)/num_dble ! t2ave

     buffer(24) = sum(TSStat%dpdx**2)/num_dble  ! dpdx2ave
     buffer(25) = sum(TSStat%dpdt**2)/num_dble  ! dpdt2ave

     buffer(34) = sum(TSStat%u**3)/num_dble   ! up3
     buffer(35) = sum(TSStat%u**4)/num_dble   ! up4
     buffer(36) = sum(TSStat%p**3)/num_dble   ! pp3
     buffer(37) = sum(TSStat%p**4)/num_dble   ! pp4

     buffer(38) = sum(TSStat%u*TSStat%v)/num_dble  ! uvave
     buffer(39) = sum(TSStat%u*TSStat%w)/num_dble  ! uwave
     buffer(40) = sum(TSStat%v*TSStat%w)/num_dble  ! vwave
     buffer(41) = sum(TSStat%u*TSStat%t)/num_dble  ! utave
     buffer(42) = sum(TSStat%v*TSStat%t)/num_dble  ! vtave
     buffer(43) = sum(TSStat%w*TSStat%t)/num_dble  ! wtave
     buffer(44) = sum(TSStat%u*TSStat%p)/num_dble  ! upave
     buffer(45) = sum(TSStat%v*TSStat%p)/num_dble  ! vpave
     buffer(46) = sum(TSStat%w*TSStat%p)/num_dble  ! wpave
     buffer(47) = sum(TSStat%t*TSStat%p)/num_dble  ! tpave

     buffer(67) = sum(TSStat%dudi)/num_dble            ! dudiave
     buffer(68) = sum(TSStat%dvdi)/num_dble            ! dvdiave
     buffer(69) = sum(TSStat%dwdi)/num_dble            ! dwdiave
     buffer(70) = sum(TSStat%dpdi)/num_dble            ! dpdiave
     buffer(71) = sum(TSStat%dtdi)/num_dble            ! dtdiave

     buffer(72) = sum(TSStat%dudj)/num_dble            ! dudjave
     buffer(73) = sum(TSStat%dvdj)/num_dble            ! dvdjave
     buffer(74) = sum(TSStat%dwdj)/num_dble            ! dwdjave
     buffer(75) = sum(TSStat%dpdj)/num_dble            ! dpdjave
     buffer(76) = sum(TSStat%dtdj)/num_dble            ! dtdjave

     buffer(77) = sum(TSStat%dudk)/num_dble            ! dudkave
     buffer(78) = sum(TSStat%dvdk)/num_dble            ! dvdkave
     buffer(79) = sum(TSStat%dwdk)/num_dble            ! dwdkave
     buffer(80) = sum(TSStat%dpdk)/num_dble            ! dpdkave
     buffer(81) = sum(TSStat%dtdk)/num_dble            ! dtdkave

     buffer(82)  = sum(TSStat%dudt)/num_dble           ! dudtave
     buffer(83)  = sum(TSStat%dvdt)/num_dble           ! dvdtave
     buffer(84)  = sum(TSStat%dwdt)/num_dble           ! dwdtave
     buffer(85)  = sum(TSStat%dtdt)/num_dble           ! dtdtave

     buffer(86) = sum(TSStat%dudx*TSStat%dudt)/num_dble ! dudxdudtave
     buffer(87) = sum(TSStat%dudt**2)/num_dble         ! dudt2
     buffer(88) = sum(TSStat%dudx**2)/num_dble         ! dudx2
     buffer(89) = sum(TSStat%dvdy**2)/num_dble         ! dvdy2
     buffer(90) = sum(TSStat%dwdz**2)/num_dble         ! dwdz2
     buffer(91) = sum(TSStat%dudx*TSStat%dvdy)/num_dble         ! dudxdvdy
     buffer(92) = sum(TSStat%dudx*TSStat%dwdz)/num_dble         ! dudxdwdz
     buffer(93) = sum(TSStat%dvdy*TSStat%dwdz)/num_dble         ! dvdydwdz

     buffer(94) = sum(TSStat%u/sqrt(gamma*rbar*TSStat%t))/num_dble      ! Mave
     buffer(95) = sum((TSStat%u/sqrt(gamma*rbar*TSStat%t))**2)/num_dble ! M2ave

     buffer(6)  = 0.d0; buffer(23) = 0.d0
     buffer(10:17) = 0.d0
     buffer(26:33) = 0.d0
     buffer(48:65) = 0.d0

     do i = 1, TSStat%nsample

         rho = TSStat%p(i)/rbar/TSStat%t(i)
         buffer(6)  = buffer(6)  + rho      ! rhoave
         buffer(23) = buffer(23) + rho**2   ! rho2ave

         div = TSStat%dudx(i) + TSStat%dvdy(i) + TSStat%dwdz(i)
         omx = TSStat%dwdy(i) - TSStat%dvdz(i)
         omy = TSStat%dudz(i) - TSStat%dwdx(i)
         omz = TSStat%dvdx(i) - TSStat%dudy(i)
         magom = sqrt(omx**2 + omy**2 + omz**2)

         buffer(10) = buffer(10) + div       ! divave
         buffer(11) = buffer(11) + omx       ! omxave
         buffer(12) = buffer(12) + omy       ! omyave
         buffer(13) = buffer(13) + omz       ! omzave
         buffer(14) = buffer(14) + magom  ! magomave

         buffer(26) = buffer(26) + div**2    ! div2ave
         buffer(27) = buffer(27) + omx**2    ! omx2ave
         buffer(28) = buffer(28) + omy**2    ! omy2ave
         buffer(29) = buffer(29) + omz**2    ! omz2ave
         buffer(30) = buffer(30) + magom**2 ! magom2ave

         s =  cp*log( TSStat%t(i)/buffer(5) ) - rbar*log( TSStat%p(i)/buffer(4) )
         t0 = TSStat%t(i) + 0.5/cp*(TSStat%u(i)**2 + TSStat%v(i)**2 + TSStat%w(i)**2)
         p0 = TSStat%p(i)*(t0/TSStat%t(i))**(gamma/(gamma-1.))
         buffer(15) = buffer(15) + s      ! save
         buffer(16) = buffer(16) + t0     ! t0ave
         buffer(17) = buffer(17) + p0     ! p0ave
         buffer(31) = buffer(31) + s**2   ! s2ave
         buffer(32) = buffer(32) + t0**2  ! t02ave
         buffer(33) = buffer(33) + p0**2  ! p02ave

         buffer(48) = buffer(48) + rho*TSStat%p(i)  ! prave
         buffer(49) = buffer(49) + s*TSStat%p(i)    ! psave
         buffer(50) = buffer(50) + s*magom            ! somave
         buffer(51) = buffer(51) + TSStat%p(i)*magom  ! pomave
         buffer(52) = buffer(52) + rho*TSStat%u(i)*t0 ! ruT0ave
         buffer(53) = buffer(53) + rho*TSStat%u(i)    ! ruave
         buffer(54) = buffer(54) + rho*TSStat%v(i)    ! rvave
         buffer(55) = buffer(55) + rho*TSStat%w(i)    ! rwave
         buffer(56) = buffer(56) + rho*TSStat%u(i)**2 ! ru2ave
         buffer(57) = buffer(57) + rho*TSStat%v(i)**2 ! rv2ave
         buffer(58) = buffer(58) + rho*TSStat%w(i)**2 ! rw2ave
         buffer(59) = buffer(59) + rho*TSStat%u(i)*TSStat%v(i) ! ruvave
         buffer(60) = buffer(60) + rho*TSStat%u(i)*TSStat%w(i) ! ruwave
         buffer(61) = buffer(61) + rho*TSStat%v(i)*TSStat%w(i) ! rvwave
         buffer(62) = buffer(62) + (rho*TSStat%u(i))**2 ! r2u2ave
         buffer(63) = buffer(63) + rho*TSStat%u(i)*TSStat%t(i)    ! rutve
         buffer(64) = buffer(64) + rho*TSStat%v(i)*TSStat%t(i)    ! rvtve
         buffer(65) = buffer(65) + rho*TSStat%w(i)*TSStat%t(i)    ! rwtve

     enddo
     buffer(6) = buffer(6)/num_dble
     buffer(23) = buffer(23)/num_dble
     buffer(10:17) = buffer(10:17)/num_dble
     buffer(26:33) = buffer(26:33)/num_dble
     buffer(48:65) = buffer(48:65)/num_dble
     buffer(66) = 1.e30

   end subroutine CalTSStat_kpStat

      subroutine CalTSStat_kpInt(buffer, buffer_grid)
        real(8), intent(out) :: buffer(nv2d_kpInt)
        real(8), intent(in) :: buffer_grid(9)
        real(8) :: num_dble, rho, buffer_rho, buffer_t, muw
        integer :: i

        do i=1, TSStat%nsample
          mu(i) = CalMu(TSStat%t(i))
        enddo

        num_dble = dble(TSStat%nsample)

        buffer = 0.d0
        !buffer(2) = 0.d0
        buffer_rho = 0.d0
        !buffer(6) = 0.d0

        do i=1, TSStat%nsample
          rho = TSStat%p(i)/rbar/TSStat%t(i)
          buffer_rho = buffer_rho + rho
          buffer(1) = buffer(1) + mu(i)*dudn(i)*buffer_grid(9)        ! tauw
          buffer(2) = buffer(2) + (mu(i)*dudn(i)*buffer_grid(9))**2   ! tauw2

          buffer(8) = buffer(8) + (mu(i)*dudn(i)*buffer_grid(9))**3   ! tauw3
          buffer(9) = buffer(9) + (mu(i)*dudn(i)*buffer_grid(9))**4   ! tauw4

          buffer(6) = buffer(6) + (mu(i)*Cp/Prandtl)*dtdn(i)*buffer_grid(9) ! qw
          buffer(7) = buffer(7) + ((mu(i)*Cp/Prandtl)*dtdn(i)*buffer_grid(9))**2 ! qw2

          buffer(10) = buffer(10) + ((mu(i)*Cp/Prandtl)*dtdn(i)*buffer_grid(9))**3 ! qw3
          buffer(11) = buffer(11) + ((mu(i)*Cp/Prandtl)*dtdn(i)*buffer_grid(9))**4 ! qw4

        enddo
        buffer_rho = buffer_rho/num_dble ! rhoave
        buffer(1) = buffer(1)/num_dble   ! tauw
        buffer(2) = buffer(2)/num_dble   ! tauw2
        buffer(3) = buffer(1)/(0.5*rhoinf*uinf*uinf) ! Cf
        buffer(4) = sqrt(abs(buffer(1)/buffer_rho))  ! utau
        buffer_t  = sum(TSStat%t)/num_dble   ! tave
        muw = CalMu(buffer_t)
        buffer(5) = muw/buffer_rho/buffer(3) ! ztau

        buffer(6) = buffer(6)/num_dble ! qw
        buffer(7) = buffer(7)/num_dble ! qw2

        buffer(8)  = buffer(8)/num_dble  ! tauw3
        buffer(9)  = buffer(9)/num_dble  ! tauw4
        buffer(10) = buffer(10)/num_dble ! qw3
        buffer(11) = buffer(11)/num_dble ! qw4

      end subroutine CalTSStat_kpInt


   subroutine WriteTSStat_kpStat(fn,kloc,nx,buffer)
      character(*), intent(in) :: fn
      integer, intent(in) :: kloc,nx
      real(8), intent(in) :: buffer(nx,nv2d_kpStat+2)
      character(4) :: fnum

       write(unit=fnum,fmt='(I04.4)') kloc

       if(trim(TSkout%fname).eq.'none') then
         TSkout%IsMultiGroup = .false.
       elseif(trim(TSkout%fname).eq.trim(fn)) then
         TSkout%IsMultiGroup = .true.
       else
         TSkout%IsMultiGroup = .false.
       endif

       TSkout%fname = trim(fn)
       TSkout%gname = "/"//fnum//'Stat'
       TSkout%rank = 1
       TSkout%dnum = nv2d_kpStat+2

       if(allocated(TSkout%dname)) then
         deallocate(TSkout%dname,TSkout%dimsf,TSkout%dimsm)
         deallocate(TSkout%block,TSkout%count,TSkout%stride)
         deallocate(TSkout%offset)
       endif

       allocate(TSkout%dname(TSkout%dnum))
       allocate(TSkout%dimsf(TSkout%rank))
       allocate(TSkout%dimsm(TSkout%rank))
       allocate(TSkout%block(TSkout%rank))
       allocate(TSkout%count(TSkout%rank))
       allocate(TSkout%stride(TSkout%rank))
       allocate(TSkout%offset(TSkout%rank))
       TSkout%dname = dname2d_kpStat
       TSkout%dimsf(1) = nx
       TSkout%dimsm  = TSkout%dimsf
       TSkout%block  = TSkout%dimsf
       TSkout%count  = 1
       TSkout%stride = 1
       TSkout%offset = 0
       TSkout%IsHSInitialized = .true.
!       TSkout%IsMultiGroup = .true.
       call WriteTSHDF5_1D(TSkout, buffer)

   end subroutine  WriteTSStat_kpStat


   subroutine WriteTSStat_kpInt(fn,kloc,nx,buffer)
      character(*), intent(in) :: fn
      integer, intent(in) :: kloc,nx
      real(8), intent(in) :: buffer(nx,nv2d_kpInt+1)
      character(4) :: fnum

       write(unit=fnum,fmt='(I04.4)') kloc

       if(trim(TSkout%fname).eq.'none') then
         TSkout%IsMultiGroup = .false.
       elseif(trim(TSkout%fname).eq.trim(fn)) then
         TSkout%IsMultiGroup = .true.
       else
         TSkout%IsMultiGroup = .false.
       endif

       TSkout%fname = trim(fn)
       TSkout%gname = "/"//fnum//'Int'
       TSkout%rank = 1
       TSkout%dnum = nv2d_kpInt+1

       if(allocated(TSkout%dname)) then
         deallocate(TSkout%dname,TSkout%dimsf,TSkout%dimsm)
         deallocate(TSkout%block,TSkout%count,TSkout%stride)
         deallocate(TSkout%offset)
       endif

       allocate(TSkout%dname(TSkout%dnum))
       allocate(TSkout%dimsf(TSkout%rank))
       allocate(TSkout%dimsm(TSkout%rank))
       allocate(TSkout%block(TSkout%rank))
       allocate(TSkout%count(TSkout%rank))
       allocate(TSkout%stride(TSkout%rank))
       allocate(TSkout%offset(TSkout%rank))
       TSkout%dname = dname2d_kpInt
       TSkout%dimsf(1) = nx
       TSkout%dimsm  = TSkout%dimsf
       TSkout%block  = TSkout%dimsf
       TSkout%count  = 1
       TSkout%stride = 1
       TSkout%offset = 0
       TSkout%IsHSInitialized = .true.

       call WriteTSHDF5_1D(TSkout, buffer)

   end subroutine  WriteTSStat_kpInt

   elemental real(8) function CalMu(tin)
     real(8), intent(in) :: tin

     CalMu = (1.458E-6)*((tin)**1.5)/(tin+110.4)

   end function CalMu


end module MTSStat_kplane

