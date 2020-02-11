

    module MTSAcousticSource
      use MTSHDF5
      implicit none

      integer, parameter :: nv_ASource = 12 + 4 + 3 + 2
      real(8), parameter :: rbar = 287.d0, gamma = 1.4, cp = 1004.5 ! (cp = 3.5*rbar)
      type tp_ASource_kplane
        integer :: nsample
        real(8), dimension(:), allocatable :: u, v, w, p, t
        real(8), dimension(:), allocatable :: dudx, dvdx, dwdx, dpdx, dtdx
        real(8), dimension(:), allocatable :: dudy, dvdy, dwdy, dpdy, dtdy
        real(8), dimension(:), allocatable :: dudz, dvdz, dwdz, dpdz, dtdz
        real(8), dimension(:), allocatable :: dudt, dvdt, dwdt, dpdt, dtdt

        real(8), dimension(:), allocatable :: dudi, dvdi, dwdi, dpdi, dtdi
        real(8), dimension(:), allocatable :: dudj, dvdj, dwdj, dpdj, dtdj
        real(8), dimension(:), allocatable :: dudk, dvdk, dwdk, dpdk, dtdk

        real(8), dimension(:), allocatable :: dsrc_Philipdx, dsrc_Philipdt

      end type tp_ASource_kplane
      type tp_ASource
        integer :: nsample
        real(8), dimension(:,:), allocatable :: u, v, w, p, t
        real(8), dimension(:,:), allocatable :: dudx, dvdx, dwdx, dpdx, dtdx
        real(8), dimension(:,:), allocatable :: dudy, dvdy, dwdy, dpdy, dtdy
        real(8), dimension(:,:), allocatable :: dudz, dvdz, dwdz, dpdz, dtdz
        real(8), dimension(:,:), allocatable :: dudt, dvdt, dwdt, dpdt, dtdt

        real(8), dimension(:,:), allocatable :: dudi, dvdi, dwdi, dpdi, dtdi
        real(8), dimension(:,:), allocatable :: dudj, dvdj, dwdj, dpdj, dtdj
        real(8), dimension(:,:), allocatable :: dudk, dvdk, dwdk, dpdk, dtdk

        real(8), dimension(:,:), allocatable :: dsrc_Philipdx, dsrc_Philipdt

      end type tp_ASource

      type(tp_ASource_kplane), private :: TS_AS_kplane
      type(tp_ASource), private :: TS_AS_iplane, TS_AS_jplane
      character(50), private :: dname_AS_kplane(nv_ASource+2), dname_AS_iplane(nv_ASource+2), dname_AS_jplane(nv_ASource+2)
      parameter( dname_AS_kplane = (/'kinetic', 'src_Philip_rms', 'linear_Philip_rms', 'upxupxrms', 'vpyvpyrms', 'wpzwpzrms', &
                                     'upyvpxrms_2x', 'upzwpxrms_2x', 'vpzwpyrms_2x', 'dis', 'Uc_src', 'nonlinear_Philips_rms', &
                                     'uave', 'wave', 'pave', 'tave', 'Klength', 'srcoKlen', 'divoKlen', 'Lambda', 'Re_Lambda', 'x', 'y'/) )
      parameter( dname_AS_iplane = (/'kinetic', 'src_Philip_rms', 'linear_Philip_rms', 'upxupxrms', 'vpyvpyrms', 'wpzwpzrms', &
                                     'upyvpxrms_2x', 'upzwpxrms_2x', 'vpzwpyrms_2x', 'dis', 'Uc_src', 'nonlinear_Philips_rms', &
                                     'uave', 'wave', 'pave', 'tave', 'Klength', 'srcoKlen', 'divoKlen', 'Lambda', 'Re_Lambda', 'y', 'z'/) )
      parameter( dname_AS_jplane = (/'kinetic', 'src_Philip_rms', 'linear_Philip_rms', 'upxupxrms', 'vpyvpyrms', 'wpzwpzrms', &
                                     'upyvpxrms_2x', 'upzwpxrms_2x', 'vpzwpyrms_2x', 'dis', 'Uc_src', 'nonlinear_Philips_rms', &
                                     'uave', 'wave', 'pave', 'tave', 'Klength', 'srcoKlen', 'divoKlen', 'Lambda', 'Re_Lambda', 'x', 'z'/) )

      type(tp_hyperslab), private :: TSAsource_kout, TSAsource_iout, TSAsource_jout

    contains


      subroutine InitTS_AS_jplane(ns, kdim)
        integer, intent(in) :: ns, kdim

        TSAsource_jout%fname = 'none'

        TS_AS_jplane%nsample = ns
        allocate(TS_AS_jplane%u(ns,kdim))
        allocate(TS_AS_jplane%v(ns,kdim))
        allocate(TS_AS_jplane%w(ns,kdim))
        allocate(TS_AS_jplane%p(ns,kdim))
        allocate(TS_AS_jplane%t(ns,kdim))

        allocate(TS_AS_jplane%dudx(ns,kdim), TS_AS_jplane%dudy(ns,kdim), TS_AS_jplane%dudz(ns,kdim))
        allocate(TS_AS_jplane%dvdx(ns,kdim), TS_AS_jplane%dvdy(ns,kdim), TS_AS_jplane%dvdz(ns,kdim))
        allocate(TS_AS_jplane%dwdx(ns,kdim), TS_AS_jplane%dwdy(ns,kdim), TS_AS_jplane%dwdz(ns,kdim))
!        allocate(TS_AS_jplane%dpdx(ns,kdim), TS_AS_jplane%dpdt(ns,kdim) )
!        allocate(TS_AS_jplane%dudt(ns,kdim), TS_AS_jplane%dvdt(ns,kdim), TS_AS_jplane%dwdt(ns,kdim))
!        allocate(TS_AS_jplane%dtdt(ns,kdim))
!        allocate(TS_AS_jplane%dudi(ns,kdim), TS_AS_jplane%dvdi(ns,kdim), TS_AS_jplane%dwdi(ns,kdim), &
!                 TS_AS_jplane%dpdi(ns,kdim), TS_AS_jplane%dtdi(ns,kdim))
!        allocate(TS_AS_jplane%dudj(ns,kdim), TS_AS_jplane%dvdj(ns,kdim), TS_AS_jplane%dwdj(ns,kdim), &
!                 TS_AS_jplane%dpdj(ns,kdim), TS_AS_jplane%dtdj(ns,kdim))
!        allocate(TS_AS_jplane%dudk(ns,kdim), TS_AS_jplane%dvdk(ns,kdim), TS_AS_jplane%dwdk(ns,kdim), &
!                 TS_AS_jplane%dpdk(ns,kdim), TS_AS_jplane%dtdk(ns,kdim))

        allocate(TS_AS_jplane%dudi(ns,kdim), TS_AS_jplane%dvdi(ns,kdim), TS_AS_jplane%dwdi(ns,kdim))
        allocate(TS_AS_jplane%dudj(ns,kdim), TS_AS_jplane%dvdj(ns,kdim), TS_AS_jplane%dwdj(ns,kdim))
        allocate(TS_AS_jplane%dudk(ns,kdim), TS_AS_jplane%dvdk(ns,kdim), TS_AS_jplane%dwdk(ns,kdim))

!        allocate(TS_AS_jplane%dsrc_Philipdx(ns,kdim), TS_AS_jplane%dsrc_Philipdt(ns,kdim))

      end subroutine InitTS_AS_jplane


      subroutine InitTS_AS_iplane(ns, kdim)
        integer, intent(in) :: ns, kdim

        TSAsource_iout%fname = 'none'

        TS_AS_iplane%nsample = ns
        allocate(TS_AS_iplane%u(ns,kdim))
        allocate(TS_AS_iplane%v(ns,kdim))
        allocate(TS_AS_iplane%w(ns,kdim))
        allocate(TS_AS_iplane%p(ns,kdim))
        allocate(TS_AS_iplane%t(ns,kdim))

        allocate(TS_AS_iplane%dudx(ns,kdim), TS_AS_iplane%dudy(ns,kdim), TS_AS_iplane%dudz(ns,kdim))
        allocate(TS_AS_iplane%dvdx(ns,kdim), TS_AS_iplane%dvdy(ns,kdim), TS_AS_iplane%dvdz(ns,kdim))
        allocate(TS_AS_iplane%dwdx(ns,kdim), TS_AS_iplane%dwdy(ns,kdim), TS_AS_iplane%dwdz(ns,kdim))
!        allocate(TS_AS_iplane%dpdx(ns,kdim), TS_AS_iplane%dpdt(ns,kdim) )
!        allocate(TS_AS_iplane%dudt(ns,kdim), TS_AS_iplane%dvdt(ns,kdim), TS_AS_iplane%dwdt(ns,kdim))
!        allocate(TS_AS_iplane%dtdt(ns,kdim))
!        allocate(TS_AS_iplane%dudi(ns,kdim), TS_AS_iplane%dvdi(ns,kdim), TS_AS_iplane%dwdi(ns,kdim), &
!                 TS_AS_iplane%dpdi(ns,kdim), TS_AS_iplane%dtdi(ns,kdim))
!        allocate(TS_AS_iplane%dudj(ns,kdim), TS_AS_iplane%dvdj(ns,kdim), TS_AS_iplane%dwdj(ns,kdim), &
!                 TS_AS_iplane%dpdj(ns,kdim), TS_AS_iplane%dtdj(ns,kdim))
!        allocate(TS_AS_iplane%dudk(ns,kdim), TS_AS_iplane%dvdk(ns,kdim), TS_AS_iplane%dwdk(ns,kdim), &
!                 TS_AS_iplane%dpdk(ns,kdim), TS_AS_iplane%dtdk(ns,kdim))

        allocate(TS_AS_iplane%dudi(ns,kdim), TS_AS_iplane%dvdi(ns,kdim), TS_AS_iplane%dwdi(ns,kdim))
        allocate(TS_AS_iplane%dudj(ns,kdim), TS_AS_iplane%dvdj(ns,kdim), TS_AS_iplane%dwdj(ns,kdim))
        allocate(TS_AS_iplane%dudk(ns,kdim), TS_AS_iplane%dvdk(ns,kdim), TS_AS_iplane%dwdk(ns,kdim))

!        allocate(TS_AS_iplane%dsrc_Philipdx(ns,kdim), TS_AS_iplane%dsrc_Philipdt(ns,kdim))

      end subroutine InitTS_AS_iplane

      subroutine InitTS_AS_kplane(ns)
        integer, intent(in) :: ns

        TSAsource_kout%fname = 'none'

        TS_AS_kplane%nsample = ns
        allocate(TS_AS_kplane%u(ns))
        allocate(TS_AS_kplane%v(ns))
        allocate(TS_AS_kplane%w(ns))
        allocate(TS_AS_kplane%p(ns))
        allocate(TS_AS_kplane%t(ns))
        allocate(TS_AS_kplane%dudx(ns), TS_AS_kplane%dudy(ns), TS_AS_kplane%dudz(ns))
        allocate(TS_AS_kplane%dvdx(ns), TS_AS_kplane%dvdy(ns), TS_AS_kplane%dvdz(ns))
        allocate(TS_AS_kplane%dwdx(ns), TS_AS_kplane%dwdy(ns), TS_AS_kplane%dwdz(ns))
        allocate(TS_AS_kplane%dpdx(ns), TS_AS_kplane%dpdt(ns) )
        allocate(TS_AS_kplane%dudt(ns), TS_AS_kplane%dvdt(ns), TS_AS_kplane%dwdt(ns))
        allocate(TS_AS_kplane%dtdt(ns))

        allocate(TS_AS_kplane%dudi(ns), TS_AS_kplane%dvdi(ns), TS_AS_kplane%dwdi(ns), TS_AS_kplane%dpdi(ns), TS_AS_kplane%dtdi(ns))
        allocate(TS_AS_kplane%dudj(ns), TS_AS_kplane%dvdj(ns), TS_AS_kplane%dwdj(ns), TS_AS_kplane%dpdj(ns), TS_AS_kplane%dtdj(ns))
        allocate(TS_AS_kplane%dudk(ns), TS_AS_kplane%dvdk(ns), TS_AS_kplane%dwdk(ns), TS_AS_kplane%dpdk(ns), TS_AS_kplane%dtdk(ns))

        allocate(TS_AS_kplane%dsrc_Philipdx(ns), TS_AS_kplane%dsrc_Philipdt(ns))

      end subroutine InitTS_AS_kplane

      subroutine CalTS_AS_Philip_jplane(nt, nx, nz, nv, buffer, buffer_grid, src_Philip_jplane)
        integer, intent(in) :: nt, nx, nz, nv
        real(8), intent(in) :: buffer(nt,nx,nz,nv), buffer_grid(nx,nz,9)
        real(8), intent(out) :: src_Philip_jplane(nt,3:nx-2,nz)
        integer :: i, j, k, n, nn
        real(8) :: didx, djdx, dkdx, didy, djdy, dkdy, didz, djdz, dkdz
        real(8) :: dudi, dvdi, dwdi, dudj, dvdj, dwdj, dudk, dvdk, dwdk
        real(8) :: dpdi, dpdj, dpdk
        real(8) :: dudx, dudy, dudz
        real(8) :: dvdx, dvdy, dvdz
        real(8) :: dwdx, dwdy, dwdz
        real(8), parameter :: c12i = 1.d0/12.d0

        src_Philip_jplane = 0.d0
        do i = 3, nx-2
          do k = 1, nz
            didx = buffer_grid(i,k,1)
            djdx = buffer_grid(i,k,2)
            dkdx = buffer_grid(i,k,3)
            didy = buffer_grid(i,k,4)
            djdy = buffer_grid(i,k,5)
            dkdy = buffer_grid(i,k,6)
            didz = buffer_grid(i,k,7)
            djdz = buffer_grid(i,k,8)
            dkdz = buffer_grid(i,k,9)
            do n=1, nt
              dudi = ( buffer(n,i-2,k,1) - 8.0*(buffer(n,i-1,k,1) - buffer(n,i+1,k,1)) - buffer(n,i+2,k,1) )*c12i
              dvdi = ( buffer(n,i-2,k,2) - 8.0*(buffer(n,i-1,k,2) - buffer(n,i+1,k,2)) - buffer(n,i+2,k,2) )*c12i
              dwdi = ( buffer(n,i-2,k,3) - 8.0*(buffer(n,i-1,k,3) - buffer(n,i+1,k,3)) - buffer(n,i+2,k,3) )*c12i

              dudj = buffer(n,i,k,6)
              dvdj = buffer(n,i,k,7)
              dwdj = buffer(n,i,k,8)

              if(k.ge.3.and.k.le.nz-2) then
                dudk = (buffer(n,i,k-2,1) - 8.0*(buffer(n,i,k-1,1) - buffer(n,i,k+1,1)) - buffer(n,i,k+2,1))*c12i
                dvdk = (buffer(n,i,k-2,2) - 8.0*(buffer(n,i,k-1,2) - buffer(n,i,k+1,2)) - buffer(n,i,k+2,2))*c12i
                dwdk = (buffer(n,i,k-2,3) - 8.0*(buffer(n,i,k-1,3) - buffer(n,i,k+1,3)) - buffer(n,i,k+2,3))*c12i
              elseif(k.eq.1) then
                dudk = -25./12.*buffer(n,i,1,1) + 4.*buffer(n,i,2,1) - 3.*buffer(n,i,3,1) + 4./3.*buffer(n,i,4,1) - 0.25*buffer(n,i,5,1)
                dvdk = -25./12.*buffer(n,i,1,2) + 4.*buffer(n,i,2,2) - 3.*buffer(n,i,3,2) + 4./3.*buffer(n,i,4,2) - 0.25*buffer(n,i,5,2)
                dwdk = -25./12.*buffer(n,i,1,3) + 4.*buffer(n,i,2,3) - 3.*buffer(n,i,3,3) + 4./3.*buffer(n,i,4,3) - 0.25*buffer(n,i,5,3)
              elseif(k.eq.2) then
                dudk = -0.25*buffer(n,i,1,1) - 5./6.*buffer(n,i,2,1) + 1.5*buffer(n,i,3,1) - 0.5*buffer(n,i,4,1) + 1./12.*buffer(n,i,5,1)
                dvdk = -0.25*buffer(n,i,1,2) - 5./6.*buffer(n,i,2,2) + 1.5*buffer(n,i,3,2) - 0.5*buffer(n,i,4,2) + 1./12.*buffer(n,i,5,2)
                dwdk = -0.25*buffer(n,i,1,3) - 5./6.*buffer(n,i,2,3) + 1.5*buffer(n,i,3,3) - 0.5*buffer(n,i,4,3) + 1./12.*buffer(n,i,5,3)
              elseif(k.eq.nz-1) then
                dudk = -1./12.*buffer(n,i,k-4,1) + 0.5*buffer(n,i,k-3,1) - 1.5*buffer(n,i,k-2,1) + 5./6.*buffer(n,i,k-1,1) + 0.25*buffer(n,i,k,1)
                dvdk = -1./12.*buffer(n,i,k-4,2) + 0.5*buffer(n,i,k-3,2) - 1.5*buffer(n,i,k-2,2) + 5./6.*buffer(n,i,k-1,2) + 0.25*buffer(n,i,k,2)
                dwdk = -1./12.*buffer(n,i,k-4,3) + 0.5*buffer(n,i,k-3,3) - 1.5*buffer(n,i,k-2,3) + 5./6.*buffer(n,i,k-1,3) + 0.25*buffer(n,i,k,3)
              elseif(k.eq.nz) then
                dudk = 0.25*buffer(n,i,k-4,1) - 4./3.*buffer(n,i,k-3,1) + 3.*buffer(n,i,k-2,1) - 4.*buffer(n,i,k-1,1) + 25./12.*buffer(n,i,k,1)
                dvdk = 0.25*buffer(n,i,k-4,2) - 4./3.*buffer(n,i,k-3,2) + 3.*buffer(n,i,k-2,2) - 4.*buffer(n,i,k-1,1) + 25./12.*buffer(n,i,k,2)
                dwdk = 0.25*buffer(n,i,k-4,3) - 4./3.*buffer(n,i,k-3,3) + 3.*buffer(n,i,k-2,3) - 4.*buffer(n,i,k-1,1) + 25./12.*buffer(n,i,k,3)
              endif
              dudx = dudi*didx + dudj*djdx + dudk*dkdx
              dvdx = dvdi*didx + dvdj*djdx + dvdk*dkdx
              dwdx = dwdi*didx + dwdj*djdx + dwdk*dkdx

              dudy = dudi*didy + dudj*djdy + dudk*dkdy
              dvdy = dvdi*didy + dvdj*djdy + dvdk*dkdy
              dwdy = dwdi*didy + dwdj*djdy + dwdk*dkdy

              dudz = dudi*didz + dudj*djdz + dudk*dkdz
              dvdz = dvdi*didz + dvdj*djdz + dvdk*dkdz
              dwdz = dwdi*didz + dwdj*djdz + dwdk*dkdz

              src_Philip_jplane(n,i,k) = gamma*(dudx**2 + dvdy**2 + dwdz**2 + 2.d0*(dudy*dvdx + dudz*dwdx + dvdz*dwdy))
            enddo ! end n loop
          enddo ! end k loop
        enddo ! end i loop

      end subroutine CalTS_AS_Philip_jplane

      subroutine CalTS_AS_Philip_iplane(nt, ny, nz, nv, buffer, buffer_grid, src_Philip_iplane)
        integer, intent(in) :: nt, ny, nz, nv
        real(8), intent(in) :: buffer(nt,ny,nz,nv), buffer_grid(ny,nz,9)
        real(8), intent(out) :: src_Philip_iplane(nt,ny,nz)
        integer :: i, j, k, n, nn
        real(8) :: didx, djdx, dkdx, didy, djdy, dkdy, didz, djdz, dkdz
        real(8) :: dudi, dvdi, dwdi, dudj, dvdj, dwdj, dudk, dvdk, dwdk
        real(8) :: dpdi, dpdj, dpdk
        real(8) :: dudx, dudy, dudz
        real(8) :: dvdx, dvdy, dvdz
        real(8) :: dwdx, dwdy, dwdz
        real(8), parameter :: c12i = 1.d0/12.d0

        src_Philip_iplane = 0.d0
        do j = 3, ny-2
          do k = 1, nz
            didx = buffer_grid(j,k,1)
            djdx = buffer_grid(j,k,2)
            dkdx = buffer_grid(j,k,3)
            didy = buffer_grid(j,k,4)
            djdy = buffer_grid(j,k,5)
            dkdy = buffer_grid(j,k,6)
            didz = buffer_grid(j,k,7)
            djdz = buffer_grid(j,k,8)
            dkdz = buffer_grid(j,k,9)
            do n=1, nt
              dudi = buffer(n,j,k,6)
              dvdi = buffer(n,j,k,7)
              dwdi = buffer(n,j,k,8)

              dudj = ( buffer(n,j-2,k,1) - 8.0*(buffer(n,j-1,k,1) - buffer(n,j+1,k,1)) - buffer(n,j+2,k,1) )*c12i
              dvdj = ( buffer(n,j-2,k,2) - 8.0*(buffer(n,j-1,k,2) - buffer(n,j+1,k,2)) - buffer(n,j+2,k,2) )*c12i
              dwdj = ( buffer(n,j-2,k,3) - 8.0*(buffer(n,j-1,k,3) - buffer(n,j+1,k,3)) - buffer(n,j+2,k,3) )*c12i

              if(k.ge.3.and.k.le.nz-2) then
                dudk = (buffer(n,j,k-2,1) - 8.0*(buffer(n,j,k-1,1) - buffer(n,j,k+1,1)) - buffer(n,j,k+2,1))*c12i
                dvdk = (buffer(n,j,k-2,2) - 8.0*(buffer(n,j,k-1,2) - buffer(n,j,k+1,2)) - buffer(n,j,k+2,2))*c12i
                dwdk = (buffer(n,j,k-2,3) - 8.0*(buffer(n,j,k-1,3) - buffer(n,j,k+1,3)) - buffer(n,j,k+2,3))*c12i
              elseif(k.eq.1) then
                dudk = -25./12.*buffer(n,j,1,1) + 4.*buffer(n,j,2,1) - 3.*buffer(n,j,3,1) + 4./3.*buffer(n,j,4,1) - 0.25*buffer(n,j,5,1)
                dvdk = -25./12.*buffer(n,j,1,2) + 4.*buffer(n,j,2,2) - 3.*buffer(n,j,3,2) + 4./3.*buffer(n,j,4,2) - 0.25*buffer(n,j,5,2)
                dwdk = -25./12.*buffer(n,j,1,3) + 4.*buffer(n,j,2,3) - 3.*buffer(n,j,3,3) + 4./3.*buffer(n,j,4,3) - 0.25*buffer(n,j,5,3)
              elseif(k.eq.2) then
                dudk = -0.25*buffer(n,j,1,1) - 5./6.*buffer(n,j,2,1) + 1.5*buffer(n,j,3,1) - 0.5*buffer(n,j,4,1) + 1./12.*buffer(n,j,5,1)
                dvdk = -0.25*buffer(n,j,1,2) - 5./6.*buffer(n,j,2,2) + 1.5*buffer(n,j,3,2) - 0.5*buffer(n,j,4,2) + 1./12.*buffer(n,j,5,2)
                dwdk = -0.25*buffer(n,j,1,3) - 5./6.*buffer(n,j,2,3) + 1.5*buffer(n,j,3,3) - 0.5*buffer(n,j,4,3) + 1./12.*buffer(n,j,5,3)
              elseif(k.eq.nz-1) then
                dudk = -1./12.*buffer(n,j,k-4,1) + 0.5*buffer(n,j,k-3,1) - 1.5*buffer(n,j,k-2,1) + 5./6.*buffer(n,j,k-1,1) + 0.25*buffer(n,j,k,1)
                dvdk = -1./12.*buffer(n,j,k-4,2) + 0.5*buffer(n,j,k-3,2) - 1.5*buffer(n,j,k-2,2) + 5./6.*buffer(n,j,k-1,2) + 0.25*buffer(n,j,k,2)
                dwdk = -1./12.*buffer(n,j,k-4,3) + 0.5*buffer(n,j,k-3,3) - 1.5*buffer(n,j,k-2,3) + 5./6.*buffer(n,j,k-1,3) + 0.25*buffer(n,j,k,3)
              elseif(k.eq.nz) then
                dudk = 0.25*buffer(n,j,k-4,1) - 4./3.*buffer(n,j,k-3,1) + 3.*buffer(n,j,k-2,1) - 4.*buffer(n,j,k-1,1) + 25./12.*buffer(n,j,k,1)
                dvdk = 0.25*buffer(n,j,k-4,2) - 4./3.*buffer(n,j,k-3,2) + 3.*buffer(n,j,k-2,2) - 4.*buffer(n,j,k-1,1) + 25./12.*buffer(n,j,k,2)
                dwdk = 0.25*buffer(n,j,k-4,3) - 4./3.*buffer(n,j,k-3,3) + 3.*buffer(n,j,k-2,3) - 4.*buffer(n,j,k-1,1) + 25./12.*buffer(n,j,k,3)
              endif
              dudx = dudi*didx + dudj*djdx + dudk*dkdx
              dvdx = dvdi*didx + dvdj*djdx + dvdk*dkdx
              dwdx = dwdi*didx + dwdj*djdx + dwdk*dkdx

              dudy = dudi*didy + dudj*djdy + dudk*dkdy
              dvdy = dvdi*didy + dvdj*djdy + dvdk*dkdy
              dwdy = dwdi*didy + dwdj*djdy + dwdk*dkdy

              dudz = dudi*didz + dudj*djdz + dudk*dkdz
              dvdz = dvdi*didz + dvdj*djdz + dvdk*dkdz
              dwdz = dwdi*didz + dwdj*djdz + dwdk*dkdz

              !src_Philip_iplane(n,j,k) = gamma*(dudx**2 + dvdy**2 + dwdz**2 + 2.d0*(dudy*dvdx + dudz*dwdx + dvdz*dwdy))
              src_Philip_iplane(n,j,k) = dudx**2 + dvdy**2 + dwdz**2 + 2.d0*(dudy*dvdx + dudz*dwdx + dvdz*dwdy)
            enddo ! end n loop
          enddo ! end k loop
        enddo ! end j loop

      end subroutine CalTS_AS_Philip_iplane


      subroutine CalTS_AS_Philip_kplane(nt, nx, ny, nv, buffer, buffer_grid, src_Philip_kplane)
        integer, intent(in) :: nt, nx, ny, nv
!     real(8), intent(in) :: buffer(nt,ny,-2:2,nv), buffer_grid(ny, 9), dt
        real(8), intent(in) :: buffer(nt, ny, nx, nv), buffer_grid(ny,nx,9)
!        real(8), intent(out), dimension(:,:,:), allocatable :: src_Philip_kplane
        real(8), intent(out) :: src_Philip_kplane(nt,ny,3:nx-2)
        integer :: i, j, n, nn
        real(8) :: didx, djdx, dkdx, didy, djdy, dkdy, didz, djdz, dkdz
        real(8) :: dudi, dvdi, dwdi, dudj, dvdj, dwdj, dudk, dvdk, dwdk
        real(8) :: dpdi, dpdj, dpdk
        real(8) :: dudx, dudy, dudz
        real(8) :: dvdx, dvdy, dvdz
        real(8) :: dwdx, dwdy, dwdz
        real(8), parameter :: c12i = 1.d0/12.d0

!        allocate(src_Philip_kplane(nt,ny,nx))
        src_Philip_kplane = 0.d0
        do i = 3, nx-2 !!!!!!!!
          do j = 3, ny-2
           ! didx, djdx, dkdx, didy, djdy, dkdy, didz, djdz, dkdz
            didx = buffer_grid(j,i,1)
            djdx = buffer_grid(j,i,2)
            dkdx = buffer_grid(j,i,3)
            didy = buffer_grid(j,i,4)
            djdy = buffer_grid(j,i,5)
            dkdy = buffer_grid(j,i,6)
            didz = buffer_grid(j,i,7)
            djdz = buffer_grid(j,i,8)
            dkdz = buffer_grid(j,i,9)
            do n = 1, nt
               dudi = ( buffer(n,j,i-2,1) - 8.0*(buffer(n,j,i-1,1) - buffer(n,j,i+1,1)) - buffer(n,j,i+2,1) )*c12i
               dvdi = ( buffer(n,j,i-2,2) - 8.0*(buffer(n,j,i-1,2) - buffer(n,j,i+1,2)) - buffer(n,j,i+2,2) )*c12i
               dwdi = ( buffer(n,j,i-2,3) - 8.0*(buffer(n,j,i-1,3) - buffer(n,j,i+1,3)) - buffer(n,j,i+2,3) )*c12i

               dudj = ( buffer(n,j-2,i,1) - 8.0*(buffer(n,j-1,i,1) - buffer(n,j+1,i,1)) - buffer(n,j+2,i,1) )*c12i
               dvdj = ( buffer(n,j-2,i,2) - 8.0*(buffer(n,j-1,i,2) - buffer(n,j+1,i,2)) - buffer(n,j+2,i,2) )*c12i
               dwdj = ( buffer(n,j-2,i,3) - 8.0*(buffer(n,j-1,i,3) - buffer(n,j+1,i,3)) - buffer(n,j+2,i,3) )*c12i

               dudk = buffer(n,j,i,6)
               dvdk = buffer(n,j,i,7)
               dwdk = buffer(n,j,i,8)

               dudx = dudi*didx + dudj*djdx + dudk*dkdx
               dvdx = dvdi*didx + dvdj*djdx + dvdk*dkdx
               dwdx = dwdi*didx + dwdj*djdx + dwdk*dkdx

               dudy = dudi*didy + dudj*djdy + dudk*dkdy
               dvdy = dvdi*didy + dvdj*djdy + dvdk*dkdy
               dwdy = dwdi*didy + dwdj*djdy + dwdk*dkdy

               dudz = dudi*didz + dudj*djdz + dudk*dkdz
               dvdz = dvdi*didz + dvdj*djdz + dvdk*dkdz
               dwdz = dwdi*didz + dwdj*djdz + dwdk*dkdz

               src_Philip_kplane(n,j,i) = gamma*(dudx**2 + dvdy**2 + dwdz**2 + 2.d0*(dudy*dvdx + dudz*dwdx + dvdz*dwdy))

            enddo ! en n loop
          enddo ! end j loop
        enddo ! end i loop

      end subroutine CalTS_AS_Philip_kplane

      subroutine CalTSDeriv_AS_jplane(nt,nz,nv,buffer,buffer_grid,dt,src_Philip_jplane)
        integer, intent(in) :: nt, nz, nv
        real(8), intent(in) :: buffer(nt,-2:2,nz,nv), buffer_grid(nz,9), dt, src_Philip_jplane(nt,-2:2,nz)
        integer :: i, j, k, n, nn
        real(8) :: didx, djdx, dkdx, didy, djdy, dkdy, didz, djdz, dkdz
        real(8) :: dsrc_Philipdi
        real(8), parameter :: c12i = 1.d0/12.d0
        real(8) :: dti

        dti = 1.d0/dt

        do k=1, nz
          TS_AS_jplane%u(:,k) = reshape(buffer(3:nt-2,0,k,1),(/TS_AS_jplane%nsample/))
          TS_AS_jplane%v(:,k) = reshape(buffer(3:nt-2,0,k,2),(/TS_AS_jplane%nsample/))
          TS_AS_jplane%w(:,k) = reshape(buffer(3:nt-2,0,k,3),(/TS_AS_jplane%nsample/))
          TS_AS_jplane%p(:,k) = reshape(buffer(3:nt-2,0,k,4),(/TS_AS_jplane%nsample/))
          TS_AS_jplane%t(:,k) = reshape(buffer(3:nt-2,0,k,5),(/TS_AS_jplane%nsample/))
        enddo

        do k=1, nz
        ! didx, djdx, dkdx, didy, djdy, dkdy, didz, djdz, dkdz
          didx = buffer_grid(k,1)
          djdx = buffer_grid(k,2)
          dkdx = buffer_grid(k,3)
          didy = buffer_grid(k,4)
          djdy = buffer_grid(k,5)
          dkdy = buffer_grid(k,6)
          didz = buffer_grid(k,7)
          djdz = buffer_grid(k,8)
          dkdz = buffer_grid(k,9)

          do n=3, nt-2
            TS_AS_jplane%dudi(n-2,k) = (buffer(n,-2,k,1) - 8.0*(buffer(n,-1,k,1) - buffer(n,1,k,1)) - buffer(n,2,k,1))*c12i
            TS_AS_jplane%dvdi(n-2,k) = (buffer(n,-2,k,2) - 8.0*(buffer(n,-1,k,2) - buffer(n,1,k,2)) - buffer(n,2,k,2))*c12i
            TS_AS_jplane%dwdi(n-2,k) = (buffer(n,-2,k,3) - 8.0*(buffer(n,-1,k,3) - buffer(n,1,k,3)) - buffer(n,2,k,3))*c12i

            TS_AS_jplane%dudj(n-2,k) = buffer(n,0,k,6)
            TS_AS_jplane%dvdj(n-2,k) = buffer(n,0,k,7)
            TS_AS_jplane%dwdj(n-2,k) = buffer(n,0,k,8)

            if(k.ge.3.and.k.le.nz-2) then
              TS_AS_jplane%dudk(n-2,k) = (buffer(n,0,k-2,1) - 8.0*(buffer(n,0,k-1,1) - buffer(n,0,k+1,1)) - buffer(n,0,k+2,1))*c12i
              TS_AS_jplane%dvdk(n-2,k) = (buffer(n,0,k-2,2) - 8.0*(buffer(n,0,k-1,2) - buffer(n,0,k+1,2)) - buffer(n,0,k+2,2))*c12i
              TS_AS_jplane%dwdk(n-2,k) = (buffer(n,0,k-2,3) - 8.0*(buffer(n,0,k-1,3) - buffer(n,0,k+1,3)) - buffer(n,0,k+2,3))*c12i
            elseif(k.eq.1) then
              TS_AS_jplane%dudk(n-2,k) = -25./12.*buffer(n,0,1,1) + 4.*buffer(n,0,2,1) - 3.*buffer(n,0,3,1) + 4./3.*buffer(n,0,4,1) - 0.25*buffer(n,0,5,1)
              TS_AS_jplane%dvdk(n-2,k) = -25./12.*buffer(n,0,1,2) + 4.*buffer(n,0,2,2) - 3.*buffer(n,0,3,2) + 4./3.*buffer(n,0,4,2) - 0.25*buffer(n,0,5,2)
              TS_AS_jplane%dwdk(n-2,k) = -25./12.*buffer(n,0,1,3) + 4.*buffer(n,0,2,3) - 3.*buffer(n,0,3,3) + 4./3.*buffer(n,0,4,3) - 0.25*buffer(n,0,5,3)
            elseif(k.eq.2) then
              TS_AS_jplane%dudk(n-2,k) = -0.25*buffer(n,0,1,1) - 5./6.*buffer(n,0,2,1) + 1.5*buffer(n,0,3,1) - 0.5*buffer(n,0,4,1) + 1./12.*buffer(n,0,5,1)
              TS_AS_jplane%dvdk(n-2,k) = -0.25*buffer(n,0,1,2) - 5./6.*buffer(n,0,2,2) + 1.5*buffer(n,0,3,2) - 0.5*buffer(n,0,4,2) + 1./12.*buffer(n,0,5,2)
              TS_AS_jplane%dwdk(n-2,k) = -0.25*buffer(n,0,1,3) - 5./6.*buffer(n,0,2,3) + 1.5*buffer(n,0,3,3) - 0.5*buffer(n,0,4,3) + 1./12.*buffer(n,0,5,3)
            elseif(k.eq.nz-1) then
              TS_AS_jplane%dudk(n-2,k) = -1./12.*buffer(n,0,k-4,1) + 0.5*buffer(n,0,k-3,1) - 1.5*buffer(n,0,k-2,1) + 5./6.*buffer(n,0,k-1,1) + 0.25*buffer(n,0,k,1)
              TS_AS_jplane%dvdk(n-2,k) = -1./12.*buffer(n,0,k-4,2) + 0.5*buffer(n,0,k-3,2) - 1.5*buffer(n,0,k-2,2) + 5./6.*buffer(n,0,k-1,2) + 0.25*buffer(n,0,k,2)
              TS_AS_jplane%dwdk(n-2,k) = -1./12.*buffer(n,0,k-4,3) + 0.5*buffer(n,0,k-3,3) - 1.5*buffer(n,0,k-2,3) + 5./6.*buffer(n,0,k-1,3) + 0.25*buffer(n,0,k,3)
            elseif(k.eq.nz) then
              TS_AS_jplane%dudk(n-2,k) = 0.25*buffer(n,0,k-4,1) - 4./3.*buffer(n,0,k-3,1) + 3.*buffer(n,0,k-2,1) - 4.*buffer(n,0,k-1,1) + 25./12.*buffer(n,0,k,1)
              TS_AS_jplane%dvdk(n-2,k) = 0.25*buffer(n,0,k-4,2) - 4./3.*buffer(n,0,k-3,2) + 3.*buffer(n,0,k-2,2) - 4.*buffer(n,0,k-1,2) + 25./12.*buffer(n,0,k,2)
              TS_AS_jplane%dwdk(n-2,k) = 0.25*buffer(n,0,k-4,3) - 4./3.*buffer(n,0,k-3,3) + 3.*buffer(n,0,k-2,3) - 4.*buffer(n,0,k-1,3) + 25./12.*buffer(n,0,k,3)
            endif

            TS_AS_jplane%dudx(n-2,k) = TS_AS_jplane%dudi(n-2,k)*didx + TS_AS_jplane%dudj(n-2,k)*djdx + TS_AS_jplane%dudk(n-2,k)*dkdx
            TS_AS_jplane%dvdx(n-2,k) = TS_AS_jplane%dvdi(n-2,k)*didx + TS_AS_jplane%dvdj(n-2,k)*djdx + TS_AS_jplane%dvdk(n-2,k)*dkdx
            TS_AS_jplane%dwdx(n-2,k) = TS_AS_jplane%dwdi(n-2,k)*didx + TS_AS_jplane%dwdj(n-2,k)*djdx + TS_AS_jplane%dwdk(n-2,k)*dkdx

            TS_AS_jplane%dudy(n-2,k) = TS_AS_jplane%dudi(n-2,k)*didy + TS_AS_jplane%dudj(n-2,k)*djdy + TS_AS_jplane%dudk(n-2,k)*dkdy
            TS_AS_jplane%dvdy(n-2,k) = TS_AS_jplane%dvdi(n-2,k)*didy + TS_AS_jplane%dvdj(n-2,k)*djdy + TS_AS_jplane%dvdk(n-2,k)*dkdy
            TS_AS_jplane%dwdy(n-2,k) = TS_AS_jplane%dwdi(n-2,k)*didy + TS_AS_jplane%dwdj(n-2,k)*djdy + TS_AS_jplane%dwdk(n-2,k)*dkdy

            TS_AS_jplane%dudz(n-2,k) = TS_AS_jplane%dudi(n-2,k)*didz + TS_AS_jplane%dudj(n-2,k)*djdz + TS_AS_jplane%dudk(n-2,k)*dkdz
            TS_AS_jplane%dvdz(n-2,k) = TS_AS_jplane%dvdi(n-2,k)*didz + TS_AS_jplane%dvdj(n-2,k)*djdz + TS_AS_jplane%dvdk(n-2,k)*dkdz
            TS_AS_jplane%dwdz(n-2,k) = TS_AS_jplane%dwdi(n-2,k)*didz + TS_AS_jplane%dwdj(n-2,k)*djdz + TS_AS_jplane%dwdk(n-2,k)*dkdz

!            dsrc_Philipdi = ( src_Philip_jplane(n,-2,k) - 8.0*(src_Philip_jplane(n,-1,k) - src_Philip_jplane(n,1,k)) &
!                            - src_Philip_jplane(n,2,k) )*c12i
!            TS_AS_jplane%dsrc_Philipdx(n-2,k) = dsrc_Philipdi*didx
!            TS_AS_jplane%dsrc_Philipdt(n-2,k) = ( src_Philip_jplane(n-2,0,k) - 8.0*(src_Philip_jplane(n-1,0,k)) &
!                                                - src_Philip_jplane(n+1,0,k) - src_Philip_jplane(n+2,0,k) )*c12i*dti

          enddo
        enddo

      end subroutine CalTSDeriv_AS_jplane


      subroutine CalTSDeriv_AS_iplane(nt,ny,nz,nv,buffer,buffer_grid,dt,src_Philip_iplane)
        integer, intent(in) :: nt, ny, nz, nv
        real(8), intent(in) :: buffer(nt,ny,nz,nv), buffer_grid(ny,nz,9), dt, src_Philip_iplane(nt,ny,nz)
        integer :: i, j, k, n, nn
        real(8) :: didx, djdx, dkdx, didy, djdy, dkdy, didz, djdz, dkdz
        real(8) :: dsrc_Philipdi
        real(8), parameter :: c12i = 1.d0/12.d0
        real(8) :: dti

        dti = 1.d0/dt

        do k=1, nz
          TS_AS_iplane%u(:,k) = reshape(buffer(3:nt-2,3:ny-2,k,1),(/TS_AS_iplane%nsample/))
          TS_AS_iplane%v(:,k) = reshape(buffer(3:nt-2,3:ny-2,k,2),(/TS_AS_iplane%nsample/))
          TS_AS_iplane%w(:,k) = reshape(buffer(3:nt-2,3:ny-2,k,3),(/TS_AS_iplane%nsample/))
          TS_AS_iplane%p(:,k) = reshape(buffer(3:nt-2,3:ny-2,k,4),(/TS_AS_iplane%nsample/))
          TS_AS_iplane%t(:,k) = reshape(buffer(3:nt-2,3:ny-2,k,5),(/TS_AS_iplane%nsample/))
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

              TS_AS_iplane%dudi(nn,k) = buffer(n,j,k,6)
              TS_AS_iplane%dvdi(nn,k) = buffer(n,j,k,7)
              TS_AS_iplane%dwdi(nn,k) = buffer(n,j,k,8)

              TS_AS_iplane%dudj(nn,k) = ( buffer(n,j-2,k,1) - 8.0*(buffer(n,j-1,k,1) - buffer(n,j+1,k,1)) - buffer(n,j+2,k,1) )*c12i
              TS_AS_iplane%dvdj(nn,k) = ( buffer(n,j-2,k,2) - 8.0*(buffer(n,j-1,k,2) - buffer(n,j+1,k,2)) - buffer(n,j+2,k,2) )*c12i
              TS_AS_iplane%dwdj(nn,k) = ( buffer(n,j-2,k,3) - 8.0*(buffer(n,j-1,k,3) - buffer(n,j+1,k,3)) - buffer(n,j+2,k,3) )*c12i

              if(k.ge.3.and.k.le.nz-2) then
                TS_AS_iplane%dudk(nn,k) = (buffer(n,j,k-2,1) - 8.0*(buffer(n,j,k-1,1) - buffer(n,j,k+1,1)) - buffer(n,j,k+2,1))*c12i
                TS_AS_iplane%dvdk(nn,k) = (buffer(n,j,k-2,2) - 8.0*(buffer(n,j,k-1,2) - buffer(n,j,k+1,2)) - buffer(n,j,k+2,2))*c12i
                TS_AS_iplane%dwdk(nn,k) = (buffer(n,j,k-2,3) - 8.0*(buffer(n,j,k-1,3) - buffer(n,j,k+1,3)) - buffer(n,j,k+2,3))*c12i
              elseif(k.eq.1) then
                TS_AS_iplane%dudk(nn,k) = -25./12.*buffer(n,j,1,1) + 4.*buffer(n,j,2,1) - 3.*buffer(n,j,3,1) + 4./3.*buffer(n,j,4,1) - 0.25*buffer(n,j,5,1)
                TS_AS_iplane%dvdk(nn,k) = -25./12.*buffer(n,j,1,2) + 4.*buffer(n,j,2,2) - 3.*buffer(n,j,3,2) + 4./3.*buffer(n,j,4,2) - 0.25*buffer(n,j,5,2)
                TS_AS_iplane%dwdk(nn,k) = -25./12.*buffer(n,j,1,3) + 4.*buffer(n,j,2,3) - 3.*buffer(n,j,3,3) + 4./3.*buffer(n,j,4,3) - 0.25*buffer(n,j,5,3)
              elseif(k.eq.2) then
                TS_AS_iplane%dudk(nn,k) = -0.25*buffer(n,j,1,1) - 5./6.*buffer(n,j,2,1) + 1.5*buffer(n,j,3,1) - 0.5*buffer(n,j,4,1) + 1./12.*buffer(n,j,5,1)
                TS_AS_iplane%dvdk(nn,k) = -0.25*buffer(n,j,1,2) - 5./6.*buffer(n,j,2,2) + 1.5*buffer(n,j,3,2) - 0.5*buffer(n,j,4,2) + 1./12.*buffer(n,j,5,2)
                TS_AS_iplane%dwdk(nn,k) = -0.25*buffer(n,j,1,3) - 5./6.*buffer(n,j,2,3) + 1.5*buffer(n,j,3,3) - 0.5*buffer(n,j,4,3) + 1./12.*buffer(n,j,5,3)

              elseif(k.eq.nz-1) then
                TS_AS_iplane%dudk(nn,k) = -1./12.*buffer(n,j,k-4,1) + 0.5*buffer(n,j,k-3,1) - 1.5*buffer(n,j,k-2,1) + 5./6.*buffer(n,j,k-1,1) + 0.25*buffer(n,j,k,1)
                TS_AS_iplane%dvdk(nn,k) = -1./12.*buffer(n,j,k-4,2) + 0.5*buffer(n,j,k-3,2) - 1.5*buffer(n,j,k-2,2) + 5./6.*buffer(n,j,k-1,2) + 0.25*buffer(n,j,k,2)
                TS_AS_iplane%dwdk(nn,k) = -1./12.*buffer(n,j,k-4,3) + 0.5*buffer(n,j,k-3,3) - 1.5*buffer(n,j,k-2,3) + 5./6.*buffer(n,j,k-1,3) + 0.25*buffer(n,j,k,3)
              elseif(k.eq.nz) then
                TS_AS_iplane%dudk(nn,k) = 0.25*buffer(n,j,k-4,1) - 4./3.*buffer(n,j,k-3,1) + 3.*buffer(n,j,k-2,1) - 4.*buffer(n,j,k-1,1) + 25./12.*buffer(n,j,k,1)
                TS_AS_iplane%dvdk(nn,k) = 0.25*buffer(n,j,k-4,2) - 4./3.*buffer(n,j,k-3,2) + 3.*buffer(n,j,k-2,2) - 4.*buffer(n,j,k-1,2) + 25./12.*buffer(n,j,k,2)
                TS_AS_iplane%dwdk(nn,k) = 0.25*buffer(n,j,k-4,3) - 4./3.*buffer(n,j,k-3,3) + 3.*buffer(n,j,k-2,3) - 4.*buffer(n,j,k-1,3) + 25./12.*buffer(n,j,k,3)
              endif
              TS_AS_iplane%dudx(nn,k) = TS_AS_iplane%dudi(nn,k)*didx + TS_AS_iplane%dudj(nn,k)*djdx + TS_AS_iplane%dudk(nn,k)*dkdx
              TS_AS_iplane%dvdx(nn,k) = TS_AS_iplane%dvdi(nn,k)*didx + TS_AS_iplane%dvdj(nn,k)*djdx + TS_AS_iplane%dvdk(nn,k)*dkdx
              TS_AS_iplane%dwdx(nn,k) = TS_AS_iplane%dwdi(nn,k)*didx + TS_AS_iplane%dwdj(nn,k)*djdx + TS_AS_iplane%dwdk(nn,k)*dkdx

              TS_AS_iplane%dudy(nn,k) = TS_AS_iplane%dudi(nn,k)*didy + TS_AS_iplane%dudj(nn,k)*djdy + TS_AS_iplane%dudk(nn,k)*dkdy
              TS_AS_iplane%dvdy(nn,k) = TS_AS_iplane%dvdi(nn,k)*didy + TS_AS_iplane%dvdj(nn,k)*djdy + TS_AS_iplane%dvdk(nn,k)*dkdy
              TS_AS_iplane%dwdy(nn,k) = TS_AS_iplane%dwdi(nn,k)*didy + TS_AS_iplane%dwdj(nn,k)*djdy + TS_AS_iplane%dwdk(nn,k)*dkdy

              TS_AS_iplane%dudz(nn,k) = TS_AS_iplane%dudi(nn,k)*didz + TS_AS_iplane%dudj(nn,k)*djdz + TS_AS_iplane%dudk(nn,k)*dkdz
              TS_AS_iplane%dvdz(nn,k) = TS_AS_iplane%dvdi(nn,k)*didz + TS_AS_iplane%dvdj(nn,k)*djdz + TS_AS_iplane%dvdk(nn,k)*dkdz
              TS_AS_iplane%dwdz(nn,k) = TS_AS_iplane%dwdi(nn,k)*didz + TS_AS_iplane%dwdj(nn,k)*djdz + TS_AS_iplane%dwdk(nn,k)*dkdz


!            dsrc_Philipdi = ( src_Philip_iplane(n,-2,k) - 8.0*(src_Philip_iplane(n,-1,k) - src_Philip_iplane(n,1,k)) &
!                            - src_Philip_iplane(n,2,k) )*c12i
!            TS_AS_iplane%dsrc_Philipdx(n-2,k) = dsrc_Philipdi*didx
!            TS_AS_iplane%dsrc_Philipdt(n-2,k) = ( src_Philip_iplane(n-2,0,k) - 8.0*(src_Philip_iplane(n-1,0,k)) &
!                                                - src_Philip_iplane(n+1,0,k) - src_Philip_iplane(n+2,0,k) )*c12i*dti
            enddo
          enddo
        enddo

      end subroutine CalTSDeriv_AS_iplane


      subroutine CalTSDeriv_AS_kplane(nt,ny,nv,buffer,buffer_grid,dt,src_Philip_kplane)
        integer, intent(in) :: nt, ny, nv
        real(8), intent(in) :: buffer(nt,ny,-2:2,nv), buffer_grid(ny, 9), dt, src_Philip_kplane(nt,ny,-2:2)
        integer :: j, n, nn
        real(8) :: didx, djdx, dkdx, didy, djdy, dkdy, didz, djdz, dkdz
        real(8) :: dsrc_Philipdi
        real(8), parameter :: c12i = 1.d0/12.d0
        real(8) :: dti

        dti = 1.d0/dt

        TS_AS_kplane%u = reshape(buffer(3:nt-2,3:ny-2,0,1),(/TS_AS_kplane%nsample/))
        TS_AS_kplane%v = reshape(buffer(3:nt-2,3:ny-2,0,2),(/TS_AS_kplane%nsample/))
        TS_AS_kplane%w = reshape(buffer(3:nt-2,3:ny-2,0,3),(/TS_AS_kplane%nsample/))
        TS_AS_kplane%p = reshape(buffer(3:nt-2,3:ny-2,0,4),(/TS_AS_kplane%nsample/))
        TS_AS_kplane%t = reshape(buffer(3:nt-2,3:ny-2,0,5),(/TS_AS_kplane%nsample/))
!        dudn = reshape(buffer(3:nt-2,3:ny-2,1,6),(/TS_AS_kplane%nsample/))
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
               TS_AS_kplane%dudi(nn) = ( buffer(n,j,-2,1) - 8.0*(buffer(n,j,-1,1) - buffer(n,j,1,1)) - buffer(n,j,2,1) )*c12i
               TS_AS_kplane%dvdi(nn) = ( buffer(n,j,-2,2) - 8.0*(buffer(n,j,-1,2) - buffer(n,j,1,2)) - buffer(n,j,2,2) )*c12i
               TS_AS_kplane%dwdi(nn) = ( buffer(n,j,-2,3) - 8.0*(buffer(n,j,-1,3) - buffer(n,j,1,3)) - buffer(n,j,2,3) )*c12i
               TS_AS_kplane%dpdi(nn) = ( buffer(n,j,-2,4) - 8.0*(buffer(n,j,-1,4) - buffer(n,j,1,4)) - buffer(n,j,2,4) )*c12i
               TS_AS_kplane%dtdi(nn) = ( buffer(n,j,-2,5) - 8.0*(buffer(n,j,-1,5) - buffer(n,j,1,5)) - buffer(n,j,2,5) )*c12i

               TS_AS_kplane%dudj(nn) = ( buffer(n,j-2,0,1) - 8.0*(buffer(n,j-1,0,1) - buffer(n,j+1,0,1)) - buffer(n,j+2,0,1) )*c12i
               TS_AS_kplane%dvdj(nn) = ( buffer(n,j-2,0,2) - 8.0*(buffer(n,j-1,0,2) - buffer(n,j+1,0,2)) - buffer(n,j+2,0,2) )*c12i
               TS_AS_kplane%dwdj(nn) = ( buffer(n,j-2,0,3) - 8.0*(buffer(n,j-1,0,3) - buffer(n,j+1,0,3)) - buffer(n,j+2,0,3) )*c12i
               TS_AS_kplane%dpdj(nn) = ( buffer(n,j-2,0,4) - 8.0*(buffer(n,j-1,0,4) - buffer(n,j+1,0,4)) - buffer(n,j+2,0,4) )*c12i
               TS_AS_kplane%dtdj(nn) = ( buffer(n,j-2,0,5) - 8.0*(buffer(n,j-1,0,5) - buffer(n,j+1,0,5)) - buffer(n,j+2,0,5) )*c12i

               TS_AS_kplane%dudk(nn) = buffer(n,j,0,6)
               TS_AS_kplane%dvdk(nn) = buffer(n,j,0,7)
               TS_AS_kplane%dwdk(nn) = buffer(n,j,0,8)
               TS_AS_kplane%dpdk(nn) = buffer(n,j,0,9)
               TS_AS_kplane%dtdk(nn) = buffer(n,j,0,10)

               TS_AS_kplane%dudx(nn) = TS_AS_kplane%dudi(nn)*didx + TS_AS_kplane%dudj(nn)*djdx + TS_AS_kplane%dudk(nn)*dkdx
               TS_AS_kplane%dvdx(nn) = TS_AS_kplane%dvdi(nn)*didx + TS_AS_kplane%dvdj(nn)*djdx + TS_AS_kplane%dvdk(nn)*dkdx
               TS_AS_kplane%dwdx(nn) = TS_AS_kplane%dwdi(nn)*didx + TS_AS_kplane%dwdj(nn)*djdx + TS_AS_kplane%dwdk(nn)*dkdx
               TS_AS_kplane%dpdx(nn) = TS_AS_kplane%dpdi(nn)*didx + TS_AS_kplane%dpdj(nn)*djdx + TS_AS_kplane%dpdk(nn)*dkdx

               TS_AS_kplane%dudy(nn) = TS_AS_kplane%dudi(nn)*didy + TS_AS_kplane%dudj(nn)*djdy + TS_AS_kplane%dudk(nn)*dkdy
               TS_AS_kplane%dvdy(nn) = TS_AS_kplane%dvdi(nn)*didy + TS_AS_kplane%dvdj(nn)*djdy + TS_AS_kplane%dvdk(nn)*dkdy
               TS_AS_kplane%dwdy(nn) = TS_AS_kplane%dwdi(nn)*didy + TS_AS_kplane%dwdj(nn)*djdy + TS_AS_kplane%dwdk(nn)*dkdy

               TS_AS_kplane%dudz(nn) = TS_AS_kplane%dudi(nn)*didz + TS_AS_kplane%dudj(nn)*djdz + TS_AS_kplane%dudk(nn)*dkdz
               TS_AS_kplane%dvdz(nn) = TS_AS_kplane%dvdi(nn)*didz + TS_AS_kplane%dvdj(nn)*djdz + TS_AS_kplane%dvdk(nn)*dkdz
               TS_AS_kplane%dwdz(nn) = TS_AS_kplane%dwdi(nn)*didz + TS_AS_kplane%dwdj(nn)*djdz + TS_AS_kplane%dwdk(nn)*dkdz

               dsrc_Philipdi = ( src_Philip_kplane(n,j,-2) - 8.0*(src_Philip_kplane(n,j,-1) - src_Philip_kplane(n,j,1)) &
                               - src_Philip_kplane(n,j,2) )*c12i
               TS_AS_kplane%dsrc_Philipdx(nn) = dsrc_Philipdi*didx
               TS_AS_kplane%dsrc_Philipdt(nn) = ( src_Philip_kplane(n-2,j,0) - 8.0*(src_Philip_kplane(n-1,j,0)) &
                                          - src_Philip_kplane(n+1,j,0) - src_Philip_kplane(n+2,j,0) )*c12i*dti
             enddo
         enddo

      end subroutine CalTSDeriv_AS_kplane

      subroutine CalTS_AS_jplane(nz,buffer)
        integer, intent(in) :: nz
        real(8), intent(out) :: buffer(nz,nv_ASource)
        real(8) :: num_dble, rho
        real(8), dimension(nz) :: src_Philip, src_Philip_ave, src_Philip2_ave
        real(8), dimension(nz) :: dwdx2ave !!!
        real(8), dimension(nz) :: upxupxave, upxupx2ave, dudxave, dudx2ave
        real(8), dimension(nz) :: vpyvpyave, vpyvpy2ave, dvdyave, dvdy2ave
        real(8), dimension(nz) :: wpzwpzave, wpzwpz2ave, dwdzave, dwdz2ave
        real(8), dimension(nz) :: upyvpxave, upyvpx2ave, dudyave, dvdxave
        real(8), dimension(nz) :: upzwpxave, upzwpx2ave, dudzave, dwdxave
        real(8), dimension(nz) :: vpzwpyave, vpzwpy2ave, dvdzave, dwdyave
        real(8), dimension(nz) :: uave, uave2, vave, vave2, wave, wave2, tave, rhoave, pave   !!!
        real(8), dimension(nz) :: nonlinear2ave !, nonlinearave
        integer :: i, j, k

        real(8), dimension(nz) :: upyupyave, upzupzave
        real(8), dimension(nz) :: vpxvpxave, vpzvpzave
        real(8), dimension(nz) :: wpxwpxave, wpywpyave
        real(8), dimension(nz) :: upxave, vpyave, wpzave
        real(8), dimension(nz) :: divpave, divp2ave, srcave
        real(8) :: divptmp


        num_dble = dble(TS_AS_jplane%nsample)

        rhoave = 0.d0
        src_Philip_ave = 0.d0; src_Philip2_ave = 0.d0

        upxupxave = 0.d0; upxupx2ave = 0.d0
        vpyvpyave = 0.d0; vpyvpy2ave = 0.d0
        wpzwpzave = 0.d0; wpzwpz2ave = 0.d0
        upyvpxave = 0.d0; upyvpx2ave = 0.d0
        upzwpxave = 0.d0; upzwpx2ave = 0.d0
        vpzwpyave = 0.d0; vpzwpy2ave = 0.d0
        nonlinear2ave = 0.d0

        upyupyave = 0.d0; upzupzave = 0.d0
        vpxvpxave = 0.d0; vpzvpzave = 0.d0
        wpxwpxave = 0.d0; wpywpyave = 0.d0
        upxave = 0.d0; vpyave = 0.d0;  wpzave = 0.d0
        divpave = 0.d0; divp2ave = 0.d0; srcave = 0.d0

        do k=1, nz
          uave(k)  = sum(TS_AS_jplane%u(:,k))/num_dble
          vave(k)  = sum(TS_AS_jplane%v(:,k))/num_dble
          wave(k)  = sum(TS_AS_jplane%w(:,k))/num_dble
          uave2(k) = sum(TS_AS_jplane%u(:,k)**2)/num_dble
          vave2(k) = sum(TS_AS_jplane%v(:,k)**2)/num_dble  !!!  change
          wave2(k) = sum(TS_AS_jplane%w(:,k)**2)/num_dble
          pave(k) = sum(TS_AS_jplane%p(:,k))/num_dble
          tave(k) =  sum(TS_AS_jplane%t(:,k))/num_dble

          buffer(k,13) = uave(k)
          buffer(k,14) = wave(k)
          buffer(k,15) = pave(k)
          buffer(k,16) = tave(k)

          dudxave(k) = sum(TS_AS_jplane%dudx(:,k))/num_dble
          dvdyave(k) = sum(TS_AS_jplane%dvdy(:,k))/num_dble
          dwdzave(k) = sum(TS_AS_jplane%dwdz(:,k))/num_dble
          dudyave(k) = sum(TS_AS_jplane%dudy(:,k))/num_dble
          dvdxave(k) = sum(TS_AS_jplane%dvdx(:,k))/num_dble
          dudzave(k) = sum(TS_AS_jplane%dudz(:,k))/num_dble
          dwdxave(k) = sum(TS_AS_jplane%dwdx(:,k))/num_dble
          dvdzave(k) = sum(TS_AS_jplane%dvdz(:,k))/num_dble
          dwdyave(k) = sum(TS_AS_jplane%dwdy(:,k))/num_dble

          buffer(k,1) = 0.5*( abs(uave2(k) - uave(k)**2) + abs(vave2(k) - vave(k)**2) + abs(wave2(k) - wave(k)**2)  ) ! change

          buffer(k,10) = 0.d0

!          buffer(k,20) = ( sqrt(abs(uave2(k) - uave(k)**2))/sqrt(abs(dudx2ave(k) - dudxave(k)**2)) + &
!                           sqrt(abs(vave2(k) - vave(k)**2))/sqrt(abs(dvdy2ave(k) - dvdyave(k)**2)) + &
!                           sqrt(abs(wave2(k) - wave(k)**2))/sqrt(abs(dwdz2ave(k) - dwdzave(k)**2)) )/3.d0 ! Lambda
          buffer(k,21) = pave(k)/tave(k)/rbar*(sqrt(abs(uave2(k) - uave(k)**2)))*buffer(k,20)/(1.458e-6*(tave(k)**1.5)/(tave(k)+110.4d0))

          do i=1, TS_AS_jplane%nsample
            rho = TS_AS_jplane%p(i,k)/rbar/TS_AS_jplane%t(i,k)
            rhoave(k) = rhoave(k) + rho

!            src_Philip(k) = gamma*(TS_AS_jplane%dudx(i,k)**2 + TS_AS_jplane%dvdy(i,k)**2 + TS_AS_jplane%dwdz(i,k)**2 &
!                           + 2.d0*(TS_AS_jplane%dudy(i,k)*TS_AS_jplane%dvdx(i,k) + TS_AS_jplane%dudz(i,k)*TS_AS_jplane%dwdx(i,k) + &
!                                   TS_AS_jplane%dvdz(i,k)*TS_AS_jplane%dwdy(i,k) ))
            src_Philip(k) = TS_AS_jplane%dudx(i,k)**2 + TS_AS_jplane%dvdy(i,k)**2 + TS_AS_jplane%dwdz(i,k)**2 &
                           + 2.d0*(TS_AS_jplane%dudy(i,k)*TS_AS_jplane%dvdx(i,k) + TS_AS_jplane%dudz(i,k)*TS_AS_jplane%dwdx(i,k) + &
                                   TS_AS_jplane%dvdz(i,k)*TS_AS_jplane%dwdy(i,k) )

            src_Philip_ave(k) = src_Philip_ave(k) + src_Philip(k)
            src_Philip2_ave(k) = src_Philip2_ave(k) + src_Philip(k)*src_Philip(k)

            upxupxave(k)  = upxupxave(k) +  (TS_AS_jplane%dudx(i,k) - dudxave(k))**2; upxupx2ave(k) = upxupx2ave(k) + (TS_AS_jplane%dudx(i,k) - dudxave(k))**4
            vpyvpyave(k)  = vpyvpyave(k) +  (TS_AS_jplane%dvdy(i,k) - dvdyave(k))**2; vpyvpy2ave(k) = vpyvpy2ave(k) + (TS_AS_jplane%dvdy(i,k) - dvdyave(k))**4
            wpzwpzave(k)  = wpzwpzave(k) +  (TS_AS_jplane%dwdz(i,k) - dwdzave(k))**2; wpzwpz2ave(k) = wpzwpz2ave(k) + (TS_AS_jplane%dwdz(i,k) - dwdzave(k))**4
            upyvpxave(k)  = upyvpxave(k) +  (TS_AS_jplane%dudy(i,k) - dudyave(k))*(TS_AS_jplane%dvdx(i,k) - dvdxave(k))
            upyvpx2ave(k) = upyvpx2ave(k)+ ((TS_AS_jplane%dudy(i,k) - dudyave(k))*(TS_AS_jplane%dvdx(i,k) - dvdxave(k)))**2
            upzwpxave(k)  = upzwpxave(k)  + (TS_AS_jplane%dudz(i,k) - dudzave(k))*(TS_AS_jplane%dwdx(i,k) - dwdxave(k))
            upzwpx2ave(k) = upzwpx2ave(k)+ ((TS_AS_jplane%dudz(i,k) - dudzave(k))*(TS_AS_jplane%dwdx(i,k) - dwdxave(k)))**2
            vpzwpyave(k)  = vpzwpyave(k)  + (TS_AS_jplane%dvdz(i,k) - dvdzave(k))*(TS_AS_jplane%dwdy(i,k) - dwdyave(k))
            vpzwpy2ave(k) = vpzwpy2ave(k)+ ((TS_AS_jplane%dvdz(i,k) - dvdzave(k))*(TS_AS_jplane%dwdy(i,k) - dwdyave(k)))**2

            upyupyave(k) = upyupyave(k) + (TS_AS_jplane%dudy(i,k) - dudyave(k))**2; upzupzave(k) = upzupzave(k) + (TS_AS_jplane%dudz(i,k) - dudzave(k))**2
            vpxvpxave(k) = vpxvpxave(k) + (TS_AS_jplane%dvdx(i,k) - dvdxave(k))**2; vpzvpzave(k) = vpzvpzave(k) + (TS_AS_jplane%dvdz(i,k) - dvdzave(k))**2
            wpxwpxave(k) = wpxwpxave(k) + (TS_AS_jplane%dwdx(i,k) - dwdxave(k))**2; wpywpyave(k) = wpywpyave(k) + (TS_AS_jplane%dwdy(i,k) - dwdyave(k))**2
            upxave(k) = upxave(K) + TS_AS_jplane%dudx(i,k) - dudxave(k)
            vpyave(k) = vpyave(K) + TS_AS_jplane%dvdy(i,k) - dvdyave(k)
            wpzave(k) = wpzave(K) + TS_AS_jplane%dwdz(i,k) - dwdzave(k)

            buffer(k,10) = buffer(k,10) + 4.0/3.0*((TS_AS_jplane%dudx(i,k) - dudxave(k))**2 + (TS_AS_jplane%dvdy(i,k) - dvdyave(k))**2 &
                        + (TS_AS_jplane%dwdz(i,k) - dwdzave(k))**2 - (TS_AS_jplane%dudx(i,k) - dudxave(k))*(TS_AS_jplane%dvdy(i,k) - dvdyave(k)) &
                        - (TS_AS_jplane%dudx(i,k) - dudxave(k))*(TS_AS_jplane%dwdz(i,k) - dwdzave(k)) -(TS_AS_jplane%dvdy(i,k) &
                        - dvdyave(k))*(TS_AS_jplane%dwdz(i,k) - dwdzave(k)) + ((TS_AS_jplane%dudy(i,k) - dudyave(k)) &
                        + (TS_AS_jplane%dvdx(i,k) - dvdxave(k)))**2 +((TS_AS_jplane%dudz(i,k) - dudzave(k)) + (TS_AS_jplane%dwdx(i,k) - dwdxave(k)))**2 &
                        + ((TS_AS_jplane%dvdz(i,k) - dvdzave(k)) + (TS_AS_jplane%dwdy(i,k) - dwdyave(k)))**2 ) ! dis

            nonlinear2ave(k) = nonlinear2ave(k) + ( (TS_AS_jplane%dudx(i,k) - dudxave(k))**2 + &
                                     (TS_AS_jplane%dvdy(i,k) - dvdyave(k))**2 + (TS_AS_jplane%dwdz(i,k) - dwdzave(k))**2 + &
                               2.d0*((TS_AS_jplane%dudy(i,k) - dudyave(k))*(TS_AS_jplane%dvdx(i,k) - dvdxave(k)) + &
                                     (TS_AS_jplane%dudz(i,k) - dudzave(k))*(TS_AS_jplane%dwdx(i,k) - dwdxave(k)) + &
                                     (TS_AS_jplane%dvdz(i,k) - dvdzave(k))*(TS_AS_jplane%dwdy(i,k) - dwdyave(k)) ) )**2

            divptmp = (TS_AS_jplane%dudx(i,k) - dudxave(k)) + (TS_AS_jplane%dvdy(i,k) - dvdyave(k)) + (TS_AS_jplane%dwdz(i,k) - dwdzave(k))
            divpave(k)  = divpave(k)  + divptmp
            divp2ave(k) = divp2ave(k) + divptmp**2
            srcave(k) = srcave(k) + (TS_AS_jplane%dudx(i,k) - dudxave(k))**2 + (TS_AS_jplane%dvdy(i,k) - dvdyave(k))**2 + (TS_AS_jplane%dwdz(i,k) - dwdzave(k))**2 + &
                                2*( (TS_AS_jplane%dudy(i,k) - dudyave(k))*(TS_AS_jplane%dvdx(i,k) - dvdxave(k)) + &
                                    (TS_AS_jplane%dudz(i,k) - dudzave(k))*(TS_AS_jplane%dwdx(i,k) - dwdxave(k)) + &
                                    (TS_AS_jplane%dvdz(i,k) - dvdzave(k))*(TS_AS_jplane%dwdy(i,k) - dwdyave(k)) )


          enddo ! end i loop
        enddo ! end k loop
        rhoave = rhoave/num_dble
        src_Philip_ave = src_Philip_ave/num_dble
        src_Philip2_ave = src_Philip2_ave/num_dble

        upxupxave  = upxupxave/num_dble; upxupx2ave = upxupx2ave/num_dble
        vpyvpyave  = vpyvpyave/num_dble; vpyvpy2ave = vpyvpy2ave/num_dble
        wpzwpzave  = wpzwpzave/num_dble; wpzwpz2ave = wpzwpz2ave/num_dble
        upyvpxave  = upyvpxave/num_dble; upyvpx2ave = upyvpx2ave/num_dble
        upzwpxave  = upzwpxave/num_dble; upzwpx2ave = upzwpx2ave/num_dble
        vpzwpyave  = vpzwpyave/num_dble; vpzwpy2ave = vpzwpy2ave/num_dble
        nonlinear2ave = nonlinear2ave/num_dble

        upyupyave = upyupyave/num_dble; upzupzave = upzupzave/num_dble
        vpxvpxave = vpxvpxave/num_dble; vpzvpzave = vpzvpzave/num_dble
        wpxwpxave = wpxwpxave/num_dble; wpywpyave = wpywpyave/num_dble
        upxave = upxave/num_dble; vpyave = vpyave/num_dble; wpzave = wpzave/num_dble
        divpave = divpave/num_dble; divp2ave = divp2ave/num_dble; srcave = srcave/num_dble


        do k=1, nz
          buffer(k,2)  = sqrt(abs(src_Philip2_ave(k) - src_Philip_ave(k)**2))  ! src_Philip_rms
!          buffer(k,3)  = 2.d0*(sum(TS_AS_jplane%dudz(:,k))/num_dble)*(sqrt(abs(sum(TS_AS_jplane%dwdx(:,k)**2)/num_dble - &
!                          (sum(TS_AS_jplane%dwdx(:,k))/num_dble)**2)))*gamma !!! linear_Philips_rms
          buffer(k,3)  = 2.d0*(sum(TS_AS_jplane%dudz(:,k))/num_dble)*(sqrt(abs(sum(TS_AS_jplane%dwdx(:,k)**2)/num_dble - &
                          (sum(TS_AS_jplane%dwdx(:,k))/num_dble)**2))) !!! linear_Philips_rms

          buffer(k,4)  = sqrt( abs(upxupx2ave(k) - upxupxave(k)**2) )      ! upxupxrms
          buffer(k,5)  = sqrt( abs(vpyvpy2ave(k) - vpyvpyave(k)**2) )      ! vpyvpyrms
          buffer(k,6)  = sqrt( abs(wpzwpz2ave(k) - wpzwpzave(k)**2) )      ! wpzwpzrms
          buffer(k,7)  = sqrt( abs(upyvpx2ave(k) - upyvpxave(k)**2) )*2.d0 ! upyvpxrms
          buffer(k,8)  = sqrt( abs(upzwpx2ave(k) - upzwpxave(k)**2) )*2.d0 ! upzwpxrms
          buffer(k,9)  = sqrt( abs(vpzwpy2ave(k) - vpzwpyave(k)**2) )*2.d0 ! vpzwpyrms

          buffer(k,10) = 1.458e-6*sqrt(tave(k)**3)/(tave(k)+110.4)/rhoave(k)*buffer(k,10)/num_dble ! dis
          buffer(k,11) = 0.d0 !-(sum(TS_AS_jplane%dsrc_Philipdt(:,k)*TS_AS_jplane%dsrc_Philipdx(:,k))/num_dble)/(1.e-30+sum(TS_AS_jplane%dsrc_Philipdx(:,k)**2)/num_dble)

!          buffer(k,12) = gamma*sqrt( abs( (upxupx2ave(k) + vpyvpy2ave(k) + wpzwpz2ave(k) + 2.d0*(upyvpx2ave(k) + upzwpx2ave(k) + vpzwpy2ave(k))) - &
!                      (upxupxave(k) + vpyvpyave(k) + wpzwpzave(k) + 2.d0*(upyvpxave(k) + upzwpxave(k) + vpzwpyave(k)))**2 )  )
          buffer(k,12) = sqrt( abs( nonlinear2ave(k) - &
                       (upxupxave(k) + vpyvpyave(k) + wpzwpzave(k) + 2.d0*(upyvpxave(k) + upzwpxave(k) + vpzwpyave(k)))**2 )  )
          buffer(k,17) = upxupxave(k) +  upyupyave(k) + upzupzave(k) + &
                         vpxvpxave(k) +  vpyvpyave(k) + vpzvpzave(k) + &
                         wpxwpxave(k) +  wpywpyave(k) + wpzwpzave(k)
          buffer(k,18) = srcave(k)
          buffer(k,19) = sqrt(abs(divp2ave(k)))

        enddo ! end k loop

      end subroutine CalTS_AS_jplane


      subroutine CalTS_AS_iplane(nz,buffer)
        integer, intent(in) :: nz
        real(8), intent(out) :: buffer(nz,nv_ASource)
        real(8) :: num_dble, rho
        real(8), dimension(nz) :: src_Philip, src_Philip_ave, src_Philip2_ave
        real(8), dimension(nz) :: dwdx2ave !!!
        real(8), dimension(nz) :: upxupxave, upxupx2ave, dudxave, dudx2ave
        real(8), dimension(nz) :: vpyvpyave, vpyvpy2ave, dvdyave, dvdy2ave
        real(8), dimension(nz) :: wpzwpzave, wpzwpz2ave, dwdzave, dwdz2ave
        real(8), dimension(nz) :: upyvpxave, upyvpx2ave, dudyave, dvdxave
        real(8), dimension(nz) :: upzwpxave, upzwpx2ave, dudzave, dwdxave
        real(8), dimension(nz) :: vpzwpyave, vpzwpy2ave, dvdzave, dwdyave
        real(8), dimension(nz) :: uave, uave2, vave, vave2, wave, wave2, tave, rhoave, pave   !!!
        real(8), dimension(nz) :: nonlinear2ave !, nonlinearave
        integer :: i, j, k

        real(8), dimension(nz) :: upyupyave, upzupzave
        real(8), dimension(nz) :: vpxvpxave, vpzvpzave
        real(8), dimension(nz) :: wpxwpxave, wpywpyave
        real(8), dimension(nz) :: upxave, vpyave, wpzave
        real(8), dimension(nz) :: divpave, divp2ave, srcave
        real(8) :: divptmp

        num_dble = dble(TS_AS_iplane%nsample)

        rhoave = 0.d0
        src_Philip_ave = 0.d0; src_Philip2_ave = 0.d0

        upxupxave = 0.d0; upxupx2ave = 0.d0
        vpyvpyave = 0.d0; vpyvpy2ave = 0.d0
        wpzwpzave = 0.d0; wpzwpz2ave = 0.d0
        upyvpxave = 0.d0; upyvpx2ave = 0.d0
        upzwpxave = 0.d0; upzwpx2ave = 0.d0
        vpzwpyave = 0.d0; vpzwpy2ave = 0.d0
        nonlinear2ave = 0.d0

        upyupyave = 0.d0; upzupzave = 0.d0
        vpxvpxave = 0.d0; vpzvpzave = 0.d0
        wpxwpxave = 0.d0; wpywpyave = 0.d0
        upxave = 0.d0; vpyave = 0.d0;  wpzave = 0.d0
        divpave = 0.d0; divp2ave = 0.d0; srcave = 0.d0

        do k=1, nz
          uave(k)  = sum(TS_AS_iplane%u(:,k))/num_dble
          vave(k)  = sum(TS_AS_iplane%v(:,k))/num_dble
          wave(k)  = sum(TS_AS_iplane%w(:,k))/num_dble
          uave2(k) = sum(TS_AS_iplane%u(:,k)**2)/num_dble
          vave2(k) = sum(TS_AS_iplane%v(:,k)**2)/num_dble
          wave2(k) = sum(TS_AS_iplane%w(:,k)**2)/num_dble
          pave(k) = sum(TS_AS_iplane%p(:,k))/num_dble
          tave(k) =  sum(TS_AS_iplane%t(:,k))/num_dble

          buffer(k,13) = uave(k)
          buffer(k,14) = wave(k)
          buffer(k,15) = pave(k)
          buffer(k,16) = tave(k)

          dudxave(k) = sum(TS_AS_iplane%dudx(:,k))/num_dble
          dudx2ave(k) = sum(TS_AS_iplane%dudx(:,k)**2)/num_dble
          dvdyave(k) = sum(TS_AS_iplane%dvdy(:,k))/num_dble
          dvdy2ave(k) = sum(TS_AS_iplane%dvdy(:,k)**2)/num_dble
          dwdzave(k) = sum(TS_AS_iplane%dwdz(:,k))/num_dble
          dwdz2ave(k) = sum(TS_AS_iplane%dwdz(:,k)**2)/num_dble
          dudyave(k) = sum(TS_AS_iplane%dudy(:,k))/num_dble
          dvdxave(k) = sum(TS_AS_iplane%dvdx(:,k))/num_dble
          dudzave(k) = sum(TS_AS_iplane%dudz(:,k))/num_dble
          dwdxave(k) = sum(TS_AS_iplane%dwdx(:,k))/num_dble
          dvdzave(k) = sum(TS_AS_iplane%dvdz(:,k))/num_dble
          dwdyave(k) = sum(TS_AS_iplane%dwdy(:,k))/num_dble

          buffer(k,1) = 0.5*( abs(uave2(k) - uave(k)**2) + abs(vave2(k) - vave(k)**2) + abs(wave2(k) - wave(k)**2)  )

          buffer(k,10) = 0.d0

          buffer(k,20) = ( sqrt(abs(uave2(k) - uave(k)**2))/sqrt(abs(dudx2ave(k) - dudxave(k)**2)) + &
                           sqrt(abs(vave2(k) - vave(k)**2))/sqrt(abs(dvdy2ave(k) - dvdyave(k)**2)) + &
                           sqrt(abs(wave2(k) - wave(k)**2))/sqrt(abs(dwdz2ave(k) - dwdzave(k)**2)) )/3.d0 ! Lambda
          buffer(k,21) = pave(k)/tave(k)/rbar*(sqrt(abs(uave2(k) - uave(k)**2)))*buffer(k,20)/(1.458e-6*(tave(k)**1.5)/(tave(k)+110.4d0))

          do i=1, TS_AS_iplane%nsample
            rho = TS_AS_iplane%p(i,k)/rbar/TS_AS_iplane%t(i,k)
            rhoave(k) = rhoave(k) + rho

            !src_Philip(k) = gamma*(TS_AS_iplane%dudx(i,k)**2 + TS_AS_iplane%dvdy(i,k)**2 + TS_AS_iplane%dwdz(i,k)**2 &
            !               + 2.d0*(TS_AS_iplane%dudy(i,k)*TS_AS_iplane%dvdx(i,k) + TS_AS_iplane%dudz(i,k)*TS_AS_iplane%dwdx(i,k) + &
            !                       TS_AS_iplane%dvdz(i,k)*TS_AS_iplane%dwdy(i,k) ))
            src_Philip(k) = TS_AS_iplane%dudx(i,k)**2 + TS_AS_iplane%dvdy(i,k)**2 + TS_AS_iplane%dwdz(i,k)**2 &
                           + 2.d0*(TS_AS_iplane%dudy(i,k)*TS_AS_iplane%dvdx(i,k) + TS_AS_iplane%dudz(i,k)*TS_AS_iplane%dwdx(i,k) + &
                                   TS_AS_iplane%dvdz(i,k)*TS_AS_iplane%dwdy(i,k) )
            src_Philip_ave(k) = src_Philip_ave(k) + src_Philip(k)
            src_Philip2_ave(k) = src_Philip2_ave(k) + src_Philip(k)*src_Philip(k)

            upxupxave(k)  = upxupxave(k) +  (TS_AS_iplane%dudx(i,k) - dudxave(k))**2; upxupx2ave(k) = upxupx2ave(k) + (TS_AS_iplane%dudx(i,k) - dudxave(k))**4
            vpyvpyave(k)  = vpyvpyave(k) +  (TS_AS_iplane%dvdy(i,k) - dvdyave(k))**2; vpyvpy2ave(k) = vpyvpy2ave(k) + (TS_AS_iplane%dvdy(i,k) - dvdyave(k))**4
            wpzwpzave(k)  = wpzwpzave(k) +  (TS_AS_iplane%dwdz(i,k) - dwdzave(k))**2; wpzwpz2ave(k) = wpzwpz2ave(k) + (TS_AS_iplane%dwdz(i,k) - dwdzave(k))**4
            upyvpxave(k)  = upyvpxave(k) +  (TS_AS_iplane%dudy(i,k) - dudyave(k))*(TS_AS_iplane%dvdx(i,k) - dvdxave(k))
            upyvpx2ave(k) = upyvpx2ave(k)+ ((TS_AS_iplane%dudy(i,k) - dudyave(k))*(TS_AS_iplane%dvdx(i,k) - dvdxave(k)))**2
            upzwpxave(k)  = upzwpxave(k)  + (TS_AS_iplane%dudz(i,k) - dudzave(k))*(TS_AS_iplane%dwdx(i,k) - dwdxave(k))
            upzwpx2ave(k) = upzwpx2ave(k)+ ((TS_AS_iplane%dudz(i,k) - dudzave(k))*(TS_AS_iplane%dwdx(i,k) - dwdxave(k)))**2
            vpzwpyave(k)  = vpzwpyave(k)  + (TS_AS_iplane%dvdz(i,k) - dvdzave(k))*(TS_AS_iplane%dwdy(i,k) - dwdyave(k))
            vpzwpy2ave(k) = vpzwpy2ave(k)+ ((TS_AS_iplane%dvdz(i,k) - dvdzave(k))*(TS_AS_iplane%dwdy(i,k) - dwdyave(k)))**2

            upyupyave(k) = upyupyave(k) + (TS_AS_iplane%dudy(i,k) - dudyave(k))**2; upzupzave(k) = upzupzave(k) + (TS_AS_iplane%dudz(i,k) - dudzave(k))**2
            vpxvpxave(k) = vpxvpxave(k) + (TS_AS_iplane%dvdx(i,k) - dvdxave(k))**2; vpzvpzave(k) = vpzvpzave(k) + (TS_AS_iplane%dvdz(i,k) - dvdzave(k))**2
            wpxwpxave(k) = wpxwpxave(k) + (TS_AS_iplane%dwdx(i,k) - dwdxave(k))**2; wpywpyave(k) = wpywpyave(k) + (TS_AS_iplane%dwdy(i,k) - dwdyave(k))**2
            upxave(k) = upxave(K) + TS_AS_iplane%dudx(i,k) - dudxave(k)
            vpyave(k) = vpyave(K) + TS_AS_iplane%dvdy(i,k) - dvdyave(k)
            wpzave(k) = wpzave(K) + TS_AS_iplane%dwdz(i,k) - dwdzave(k)

            buffer(k,10) = buffer(k,10) + 4.0/3.0*((TS_AS_iplane%dudx(i,k) - dudxave(k))**2 + (TS_AS_iplane%dvdy(i,k) - dvdyave(k))**2 &
                        + (TS_AS_iplane%dwdz(i,k) - dwdzave(k))**2 - (TS_AS_iplane%dudx(i,k) - dudxave(k))*(TS_AS_iplane%dvdy(i,k) - dvdyave(k)) &
                        - (TS_AS_iplane%dudx(i,k) - dudxave(k))*(TS_AS_iplane%dwdz(i,k) - dwdzave(k)) -(TS_AS_iplane%dvdy(i,k) &
                        - dvdyave(k))*(TS_AS_iplane%dwdz(i,k) - dwdzave(k)) + ((TS_AS_iplane%dudy(i,k) - dudyave(k)) &
                        + (TS_AS_iplane%dvdx(i,k) - dvdxave(k)))**2 +((TS_AS_iplane%dudz(i,k) - dudzave(k)) + (TS_AS_iplane%dwdx(i,k) - dwdxave(k)))**2 &
                        + ((TS_AS_iplane%dvdz(i,k) - dvdzave(k)) + (TS_AS_iplane%dwdy(i,k) - dwdyave(k)))**2 ) ! dis

            nonlinear2ave(k) = nonlinear2ave(k) + ( (TS_AS_iplane%dudx(i,k) - dudxave(k))**2 + &
                                     (TS_AS_iplane%dvdy(i,k) - dvdyave(k))**2 + (TS_AS_iplane%dwdz(i,k) - dwdzave(k))**2 + &
                               2.d0*((TS_AS_iplane%dudy(i,k) - dudyave(k))*(TS_AS_iplane%dvdx(i,k) - dvdxave(k)) + &
                                     (TS_AS_iplane%dudz(i,k) - dudzave(k))*(TS_AS_iplane%dwdx(i,k) - dwdxave(k)) + &
                                     (TS_AS_iplane%dvdz(i,k) - dvdzave(k))*(TS_AS_iplane%dwdy(i,k) - dwdyave(k)) ) )**2

            divptmp = (TS_AS_iplane%dudx(i,k) - dudxave(k)) + (TS_AS_iplane%dvdy(i,k) - dvdyave(k)) + (TS_AS_iplane%dwdz(i,k) - dwdzave(k))
            divpave(k)  = divpave(k)  + divptmp
            divp2ave(k) = divp2ave(k) + divptmp**2
            srcave(k) = srcave(k) + (TS_AS_iplane%dudx(i,k) - dudxave(k))**2 + (TS_AS_iplane%dvdy(i,k) - dvdyave(k))**2 + (TS_AS_iplane%dwdz(i,k) - dwdzave(k))**2 + &
                                    2*( (TS_AS_iplane%dudy(i,k) - dudyave(k))*(TS_AS_iplane%dvdx(i,k) - dvdxave(k)) + &
                                        (TS_AS_iplane%dudz(i,k) - dudzave(k))*(TS_AS_iplane%dwdx(i,k) - dwdxave(k)) + &
                                        (TS_AS_iplane%dvdz(i,k) - dvdzave(k))*(TS_AS_iplane%dwdy(i,k) - dwdyave(k)) )

          enddo ! end i loop
        enddo ! end k loop
        rhoave = rhoave/num_dble
        src_Philip_ave = src_Philip_ave/num_dble
        src_Philip2_ave = src_Philip2_ave/num_dble

        upxupxave  = upxupxave/num_dble; upxupx2ave = upxupx2ave/num_dble
        vpyvpyave  = vpyvpyave/num_dble; vpyvpy2ave = vpyvpy2ave/num_dble
        wpzwpzave  = wpzwpzave/num_dble; wpzwpz2ave = wpzwpz2ave/num_dble
        upyvpxave  = upyvpxave/num_dble; upyvpx2ave = upyvpx2ave/num_dble
        upzwpxave  = upzwpxave/num_dble; upzwpx2ave = upzwpx2ave/num_dble
        vpzwpyave  = vpzwpyave/num_dble; vpzwpy2ave = vpzwpy2ave/num_dble
        nonlinear2ave = nonlinear2ave/num_dble

        upyupyave = upyupyave/num_dble; upzupzave = upzupzave/num_dble
        vpxvpxave = vpxvpxave/num_dble; vpzvpzave = vpzvpzave/num_dble
        wpxwpxave = wpxwpxave/num_dble; wpywpyave = wpywpyave/num_dble
        upxave = upxave/num_dble; vpyave = vpyave/num_dble; wpzave = wpzave/num_dble
        divpave = divpave/num_dble; divp2ave = divp2ave/num_dble; srcave = srcave/num_dble

        do k=1, nz
          buffer(k,2)  = sqrt(abs(src_Philip2_ave(k) - src_Philip_ave(k)**2))  ! src_Philip_rms
          !buffer(k,3)  = 2.d0*(sum(TS_AS_iplane%dudz(:,k))/num_dble)*(sqrt(abs(sum(TS_AS_iplane%dwdx(:,k)**2)/num_dble - &
          !                (sum(TS_AS_iplane%dwdx(:,k))/num_dble)**2)))*gamma !!! linear_Philips_rms
          buffer(k,3)  = 2.d0*(sum(TS_AS_iplane%dudz(:,k))/num_dble)*(sqrt(abs(sum(TS_AS_iplane%dwdx(:,k)**2)/num_dble - &
                          (sum(TS_AS_iplane%dwdx(:,k))/num_dble)**2))) !!! linear_Philips_rms

          buffer(k,4)  = sqrt( abs(upxupx2ave(k) - upxupxave(k)**2) )       ! upxupxrms
          buffer(k,5)  = sqrt( abs(vpyvpy2ave(k) - vpyvpyave(k)**2) )       ! vpyvpyrms
          buffer(k,6)  = sqrt( abs(wpzwpz2ave(k) - wpzwpzave(k)**2) )       ! wpzwpzrms
          buffer(k,7)  = sqrt( abs(upyvpx2ave(k) - upyvpxave(k)**2) )*2.d0  ! upyvpxrms
          buffer(k,8)  = sqrt( abs(upzwpx2ave(k) - upzwpxave(k)**2) )*2.d0  ! upzwpxrms
          buffer(k,9)  = sqrt( abs(vpzwpy2ave(k) - vpzwpyave(k)**2) )*2.d0  ! vpzwpyrms

          buffer(k,10) = 1.458e-6*sqrt(tave(k)**3)/(tave(k)+110.4)/rhoave(k)*buffer(k,10)/num_dble ! dis
          buffer(k,11) = 0.d0

          !buffer(k,11) = -(sum(TS_AS_iplane%dsrc_Philipdt(:,k)*TS_AS_iplane%dsrc_Philipdx(:,k))/num_dble)/(1.e-30+sum(TS_AS_iplane%dsrc_Philipdx(:,k)**2)/num_dble)

          !buffer(k,12) = gamma*sqrt( abs( (upxupx2ave(k) + vpyvpy2ave(k) + wpzwpz2ave(k) + 2.d0*(upyvpx2ave(k) + upzwpx2ave(k) + vpzwpy2ave(k))) - &
          !            (upxupxave(k) + vpyvpyave(k) + wpzwpzave(k) + 2.d0*(upyvpxave(k) + upzwpxave(k) + vpzwpyave(k)))**2 )  )
          buffer(k,12) =  sqrt( abs( nonlinear2ave(k) - &
                       (upxupxave(k) + vpyvpyave(k) + wpzwpzave(k) + 2.d0*(upyvpxave(k) + upzwpxave(k) + vpzwpyave(k)))**2 )  )

          buffer(k,17) = upxupxave(k) +  upyupyave(k) + upzupzave(k) + &
                         vpxvpxave(k) +  vpyvpyave(k) + vpzvpzave(k) + &
                         wpxwpxave(k) +  wpywpyave(k) + wpzwpzave(k)
          !buffer(k,18) = upxupxave(k) + vpyvpyave(k) + wpzwpzave(k) + 2.d0*(upyvpxave(k) + upzwpxave(k) + vpzwpyave(k))
          !buffer(k,19) = upxave(k) + vpyave(k) + wpzave(k)
          buffer(k,18) = srcave(k)
          buffer(k,19) = sqrt(abs(divp2ave(k)))

!          buffer(k,20) = sqrt(abs(dudx2ave(k) - dudxave(k)**2 ) )  ! dudx_rms
!          buffer(k,21) = sqrt(abs(dvdy2ave(k) - dvdyave(k)**2 ) )  ! dvdy_rms
!          buffer(k,22) = sqrt(abs(dwdz2ave(k) - dwdzave(k)**2 ) )  ! dwdz_rms

        enddo ! end k loop

      end subroutine CalTS_AS_iplane




      subroutine CalTS_AS_kplane(buffer)
        real(8), intent(out) :: buffer(nv_ASource)
        real(8) :: num_dble, rho
        real(8) :: src_Philip, src_Philip_ave, src_Philip2_ave
        real(8) :: dwdx2ave !!!
        real(8) :: upxupxave, upxupx2ave, dudxave
        real(8) :: vpyvpyave, vpyvpy2ave, dvdyave
        real(8) :: wpzwpzave, wpzwpz2ave, dwdzave
        real(8) :: upyvpxave, upyvpx2ave, dudyave, dvdxave
        real(8) :: upzwpxave, upzwpx2ave, dudzave, dwdxave
        real(8) :: vpzwpyave, vpzwpy2ave, dvdzave, dwdyave
        real(8) :: uave, uave2, vave, vave2, wave, wave2, tave, rhoave, pave   !!!
        integer :: i

        num_dble = dble(TS_AS_kplane%nsample)
        uave   = sum(TS_AS_kplane%u)/num_dble     ! uave
        vave   = sum(TS_AS_kplane%v)/num_dble     ! vave
        wave   = sum(TS_AS_kplane%w)/num_dble     ! wave
        uave2  = sum(TS_AS_kplane%u**2)/num_dble  ! u2
        vave2  = sum(TS_AS_kplane%v**2)/num_dble  ! v2
        wave2  = sum(TS_AS_kplane%w**2)/num_dble  ! w2
        pave   = sum(TS_AS_kplane%p)/num_dble     ! pave
        tave   = sum(TS_AS_kplane%t)/num_dble     ! tave

        buffer(13) = uave
        buffer(14) = wave
        buffer(15) = pave
        buffer(16) = tave

        rhoave = 0.d0
        src_Philip_ave = 0.d0; src_Philip2_ave = 0.d0

        upxupxave = 0.d0; upxupx2ave = 0.d0
        vpyvpyave = 0.d0; vpyvpy2ave = 0.d0
        wpzwpzave = 0.d0; wpzwpz2ave = 0.d0
        upyvpxave = 0.d0; upyvpx2ave = 0.d0
        upzwpxave = 0.d0; upzwpx2ave = 0.d0
        vpzwpyave = 0.d0; vpzwpy2ave = 0.d0

        dudxave = sum(TS_AS_kplane%dudx)/num_dble
        dvdyave = sum(TS_AS_kplane%dvdy)/num_dble
        dwdzave = sum(TS_AS_kplane%dwdz)/num_dble
        dudyave = sum(TS_AS_kplane%dudy)/num_dble
        dvdxave = sum(TS_AS_kplane%dvdx)/num_dble
        dudzave = sum(TS_AS_kplane%dudz)/num_dble
        dwdxave = sum(TS_AS_kplane%dwdx)/num_dble
        dvdzave = sum(TS_AS_kplane%dvdz)/num_dble
        dwdyave = sum(TS_AS_kplane%dwdy)/num_dble

        buffer(1) = 0.5*( abs(uave2-uave**2) + abs(vave2-vave**2) + abs(wave2-wave**2) )

        buffer(10) = 0.d0 ! dis

        do i=1, TS_AS_kplane%nsample

          rho = TS_AS_kplane%p(i)/rbar/TS_AS_kplane%t(i)
          rhoave = rhoave + rho

          src_Philip = gamma*(TS_AS_kplane%dudx(i)**2 + TS_AS_kplane%dvdy(i)**2 + TS_AS_kplane%dwdz(i)**2 &
                      + 2.d0*(TS_AS_kplane%dudy(i)*TS_AS_kplane%dvdx(i) + TS_AS_kplane%dudz(i)*TS_AS_kplane%dwdx(i) + &
                              TS_AS_kplane%dvdz(i)*TS_AS_kplane%dwdy(i) ))
          src_Philip_ave = src_Philip_ave + src_Philip
          src_Philip2_ave = src_Philip2_ave + src_Philip*src_Philip

          upxupxave = upxupxave + (TS_AS_kplane%dudx(i) - dudxave)**2; upxupx2ave = upxupx2ave + (TS_AS_kplane%dudx(i) - dudxave)**4
          vpyvpyave = vpyvpyave + (TS_AS_kplane%dvdy(i) - dvdyave)**2; vpyvpy2ave = vpyvpy2ave + (TS_AS_kplane%dvdy(i) - dvdyave)**4
          wpzwpzave = wpzwpzave + (TS_AS_kplane%dwdz(i) - dwdzave)**2; wpzwpz2ave = wpzwpz2ave + (TS_AS_kplane%dwdz(i) - dwdzave)**4
          upyvpxave = upyvpxave + (TS_AS_kplane%dudy(i) - dudyave)*(TS_AS_kplane%dvdx(i) - dvdxave)
          upyvpx2ave = upyvpx2ave + ((TS_AS_kplane%dudy(i) - dudyave)*(TS_AS_kplane%dvdx(i) - dvdxave))**2
          upzwpxave = upzwpxave + (TS_AS_kplane%dudz(i) - dudzave)*(TS_AS_kplane%dwdx(i) - dwdxave)
          upzwpx2ave = upzwpx2ave + ((TS_AS_kplane%dudz(i) - dudzave)*(TS_AS_kplane%dwdx(i) - dwdxave))**2
          vpzwpyave = vpzwpyave + (TS_AS_kplane%dvdz(i) - dvdzave)*(TS_AS_kplane%dwdy(i) - dwdyave)
          vpzwpy2ave = vpzwpy2ave + ((TS_AS_kplane%dvdz(i) - dvdzave)*(TS_AS_kplane%dwdy(i) - dwdyave))**2

          buffer(10) = buffer(10) + 4.0/3.0*((TS_AS_kplane%dudx(i) - dudxave)**2 + (TS_AS_kplane%dvdy(i) - dvdyave)**2 &
                      + (TS_AS_kplane%dwdz(i) - dwdzave)**2 - (TS_AS_kplane%dudx(i) - dudxave)*(TS_AS_kplane%dvdy(i) - dvdyave) &
                      - (TS_AS_kplane%dudx(i) - dudxave)*(TS_AS_kplane%dwdz(i) - dwdzave) -(TS_AS_kplane%dvdy(i) &
                      - dvdyave)*(TS_AS_kplane%dwdz(i) - dwdzave) + ((TS_AS_kplane%dudy(i) - dudyave) &
                      + (TS_AS_kplane%dvdx(i) - dvdxave))**2 +((TS_AS_kplane%dudz(i) - dudzave) + (TS_AS_kplane%dwdx(i) - dwdxave))**2 &
                      + ((TS_AS_kplane%dvdz(i) - dvdzave) + (TS_AS_kplane%dwdy(i) - dwdyave))**2 ) ! dis

        enddo

        rhoave = rhoave/num_dble
        src_Philip_ave = src_Philip_ave/num_dble
        src_Philip2_ave = src_Philip2_ave/num_dble
        buffer(2)  = sqrt(abs(src_Philip2_ave - src_Philip_ave**2))  ! src_Philip_rms
        buffer(3)  = 2.d0*(sum(TS_AS_kplane%dudz)/num_dble)*(sqrt(abs(sum(TS_AS_kplane%dwdx**2)/num_dble - &
                          (sum(TS_AS_kplane%dwdx)/num_dble)**2)))*gamma !!! linear_Philips_rms
        upxupxave  = upxupxave/num_dble; upxupx2ave = upxupx2ave/num_dble
        vpyvpyave  = vpyvpyave/num_dble; vpyvpy2ave = vpyvpy2ave/num_dble
        wpzwpzave  = wpzwpzave/num_dble; wpzwpz2ave = wpzwpz2ave/num_dble
        upyvpxave  = upyvpxave/num_dble; upyvpx2ave = upyvpx2ave/num_dble
        upzwpxave  = upzwpxave/num_dble; upzwpx2ave = upzwpx2ave/num_dble
        vpzwpyave  = vpzwpyave/num_dble; vpzwpy2ave = vpzwpy2ave/num_dble
        buffer(4)  = sqrt( abs(upxupx2ave - upxupxave**2) )      ! upxupxrms
        buffer(5)  = sqrt( abs(vpyvpy2ave - vpyvpyave**2) )      ! vpyvpyrms
        buffer(6)  = sqrt( abs(wpzwpz2ave - wpzwpzave**2) )      ! wpzwpzrms
        buffer(7)  = sqrt( abs(upyvpx2ave - upyvpxave**2) )*2.d0 ! upyvpxrms
        buffer(8)  = sqrt( abs(upzwpx2ave - upzwpxave**2) )*2.d0 ! upzwpxrms
        buffer(9)  = sqrt( abs(vpzwpy2ave - vpzwpyave**2) )*2.d0 ! vpzwpyrms

        buffer(10) = 1.458e-6*sqrt(tave**3)/(tave+110.4)/rhoave*buffer(10)/num_dble ! dis
        buffer(11) = -(sum(TS_AS_kplane%dsrc_Philipdt*TS_AS_kplane%dsrc_Philipdx)/num_dble)/(1.e-30+sum(TS_AS_kplane%dsrc_Philipdx**2)/num_dble)

        buffer(12) = gamma*sqrt( abs( (upxupx2ave + vpyvpy2ave + wpzwpz2ave + 2.d0*(upyvpx2ave + upzwpx2ave + vpzwpy2ave)) - &
                      (upxupxave + vpyvpyave + wpzwpzave + 2.d0*(upyvpxave + upzwpxave + vpzwpyave))**2 )  )

      end subroutine CalTS_AS_kplane

!parameter( dname_AS_kplane = (/'kinetic', 'src_Philip_rms', 'linear_Philip_rms', 'upxupxrms', 'vpyvpyrms', 'wpzwpzrms', &
!                                     'upyvpxrms', 'upzwpxrms', 'vpzwpyrms', 'dis', 'U_src', 'x', 'z'/) )


      subroutine WriteTSASource_jplane(fn,jloc,nx,nz,buffer)
            character(*), intent(in) :: fn
      integer, intent(in) :: jloc, nx, nz
      real(8), intent(in) :: buffer(nz,nx,nv_ASource+2)
      character(4) :: fnum

       write(unit=fnum,fmt='(I04.4)') jloc

       if(trim(TSAsource_jout%fname).eq.'none') then
         TSAsource_jout%IsMultiGroup = .false.
       elseif(trim(TSAsource_jout%fname).eq.trim(fn)) then
         TSAsource_jout%IsMultiGroup = .true.
       else
         TSAsource_jout%IsMultiGroup = .false.
       endif

       TSAsource_jout%fname = trim(fn)
       TSAsource_jout%gname = "/"//fnum//'ASource'
       TSAsource_jout%rank = 2
       TSAsource_jout%dnum = nv_ASource+2

       if(allocated(TSAsource_jout%dname)) then
         deallocate(TSAsource_jout%dname,TSAsource_jout%dimsf,TSAsource_jout%dimsm)
         deallocate(TSAsource_jout%block,TSAsource_jout%count,TSAsource_jout%stride)
         deallocate(TSAsource_jout%offset)
       endif

       allocate(TSAsource_jout%dname(TSAsource_jout%dnum))
       allocate(TSAsource_jout%dimsf(TSAsource_jout%rank))
       allocate(TSAsource_jout%dimsm(TSAsource_jout%rank))
       allocate(TSAsource_jout%block(TSAsource_jout%rank))
       allocate(TSAsource_jout%count(TSAsource_jout%rank))
       allocate(TSAsource_jout%stride(TSAsource_jout%rank))
       allocate(TSAsource_jout%offset(TSAsource_jout%rank))
       TSAsource_jout%dname = dname_AS_jplane
       TSAsource_jout%dimsf(1) = nz
       TSAsource_jout%dimsf(2) = nx
       TSAsource_jout%dimsm  = TSAsource_jout%dimsf
       TSAsource_jout%block  = TSAsource_jout%dimsf
       TSAsource_jout%count  = 1
       TSAsource_jout%stride = 1
       TSAsource_jout%offset = 0
       TSAsource_jout%IsHSInitialized = .true.
!       TSAsource_kout%IsMultiGroup = .true.
       call WriteTSHDF5_2D(TSAsource_jout, buffer)


      end subroutine WriteTSASource_jplane


      subroutine WriteTSASource_iplane(fn,iloc,nz,buffer)
            character(*), intent(in) :: fn
      integer, intent(in) :: iloc, nz
      real(8), intent(in) :: buffer(1,nz,nv_ASource+2)
      character(4) :: fnum

       write(unit=fnum,fmt='(I04.4)') iloc

       if(trim(TSAsource_iout%fname).eq.'none') then
         TSAsource_iout%IsMultiGroup = .false.
       elseif(trim(TSAsource_iout%fname).eq.trim(fn)) then
         TSAsource_iout%IsMultiGroup = .true.
       else
         TSAsource_iout%IsMultiGroup = .false.
       endif

       TSAsource_iout%fname = trim(fn)
       TSAsource_iout%gname = "/"//fnum//'ASource'
       TSAsource_iout%rank = 2
       TSAsource_iout%dnum = nv_ASource+2

       if(allocated(TSAsource_iout%dname)) then
         deallocate(TSAsource_iout%dname,TSAsource_iout%dimsf,TSAsource_iout%dimsm)
         deallocate(TSAsource_iout%block,TSAsource_iout%count,TSAsource_iout%stride)
         deallocate(TSAsource_iout%offset)
       endif

       allocate(TSAsource_iout%dname(TSAsource_iout%dnum))
       allocate(TSAsource_iout%dimsf(TSAsource_iout%rank))
       allocate(TSAsource_iout%dimsm(TSAsource_iout%rank))
       allocate(TSAsource_iout%block(TSAsource_iout%rank))
       allocate(TSAsource_iout%count(TSAsource_iout%rank))
       allocate(TSAsource_iout%stride(TSAsource_iout%rank))
       allocate(TSAsource_iout%offset(TSAsource_iout%rank))
       TSAsource_iout%dname = dname_AS_iplane
       TSAsource_iout%dimsf(1) = 1
       TSAsource_iout%dimsf(2) = nz
       TSAsource_iout%dimsm  = TSAsource_iout%dimsf
       TSAsource_iout%block  = TSAsource_iout%dimsf
       TSAsource_iout%count  = 1
       TSAsource_iout%stride = 1
       TSAsource_iout%offset = 0
       TSAsource_iout%IsHSInitialized = .true.
!       TSAsource_kout%IsMultiGroup = .true.
       call WriteTSHDF5_2D(TSAsource_iout, buffer)


      end subroutine WriteTSASource_iplane



      subroutine WriteTSASource_kplane(fn,kloc,nx,buffer)
            character(*), intent(in) :: fn
      integer, intent(in) :: kloc,nx
      real(8), intent(in) :: buffer(nx,nv_ASource+2)
      character(4) :: fnum

       write(unit=fnum,fmt='(I04.4)') kloc

       if(trim(TSAsource_kout%fname).eq.'none') then
         TSAsource_kout%IsMultiGroup = .false.
       elseif(trim(TSAsource_kout%fname).eq.trim(fn)) then
         TSAsource_kout%IsMultiGroup = .true.
       else
         TSAsource_kout%IsMultiGroup = .false.
       endif

       TSAsource_kout%fname = trim(fn)
       TSAsource_kout%gname = "/"//fnum//'ASource'
       TSAsource_kout%rank = 1
       TSAsource_kout%dnum = nv_ASource+2

       if(allocated(TSAsource_kout%dname)) then
         deallocate(TSAsource_kout%dname,TSAsource_kout%dimsf,TSAsource_kout%dimsm)
         deallocate(TSAsource_kout%block,TSAsource_kout%count,TSAsource_kout%stride)
         deallocate(TSAsource_kout%offset)
       endif

       allocate(TSAsource_kout%dname(TSAsource_kout%dnum))
       allocate(TSAsource_kout%dimsf(TSAsource_kout%rank))
       allocate(TSAsource_kout%dimsm(TSAsource_kout%rank))
       allocate(TSAsource_kout%block(TSAsource_kout%rank))
       allocate(TSAsource_kout%count(TSAsource_kout%rank))
       allocate(TSAsource_kout%stride(TSAsource_kout%rank))
       allocate(TSAsource_kout%offset(TSAsource_kout%rank))
       TSAsource_kout%dname = dname_AS_kplane
       TSAsource_kout%dimsf(1) = nx
       TSAsource_kout%dimsm  = TSAsource_kout%dimsf
       TSAsource_kout%block  = TSAsource_kout%dimsf
       TSAsource_kout%count  = 1
       TSAsource_kout%stride = 1
       TSAsource_kout%offset = 0
       TSAsource_kout%IsHSInitialized = .true.
!       TSAsource_kout%IsMultiGroup = .true.
       call WriteTSHDF5_1D(TSAsource_kout, buffer)


      end subroutine WriteTSASource_kplane


    end module MTSAcousticSource
