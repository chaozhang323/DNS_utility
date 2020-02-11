

    module MTSAngle
      use MTSHDF5
      implicit none

      integer, parameter :: nv_Angle = 5
      real(8), parameter :: rbar = 287.d0, gamma = 1.4, cp = 1004.5 ! (cp = 3.5*rbar)
      integer, parameter :: num_i = 181 ! number of angle

      type tp_Angle
        integer :: nsample
        real(8), dimension(:,:), allocatable :: u, w, p

      end type tp_Angle

      type(tp_Angle), private :: TS_AN_iplane
      character(50), private :: dname_Angle_iplane(nv_Angle+1)
      parameter( dname_Angle_iplane = (/'uave', 'wave', 'pave', 'angle_ave', 'angle_ave2', 'z'/) )

      type(tp_hyperslab), private :: TSAsource_iout

      real(8), dimension(:), allocatable :: theta, theta_rad

    contains


      subroutine InitTS_AN_iplane(ns, kdim)
        integer, intent(in) :: ns, kdim
        integer :: i

        TSAsource_iout%fname = 'none'

        TS_AN_iplane%nsample = ns
        allocate(TS_AN_iplane%u(ns,kdim))
        allocate(TS_AN_iplane%w(ns,kdim))
        allocate(TS_AN_iplane%p(ns,kdim))

        allocate(theta(num_i),theta_rad(num_i))

        do i=1, num_i
          theta(i) = (i-1)*180.d0/(num_i-1)
          theta_rad(i) = theta(i)*4.d0*atan(1.d0)/180.d0
        enddo

      end subroutine InitTS_AN_iplane

      subroutine CalTS_AN_iplane(nz,bufferin,buffer)
        integer, intent(in) :: nz
        real(8), intent(in) :: bufferin(:,:,:,:)
        real(8), intent(out) :: buffer(nz,nv_Angle)
        real(8) :: num_dble, rho, M, pinf
        real(8), dimension(nz) :: uave, uave2, wave, wave2, pave, angleave, angleave2   !!!
        integer :: i, j, k, nt, ny

        pinf = rhoinf*rbar*Tinf
        nt = size(bufferin,dim=1)
        ny = size(bufferin,dim=2)

print *, 'nt = ', nt, 'ny = ', ny
print *, 'Uinf = ', uinf, 'rhoinf = ', rhoinf, 'Tinf = ', Tinf, 'pinf = ', pinf

        num_dble = dble(TS_AN_iplane%nsample)
        M = Uinf/sqrt(gamma*rbar*Tinf)

        do k=1, nz
          TS_AN_iplane%u(:,k) = reshape(bufferin(1:nt,1:ny,k,1),(/TS_AN_iplane%nsample/))
          TS_AN_iplane%w(:,k) = reshape(bufferin(1:nt,1:ny,k,2),(/TS_AN_iplane%nsample/))
          TS_AN_iplane%p(:,k) = reshape(bufferin(1:nt,1:ny,k,3),(/TS_AN_iplane%nsample/))
        enddo

        angleave = 0.d0; angleave2 = 0.d0
        do k=1, nz
          uave(k)  = sum(TS_AN_iplane%u(:,k))/num_dble
          wave(k)  = sum(TS_AN_iplane%w(:,k))/num_dble
          pave(k)  = sum(TS_AN_iplane%p(:,k))/num_dble

          buffer(k,1) = uave(k)
          buffer(k,2) = wave(k)
          buffer(k,3) = pave(k)

          do i=1, TS_AN_iplane%nsample
            TS_AN_iplane%u(i,k)  = TS_AN_iplane%u(i,k)  - uave(k)
            TS_AN_iplane%w(i,k)  = TS_AN_iplane%w(i,k)  - wave(k)
            TS_AN_iplane%p(i,k)  = TS_AN_iplane%p(i,k)  - pave(k)
            angleave(k)  = angleave(k)  +   Cal_theta(Uinf, pinf, M, TS_AN_iplane%u(i,k), TS_AN_iplane%w(i,k), TS_AN_iplane%p(i,k))
            angleave2(k) = angleave2(k) + ( Cal_theta(Uinf, pinf, M, TS_AN_iplane%u(i,k), TS_AN_iplane%w(i,k), TS_AN_iplane%p(i,k)) )**2
          enddo ! end i loop
          buffer(k,4) = angleave(k)/num_dble
          buffer(k,5) = angleave2(k)/num_dble
        enddo ! end k loop


      end subroutine CalTS_AN_iplane

      real function Cal_theta(Uinf, pinf, M, up, wp, pp)
        real(8), intent(in) :: Uinf, pinf, M, up, wp, pp
        integer :: i, j !, num_i
        !integer :: max_i, min_i
        !real(8), dimension(:), allocatable :: theta, theta_rad
        !real(8), dimension(:), allocatable :: unp
        !real(8), dimension(:), allocatable :: u_other
        real(8) :: A
        real(8) :: tmp1, tmp2

        A = Uinf/(gamma*M*pinf)

        tmp1 = ( 2.d0*A*pp*wp + sqrt( (2.d0*pp*wp)**2 - 2.d0*(up**2+wp**2)*(A*A*pp*pp-up**2) ) )/(2.d0*(up**2+wp**2)+1.e-30)
        tmp2 = ( 2.d0*A*pp*wp - sqrt( (2.d0*pp*wp)**2 - 2.d0*(up**2+wp**2)*(A*A*pp*pp-up**2) ) )/(2.d0*(up**2+wp**2)+1.e-30)

        if((tmp1.gt.0.and.tmp2.gt.0).or.(tmp1.lt.0.and.tmp2.lt.0)) then
          print *, 'error in find angle'
          print *, 'tmp1 = ', tmp1 , 'tmp2 = ', tmp2
          stop
        endif

        if(tmp1.gt.0) then
          print *, 'tmp1 = ', tmp1
          Cal_theta = asin(tmp1)
        else
          Cal_theta = asin(tmp2)
        endif
        print *, 'Cal_theta = ', Cal_theta
!        !num_i = 2001
!        !allocate(theta(num_i),theta_rad(num_i))
!        allocate(unp(num_i))
!        allocate(u_other(num_i))
!
!        do i=1, num_i
!          !theta(i) = (i-1)*180.d0/(num_i-1)
!          !theta_rad(i) = theta(i)*4.d0*atan(1.d0)/180.d0
!          unp(i) = up*cos(theta_rad(i)) + wp*sin(theta_rad(i))
!        enddo
!
!        u_other(:) = abs( unp(:)/Uinf - (pp/pinf)/(gamma*M) )
!
!        min_i = minloc(u_other,1)
!        ! print *, 'min_i = ', min_i
!        ! print *, 'theta(min_i) = ', theta(min_i)
!        Cal_theta = theta(min_i)

      end function Cal_theta


      subroutine WriteTSAngle_iplane(fn,iloc,nz,buffer)
            character(*), intent(in) :: fn
      integer, intent(in) :: iloc, nz
      real(8), intent(in) :: buffer(1,nz,nv_Angle+1)
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
       TSAsource_iout%gname = "/"//fnum//'Angle'
       TSAsource_iout%rank = 2
       TSAsource_iout%dnum = nv_Angle+1

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
       TSAsource_iout%dname = dname_Angle_iplane
       TSAsource_iout%dimsf(1) = 1
       TSAsource_iout%dimsf(2) = nz
       TSAsource_iout%dimsm  = TSAsource_iout%dimsf
       TSAsource_iout%block  = TSAsource_iout%dimsf
       TSAsource_iout%count  = 1
       TSAsource_iout%stride = 1
       TSAsource_iout%offset = 0
       TSAsource_iout%IsHSInitialized = .true.
       call WriteTSHDF5_2D(TSAsource_iout, buffer)


      end subroutine WriteTSAngle_iplane


    end module MTSAngle
