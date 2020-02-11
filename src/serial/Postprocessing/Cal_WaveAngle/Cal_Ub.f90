

program Cal_Ub

  use modRWHDF5
  implicit none

  integer :: file_be, file_end, file_skip, num_file
  character(400) :: fname
  character(400) :: filepath, filename, gridname
  real(8) :: dpdx
  real(8) :: dpdt
  real(8), dimension(:), allocatable :: dpdx2ave, dpdtdpdxave, uave, pave, tave, wave  !!!!
  real(8), dimension(:), allocatable :: upave, wpave, uwave, urms, wrms, prms, u2ave, w2ave, p2ave
  real(8), dimension(:), allocatable :: uconv_p
  real(8), dimension(:,:), allocatable :: urms_2d, wrms_2d, prms_2d
  real(8) :: dt_sample, didx_ip
  integer :: nt, nt_total !!!!!!!!!!!!!
  integer :: num_iplane, jbe_i, jend_i, jsp_i, kbe_i, kend_i, ksp_i
  integer :: ibe, iend, iskip, kbe, kend, kskip
  integer :: ilen, jlen, klen
  real(8), dimension(:,:,:,:,:), allocatable :: buffer, buffer_tmp
  real(8), dimension(:,:,:,:), allocatable :: buffer_grid
  real(8), dimension(:,:,:), allocatable :: buffer_Acoustic, buffer_grid_Acoustic, angle_2d
  character(10) :: gname
  integer :: num_dset
  character(10), dimension(:), allocatable :: dname
  character(8) :: fnum_8
  integer :: i, j, k, n, m
  real(8) :: tmp(-2:2)
  integer :: dims(1)
  integer :: numpt, iFormat
  real(8), dimension(:), allocatable :: theta_output
  real(8), dimension(:,:), allocatable :: theta_output_2d, theta_output_2d_2

  type(tp_rdwt_hdf5) :: fsol, grd

  call Input()
  call InitHDF5()

if(iFormat.eq.0) then

  num_file = (file_end-file_be)/file_skip + 1

  print *, 'Number of Files = ',  num_file
  print *, 'File begin number: ', file_be
  print *, 'File end number: ',   file_end
  print *, 'File skip number: ',  file_skip
  print *, '*********************************************'

  jlen = (jend_i-jbe_i)/jsp_i + 1
  klen = (kend_i-kbe_i)/ksp_i + 1
  nt_total = nt*num_file

  num_dset = 5

  allocate(buffer_tmp(nt,jlen,klen,num_iplane,num_dset))    !!!
  allocate(buffer(nt_total,jlen,klen,num_iplane,num_dset))  !!!!
  allocate(dname(num_dset))

  fsol%gname = "/iplane"
  fsol%dnum = num_dset
  fsol%rank = 4
  allocate(fsol%dname(fsol%dnum),fsol%dimsf(fsol%rank))
  fsol%dname(1) = "p"
  fsol%dname(2) = "pi"
  fsol%dname(3) = "u"
  fsol%dname(4) = "w"
  fsol%dname(5) = "T"
  fsol%IsHSInitialized = .true.
  fsol%dimsf = (/nt,jlen,klen,1/) !!!!

  m = 0
  do n=1, num_file
    write(unit=fnum_8,fmt='(I08.8)') (file_be + (n-1)*file_skip)
    fsol%fname = trim(filepath)//'timeseries_'//fnum_8//'.h5'
    call ReadHDF5_4D(fsol,buffer_tmp)
    ! call Readtsflow_partial(fname,gname,nt,jlen,klen,num_iplane,num_dset,dname,buffer_tmp)!!!
    m = m + 1
    buffer(((m-1)*nt+1):m*nt,:,:,:,:) = buffer_tmp(1:nt,:,:,:,:)
  enddo

  allocate(dpdtdpdxave(klen))
  allocate(dpdx2ave(klen))
  allocate(uave(klen))
  allocate(uconv_p(klen))
  allocate(buffer_grid(jlen,klen,1,12)) !!!!
  allocate(pave(klen))
  allocate(tave(klen))
  allocate(wave(klen))
  allocate(upave(klen), wpave(klen), uwave(klen), urms(klen),wrms(klen), prms(klen) )
  allocate(u2ave(klen), w2ave(klen), p2ave(klen))
  !  pave, tave, M
  !  upave, wpave, uwave, urms, wrms, prms
  do k=1, klen
    numpt = 0
    dpdtdpdxave(k) = 0.
    dpdx2ave(k) = 0.
    uave(k) = 0.

    pave(k) = 0.; tave(k) = 0.;  wave(k) = 0.
    upave(k) = 0.; wpave(k) = 0.; uwave(k) = 0.
    u2ave(k) = 0.; w2ave(k) = 0.; p2ave(k) = 0.

    do n=3, nt_total-2
      do j=3, jlen-2
        tmp(-2:2) = buffer(n-2:n+2,j,k,1,1)
        dpdt = ( tmp(-2)-8.d0*(tmp(-1)-tmp(1))-tmp(2) )/(12.d0*dt_sample)
        dpdx = buffer(n,j,k,1,2)*didx_ip

        uave(k) = uave(k) + buffer(n,j,k,1,3)
        dpdtdpdxave(k) = dpdtdpdxave(k) + dpdt*dpdx
        dpdx2ave(k) = dpdx2ave(k) + dpdx**2

        pave(k) = pave(k) + buffer(n,j,k,1,1)
        tave(k) = tave(k) + buffer(n,j,k,1,5)
        wave(k) = wave(k) + buffer(n,j,k,1,4)
        upave(k) = upave(k) + buffer(n,j,k,1,3)*buffer(n,j,k,1,1)
        wpave(k) = wpave(k) + buffer(n,j,k,1,4)*buffer(n,j,k,1,1)
        uwave(k) = uwave(k) + buffer(n,j,k,1,3)*buffer(n,j,k,1,4)
        u2ave(k) = u2ave(k) + buffer(n,j,k,1,3)*buffer(n,j,k,1,3)
        w2ave(k) = w2ave(k) + buffer(n,j,k,1,4)*buffer(n,j,k,1,4)
        p2ave(k) = p2ave(k) + buffer(n,j,k,1,1)*buffer(n,j,k,1,1)

        numpt = numpt + 1
      enddo
    enddo

    dpdtdpdxave(k) = dpdtdpdxave(k)/numpt
    dpdx2ave(k) = dpdx2ave(k)/numpt
    uave(k) = uave(k)/numpt
    uconv_p(k) = -dpdtdpdxave(k)/dpdx2ave(k)

    pave(k) = pave(k)/numpt
    tave(k) = tave(k)/numpt
    wave(k) = wave(k)/numpt
    upave(k) = upave(k)/numpt
    wpave(k) = wpave(k)/numpt
    uwave(k) = uwave(k)/numpt
    u2ave(k) = u2ave(k)/numpt
    w2ave(k) = w2ave(k)/numpt
    p2ave(k) = p2ave(k)/numpt

    urms(k) = sqrt(abs(u2ave(k)-uave(k)**2))
    wrms(k) = sqrt(abs(w2ave(k)-wave(k)**2))
    prms(k) = sqrt(abs(p2ave(k)-pave(k)**2))

  enddo

  fname = trim(filepath)//'timeseries_GridMetrics.h5'
  gname = "/iplane"
  call ReadtsGridMetrics(fname,gname,jlen,klen,num_iplane,buffer_grid)

  dims(1) = klen
  call Write1dHDF5("ave.h5",dims,1,"z",buffer_grid(1,1:klen,1,3),0)
  call Write1dHDF5("ave.h5",dims,1,"u",uave,1)
  call Write1dHDF5("ave.h5",dims,1,"uconv_p",uconv_p,1)
  call Write1dHDF5("ave.h5",dims,1,"p",pave,1)
  call Write1dHDF5("ave.h5",dims,1,"w",wave,1)
  call Write1dHDF5("ave.h5",dims,1,"up",upave,1)
  call Write1dHDF5("ave.h5",dims,1,"wp",wpave,1)
  call Write1dHDF5("ave.h5",dims,1,"urms",urms,1)
  call Write1dHDF5("ave.h5",dims,1,"wrms",wrms,1)
  call Write1dHDF5("ave.h5",dims,1,"prms",prms,1)
  call Write1dHDF5("ave.h5",dims,1,"uw",uwave,1)
  call Write1dHDF5("ave.h5",dims,1,"t",tave,1)

  allocate(theta_output(klen))
  theta_output = 0.d0
  do k=2, klen
    theta_output(k) = Cal_theta(uave(k),pave(k),uave(k)/(sqrt(1.4*287.d0*tave(k))),upave(k)-uave(k)*pave(k), &
                                wpave(k)-wave(k)*pave(k),urms(k),wrms(k),prms(k),uwave(k)-uave(k)*wave(k),tave(k) )
  enddo ! end k loop
  call Write1dHDF5("ave.h5",dims,1,'theta1',theta_output,1)
  do k=2, klen
    theta_output(k) = Cal_theta2(uave(k),pave(k),uave(k)/(sqrt(1.4*287.d0*tave(k))),upave(k)-uave(k)*pave(k), &
                                wpave(k)-wave(k)*pave(k),urms(k),wrms(k),prms(k),uwave(k)-uave(k)*wave(k),tave(k) )
  enddo ! end k loop
  call Write1dHDF5("ave.h5",dims,1,'theta2',theta_output,1)

else ! iFormat.eq.1

  num_dset =  10  !!!!!!!!!!!!!!!!!!!!!!
  ilen = (iend - ibe)/iskip + 1
  klen = (kend - kbe)/kskip + 1
  allocate(buffer_Acoustic(klen,ilen,num_dset))
  allocate(urms_2d(klen,ilen),wrms_2d(klen,ilen),prms_2d(klen,ilen))

  fsol%gname = 'Stat2d'
  fsol%dnum  = num_dset
  fsol%rank  = 2
  allocate(fsol%dname(fsol%dnum),fsol%dimsf(fsol%rank))

  fsol%dname(1) = 'uave'
  fsol%dname(2) = 'wave'
  fsol%dname(3) = 'pave'
  fsol%dname(4) = 'tave'
  fsol%dname(5) = 'u2'
  fsol%dname(6) = 'w2'
  fsol%dname(7) = 'p2'
  fsol%dname(8) = 'up'
  fsol%dname(9) = 'uw'
  fsol%dname(10) = 'wp'

  fsol%IsHSInitialized = .true.
  fsol%dimsf = (/klen,ilen/)

  fsol%fname = trim(filename)
  print *, 'Reading file: ', trim(fsol%fname)
  call ReadHDF5_2D(fsol,buffer_Acoustic)
  print *, 'finish reading file'
  do k=1, klen
    do i=1, ilen
      urms_2d(k,i) = sqrt(abs(buffer_Acoustic(k,i,5) - buffer_Acoustic(k,i,1)**2 ))
      wrms_2d(k,i) = sqrt(abs(buffer_Acoustic(k,i,6) - buffer_Acoustic(k,i,2)**2 ))
      prms_2d(k,i) = sqrt(abs(buffer_Acoustic(k,i,7) - buffer_Acoustic(k,i,3)**2 ))
    enddo
  enddo

  allocate(theta_output_2d(klen,ilen),theta_output_2d_2(klen,ilen))
  theta_output_2d = 0.d0

  do k=2, klen
    do i=1, ilen
!      print *, 'k = ', k, 'i = ', i
      theta_output_2d(k,i) = Cal_theta( buffer_Acoustic(k,i,1), buffer_Acoustic(k,i,3), buffer_Acoustic(k,i,1)/(sqrt(1.4*287.d0*buffer_Acoustic(k,i,4))), &
                                        buffer_Acoustic(k,i,8) - buffer_Acoustic(k,i,1)*buffer_Acoustic(k,i,3), &
                                        buffer_Acoustic(k,i,10) - buffer_Acoustic(k,i,2)*buffer_Acoustic(k,i,3), &
                                        urms_2d(k,i), wrms_2d(k,i), prms_2d(k,i), &
                                        buffer_Acoustic(k,i,9) - buffer_Acoustic(k,i,1)*buffer_Acoustic(k,i,2), buffer_Acoustic(k,i,4)  )
    enddo
  enddo
  do k=2, klen
    do i=1, ilen
      theta_output_2d_2(k,i) = Cal_theta2( buffer_Acoustic(k,i,1), buffer_Acoustic(k,i,3), buffer_Acoustic(k,i,1)/(sqrt(1.4*287.d0*buffer_Acoustic(k,i,4))), &
                                        buffer_Acoustic(k,i,8) - buffer_Acoustic(k,i,1)*buffer_Acoustic(k,i,3), &
                                        buffer_Acoustic(k,i,10) - buffer_Acoustic(k,i,2)*buffer_Acoustic(k,i,3), &
                                        urms_2d(k,i), wrms_2d(k,i), prms_2d(k,i), &
                                        buffer_Acoustic(k,i,9) - buffer_Acoustic(k,i,1)*buffer_Acoustic(k,i,2), buffer_Acoustic(k,i,4)  )
    enddo
  enddo

  print *, 'reading grid file'
  call ReadGrid()
  print *, 'finish reading grid file'

  deallocate(fsol%dname)
  fsol%gname = 'Angle'
  fsol%dnum = 4
  fsol%rank = 2
  allocate(fsol%dname(fsol%dnum))
  fsol%dname(1) = 'x'
  fsol%dname(2) = 'z'
  fsol%dname(3) = 'angle1'
  fsol%dname(4) = 'angle2'

  fsol%IsHSInitialized = .true.
  fsol%dimsf = (/klen,ilen/)
  fsol%fname = 'wave-angle.h5'

  allocate(angle_2d(klen,ilen,4))
  do k=1, klen
    do i=1,ilen
      angle_2d(k,i,1) = buffer_grid_Acoustic(k,i,1)
      angle_2d(k,i,2) = buffer_grid_Acoustic(k,i,2)
      angle_2d(k,i,3) = theta_output_2d(k,i)
      angle_2d(k,i,4) = theta_output_2d_2(k,i)
    enddo
  enddo

  print *, 'writing file: ', trim(fsol%fname)
  call WriteHDF5_2D(fsol,angle_2d)



endif



  call FinalizeHDF5()


contains

  subroutine ReadGrid()
    implicit none

    grd%fname = trim(gridname)
    grd%gname = '/'
    call DetectHDF5(grd)

    allocate(buffer_grid_Acoustic(grd%dimsf(1),grd%dimsf(2),2))
    grd%dnum = 2
    grd%dname(1) = 'x'
    grd%dname(2) = 'z'

    print *, 'Reading grid file: ', trim(gridname)
    call ReadHDF5_2D(grd, buffer_grid_Acoustic)

  end subroutine ReadGrid

  subroutine Input()
    implicit none

    read(*,*)
    read(*,*) iFormat

    if(iFormat.eq.0) then

      read(*,*)
      read(*,*)
      read(*,'(a)') filepath
      read(*,*)
      read(*,*) file_be, file_end, file_skip
      read(*,*)
      read(*,*) dt_sample, didx_ip
      read(*,*)
      read(*,*) num_iplane, nt, jbe_i, jend_i, jsp_i, kbe_i, kend_i, ksp_i

    elseif(iFormat.eq.1) then

      do i=1, 11
        read(*,*)
      enddo
      read(*,'(a)') filename
      read(*,*)
      read(*,'(a)') gridname
      read(*,*)
      read(*,*) ibe, iend, iskip, kbe, kend, kskip

    endif



  end subroutine Input


  real function Cal_theta(U, p, M, up, wp, urms, wrms, prms, uw, T)
    implicit none

    real(8), intent(in) :: U, p, M, up, wp, urms, wrms, prms, uw, T
    integer :: i, j, num_i
    integer :: max_i, min_i
    real(8), dimension(:), allocatable :: theta, theta_rad
    real(8), dimension(:), allocatable :: unp, unrms, C_unp
    real(8), dimension(:), allocatable :: u_other
    real(8) :: gamma = 1.4

    num_i = 2001
    allocate(theta(num_i),theta_rad(num_i))
    allocate(unp(num_i),unrms(num_i),C_unp(num_i))
    allocate(u_other(num_i))

    do i=1, num_i
      theta(i) = (i-1)*180.d0/(num_i-1)
      theta_rad(i) = theta(i)*4.d0*atan(1.d0)/180.d0
      unp(i) = up*cos(theta_rad(i)) + wp*sin(theta_rad(i))
      unrms(i) = sqrt(abs( (urms*cos(theta_rad(i)))**2 + &
                 (wrms*sin(theta_rad(i)))**2 + 2.d0*uw*sin(theta_rad(i))*cos(theta_rad(i)) ))
      C_unp(i) = unp(i)/(prms*unrms(i)+1.e-30)
    enddo

    max_i = maxloc(C_unp(:),1)
!    print *, 'max_i = ', max_i
!    print *, 'theta(max_i) = ', theta(max_i)

    Cal_theta = theta(max_i)

  end function Cal_theta


  real function Cal_theta2(U, p, M, up, wp, urms, wrms, prms, uw, T)
    implicit none

    real(8), intent(in) :: U, p, M, up, wp, urms, wrms, prms, uw, T
    integer :: i, j, num_i
    integer :: max_i, min_i
    real(8), dimension(:), allocatable :: theta, theta_rad
    real(8), dimension(:), allocatable :: unp, unrms, C_unp
    real(8), dimension(:), allocatable :: u_other
    real(8) :: gamma = 1.4

    num_i = 2001
    allocate(theta(num_i),theta_rad(num_i))
    allocate(unp(num_i),unrms(num_i),C_unp(num_i))
    allocate(u_other(num_i))

    do i=1, num_i
      theta(i) = (i-1)*180.d0/(num_i-1)
      theta_rad(i) = theta(i)*4.d0*atan(1.d0)/180.d0
      unp(i) = up*cos(theta_rad(i)) + wp*sin(theta_rad(i))
      unrms(i) = sqrt(abs( (urms*cos(theta_rad(i)))**2 + &
                 (wrms*sin(theta_rad(i)))**2 + 2.d0*uw*sin(theta_rad(i))*cos(theta_rad(i)) ))
      C_unp(i) = unp(i)/(prms*unrms(i)+1.e-30)
    enddo

    u_other(:) = abs( unrms(:)/U - (prms/p)/(gamma*M) )

    min_i = minloc(u_other,1)
!    print *, 'min_i = ', min_i
!    print *, 'theta(min_i) = ', theta(min_i)
    Cal_theta2 = theta(min_i)

  end function Cal_theta2













end program Cal_Ub





