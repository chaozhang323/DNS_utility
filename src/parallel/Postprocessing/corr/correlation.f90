program correlation
  use decomp_2d
  use MFileIO
  use MPRWHDF5
  use decomp2d_fftw
  use MFFTWindow
  implicit none

!   include 'mpif.h'
  real(8), parameter :: R=8314.3D0

  character(400) :: fname, datapath, gridfilename, groupname
  character(8) :: fnum, fnum3, fnum4
  character(4) :: fnum1, fnum2, fnum5, fnum6
  integer :: imax, jmax, kmax
  integer :: iil, iih, stride, itimeave
  real(8) :: rbar, Rm

  integer :: icalcorrxz, icalcorrxy, icalcorryz, icalcorr3D, icalWNSD_ky, icalWNSD_kx
  type tp_corr
     integer :: ibe, iend, jbe, jend, kbe, kend
     integer :: num_kref
     integer, dimension(:), allocatable :: kref
     integer :: iwinl, iwinr
     integer :: p_row, p_col ! processor grid
     type(DECOMP_INFO) :: decomp
     integer :: iinterp_x, imax_new, iinterp_y, jmax_new, iwindow
  end type tp_corr
  type(tp_corr) :: corrxz, corrxy, corryz, corr3D, WNSD_ky, WNSD_kx

  integer, parameter :: nvarcorr_xz = 1, nvarcorr_xy = 3, nvarave = 1, nvarcorr_yz = 3, nvarcorr_3D = 3
  integer, dimension(:), allocatable :: varindex
  integer :: num_output, varidx
  character(10) :: dname(14)
  !>                   1   2   3   4   5    6     7   8     9   10   11   12   13   14
  parameter(dname = (/'u','v','w','p','T','rho','P0','T0','ru','rv','rw','uv','uw','vw'/))
  real(8), dimension(:,:,:,:), allocatable :: buffer_flow, buffer_corr, buffer_rms1, buffer_flow_ref, buffer_cal, buffer_cal_ref
  real(8), dimension(:,:,:), allocatable :: buffer_ave1, buffer_ave2, avetmp1, avetmp2
  real(8), dimension(:,:,:), allocatable :: buffer_var1, buffer_var2, vartmp1, vartmp2, buffer_rms2
  real(8), dimension(:,:,:,:), allocatable :: corrtmp
  real(8), dimension(:,:,:), allocatable :: xx, yy, zz, xx_tmp
  real(8), dimension(:,:,:), allocatable :: buffer_2D

  integer :: ierr, errcode
  integer :: myid, numprocs, kloc
  integer :: kmax_rd, imax_rd, jmax_rd, n, i, ii, j, k, kk, m, nn
  integer :: kshift, ishift, jshift
  integer :: icorrxz, kcorrxz, iave_st, iave_end, iref, nsamples_pervol, numvol
  logical, dimension(3) :: periodic_bc
  integer :: kstart, klen
  integer :: icorrxy, jcorrxy, kcorrxy
  integer :: ist_rd, iend_rd, ilen_rd
  integer :: iWNSD_ky, jWNSD_ky, kWNSD_ky, iWNSD_kx, jWNSD_kx, kWNSD_kx

  integer :: icorryz, jcorryz, kcorryz, jref
  complex, dimension(:,:,:), allocatable :: specttmp1, spect1dave, spect1dtmp
  real(8), dimension(:,:,:), allocatable :: corrtmp1

  complex, dimension(:,:,:), allocatable :: spect11tmp, spect22tmp
  real(8), dimension(:,:,:), allocatable :: corr11tmp, corr22tmp, corr11, corr22
  real(8) :: xlen, ddx_min, ddx_max, pi

  integer :: icorr3D, jcorr3D, kcorr3D
  real(8), dimension(:), allocatable :: ddx, ddy

  integer :: dim1_line, dim2_line
  integer :: coords1(2), coords2(2)
  integer :: sumtmp

  ! initialize MPI
  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
  if(myid.eq.0) print *, 'Started reading parameters'
  call Input()
  if(myid.eq.0) print *, 'Finished reading parameters'
  call InitHDF5()


  if(icalcorr3D.gt.0) then
    ! do automatic domain decomposition using 2DECOMP&FFT

    kmax_rd = corr3D%kend - corr3D%kbe + 1
    imax_rd = corr3D%iend - corr3D%ibe + 1
    jmax_rd = corr3D%jend - corr3D%jbe + 1
    periodic_bc = (/.false.,.false.,.true./)
    call decomp_2d_init(kmax_rd, imax_rd, jmax_rd, corr3D%p_row, corr3D%p_col)
    call decomp_info_init(kmax_rd, imax_rd, jmax_rd, corr3D%decomp)
    call fftw_init_zpencil(kmax_rd, imax_rd, jmax_rd) ! PHYSICAL_IN_Z = 3
    call MPI_Cart_Sub(DECOMP_2D_COMM_CART_Z,(/.false.,.true./),dim2_line,ierr)  !!!!

    ! shift relative in x & y directions to DNS volume (for reading desired data within the full DNS volume)
    kshift = corr3D%kbe - 1
    ishift = corr3D%ibe - 1
    jshift = corr3D%jbe - 1

    ! Indexes ranges for reading DNS data
    ist_rd  = ishift+corr3D%decomp%zst(2)-corr3D%iwinl
    iend_rd = ishift+corr3D%decomp%zen(2)+corr3D%iwinr
    ilen_rd = iend_rd - ist_rd + 1

    ! Local length of k array size
    klen = corrxy%decomp%zsz(1)
    ! Sanity Check of domain decomposition
!    if( .not.( (klen.eq.1).or.(klen.eq.kmax_rd/numprocs) ) ) then
!        if(myid.eq.0) print *, 'local number of k planes klen =', klen
!        if(myid.eq.0) print *, 'Domain decomposition error, klen must be equal to 1 or kmax_rd/numprocs'
!        errcode = 227
!        call MPI_Abort(MPI_comm_world,errcode,ierr)
!    endif

    allocate( xx(kshift+corr3D%decomp%zst(1):kshift+corr3D%decomp%zen(1), &
                 ist_rd:iend_rd, 1:jmax), &
              yy(kshift+corr3D%decomp%zst(1):kshift+corr3D%decomp%zen(1), &
                 ist_rd:iend_rd, 1:jmax), &
              zz(kshift+corr3D%decomp%zst(1):kshift+corr3D%decomp%zen(1), &
                 ist_rd:iend_rd, 1:jmax) )

    fname = trim(gridfilename)
    if(myid.eq.0)  print *, 'Reading grid: ', trim(fname)
    !   kstart = corrxy%kref(k)
    call ReadHDF5grid_P(trim(fname),corr3D%decomp%zsz(1),       ilen_rd, jmax, &
                                    kshift+corr3D%decomp%zst(1), ist_rd, 1,    &
                                    xx, yy, zz)
    allocate( buffer_flow(     kshift+corr3D%decomp%zst(1):kshift+corr3D%decomp%zen(1),  &
                               ist_rd:iend_rd, jmax, 5) )
    allocate( buffer_flow_ref( kshift+corr3D%decomp%zst(1):kshift+corr3D%decomp%zen(1),  &
                               ist_rd:iend_rd, jmax, 5) )
    allocate( buffer_cal(      kshift+corr3D%decomp%zst(1):kshift+corr3D%decomp%zen(1),  &
                               ist_rd:iend_rd, jmax, 1) )
    allocate( buffer_cal_ref(  kshift+corr3D%decomp%zst(1):kshift+corr3D%decomp%zen(1),  &
                               ist_rd:iend_rd, jmax, 1) )
    ! Size of correlation array per core
    icorr3D = corr3D%iwinl + corr3D%iwinr + 1
    jcorr3D = jmax
    kcorr3D = corr3D%decomp%zsz(1)

    iave_st  = ishift+corr3D%decomp%zst(2)
    iave_end = ishift+corr3D%decomp%zen(2)
    iref = (iave_st + iave_end)/2
    !nsamples_pervol = (iave_end-iave_st+1)*jmax  ! # of samples for averaging in a single DNS volume
    nsamples_pervol = imax_rd*jmax  ! # of samples for averaging in a single DNS volume
    numvol = (iih-iil)/stride+1  ! Total # DNS volumes for averaging
    if(myid.eq.0) print *, 'iave_st =', iave_st, 'iave_end =', iave_end, 'icorr3D =', icorr3D

    allocate(     corrtmp(-corr3D%iwinl:corr3D%iwinr, 1:jmax, nvarcorr_3D, corr3D%decomp%zst(1):corr3D%decomp%zen(1)) )   !!!!!
    allocate( buffer_corr(-corr3D%iwinl:corr3D%iwinr, 1:jmax, nvarcorr_3D, corr3D%decomp%zst(1):corr3D%decomp%zen(1)) )

    allocate(    specttmp1(corr3D%decomp%zst(1):corr3D%decomp%zen(1), iave_st:iave_end, 1:jmax/2+1) )
    allocate(     corrtmp1(corr3D%decomp%zst(1):corr3D%decomp%zen(1), iave_st:iave_end, 1:jmax) )

    allocate(    spect11tmp(corr3D%decomp%zst(1):corr3D%decomp%zen(1), iave_st:iave_end, 1:jmax/2+1) )
    allocate(    spect22tmp(corr3D%decomp%zst(1):corr3D%decomp%zen(1), iave_st:iave_end, 1:jmax/2+1) )
    allocate(     corr11tmp(corr3D%decomp%zst(1):corr3D%decomp%zen(1), iave_st:iave_end, 1:jmax) )
    allocate(     corr22tmp(corr3D%decomp%zst(1):corr3D%decomp%zen(1), iave_st:iave_end, 1:jmax) )

    do nn=1, num_output
      varidx = varindex(nn)
      if(myid.eq.0) then
        print *, '###############################################################'
        print *, 'Calculating corr3D correlation for variable: ', trim(dname(varidx))
        print *, '###############################################################'
      endif

      do kk=1, corr3D%num_kref
        corrtmp = 0.d0 ! correlation
        buffer_corr = 0.d0
        do n=iil, iih, stride
          write(unit=fnum,fmt='(I08.8)') n
          fname = trim(datapath)//fnum//'.h5'
          if(myid.eq.0)  print *, 'Reading file: ', trim(fname)
          ! Start to read required data from DNS volume
          kstart = corr3D%kref(kk)
          call ReadHDF5sol_P( trim(fname),trim(groupname), corr3D%decomp%zsz(1), ilen_rd, jmax, &
                                           kshift+corr3D%decomp%zst(1), ist_rd, 1, &
                                           buffer_flow(:,:,:,1:5) )
          if(myid.eq.0)  print *, 'Reading reference kplane data: ', trim(fname)
          call ReadHDF5sol_P( trim(fname),trim(groupname), 1, ilen_rd, jmax, &
                                           kstart, ist_rd, 1, &
                     buffer_flow_ref(kshift+corr3D%decomp%zst(1):kshift+corr3D%decomp%zst(1),:,:,1:5) )
          do k = kshift+corr3D%decomp%zst(1)+1, kshift+corr3D%decomp%zen(1)
            buffer_flow_ref(k,:,:,1:5) = buffer_flow_ref(kshift+corr3D%decomp%zst(1),:,:,1:5)
          enddo

          if(varidx.eq.6) then
            buffer_cal(:,:,:,1)     = buffer_flow(:,:,:,4)/(rbar*buffer_flow(:,:,:,5))          ! rho
            buffer_cal_ref(:,:,:,1) = buffer_flow_ref(:,:,:,4)/(rbar*buffer_flow_ref(:,:,:,5))
          elseif(varidx.eq.7) then
            buffer_cal(:,:,:,1)     = buffer_flow(:,:,:,4)*( ( buffer_flow(:,:,:,5)+ 0.5*(1.4-1.0)/1.4/rbar*(buffer_flow(:,:,:,1)**2 + &
                                      buffer_flow(:,:,:,2)**2 + buffer_flow(:,:,:,3)**2 ) )/buffer_flow(:,:,:,5) )**(1.4/(1.4-1.d0))  ! P0
            buffer_cal_ref(:,:,:,1) = buffer_flow_ref(:,:,:,4)*( ( buffer_flow_ref(:,:,:,5)+ 0.5*(1.4-1.0)/1.4/rbar*(buffer_flow_ref(:,:,:,1)**2 + &
                                      buffer_flow_ref(:,:,:,2)**2 + buffer_flow_ref(:,:,:,3)**2 ) )/buffer_flow_ref(:,:,:,5) )**(1.4/(1.4-1.d0))  ! P0
          elseif(varidx.eq.8) then
            buffer_cal(:,:,:,1)     = buffer_flow(:,:,:,5)+ 0.5*(1.4-1.0)/1.4/rbar*(buffer_flow(:,:,:,1)**2 + &
                                      buffer_flow(:,:,:,2)**2 + buffer_flow(:,:,:,3)**2 )   ! T0
            buffer_cal_ref(:,:,:,1) = buffer_flow_ref(:,:,:,5)+ 0.5*(1.4-1.0)/1.4/rbar*(buffer_flow_ref(:,:,:,1)**2 + &
                                      buffer_flow_ref(:,:,:,2)**2 + buffer_flow_ref(:,:,:,3)**2 )   ! T0
          elseif(varidx.eq.9) then
            buffer_cal(:,:,:,1)     = buffer_flow(:,:,:,4)/(rbar*buffer_flow(:,:,:,5))*buffer_flow(:,:,:,1)              ! rho*u
            buffer_cal_ref(:,:,:,1) = buffer_flow_ref(:,:,:,4)/(rbar*buffer_flow_ref(:,:,:,5))*buffer_flow_ref(:,:,:,1)  ! rho*u
          elseif(varidx.eq.10) then
            buffer_cal(:,:,:,1)     = buffer_flow(:,:,:,4)/(rbar*buffer_flow(:,:,:,5))*buffer_flow(:,:,:,2)              ! rho*v
            buffer_cal_ref(:,:,:,1) = buffer_flow_ref(:,:,:,4)/(rbar*buffer_flow_ref(:,:,:,5))*buffer_flow_ref(:,:,:,2)  ! rho*v
          elseif(varidx.eq.11) then
            buffer_cal(:,:,:,1)     = buffer_flow(:,:,:,4)/(rbar*buffer_flow(:,:,:,5))*buffer_flow(:,:,:,3)              ! rho*w
            buffer_cal_ref(:,:,:,1) = buffer_flow_ref(:,:,:,4)/(rbar*buffer_flow_ref(:,:,:,5))*buffer_flow_ref(:,:,:,3)  ! rho*w
          elseif(varidx.eq.12) then
            buffer_cal(:,:,:,1)     = buffer_flow(:,:,:,1)*buffer_flow(:,:,:,2)          ! u*v
            buffer_cal_ref(:,:,:,1) = buffer_flow_ref(:,:,:,1)*buffer_flow_ref(:,:,:,2)  ! u*v
          elseif(varidx.eq.13) then
            buffer_cal(:,:,:,1)     = buffer_flow(:,:,:,1)*buffer_flow(:,:,:,3)          ! u*w
            buffer_cal_ref(:,:,:,1) = buffer_flow_ref(:,:,:,1)*buffer_flow_ref(:,:,:,3)  ! u*w
          elseif(varidx.eq.14) then
            buffer_cal(:,:,:,1)     = buffer_flow(:,:,:,2)*buffer_flow(:,:,:,3)          ! v*w
            buffer_cal_ref(:,:,:,1) = buffer_flow_ref(:,:,:,2)*buffer_flow_ref(:,:,:,3)  ! v*w
          else
            buffer_cal(:,:,:,1)     = buffer_flow(:,:,:,varidx)
            buffer_cal_ref(:,:,:,1) = buffer_flow_ref(:,:,:,varidx)
          endif
          ! Finished reading the required data

          ! Two-point correlation by Wiener-Khinchin theorem (make use of FFT)
          ! Ref: https://en.wikipedia.org/wiki/Autocorrelation
          ! Note: FFTW package is not normalized
          !       (in other words, applying the forward and then the backward transform will multiply the input by n = jmax)

          ! FFT(p)*congj(FFT(p))
          do ii = -corr3D%iwinl,corr3D%iwinr
            call fftw_crossspect_1m_z(buffer_cal(kshift+corr3D%decomp%zst(1):kshift+corr3D%decomp%zen(1), iave_st+ii:iave_end+ii,1:jmax,1), &
                                  buffer_cal_ref(kshift+corr3D%decomp%zst(1):kshift+corr3D%decomp%zen(1), iave_st:iave_end,      1:jmax,1), &
                                       specttmp1(corr3D%decomp%zst(1):corr3D%decomp%zen(1), iave_st:iave_end, 1:jmax/2+1) )
            call fftw_c2r_1m_z(specttmp1, corrtmp1(corr3D%decomp%zst(1):corr3D%decomp%zen(1),iave_st:iave_end,1:jmax))
            corrtmp1 = corrtmp1/dble(jmax) ! FFT normalization
            do k = corr3D%decomp%zst(1), corr3D%decomp%zen(1)
              corrtmp(ii,1:jmax,1,k) = corrtmp(ii,1:jmax,1,k) + sum(corrtmp1(k,iave_st:iave_end,1:jmax), dim=1)
            enddo

            ! ref
            call fftw_autospect_1m_z(  buffer_cal_ref(kshift+corr3D%decomp%zst(1):kshift+corr3D%decomp%zen(1), iave_st:iave_end, 1:jmax,1), &
                                           spect11tmp(corr3D%decomp%zst(1):corr3D%decomp%zen(1), iave_st:iave_end, 1:jmax/2+1) )
            call fftw_c2r_1m_z(spect11tmp, corr11tmp(corr3D%decomp%zst(1):corr3D%decomp%zen(1),iave_st:iave_end,1:jmax))
            corr11tmp = corr11tmp/dble(jmax)
            do k = corr3D%decomp%zst(1), corr3D%decomp%zen(1)
              corrtmp(ii,1:jmax,2,k) = corrtmp(ii,1:jmax,2,k) + sum(corr11tmp(k,iave_st:iave_end,1:jmax), dim=1)
            enddo

            call fftw_autospect_1m_z(  buffer_cal(kshift+corr3D%decomp%zst(1):kshift+corr3D%decomp%zen(1), iave_st+ii:iave_end+ii, 1:jmax,1), &
                                       spect22tmp(corr3D%decomp%zst(1):corr3D%decomp%zen(1), iave_st:iave_end, 1:jmax/2+1) )
            call fftw_c2r_1m_z(spect22tmp, corr22tmp(corr3D%decomp%zst(1):corr3D%decomp%zen(1),iave_st:iave_end,1:jmax))
            corr22tmp = corr22tmp/dble(jmax)
            do k = corr3D%decomp%zst(1), corr3D%decomp%zen(1)
              corrtmp(ii,1:jmax,3,k) = corrtmp(ii,1:jmax,3,k) + sum(corr22tmp(k,iave_st:iave_end,1:jmax), dim=1)
            enddo
           enddo ! end ii loop

         enddo ! end n loop

         if(itimeave.ne.0) then
           ! Average over all the DNS volumes

           call MPI_Allreduce(corrtmp,buffer_corr,icorr3D*jcorr3D*kcorr3D*nvarcorr_3D,&
                              MPI_DOUBLE_PRECISION,MPI_SUM,dim2_line,ierr)

           buffer_corr = buffer_corr/dble(nsamples_pervol*numvol)

           ! Write (time) volume-averaged 3D correlation
           write(unit=fnum1,fmt='(I04.4)') corr3D%ibe
           write(unit=fnum2,fmt='(I04.4)') corr3D%iend
           write(unit=fnum3,fmt='(I08.8)') iil
           write(unit=fnum4,fmt='(I08.8)') iih
           write(unit=fnum5,fmt='(I04.4)') corr3D%kref(1)
           write(unit=fnum6,fmt='(I04.4)') corr3D%kref(corr3D%num_kref)

           fname = 'corr3D_i'//fnum1//'-'//fnum2//'_timeave'//fnum3//'-'//fnum4//'_kplane'//fnum5//'-'//fnum6//'_'//trim(dname(varidx))//'.dat'
           if(myid.eq.0) print *, 'Writing corr3D file: ', trim(fname)
           call SaveCorr3D_Ascii(trim(fname))
         endif
      enddo ! end kk loop
    enddo ! end nn loop

    deallocate( xx, yy, zz, buffer_flow, buffer_flow_ref )
    deallocate( corrtmp, corrtmp1, buffer_corr, specttmp1)
    call fftw_finalize
    call decomp_info_finalize(corr3D%decomp)
    call decomp_2d_finalize
  endif ! end icalcorr3D.gt.0



  if(icalcorrxz.gt.0) then
    ! do automatic domain decomposition using 2DECOMP&FFT
    kmax_rd = corrxz%kend - corrxz%kbe + 1
    imax_rd = (corrxz%iend+corrxz%iwinr) - (corrxz%ibe-corrxz%iwinl) + 1
    jmax_rd = corrxz%jend - corrxz%jbe + 1
    periodic_bc = (/.false.,.false.,.true./)
    call decomp_2d_init(kmax_rd, imax_rd, jmax_rd, corrxz%p_row, corrxz%p_col)
    call decomp_info_init(kmax_rd, imax_rd, jmax_rd, corrxz%decomp)
    ! shift relative to DNS volume (for reading desired data within the full DNS volume)
    kshift = corrxz%kbe - 1
    ishift = corrxz%ibe - corrxz%iwinl - 1 ! ibe >= iwinl + 1 required
    jshift = corrxz%jbe - 1

    allocate( xx(kshift+corrxz%decomp%xst(1):kshift+corrxz%decomp%xen(1), &
                 ishift+corrxz%decomp%xst(2):ishift+corrxz%decomp%xen(2), &
                 jshift+corrxz%decomp%xst(3):jshift+corrxz%decomp%xen(3)), &
              yy(kshift+corrxz%decomp%xst(1):kshift+corrxz%decomp%xen(1), &
                 ishift+corrxz%decomp%xst(2):ishift+corrxz%decomp%xen(2), &
                 jshift+corrxz%decomp%xst(3):jshift+corrxz%decomp%xen(3)), &
              zz(kshift+corrxz%decomp%xst(1):kshift+corrxz%decomp%xen(1), &
                 ishift+corrxz%decomp%xst(2):ishift+corrxz%decomp%xen(2), &
                 jshift+corrxz%decomp%xst(3):jshift+corrxz%decomp%xen(3)) )

    fname = trim(gridfilename)
    if(myid.eq.0)  print *, 'Reading grid: ', trim(fname)
    call ReadHDF5grid_P(trim(fname),corrxz%decomp%xsz(1),&
                                    corrxz%decomp%xsz(2),&
                                    corrxz%decomp%xsz(3),&
                                    kshift + corrxz%decomp%xst(1),&
                                    ishift + corrxz%decomp%xst(2),&
                                    jshift + corrxz%decomp%xst(3),&
                                    xx, yy, zz)

    !>  1  2  3  4  5   6    7   8  9   10  11
    !> (u, v, w, p, T, rho, P0, T0, ru, rv, rw)
    allocate( buffer_flow(kshift+corrxz%decomp%xst(1):kshift+corrxz%decomp%xen(1), &
                          ishift+corrxz%decomp%xst(2):ishift+corrxz%decomp%xen(2), &
                          jshift+corrxz%decomp%xst(3):jshift+corrxz%decomp%xen(3), 5) )
    allocate(  buffer_cal(kshift+corrxz%decomp%xst(1):kshift+corrxz%decomp%xen(1), &
                          ishift+corrxz%decomp%xst(2):ishift+corrxz%decomp%xen(2), &
                          jshift+corrxz%decomp%xst(3):jshift+corrxz%decomp%xen(3), 1) )


    ! Size of correlation array for output
    icorrxz = corrxz%iwinl + corrxz%iwinr + 1
    kcorrxz = kmax_rd ! (N1 = corrxz%kend - corrxz%kbe + 1)

    iave_st = corrxz%ibe
    iave_end = corrxz%iend
    iref = (iave_st + iave_end)/2
    nsamples_pervol = (iave_end-iave_st+1)*jmax  ! # of samples for averaging in a single DNS volume
    numvol = (iih-iil)/stride+1  ! Total # DNS volumes for averaging 
    if(myid.eq.0) print *, 'iave_st =', iave_st, 'iave_end =', iave_end, 'icorrxz =', icorrxz

    allocate(      avetmp1(-corrxz%iwinl:corrxz%iwinr, corrxz%kbe:corrxz%kend, nvarave) )
    allocate(      avetmp2(-corrxz%iwinl:corrxz%iwinr, corrxz%kbe:corrxz%kend, nvarave) )
    allocate(  buffer_ave1(-corrxz%iwinl:corrxz%iwinr, corrxz%kbe:corrxz%kend, nvarave) )
    allocate(  buffer_ave2(-corrxz%iwinl:corrxz%iwinr, corrxz%kbe:corrxz%kend, nvarave) )
    allocate(      vartmp1(-corrxz%iwinl:corrxz%iwinr, corrxz%kbe:corrxz%kend, nvarave) )
    allocate(      vartmp2(-corrxz%iwinl:corrxz%iwinr, corrxz%kbe:corrxz%kend, nvarave) )
    allocate(  buffer_var1(-corrxz%iwinl:corrxz%iwinr, corrxz%kbe:corrxz%kend, nvarave) )
    allocate(  buffer_var2(-corrxz%iwinl:corrxz%iwinr, corrxz%kbe:corrxz%kend, nvarave) )
    allocate(     corrtmp(-corrxz%iwinl:corrxz%iwinr, corrxz%kbe:corrxz%kend, nvarcorr_xz, corrxz%num_kref) )
    allocate( buffer_corr(-corrxz%iwinl:corrxz%iwinr, corrxz%kbe:corrxz%kend, nvarcorr_xz, corrxz%num_kref) )

    do nn=1, num_output
      varidx = varindex(nn)
      if(myid.eq.0) then
        print *, '###############################################################'
        print *, 'Calculating x-z correlation for variable: ', trim(dname(varidx))
        print *, '###############################################################'
      endif
      avetmp1 = 0.d0  ! mean
      avetmp2 = 0.d0  ! mean
      vartmp1 = 0.d0  ! variance
      vartmp2 = 0.d0  ! variance
      corrtmp = 0.d0  ! correlation
      do n=iil, iih, stride
        write(unit=fnum,fmt='(I08.8)') n
        fname = trim(datapath)//fnum//'.h5'
        if(myid.eq.0)  print *, 'Reading file: ', trim(fname)
        ! Start to read required data from DNS volume
        call ReadHDF5sol_P( trim(fname),trim(groupname),corrxz%decomp%xsz(1),&
                                        corrxz%decomp%xsz(2),&
                                        corrxz%decomp%xsz(3),&
                                        kshift + corrxz%decomp%xst(1),&
                                        ishift + corrxz%decomp%xst(2),&
                                        jshift + corrxz%decomp%xst(3),&
                                        buffer_flow(:,:,:,1:5) )
        if(varidx.eq.6) then
          buffer_cal(:,:,:,1) = buffer_flow(:,:,:,4)/(rbar*buffer_flow(:,:,:,5)) ! rho
        elseif(varidx.eq.7) then
          buffer_cal(:,:,:,1) = buffer_flow(:,:,:,4)*( ( buffer_flow(:,:,:,5)+ 0.5*(1.4-1.0)/1.4/rbar*(buffer_flow(:,:,:,1)**2 + &
                                buffer_flow(:,:,:,2)**2 + buffer_flow(:,:,:,3)**2 ) )/buffer_flow(:,:,:,5) )**(1.4/(1.4-1.d0))  ! P0
        elseif(varidx.eq.8) then
          buffer_cal(:,:,:,1) = buffer_flow(:,:,:,5)+ 0.5*(1.4-1.0)/1.4/rbar*(buffer_flow(:,:,:,1)**2 + &
                                buffer_flow(:,:,:,2)**2 + buffer_flow(:,:,:,3)**2 )   ! T0
        elseif(varidx.eq.9) then
          buffer_cal(:,:,:,1) = buffer_flow(:,:,:,4)/(rbar*buffer_flow(:,:,:,5))*buffer_flow(:,:,:,1)  ! rho*u
        elseif(varidx.eq.10) then
          buffer_cal(:,:,:,1) = buffer_flow(:,:,:,4)/(rbar*buffer_flow(:,:,:,5))*buffer_flow(:,:,:,2)  ! rho*v
        elseif(varidx.eq.11) then
          buffer_cal(:,:,:,1) = buffer_flow(:,:,:,4)/(rbar*buffer_flow(:,:,:,5))*buffer_flow(:,:,:,3)  ! rho*w
        elseif(varidx.eq.12) then
          buffer_cal(:,:,:,1) = buffer_flow(:,:,:,1)*buffer_flow(:,:,:,2)          ! u*v
        elseif(varidx.eq.13) then
          buffer_cal(:,:,:,1) = buffer_flow(:,:,:,1)*buffer_flow(:,:,:,3)          ! u*w
        elseif(varidx.eq.14) then
          buffer_cal(:,:,:,1) = buffer_flow(:,:,:,2)*buffer_flow(:,:,:,3)          ! v*w
        else
          buffer_cal(:,:,:,1) = buffer_flow(:,:,:,varidx)
        endif
        ! Finished reading the required data

        do m = 1, corrxz%num_kref
          kloc = corrxz%kref(m)
          ! Calculate average qunatities (needed for deriving correlation coefficient)
          do i = iave_st, iave_end
            do ii = -corrxz%iwinl,corrxz%iwinr
              avetmp1(ii,kloc,1) = avetmp1(ii,kloc,1) + sum(buffer_cal(kloc,i,:,1))     ! varave
              vartmp1(ii,kloc,1) = vartmp1(ii,kloc,1) + sum(buffer_cal(kloc,i,:,1)**2)  ! var2ave
            enddo
          enddo
        enddo ! end m loop

        do k = corrxz%kbe, corrxz%kend
          do i = iave_st, iave_end
            do ii = -corrxz%iwinl,corrxz%iwinr
              avetmp2(ii,k,1) = avetmp2(ii,k,1) + sum(buffer_cal(k,i+ii,:,1))     ! varave
              vartmp2(ii,k,1) = vartmp2(ii,k,1) + sum(buffer_cal(k,i+ii,:,1)**2)  ! var2ave
            enddo
          enddo
        enddo

        do m = 1, corrxz%num_kref
          kloc = corrxz%kref(m)
          ! Calculate xz correlations (not yet normalized by rms values)
          do k = corrxz%kbe, corrxz%kend
            do ii = -corrxz%iwinl,corrxz%iwinr
              do i = iave_st, iave_end
                corrtmp(ii,k,1,m) = corrtmp(ii,k,1,m) + sum(buffer_cal(kloc,i,:,1)*buffer_cal(k,i+ii,:,1)) ! \bar(varref*var)
              enddo
            enddo
          enddo
        enddo ! end loop m

        if(itimeave.eq.0) then
          call MPI_Allreduce(avetmp1,buffer_ave1,icorrxz*kcorrxz*nvarave,&
                          MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
          call MPI_Allreduce(avetmp2,buffer_ave2,icorrxz*kcorrxz*nvarave,&
                          MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
          call MPI_Allreduce(vartmp1,buffer_var1,icorrxz*kcorrxz*nvarave,&
                          MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
          call MPI_Allreduce(vartmp2,buffer_var2,icorrxz*kcorrxz*nvarave,&
                          MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
          call MPI_Allreduce(corrtmp,buffer_corr,icorrxz*kcorrxz*nvarcorr_xz*corrxz%num_kref,&
                          MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

          buffer_ave1 = buffer_ave1/dble(nsamples_pervol)
          buffer_ave2 = buffer_ave2/dble(nsamples_pervol)
          buffer_var1 = buffer_var1/dble(nsamples_pervol)
          buffer_var2 = buffer_var2/dble(nsamples_pervol)
          buffer_corr = buffer_corr/dble(nsamples_pervol)

          do m = 1, corrxz%num_kref
            kloc = corrxz%kref(m)
            do k = corrxz%kbe, corrxz%kend
              do ii = -corrxz%iwinl,corrxz%iwinr
                buffer_corr(ii,k,1,m) = buffer_corr(ii,k,1,m) - buffer_ave1(ii,kloc,1)*buffer_ave2(ii,k,1)
              enddo
            enddo
          enddo ! end m loop
          ! Write xz correlation for the particular volume
          if(myid.eq.0) then
            write(unit=fnum1,fmt='(I04.4)') corrxz%ibe
            write(unit=fnum2,fmt='(I04.4)') corrxz%iend
            fname = 'corrxz_i'//fnum1//'-'//fnum2//'_'//fnum//'.dat'
            print *, 'Writing corrxz file: ', trim(fname)
            call SaveCorrxz_Ascii(trim(fname))
          endif
          avetmp1 = 0.d0
          avetmp2 = 0.d0
          vartmp1 = 0.d0
          vartmp2 = 0.d0
          corrtmp = 0.d0
        endif
 
      enddo ! end n loop

      if(itimeave.ne.0) then
       ! Average over all the DNS volumes
        call MPI_Allreduce(avetmp1,buffer_ave1,icorrxz*kcorrxz*nvarave,&
                          MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_Allreduce(avetmp2,buffer_ave2,icorrxz*kcorrxz*nvarave,&
                          MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_Allreduce(vartmp1,buffer_var1,icorrxz*kcorrxz*nvarave,&
                          MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_Allreduce(vartmp2,buffer_var2,icorrxz*kcorrxz*nvarave,&
                          MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_Allreduce(corrtmp,buffer_corr,icorrxz*kcorrxz*nvarcorr_xz*corrxz%num_kref,&
                          MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

        buffer_ave1 = buffer_ave1/dble(nsamples_pervol*numvol)
        buffer_ave2 = buffer_ave2/dble(nsamples_pervol*numvol)
        buffer_var1 = buffer_var1/dble(nsamples_pervol*numvol)
        buffer_var2 = buffer_var2/dble(nsamples_pervol*numvol)
        buffer_corr = buffer_corr/dble(nsamples_pervol*numvol)

        do m = 1, corrxz%num_kref
          kloc = corrxz%kref(m)
          do k = corrxz%kbe, corrxz%kend
            do ii = -corrxz%iwinl,corrxz%iwinr
              buffer_corr(ii,k,1,m) = buffer_corr(ii,k,1,m) - buffer_ave1(ii,kloc,1)*buffer_ave2(ii,k,1)
            enddo
          enddo
        enddo ! end m loop

        ! Write (time) volume-averaged xz correlation
        if(myid.eq.0) then
           write(unit=fnum1,fmt='(I04.4)') corrxz%ibe
           write(unit=fnum2,fmt='(I04.4)') corrxz%iend
           write(unit=fnum3,fmt='(I08.8)') iil
           write(unit=fnum4,fmt='(I08.8)') iih
           write(unit=fnum5,fmt='(I04.4)') corrxz%kref(1)
           write(unit=fnum6,fmt='(I04.4)') corrxz%kref(corrxz%num_kref)
           fname = 'corrxz_i'//fnum1//'-'//fnum2//'_timeave'//fnum3//'-'//fnum4//'_kref'//fnum5//'-'//fnum6//'_'//trim(dname(varindex(nn)))//'.dat'
           print *, 'Writing corrxz file: ', trim(fname)
           call SaveCorrxz_Ascii(trim(fname))
        endif
      endif
    
    enddo ! end nn loop

    deallocate( xx, yy, zz, buffer_flow )
    deallocate( avetmp1, buffer_ave1, vartmp1, buffer_var1, corrtmp, buffer_corr )
    deallocate( avetmp2, buffer_ave2, vartmp2, buffer_var2 )
    call decomp_info_finalize(corrxz%decomp)
    call decomp_2d_finalize
  endif ! end icalcorrxz.gt.0

  if(icalcorrxy.gt.0) then
    ! do automatic domain decomposition using 2DECOMP&FFT

    kmax_rd = corrxy%num_kref
    imax_rd = corrxy%iend - corrxy%ibe + 1
    jmax_rd = jmax
    periodic_bc = (/.false.,.false.,.true./)
    call decomp_2d_init(kmax_rd, imax_rd, jmax_rd, corrxy%p_row, corrxy%p_col)
    call decomp_info_init(kmax_rd, imax_rd, jmax_rd, corrxy%decomp)
    call fftw_init_zpencil(kmax_rd, imax_rd, jmax_rd) ! PHYSICAL_IN_Z = 3
    call MPI_Cart_Sub(DECOMP_2D_COMM_CART_Z,(/.false.,.true./),dim2_line,ierr)

    ! shift relative in x & y directions to DNS volume (for reading desired data within the full DNS volume)
    ishift = corrxy%ibe - 1
    ! Indexes ranges for reading DNS data
    ist_rd  = ishift+corrxy%decomp%zst(2)-corrxy%iwinl
    iend_rd = ishift+corrxy%decomp%zen(2)+corrxy%iwinr
    ilen_rd = iend_rd - ist_rd + 1

    ! Local length of k array size
    klen = corrxy%decomp%zsz(1)
    ! Sanity Check of domain decomposition
    if( .not.( (klen.eq.1).or.(klen.eq.kmax_rd/numprocs) ) ) then
        if(myid.eq.0) print *, 'local number of k planes klen =', klen
        if(myid.eq.0) print *, 'Domain decomposition error, klen must be equal to 1 or kmax_rd/numprocs'
        errcode = 227
        call MPI_Abort(MPI_comm_world,errcode,ierr)
    endif

    allocate( xx(corrxy%decomp%zst(1):corrxy%decomp%zen(1), ist_rd:iend_rd, 1:jmax), &
              yy(corrxy%decomp%zst(1):corrxy%decomp%zen(1), ist_rd:iend_rd, 1:jmax), &
              zz(corrxy%decomp%zst(1):corrxy%decomp%zen(1), ist_rd:iend_rd, 1:jmax) )

    fname = trim(gridfilename)
    if(myid.eq.0)  print *, 'Reading grid: ', trim(fname)
    do k = corrxy%decomp%zst(1), corrxy%decomp%zen(1)
       kstart = corrxy%kref(k)
       call ReadHDF5grid_P(trim(fname),1, ilen_rd, jmax,&
                                       kstart, ist_rd, 1, &
                                       xx(k:k,:,:), yy(k:k,:,:), zz(k:k,:,:))
    enddo

    allocate( buffer_flow(corrxy%decomp%zst(1):corrxy%decomp%zen(1),  &
                          ist_rd:iend_rd, jmax, 5) )
    allocate( buffer_cal(corrxy%decomp%zst(1):corrxy%decomp%zen(1),  &
                         ist_rd:iend_rd, jmax, 1) )

    ! Size of correlation array per core
    icorrxy = corrxy%iwinl + corrxy%iwinr + 1
    jcorrxy = jmax
    kcorrxy = corrxy%decomp%zsz(1)

    iave_st = ishift+corrxy%decomp%zst(2)
    iave_end = ishift+corrxy%decomp%zen(2)
    iref = (iave_st + iave_end)/2
    !nsamples_pervol = (iave_end-iave_st+1)*jmax  ! # of samples for averaging in a single DNS volume
    nsamples_pervol = imax_rd*jmax  ! # of samples for averaging in a single DNS volume
    numvol = (iih-iil)/stride+1  ! Total # DNS volumes for averaging
    if(myid.eq.0) print *, 'iave_st =', iave_st, 'iave_end =', iave_end, 'icorrxy =', icorrxy

    allocate(     corrtmp(-corrxy%iwinl:corrxy%iwinr, 1:jmax, nvarcorr_xy, corrxy%decomp%zst(1):corrxy%decomp%zen(1)) )
    allocate( buffer_corr(-corrxy%iwinl:corrxy%iwinr, 1:jmax, nvarcorr_xy, corrxy%decomp%zst(1):corrxy%decomp%zen(1)) )

    allocate(    specttmp1(corrxy%decomp%zst(1):corrxy%decomp%zen(1), iave_st:iave_end, 1:jmax/2+1) )
    allocate(     corrtmp1(corrxy%decomp%zst(1):corrxy%decomp%zen(1), iave_st:iave_end, 1:jmax) )
    allocate(    spect11tmp(corrxy%decomp%zst(1):corrxy%decomp%zen(1), iave_st:iave_end, 1:jmax/2+1) )
    allocate(    spect22tmp(corrxy%decomp%zst(1):corrxy%decomp%zen(1), iave_st:iave_end, 1:jmax/2+1) )
    allocate(     corr11tmp(corrxy%decomp%zst(1):corrxy%decomp%zen(1), iave_st:iave_end, 1:jmax) )
    allocate(     corr22tmp(corrxy%decomp%zst(1):corrxy%decomp%zen(1), iave_st:iave_end, 1:jmax) )

    do nn=1, num_output
      varidx = varindex(nn)
      if(myid.eq.0) then
        print *, '###############################################################'
        print *, 'Calculating x-y correlation for variable: ', trim(dname(varidx))
        print *, '###############################################################'
      endif
      corrtmp = 0.d0 ! correlation
      do n=iil, iih, stride
        write(unit=fnum,fmt='(I08.8)') n
        fname = trim(datapath)//fnum//'.h5'
        if(myid.eq.0)  print *, 'Reading file: ', trim(fname)
        ! Start to read required data from DNS volume
        do k = corrxy%decomp%zst(1), corrxy%decomp%zen(1)
          kstart = corrxy%kref(k)
          call ReadHDF5sol_P( trim(fname),trim(groupname), 1, ilen_rd, jmax, kstart, ist_rd, 1, &
                              buffer_flow(k:k,:,:,1:5) )

          if(varidx.eq.6) then
            buffer_cal(k,:,:,1) = buffer_flow(k,:,:,4)/(rbar*buffer_flow(k,:,:,5))  ! rho
          elseif(varidx.eq.7) then
            buffer_cal(k,:,:,1) = buffer_flow(k,:,:,4)*( ( buffer_flow(k,:,:,5)+ 0.5*(1.4-1.0)/1.4/rbar*(buffer_flow(k,:,:,1)**2 + &
                                buffer_flow(k,:,:,2)**2 + buffer_flow(k,:,:,3)**2 ) )/buffer_flow(k,:,:,5) )**(1.4/(1.4-1.d0))  ! P0
          elseif(varidx.eq.8) then
            buffer_cal(k,:,:,1) = buffer_flow(k,:,:,5)+ 0.5*(1.4-1.0)/1.4/rbar*(buffer_flow(k,:,:,1)**2 + &
                              buffer_flow(k,:,:,2)**2 + buffer_flow(k,:,:,3)**2 )   ! T0
          elseif(varidx.eq.9) then
            buffer_cal(k,:,:,1) = buffer_flow(k,:,:,4)/(rbar*buffer_flow(k,:,:,5))*buffer_flow(k,:,:,1)  ! rho*u
          elseif(varidx.eq.10) then
            buffer_cal(k,:,:,1) = buffer_flow(k,:,:,4)/(rbar*buffer_flow(k,:,:,5))*buffer_flow(k,:,:,2)  ! rho*v
          elseif(varidx.eq.11) then
            buffer_cal(k,:,:,1) = buffer_flow(k,:,:,4)/(rbar*buffer_flow(k,:,:,5))*buffer_flow(k,:,:,3)  ! rho*w
          elseif(varidx.eq.12) then
            buffer_cal(k,:,:,1) = buffer_flow(k,:,:,1)*buffer_flow(k,:,:,2)          ! u*v
          elseif(varidx.eq.13) then
            buffer_cal(k,:,:,1) = buffer_flow(k,:,:,1)*buffer_flow(k,:,:,3)          ! u*w
          elseif(varidx.eq.14) then
            buffer_cal(k,:,:,1) = buffer_flow(k,:,:,2)*buffer_flow(k,:,:,3)          ! v*w
          else
            buffer_cal(k,:,:,1) = buffer_flow(k,:,:,varidx)
          endif
        enddo ! end k loop
        ! Finished reading the required data

        ! Two-point correlation by Wiener-Khinchin theorem (make use of FFT)
        ! Ref: https://en.wikipedia.org/wiki/Autocorrelation
        ! Note: FFTW package is not normalized
        !       (in other words, applying the forward and then the backward transform will multiply the input by n = jmax)

        ! FFT(var)*congj(FFT(var))
        do ii = -corrxy%iwinl,corrxy%iwinr
          call fftw_crossspect_1m_z(buffer_cal(corrxy%decomp%zst(1):corrxy%decomp%zen(1), iave_st+ii:iave_end+ii,1:jmax,1), &
                                    buffer_cal(corrxy%decomp%zst(1):corrxy%decomp%zen(1), iave_st:iave_end,      1:jmax,1), &
                                    specttmp1(corrxy%decomp%zst(1):corrxy%decomp%zen(1), iave_st:iave_end, 1:jmax/2+1) )
          call fftw_c2r_1m_z(specttmp1, corrtmp1(corrxy%decomp%zst(1):corrxy%decomp%zen(1),iave_st:iave_end,1:jmax))
          corrtmp1 = corrtmp1/dble(jmax) ! FFT normalization  !!!!!
          do k = corrxy%decomp%zst(1), corrxy%decomp%zen(1)
            corrtmp(ii,1:jmax,1,k) = corrtmp(ii,1:jmax,1,k) + sum(corrtmp1(k,iave_st:iave_end,1:jmax), dim=1)
          enddo

          ! ref
          call fftw_autospect_1m_z(  buffer_cal(corrxy%decomp%zst(1):corrxy%decomp%zen(1), iave_st:iave_end, 1:jmax,1), &
                                     spect11tmp(corrxy%decomp%zst(1):corrxy%decomp%zen(1), iave_st:iave_end, 1:jmax/2+1) )
          call fftw_c2r_1m_z(spect11tmp, corr11tmp(corrxy%decomp%zst(1):corrxy%decomp%zen(1),iave_st:iave_end,1:jmax))
          corr11tmp = corr11tmp/dble(jmax)
          do k = corrxy%decomp%zst(1), corrxy%decomp%zen(1)
            corrtmp(ii,1:jmax,2,k) = corrtmp(ii,1:jmax,2,k) + sum(corr11tmp(k,iave_st:iave_end,1:jmax), dim=1)
          enddo

          call fftw_autospect_1m_z(  buffer_cal(corrxy%decomp%zst(1):corrxy%decomp%zen(1), iave_st+ii:iave_end+ii, 1:jmax,1), &
                                     spect22tmp(corrxy%decomp%zst(1):corrxy%decomp%zen(1), iave_st:iave_end, 1:jmax/2+1) )
          call fftw_c2r_1m_z(spect22tmp, corr22tmp(corrxy%decomp%zst(1):corrxy%decomp%zen(1),iave_st:iave_end,1:jmax))
          corr22tmp = corr22tmp/dble(jmax)
          do k = corrxy%decomp%zst(1), corrxy%decomp%zen(1)
            corrtmp(ii,1:jmax,3,k) = corrtmp(ii,1:jmax,3,k) + sum(corr22tmp(k,iave_st:iave_end,1:jmax), dim=1)
          enddo
        enddo ! end ii loop


        if(itimeave.eq.0) then
          call MPI_Allreduce(corrtmp,buffer_corr,icorrxy*jcorrxy*kcorrxy*nvarcorr_xy,&
                             MPI_DOUBLE_PRECISION,MPI_SUM,dim2_line,ierr)
          buffer_corr = buffer_corr/dble(nsamples_pervol)

          ! Write xy correlation for the particular volume
          write(unit=fnum1,fmt='(I04.4)') corrxy%ibe
          write(unit=fnum2,fmt='(I04.4)') corrxy%iend
          fname = 'corrxy_i'//fnum1//'-'//fnum2//'_'//fnum//'.dat'
          if(myid.eq.0)  print *, 'Writing corrxy file: ', trim(fname)
          call SaveCorrxy_Ascii(trim(fname))

         !! writing zgrid file
!         fname = 'corrxy_i'//fnum1//'-'//fnum2//'_'//fnum//'_zgrid.dat'
!         if(myid.eq.0) then
!           call Writezgrid(fname)
!         endif

          avetmp1 = 0.d0
          avetmp2 = 0.d0
          vartmp1 = 0.d0
          vartmp2 = 0.d0
          corrtmp = 0.d0
        endif

      enddo ! end n loop

      if(itimeave.ne.0) then
        ! Average over all the DNS volumes
        call MPI_Allreduce(corrtmp,buffer_corr,icorrxy*jcorrxy*kcorrxy*nvarcorr_xy,&
                           MPI_DOUBLE_PRECISION,MPI_SUM,dim2_line,ierr)
        buffer_corr = buffer_corr/dble(nsamples_pervol*numvol)

        ! Write (time) volume-averaged xy correlation
        write(unit=fnum1,fmt='(I04.4)') corrxy%ibe
        write(unit=fnum2,fmt='(I04.4)') corrxy%iend
        write(unit=fnum3,fmt='(I08.8)') iil
        write(unit=fnum4,fmt='(I08.8)') iih
        write(unit=fnum5,fmt='(I04.4)') corrxy%kref(1)
        write(unit=fnum6,fmt='(I04.4)') corrxy%kref(corrxy%num_kref)

        fname = 'corrxy_i'//fnum1//'-'//fnum2//'_timeave'//fnum3//'-'//fnum4//'_kplane'//fnum5//'-'//fnum6//'_'//trim(dname(varidx))//'.dat'
        if(myid.eq.0) print *, 'Writing corrxy file: ', trim(fname)
        call SaveCorrxy_Ascii(trim(fname))

!        fname = 'corrxy_i'//fnum1//'-'//fnum2//'_timeave'//fnum3//'-'//fnum4//'_zgrid.dat'
!        if(myid.eq.0) then
!          call Writezgrid(fname)
!        endif

      endif

    enddo ! end nn loop

    deallocate( xx, yy, zz, buffer_flow )
    deallocate( corrtmp, corrtmp1, buffer_corr, specttmp1)
    call fftw_finalize
    call decomp_info_finalize(corrxy%decomp)
    call decomp_2d_finalize
  endif ! end icalcorrxy.gt.0


  !> corryz correlation
  if(icalcorryz.gt.0) then
    kmax_rd = corryz%kend - corryz%kbe + 1
    imax_rd = corryz%iend - corryz%ibe + 1
    jmax_rd = corryz%jend - corryz%jbe + 1
    periodic_bc = (/.false.,.false.,.true./)
    call decomp_2d_init(kmax_rd, imax_rd, jmax_rd, corryz%p_row, corryz%p_col)
    call decomp_info_init(kmax_rd, imax_rd, jmax_rd, corryz%decomp)
    !call fftw_init_zpencil(kmax_rd, imax_rd, jmax_rd) ! PHYSICAL_IN_Z = 3
    call fftw_init_zpencil(1, imax_rd, jmax_rd) ! PHYSICAL_IN_Z = 3

    kshift = corryz%kbe - 1
    ishift = corryz%ibe - 1
    jshift = corryz%jbe - 1

    allocate( xx(kshift+corryz%decomp%zst(1):kshift+corryz%decomp%zen(1), &
                 ishift+corryz%decomp%zst(2):ishift+corryz%decomp%zen(2), &
                 jshift+corryz%decomp%zst(3):jshift+corryz%decomp%zen(3)), &
              yy(kshift+corryz%decomp%zst(1):kshift+corryz%decomp%zen(1), &
                 ishift+corryz%decomp%zst(2):ishift+corryz%decomp%zen(2), &
                 jshift+corryz%decomp%zst(3):jshift+corryz%decomp%zen(3)), &
              zz(kshift+corryz%decomp%zst(1):kshift+corryz%decomp%zen(1), &
                 ishift+corryz%decomp%zst(2):ishift+corryz%decomp%zen(2), &
                 jshift+corryz%decomp%zst(3):jshift+corryz%decomp%zen(3)) )


    fname = trim(gridfilename)
    if(myid.eq.0)  print *, 'Reading grid: ', trim(fname)
    call ReadHDF5grid_P(trim(fname),corryz%decomp%zsz(1),&
                                    corryz%decomp%zsz(2),&
                                    corryz%decomp%zsz(3),&
                                    kshift + corryz%decomp%zst(1),&
                                    ishift + corryz%decomp%zst(2),&
                                    jshift + corryz%decomp%zst(3),&
                                    xx, yy, zz)

    allocate( buffer_flow(kshift+corryz%decomp%zst(1):kshift+corryz%decomp%zen(1), &
                          ishift+corryz%decomp%zst(2):ishift+corryz%decomp%zen(2), &
                          jshift+corryz%decomp%zst(3):jshift+corryz%decomp%zen(3), 5) )
    allocate( buffer_cal(kshift+corryz%decomp%zst(1):kshift+corryz%decomp%zen(1), &
                         ishift+corryz%decomp%zst(2):ishift+corryz%decomp%zen(2), &
                         jshift+corryz%decomp%zst(3):jshift+corryz%decomp%zen(3), 1) )

    ! Size of correlation array for output
    jcorryz = jmax_rd
    kcorryz = kmax_rd

    jref = (corryz%jend + corryz%jbe)/2

    nsamples_pervol = imax_rd*jmax  ! # of samples for averaging in a single DNS volume
    numvol = (iih-iil)/stride+1  ! Total # DNS volumes for averaging


    allocate(     corrtmp(1:jmax, corryz%kbe:corryz%kend, nvarcorr_yz, corryz%num_kref) )
    allocate( buffer_corr(1:jmax, corryz%kbe:corryz%kend, nvarcorr_yz, corryz%num_kref) )

    allocate(    specttmp1(1:1, ishift+corryz%decomp%zst(2):ishift+corryz%decomp%zen(2), 1:jmax/2+1) )
    allocate(     corrtmp1(1:1, ishift+corryz%decomp%zst(2):ishift+corryz%decomp%zen(2), 1:jmax) )
    allocate(    spect11tmp(1:1, ishift+corryz%decomp%zst(2):ishift+corryz%decomp%zen(2), 1:jmax/2+1) )
    allocate(    spect22tmp(1:1, ishift+corryz%decomp%zst(2):ishift+corryz%decomp%zen(2), 1:jmax/2+1) )
    allocate(     corr11tmp(1:1, ishift+corryz%decomp%zst(2):ishift+corryz%decomp%zen(2), 1:jmax ))
    allocate(     corr22tmp(1:1, ishift+corryz%decomp%zst(2):ishift+corryz%decomp%zen(2), 1:jmax ))

    do nn=1, num_output
      varidx = varindex(nn)
      if(myid.eq.0) then
        print *, '###############################################################'
        print *, 'Calculating y-z correlation for variable: ', trim(dname(varidx))
        print *, '###############################################################'
      endif
      corrtmp = 0.d0 ! correlation
      buffer_corr = 0.d0
      do n=iil, iih, stride
        write(unit=fnum,fmt='(I08.8)') n
        fname = trim(datapath)//fnum//'.h5'
        if(myid.eq.0)  print *, 'Reading file: ', trim(fname)
        ! Start to read required data from DNS volume
        call ReadHDF5sol_P( trim(fname),trim(groupname),corryz%decomp%zsz(1),&
                                        corryz%decomp%zsz(2),&
                                        corryz%decomp%zsz(3),&
                                        kshift + corryz%decomp%zst(1),&
                                        ishift + corryz%decomp%zst(2),&
                                        jshift + corryz%decomp%zst(3),&
                                        buffer_flow(:,:,:,1:5) )
        if(varidx.eq.6) then
          buffer_cal(:,:,:,1) = buffer_flow(:,:,:,4)/(rbar*buffer_flow(:,:,:,5)) ! rho
        elseif(varidx.eq.7) then
          buffer_cal(:,:,:,1) = buffer_flow(:,:,:,4)*( ( buffer_flow(:,:,:,5)+ 0.5*(1.4-1.0)/1.4/rbar*(buffer_flow(:,:,:,1)**2 + &
                                buffer_flow(:,:,:,2)**2 + buffer_flow(:,:,:,3)**2 ) )/buffer_flow(:,:,:,5) )**(1.4/(1.4-1.d0))  ! P0
        elseif(varidx.eq.8) then
          buffer_cal(:,:,:,1) = buffer_flow(:,:,:,5)+ 0.5*(1.4-1.0)/1.4/rbar*(buffer_flow(:,:,:,1)**2 + &
                                buffer_flow(:,:,:,2)**2 + buffer_flow(:,:,:,3)**2 )   ! T0
        elseif(varidx.eq.9) then
          buffer_cal(:,:,:,1) = buffer_flow(:,:,:,4)/(rbar*buffer_flow(:,:,:,5))*buffer_flow(:,:,:,1)  ! rho*u
        elseif(varidx.eq.10) then
          buffer_cal(:,:,:,1) = buffer_flow(:,:,:,4)/(rbar*buffer_flow(:,:,:,5))*buffer_flow(:,:,:,2)  ! rho*v
        elseif(varidx.eq.11) then
          buffer_cal(:,:,:,1) = buffer_flow(:,:,:,4)/(rbar*buffer_flow(:,:,:,5))*buffer_flow(:,:,:,3)  ! rho*w
        elseif(varidx.eq.12) then
          buffer_cal(:,:,:,1) = buffer_flow(:,:,:,1)*buffer_flow(:,:,:,2)          ! u*v
        elseif(varidx.eq.13) then
          buffer_cal(:,:,:,1) = buffer_flow(:,:,:,1)*buffer_flow(:,:,:,3)          ! u*w
        elseif(varidx.eq.14) then
          buffer_cal(:,:,:,1) = buffer_flow(:,:,:,2)*buffer_flow(:,:,:,3)          ! v*w
        else
          buffer_cal(:,:,:,1) = buffer_flow(:,:,:,varidx)
        endif
        ! Finished reading the required data


        ! FFT(var)*congj(FFT(var))
        do m=1, corryz%num_kref
          kloc = corryz%kref(m)
          do kk = corryz%kbe, corryz%kend
            !print *, 'kk = ', kk
            call fftw_crossspect_1m_z(buffer_cal(kk:kk,     ishift+corryz%decomp%zst(2):ishift+corryz%decomp%zen(2),  1:jmax,1), &
                                      buffer_cal(kloc:kloc, ishift+corryz%decomp%zst(2):ishift+corryz%decomp%zen(2),  1:jmax,1), &
                                      specttmp1(  1:1,     ishift+corryz%decomp%zst(2):ishift+corryz%decomp%zen(2),  1:jmax/2+1) )

            call fftw_c2r_1m_z(specttmp1(1:1,ishift+corryz%decomp%zst(2):ishift+corryz%decomp%zen(2),1:jmax/2+1), corrtmp1(1:1,ishift+corryz%decomp%zst(2):ishift+corryz%decomp%zen(2),1:jmax))
            corrtmp1 = corrtmp1/dble(jmax) ! FFT normalization
            corrtmp(1:jmax,kk,1,m) = corrtmp(1:jmax,kk,1,m) + sum(corrtmp1(1,ishift+corryz%decomp%zst(2):ishift+corryz%decomp%zen(2),1:jmax),dim=1)

            ! ref
            call fftw_autospect_1m_z( buffer_cal(kloc:kloc, ishift+corryz%decomp%zst(2):ishift+corryz%decomp%zen(2),  1:jmax,1), &
                                      spect11tmp(  1:1,     ishift+corryz%decomp%zst(2):ishift+corryz%decomp%zen(2),  1:jmax/2+1) )

            call fftw_c2r_1m_z(spect11tmp(1:1,ishift+corryz%decomp%zst(2):ishift+corryz%decomp%zen(2),1:jmax/2+1), corr11tmp(1:1,ishift+corryz%decomp%zst(2):ishift+corryz%decomp%zen(2),1:jmax))
            corr11tmp = corr11tmp/dble(jmax) ! FFT normalization
            corrtmp(1:jmax,kk,2,m) = corrtmp(1:jmax,kk,2,m) + sum(corr11tmp(1,ishift+corryz%decomp%zst(2):ishift+corryz%decomp%zen(2),1:jmax),dim=1)

            call fftw_autospect_1m_z( buffer_cal(kk:kk, ishift+corryz%decomp%zst(2):ishift+corryz%decomp%zen(2),  1:jmax,1), &
                                      spect22tmp(  1:1,     ishift+corryz%decomp%zst(2):ishift+corryz%decomp%zen(2),  1:jmax/2+1) )
            call fftw_c2r_1m_z(spect22tmp(1:1,ishift+corryz%decomp%zst(2):ishift+corryz%decomp%zen(2),1:jmax/2+1), corr22tmp(1:1,ishift+corryz%decomp%zst(2):ishift+corryz%decomp%zen(2),1:jmax))
            corr22tmp = corr22tmp/dble(jmax) ! FFT normalization
            corrtmp(1:jmax,kk,3,m) = corrtmp(1:jmax,kk,3,m) + sum(corr22tmp(1,ishift+corryz%decomp%zst(2):ishift+corryz%decomp%zen(2),1:jmax),dim=1)
          enddo ! end kk loop
        enddo ! end m loop


        if(itimeave.eq.0) then
          call MPI_Allreduce(corrtmp,buffer_corr,jcorryz*kcorryz*nvarcorr_yz*corryz%num_kref,&
                             MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
          buffer_corr = buffer_corr/dble(nsamples_pervol)

          if(myid.eq.0) then
            ! Write yz correlation for the particular volume
            write(unit=fnum1,fmt='(I04.4)') corryz%kbe
            write(unit=fnum2,fmt='(I04.4)') corryz%kend
            fname = 'corryz_k'//fnum1//'-'//fnum2//'_'//fnum//'_'//trim(dname(varidx))//'.dat'
            if(myid.eq.0)  print *, 'Writing corryz file: ', trim(fname)
            call SaveCorryz_Ascii(trim(fname))
          endif
        endif ! end itimeave.eq.0
      enddo ! end n loop


      if(itimeave.ne.0) then
        call MPI_Allreduce(corrtmp,buffer_corr,jcorryz*kcorryz*nvarcorr_yz*corryz%num_kref,&
                           MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        buffer_corr = buffer_corr/dble(nsamples_pervol*numvol)
        if(myid.eq.0) then
          ! Write yz correlation for the particular volume
          write(unit=fnum1,fmt='(I04.4)') corryz%kbe
          write(unit=fnum2,fmt='(I04.4)') corryz%kend
          write(unit=fnum3,fmt='(I08.8)') iil
          write(unit=fnum4,fmt='(I08.8)') iih
          write(unit=fnum5,fmt='(I04.4)') corryz%kref(1)
          write(unit=fnum6,fmt='(I04.4)') corryz%kref(corryz%num_kref)
          fname = 'corryz_k'//fnum1//'-'//fnum2//'_timeave'//fnum3//'-'//fnum4//'_kref'//fnum5//'-'//fnum6//'_'//trim(dname(varidx))//'.dat'
          if(myid.eq.0)  print *, 'Writing corryz file: ', trim(fname)
          call SaveCorryz_Ascii(trim(fname))
        endif
      endif ! end itimeave.eq.0

    enddo ! end nn loop

    deallocate( xx, yy, zz, buffer_flow )
    deallocate( corrtmp, corrtmp1, buffer_corr, specttmp1)
    call fftw_finalize
    call decomp_info_finalize(corryz%decomp)
    call decomp_2d_finalize
  endif ! end icalcorryz.gt.0



  if(icalWNSD_ky.gt.0) then
    ! do automatic domain decomposition using 2DECOMP&FFT

    kmax_rd = WNSD_ky%num_kref
    imax_rd = WNSD_ky%iend - WNSD_ky%ibe + 1
    jmax_rd = jmax
    periodic_bc = (/.false.,.false.,.true./)
    call decomp_2d_init(kmax_rd, imax_rd, jmax_rd, WNSD_ky%p_row, WNSD_ky%p_col)
    call decomp_info_init(kmax_rd, imax_rd, jmax_rd, WNSD_ky%decomp)
    call fftw_init_zpencil(kmax_rd, imax_rd, jmax_rd) ! PHYSICAL_IN_Z = 3
    call MPI_Cart_Sub(DECOMP_2D_COMM_CART_Z,(/.false.,.true./),dim2_line,ierr)

    ! shift relative in x & y directions to DNS volume (for reading desired data within the full DNS volume)
    ishift = WNSD_ky%ibe - 1
    ! Indexes ranges for reading DNS data
    ist_rd  = ishift+WNSD_ky%decomp%zst(2)
    iend_rd = ishift+WNSD_ky%decomp%zen(2)
    ilen_rd = iend_rd - ist_rd + 1

    ! Local length of k array size
    klen = WNSD_ky%decomp%zsz(1)
    ! Sanity Check of domain decomposition
    if( .not.( (klen.eq.1).or.(klen.eq.kmax_rd/numprocs) ) ) then
        if(myid.eq.0) print *, 'local number of k planes klen =', klen
        if(myid.eq.0) print *, 'Domain decomposition error, klen must be equal to 1 or kmax_rd/numprocs'
        errcode = 227
        call MPI_Abort(MPI_comm_world,errcode,ierr)
    endif

    allocate( xx(WNSD_ky%decomp%zst(1):WNSD_ky%decomp%zen(1), ist_rd:iend_rd, 1:jmax), &
              yy(WNSD_ky%decomp%zst(1):WNSD_ky%decomp%zen(1), ist_rd:iend_rd, 1:jmax), &
              zz(WNSD_ky%decomp%zst(1):WNSD_ky%decomp%zen(1), ist_rd:iend_rd, 1:jmax) )

    fname = trim(gridfilename)
    if(myid.eq.0)  print *, 'Reading grid: ', trim(fname)
    do k = WNSD_ky%decomp%zst(1), WNSD_ky%decomp%zen(1)
       kstart = WNSD_ky%kref(k)
       call ReadHDF5grid_P(trim(fname),1, ilen_rd, jmax,&
                                       kstart, ist_rd, 1, &
                                       xx(k:k,:,:), yy(k:k,:,:), zz(k:k,:,:))
    enddo

    allocate( buffer_flow(WNSD_ky%decomp%zst(1):WNSD_ky%decomp%zen(1),  &
                          ist_rd:iend_rd, jmax, 5) )
    allocate( buffer_cal(WNSD_ky%decomp%zst(1):WNSD_ky%decomp%zen(1),  &
                         ist_rd:iend_rd, jmax, 1) )

    ! Size of correlation array per core
    iWNSD_ky = 1
    jWNSD_ky = jmax
    kWNSD_ky = WNSD_ky%decomp%zsz(1)

    iave_st = ishift+WNSD_ky%decomp%zst(2)
    iave_end = ishift+WNSD_ky%decomp%zen(2)
    iref = (iave_st + iave_end)/2
    !nsamples_pervol = (iave_end-iave_st+1)*jmax  ! # of samples for averaging in a single DNS volume
    nsamples_pervol = imax_rd*jmax  ! # of samples for averaging in a single DNS volume
    numvol = (iih-iil)/stride+1  ! Total # DNS volumes for averaging
    if(myid.eq.0) print *, 'iave_st =', iave_st, 'iave_end =', iave_end, 'iWNSD_ky =', iWNSD_ky

    allocate(    specttmp1(WNSD_ky%decomp%zst(1):WNSD_ky%decomp%zen(1), iave_st:iave_end, 1:jmax/2+1) )
    allocate(   spect1dtmp(WNSD_ky%decomp%zst(1):WNSD_ky%decomp%zen(1),                1, 1:jmax/2+1) )
    allocate(   spect1dave(WNSD_ky%decomp%zst(1):WNSD_ky%decomp%zen(1),                1, 1:jmax/2+1) )

    allocate(  buffer_ave1(WNSD_ky%decomp%zst(1):WNSD_ky%decomp%zen(1), 1, 1) )
    allocate(  buffer_var1(WNSD_ky%decomp%zst(1):WNSD_ky%decomp%zen(1), 1, 1) )

    do nn=1, num_output
      varidx = varindex(nn)
      if(myid.eq.0) then
        print *, '###############################################################'
        print *, 'Calculating wavenumber spectrum in j-direction for variable: ', trim(dname(varidx))
        print *, '###############################################################'
      endif
      spect1dtmp = 0.d0 ! correlation
      do n=iil, iih, stride
        write(unit=fnum,fmt='(I08.8)') n
        fname = trim(datapath)//fnum//'.h5'
        if(myid.eq.0)  print *, 'Reading file: ', trim(fname)
        ! Start to read required data from DNS volume
        do k = WNSD_ky%decomp%zst(1), WNSD_ky%decomp%zen(1)
          kstart = WNSD_ky%kref(k)
          call ReadHDF5sol_P( trim(fname),trim(groupname), 1, ilen_rd, jmax, kstart, ist_rd, 1, &
                              buffer_flow(k:k,:,:,1:5) )

          if(varidx.eq.6) then
            buffer_cal(k,:,:,1) = buffer_flow(k,:,:,4)/(rbar*buffer_flow(k,:,:,5))  ! rho
          elseif(varidx.eq.7) then
            buffer_cal(k,:,:,1) = buffer_flow(k,:,:,4)*( ( buffer_flow(k,:,:,5)+ 0.5*(1.4-1.0)/1.4/rbar*(buffer_flow(k,:,:,1)**2 + &
                                buffer_flow(k,:,:,2)**2 + buffer_flow(k,:,:,3)**2 ) )/buffer_flow(k,:,:,5) )**(1.4/(1.4-1.d0))  ! P0
          elseif(varidx.eq.8) then
            buffer_cal(k,:,:,1) = buffer_flow(k,:,:,5)+ 0.5*(1.4-1.0)/1.4/rbar*(buffer_flow(k,:,:,1)**2 + &
                              buffer_flow(k,:,:,2)**2 + buffer_flow(k,:,:,3)**2 )   ! T0
          elseif(varidx.eq.9) then
            buffer_cal(k,:,:,1) = buffer_flow(k,:,:,4)/(rbar*buffer_flow(k,:,:,5))*buffer_flow(k,:,:,1)  ! rho*u
          elseif(varidx.eq.10) then
            buffer_cal(k,:,:,1) = buffer_flow(k,:,:,4)/(rbar*buffer_flow(k,:,:,5))*buffer_flow(k,:,:,2)  ! rho*v
          elseif(varidx.eq.11) then
            buffer_cal(k,:,:,1) = buffer_flow(k,:,:,4)/(rbar*buffer_flow(k,:,:,5))*buffer_flow(k,:,:,3)  ! rho*w
          elseif(varidx.eq.12) then
            buffer_cal(k,:,:,1) = buffer_flow(k,:,:,1)*buffer_flow(k,:,:,2)          ! u*v
          elseif(varidx.eq.13) then
            buffer_cal(k,:,:,1) = buffer_flow(k,:,:,1)*buffer_flow(k,:,:,3)          ! u*w
          elseif(varidx.eq.14) then
            buffer_cal(k,:,:,1) = buffer_flow(k,:,:,2)*buffer_flow(k,:,:,3)          ! v*w
          else
            buffer_cal(k,:,:,1) = buffer_flow(k,:,:,varidx)
          endif
        enddo ! end k loop
        ! Finished reading the required data

        do k=WNSD_ky%decomp%zst(1), WNSD_ky%decomp%zen(1)
          buffer_ave1(k,1,1) = sum( buffer_cal(k,iave_st:iave_end,1:jmax_rd,1) )
          buffer_var1(k,1,1) = sum( buffer_cal(k,iave_st:iave_end,1:jmax_rd,1)**2 )
        enddo

        ! Two-point correlation by Wiener-Khinchin theorem (make use of FFT)
        ! Ref: https://en.wikipedia.org/wiki/Autocorrelation
        ! Note: FFTW package is not normalized
        !       (in other words, applying the forward and then the backward transform will multiply the input by n = jmax)

        call fftw_autospect_1m_z(  buffer_cal(WNSD_ky%decomp%zst(1):WNSD_ky%decomp%zen(1), iave_st:iave_end, 1:jmax,1), &
                                    specttmp1(WNSD_ky%decomp%zst(1):WNSD_ky%decomp%zen(1), iave_st:iave_end, 1:jmax/2+1) )
        specttmp1 = specttmp1/dble(jmax)
        do k = WNSD_ky%decomp%zst(1), WNSD_ky%decomp%zen(1)
          spect1dtmp(k,1,1:jmax/2+1) = spect1dtmp(k,1,1:jmax/2+1) + sum(specttmp1(k,iave_st:iave_end,1:jmax/2+1), dim=1)
        enddo

      enddo ! end n loop

      ! Average over all the DNS volumes
      !call MPI_Allreduce(spect1dtmp,spect1dave,iWNSD_ky*jWNSD_ky*kWNSD_ky,&
      !                   MPI_DOUBLE_PRECISION,MPI_SUM,dim2_line,ierr)
      spect1dave = spect1dtmp
      spect1dave = spect1dave/dble(nsamples_pervol*numvol)

      buffer_ave1 = buffer_ave1/dble(ilen_rd*jmax_rd*numvol)
      buffer_var1 = buffer_var1/dble(ilen_rd*jmax_rd*numvol)
      buffer_var1 = abs(buffer_var1 - buffer_ave1**2)

      ! Write (time) volume-averaged xy correlation
      write(unit=fnum1,fmt='(I04.4)') WNSD_ky%ibe
      write(unit=fnum2,fmt='(I04.4)') WNSD_ky%iend
      write(unit=fnum3,fmt='(I08.8)') iil
      write(unit=fnum4,fmt='(I08.8)') iih
      write(unit=fnum5,fmt='(I04.4)') WNSD_ky%kref(1)
      write(unit=fnum6,fmt='(I04.4)') WNSD_ky%kref(WNSD_ky%num_kref)

      fname = 'WNSD_ky_i'//fnum1//'-'//fnum2//'_timeave'//fnum3//'-'//fnum4//'_kplane'//fnum5//'-'//fnum6//'_'//trim(dname(varidx))//'.dat'
      if(myid.eq.0) print *, 'Writing WNSD_ky file: ', trim(fname)
      call SaveWNSD_ky_Ascii(trim(fname))

    enddo ! end nn loop

    deallocate( xx, yy, zz, buffer_flow )
    deallocate(specttmp1,spect1dtmp,spect1dave)
    call fftw_finalize
    call decomp_info_finalize(WNSD_ky%decomp)
    call decomp_2d_finalize
  endif ! end icalWNSD_ky.gt.0


  if(icalWNSD_kx.gt.0) then
    ! do automatic domain decomposition using 2DECOMP&FFT

    kmax_rd = WNSD_kx%num_kref
    imax_rd = WNSD_kx%iend - WNSD_kx%ibe + 1
    jmax_rd = WNSD_kx%jend - WNSD_kx%jbe + 1
    !periodic_bc = (/.false.,.false.,.true./)
    periodic_bc = (/.false.,.true.,.false./)
    call decomp_2d_init(kmax_rd, imax_rd, jmax_rd, WNSD_kx%p_row, WNSD_kx%p_col)
    call decomp_info_init(kmax_rd, imax_rd, jmax_rd, WNSD_kx%decomp)
    call fftw_init_ypencil(kmax_rd, imax_rd, jmax_rd)
    !call MPI_Cart_Sub(DECOMP_2D_COMM_CART_Y,(/.false.,.true./),dim2_line,ierr)
    call InitFFTWindow(imax_rd, WNSD_kx%iwindow)

    ! shift relative in x & y directions to DNS volume (for reading desired data within the full DNS volume)
    ishift = WNSD_kx%ibe - 1
    ! Indexes ranges for reading DNS data
    ist_rd  = WNSD_kx%ibe
    iend_rd = WNSD_kx%iend
    ilen_rd = iend_rd - ist_rd + 1

    ! Local length of k array size
    klen = WNSD_kx%decomp%ysz(1)

!print *, 'klen = ', klen, '', WNSD_kx%decomp%ysz(2), WNSD_kx%decomp%ysz(3)


!    ! Sanity Check of domain decomposition
!    if( .not.( (klen.eq.1).or.(klen.eq.kmax_rd/numprocs) ) ) then
!        if(myid.eq.0) print *, 'local number of k planes klen =', klen
!        if(myid.eq.0) print *, 'Domain decomposition error, klen must be equal to 1 or kmax_rd/numprocs'
!        errcode = 227
!        call MPI_Abort(MPI_comm_world,errcode,ierr)
!    endif
!
    allocate( xx(WNSD_kx%decomp%yst(1):WNSD_kx%decomp%yen(1), ist_rd:iend_rd, 1:jmax_rd), &
              yy(WNSD_kx%decomp%yst(1):WNSD_kx%decomp%yen(1), ist_rd:iend_rd, 1:jmax_rd), &
              zz(WNSD_kx%decomp%yst(1):WNSD_kx%decomp%yen(1), ist_rd:iend_rd, 1:jmax_rd) )
    allocate( xx_tmp(WNSD_kx%decomp%yst(1):WNSD_kx%decomp%yen(1), ist_rd:iend_rd, 1:jmax_rd) )

    fname = trim(gridfilename)
    if(myid.eq.0)  print *, 'Reading grid: ', trim(fname)
    do k = WNSD_kx%decomp%yst(1), WNSD_kx%decomp%yen(1)
       kstart = WNSD_kx%kref(k)
       call ReadHDF5grid_P(trim(fname),1, ilen_rd, jmax_rd,&
                                       kstart, ist_rd, 1, &
                                       xx_tmp(k:k,:,:), yy(k:k,:,:), zz(k:k,:,:))
    enddo

    deallocate(yy,zz)
    ! check x-grid distribution
    allocate(ddx(ist_rd:iend_rd-1))
    do i=ist_rd, iend_rd-1
      ddx(i) = xx_tmp(WNSD_kx%decomp%yst(1),i+1,1) - xx_tmp(WNSD_kx%decomp%yst(1),i,1)
    enddo
    ddx_min = minval(ddx)
    ddx_max = maxval(ddx)
    if(abs(ddx_max-ddx_min).gt.1.e-10) then
      WNSD_kx%iinterp_x = 1
    else
      WNSD_kx%iinterp_x = 0
    endif

    allocate( buffer_flow(WNSD_kx%decomp%yst(1):WNSD_kx%decomp%yen(1),  &
                          ist_rd:iend_rd, 1:jmax_rd, 5) )
    allocate(  buffer_cal(WNSD_kx%decomp%yst(1):WNSD_kx%decomp%yen(1),  &
                          1:ilen_rd, 1:jmax_rd, 2) )
    allocate(buffer_2D(WNSD_kx%decomp%yst(1):WNSD_kx%decomp%yen(1),1:ilen_rd,1))

    ! Size of correlation array per core
    iWNSD_kx = WNSD_kx%iend - WNSD_kx%ibe + 1
    jWNSD_kx = 1
    kWNSD_kx = WNSD_kx%decomp%ysz(1)

    iave_st = WNSD_kx%ibe
    iave_end = WNSD_kx%iend
    iref = (iave_st + iave_end)/2
    nsamples_pervol = jmax_rd*iWNSD_kx  ! # of samples for averaging in a single DNS volume
    numvol = (iih-iil)/stride+1  ! Total # DNS volumes for averaging
    if(myid.eq.0) print *, 'ist =', iave_st, 'iend =', iave_end, 'iWNSD_kx =', iWNSD_kx

    allocate(    specttmp1(WNSD_kx%decomp%yst(1):WNSD_kx%decomp%yen(1), 1:iWNSD_kx/2+1, 1:jmax_rd) )
    allocate(   spect1dtmp(WNSD_kx%decomp%yst(1):WNSD_kx%decomp%yen(1), 1:iWNSD_kx/2+1, 1) )
    allocate(   spect1dave(WNSD_kx%decomp%yst(1):WNSD_kx%decomp%yen(1), 1:iWNSD_kx/2+1, 1) )

    allocate( buffer_ave1(WNSD_kx%decomp%yst(1):WNSD_kx%decomp%yen(1), 1, 1) )
    allocate( buffer_var1(WNSD_kx%decomp%yst(1):WNSD_kx%decomp%yen(1), 1, 1) )

    do nn=1, num_output
      varidx = varindex(nn)
      if(myid.eq.0) then
        print *, '###############################################################'
        print *, 'Calculating wavenumber spectrum in i-direction for variable: ', trim(dname(varidx))
        print *, '###############################################################'
      endif
      spect1dtmp = 0.d0 ! correlation
      buffer_ave1 = 0.d0
      buffer_var1 = 0.d0
      do n=iil, iih, stride
        write(unit=fnum,fmt='(I08.8)') n
        fname = trim(datapath)//fnum//'.h5'
        if(myid.eq.0)  print *, 'Reading file: ', trim(fname)
        ! Start to read required data from DNS volume
        do k = WNSD_kx%decomp%yst(1), WNSD_kx%decomp%yen(1)
          kstart = WNSD_kx%kref(k)
          call ReadHDF5sol_P( trim(fname),trim(groupname), 1, ilen_rd, jmax_rd, kstart, ist_rd, 1, &
                              buffer_flow(k:k,:,:,1:5) )

          if(varidx.eq.6) then
            buffer_cal(k,1:ilen_rd,:,2) = buffer_flow(k,:,:,4)/(rbar*buffer_flow(k,:,:,5))  ! rho
          elseif(varidx.eq.7) then
            buffer_cal(k,1:ilen_rd,:,2) = buffer_flow(k,:,:,4)*( ( buffer_flow(k,:,:,5)+ 0.5*(1.4-1.0)/1.4/rbar*(buffer_flow(k,:,:,1)**2 + &
                                buffer_flow(k,:,:,2)**2 + buffer_flow(k,:,:,3)**2 ) )/buffer_flow(k,:,:,5) )**(1.4/(1.4-1.d0))  ! P0
          elseif(varidx.eq.8) then
            buffer_cal(k,1:ilen_rd,:,2) = buffer_flow(k,:,:,5)+ 0.5*(1.4-1.0)/1.4/rbar*(buffer_flow(k,:,:,1)**2 + &
                              buffer_flow(k,:,:,2)**2 + buffer_flow(k,:,:,3)**2 )   ! T0
          elseif(varidx.eq.9) then
            buffer_cal(k,1:ilen_rd,:,2) = buffer_flow(k,:,:,4)/(rbar*buffer_flow(k,:,:,5))*buffer_flow(k,:,:,1)  ! rho*u
          elseif(varidx.eq.10) then
            buffer_cal(k,1:ilen_rd,:,2) = buffer_flow(k,:,:,4)/(rbar*buffer_flow(k,:,:,5))*buffer_flow(k,:,:,2)  ! rho*v
          elseif(varidx.eq.11) then
            buffer_cal(k,1:ilen_rd,:,2) = buffer_flow(k,:,:,4)/(rbar*buffer_flow(k,:,:,5))*buffer_flow(k,:,:,3)  ! rho*w
          elseif(varidx.eq.12) then
            buffer_cal(k,1:ilen_rd,:,2) = buffer_flow(k,:,:,1)*buffer_flow(k,:,:,2)          ! u*v
          elseif(varidx.eq.13) then
            buffer_cal(k,1:ilen_rd,:,2) = buffer_flow(k,:,:,1)*buffer_flow(k,:,:,3)          ! u*w
          elseif(varidx.eq.14) then
            buffer_cal(k,1:ilen_rd,:,2) = buffer_flow(k,:,:,2)*buffer_flow(k,:,:,3)          ! v*w
          else
            buffer_cal(k,1:ilen_rd,:,2) = buffer_flow(k,:,:,varidx)
          endif
        enddo ! end k loop


        if(WNSD_kx%iinterp_x.gt.0) then
          xlen = xx_tmp(WNSD_kx%decomp%yst(1),iend_rd,1) - xx_tmp(WNSD_kx%decomp%yst(1),ist_rd,1)
          do i=ist_rd, iend_rd
            xx(WNSD_kx%decomp%yst(1):WNSD_kx%decomp%yen(1), i, 1:jmax_rd) =  xx_tmp(WNSD_kx%decomp%yst(1):WNSD_kx%decomp%yen(1),ist_rd,1:jmax_rd) + dble(i-ist_rd)*xlen/dble(ilen_rd-1)
          enddo
          if(nrank.eq.0) then
            print *, 'interpolating in X direction ... '
          endif
          do k=WNSD_kx%decomp%yst(1), WNSD_kx%decomp%yen(1)
            do j=1, jmax_rd
              call SplineInterp(xx_tmp(k,:,j),buffer_cal(k,:,j,2),ilen_rd, &
                                    xx(k,:,j),buffer_cal(k,:,j,1),ilen_rd,2)
            enddo
          enddo
          if(nrank.eq.0) then
            print *, 'writing interpolation debug file: check_interp.dat' 
            open(7,file='check_interp.dat',status='unknown')
              write(7,'(2a)') 'variables=xold,varold,xnew,varnew'
              do i=ist_rd, iend_rd
                write(7,'(4E20.12)') xx_tmp(WNSD_kx%decomp%yst(1),i,1), buffer_cal(WNSD_kx%decomp%yst(1),i-ist_rd+1,1,2), &
                                         xx(WNSD_kx%decomp%yst(1),i,1), buffer_cal(WNSD_kx%decomp%yst(1),i-ist_rd+1,1,1)
              enddo
            close(7)
          endif

        else
          buffer_cal(:,:,:,1) = buffer_cal(:,:,:,2)
          xx(WNSD_kx%decomp%yst(1):WNSD_kx%decomp%yen(1), ist_rd:iend_rd, 1:jmax_rd) = xx_tmp(WNSD_kx%decomp%yst(1):WNSD_kx%decomp%yen(1), ist_rd:iend_rd, 1:jmax_rd)
        endif

!        do k=WNSD_kx%decomp%yst(1), WNSD_kx%decomp%yen(1)
!          buffer_ave1(k,1,1) = sum( buffer_cal(k,1:ilen_rd,1:jmax_rd,1) )
!          buffer_var1(k,1,1) = sum( buffer_cal(k,1:ilen_rd,1:jmax_rd,1)**2 )
!        enddo

        ! Finished reading the required data
        buffer_2D(WNSD_kx%decomp%yst(1):WNSD_kx%decomp%yen(1),1:ilen_rd,1) = sum(buffer_cal(WNSD_kx%decomp%yst(1):WNSD_kx%decomp%yen(1),  &
                         1:ilen_rd, 1:jmax_rd, 1),dim=3)/dble(jmax_rd)
        do k=WNSD_kx%decomp%yst(1), WNSD_kx%decomp%yen(1)
          do j=1, jmax_rd
            buffer_cal(k,1:ilen_rd,j,1) = (buffer_cal(k,1:ilen_rd,j,1) - buffer_2D(k,1:ilen_rd,1))*wcoeff(1:ilen_rd)
          enddo
        enddo

        do k=WNSD_kx%decomp%yst(1), WNSD_kx%decomp%yen(1)
          buffer_ave1(k,1,1) = sum( buffer_cal(k,1:ilen_rd,1:jmax_rd,1) )
          buffer_var1(k,1,1) = sum( buffer_cal(k,1:ilen_rd,1:jmax_rd,1)**2 )
        enddo

        ! Two-point correlation by Wiener-Khinchin theorem (make use of FFT)
        ! Ref: https://en.wikipedia.org/wiki/Autocorrelation
        ! Note: FFTW package is not normalized
        !       (in other words, applying the forward and then the backward transform will multiply the input by n = imax_rd)
        do j=1, jmax_rd
          call fftw_autospect_1m_y(  buffer_cal(WNSD_kx%decomp%yst(1):WNSD_kx%decomp%yen(1),      1:ilen_rd, j:j,1), &
                                      specttmp1(WNSD_kx%decomp%yst(1):WNSD_kx%decomp%yen(1),  1:ilen_rd/2+1, j:j) )
        enddo
        specttmp1 = specttmp1/dble(iWNSD_kx)
        do k = WNSD_kx%decomp%yst(1), WNSD_kx%decomp%yen(1)
          spect1dtmp(k,1:iWNSD_kx/2+1,1) = spect1dtmp(k,1:iWNSD_kx/2+1,1) + sum(specttmp1(k,1:iWNSD_kx/2+1,1:jmax_rd), dim=2)
        enddo

      enddo ! end n loop

      ! Average over all the DNS volumes
      spect1dave = spect1dtmp/dble(nsamples_pervol*numvol)
      buffer_ave1 = buffer_ave1/dble(ilen_rd*jmax_rd*numvol)
      buffer_var1 = buffer_var1/dble(ilen_rd*jmax_rd*numvol)
      buffer_var1 = abs(buffer_var1 - buffer_ave1**2)

      ! Write (time) volume-averaged wavenumber spectrum kx
      write(unit=fnum1,fmt='(I04.4)') WNSD_kx%ibe
      write(unit=fnum2,fmt='(I04.4)') WNSD_kx%iend
      write(unit=fnum3,fmt='(I08.8)') iil
      write(unit=fnum4,fmt='(I08.8)') iih
      write(unit=fnum5,fmt='(I04.4)') WNSD_kx%kref(1)
      write(unit=fnum6,fmt='(I04.4)') WNSD_kx%kref(WNSD_kx%num_kref)

      fname = 'WNSD_kx_i'//fnum1//'-'//fnum2//'_timeave'//fnum3//'-'//fnum4//'_kplane'//fnum5//'-'//fnum6//'_'//trim(dname(varidx))//'.dat'
      if(myid.eq.0) print *, 'Writing WNSD_kx file: ', trim(fname)
      call SaveWNSD_kx_Ascii(trim(fname))

    enddo ! end nn loop

    deallocate( xx, buffer_flow, buffer_cal, buffer_ave1, buffer_var1)
    deallocate(specttmp1,spect1dtmp,spect1dave)
    call fftw_finalize
    call decomp_info_finalize(WNSD_kx%decomp)
    call decomp_2d_finalize
  endif ! end icalWNSD_kx.gt.0

  call FinalizeHDF5()
  call MPI_FINALIZE(ierr)

  contains

    subroutine Writezgrid(fn)
      implicit none
      character(*), intent(in) :: fn

      open(7,file=trim(fn),status='unknown')
        write(7,'(a)') 'variables=z'
        do k=1, corrxy%num_kref
          write(7,*) zz(corrxy%kref(k),1,1)
        enddo
      close(7)
    end subroutine Writezgrid

    subroutine Input()
      integer :: i

      if(myid.eq.0) then
        read(*,*)
        read(*,'(A)')gridfilename
        read(*,*)
        read(*,'(A)')datapath  !path where data files are stored
        read(*,*)
        read(*,'(A)')groupname
        read(*,*)
        read(*,*)
        read(*,*) imax, jmax, kmax
        read(*,*)
        read(*,*) iil, iih, stride, itimeave
        read(*,*)
        read(*,*) Rm
        read(*,*)
        read(*,*) num_output
        allocate(varindex(num_output))
        read(*,*) (varindex(i),i=1,num_output)
        read(*,*)
        read(*,*)
        read(*,*) icalcorrxz 
        if(icalcorrxz.gt.0) then
           read(*,*)
           read(*,*) corrxz%num_kref
           allocate(corrxz%kref(corrxz%num_kref))
           read(*,*) (corrxz%kref(i), i = 1, corrxz%num_kref)
           read(*,*)
           read(*,*) corrxz%ibe, corrxz%iend, corrxz%kbe, corrxz%kend, corrxz%iwinl, corrxz%iwinr
           read(*,*)
           corrxz%jbe = 1
           corrxz%jend = jmax
           corrxz%p_row = 1            ! 1 core in x direction (x pencil)
           corrxz%p_col = numprocs     ! all the cores in y direction
        else
           read(*,*)
           read(*,*)
           read(*,*)
           read(*,*)
           read(*,*)
           read(*,*)
        endif

        read(*,*)
        read(*,*) icalcorrxy
        if(icalcorrxy.gt.0) then
           read(*,*)
           read(*,*) corrxy%num_kref
           allocate(corrxy%kref(corrxy%num_kref))
           read(*,*) (corrxy%kref(i), i = 1, corrxy%num_kref)
           read(*,*)
           read(*,*) corrxy%ibe, corrxy%iend, corrxy%iwinl, corrxy%iwinr
           read(*,*)
           corrxy%jbe = 1
           corrxy%jend = jmax
           corrxy%kbe = 1
           corrxy%kend = kmax
           if(numprocs.ge.corrxy%num_kref) then
              corrxy%p_row = corrxy%num_kref              ! p_row = klane
              corrxy%p_col = numprocs/corrxy%num_kref     ! Z-Pencil
              print *, '2D Domain decomposition (grid processor) for corrxy ...'
              print *, 'Make sure that min(nx,ny) >= p_row and min(ny,nz) >= p_col'
              print *, 'p_row is equal to num_kref: p_row =', corrxy%p_row
              print *, 'p_col is equal to numprocs/num_kref: p_col =', corrxy%p_col
           else
              corrxy%p_row = numprocs
              corrxy%p_col = 1
           endif
        else
           read(*,*)
           read(*,*)
           read(*,*)
           read(*,*)
           read(*,*)
           read(*,*)
        endif

        read(*,*)
        read(*,*)  icalcorryz
        if(icalcorryz.gt.0) then
          read(*,*)
          read(*,*) corryz%num_kref
          allocate(corryz%kref(corryz%num_kref))
          read(*,*) (corryz%kref(i),i=1,corryz%num_kref)
          read(*,*)
          read(*,*) corryz%ibe, corryz%iend, corryz%kbe, corryz%kend
          read(*,*)
          corryz%jbe = 1
          corryz%jend = jmax
          corryz%p_row = 1                    ! 1 core in z direction
          corryz%p_col = numprocs             ! all cores in x direction (z pencil)
        else
          read(*,*)
          read(*,*)
          read(*,*)
          read(*,*)
          read(*,*)
          read(*,*)
        endif

        read(*,*)
        read(*,*) icalcorr3D
        if(icalcorr3D.gt.0) then
          read(*,*)
          read(*,*) corr3D%num_kref
          allocate(corr3D%kref(corr3D%num_kref))
          read(*,*) (corr3D%kref(i),i=1,corr3D%num_kref)
          read(*,*)
          read(*,*) corr3D%ibe, corr3D%iend, corr3D%iwinl, corr3D%iwinr, corr3D%kbe, corr3D%kend
          corr3D%jbe = 1
          corr3D%jend = jmax
!          if(numprocs.gt.(corr3D%iend-corr3D%ibe+1)) then
!            corr3D%p_row = numprocs/(corr3D%iend-corr3D%ibe+1)
!            corr3D%p_col = corr3D%iend-corr3D%ibe+1
!          else
            corr3D%p_row = numprocs
            corr3D%p_col = 1
!          endif
        else
          read(*,*)
          read(*,*)
          read(*,*)
          read(*,*)
          read(*,*)
          read(*,*)
        endif ! end icalcorr3D.gt.0

        read(*,*)
        read(*,*) icalWNSD_ky
        if(icalWNSD_ky.gt.0) then
          read(*,*)
          read(*,*) WNSD_ky%num_kref
          allocate(WNSD_ky%kref(WNSD_ky%num_kref))
          read(*,*) (WNSD_ky%kref(i),i=1,WNSD_ky%num_kref)
          read(*,*)
          read(*,*)  WNSD_ky%ibe, WNSD_ky%iend
          read(*,*)
          WNSD_ky%jbe  = 1
          WNSD_ky%jend = jmax
          WNSD_ky%kbe = 1
          WNSD_ky%kend = kmax
          WNSD_ky%p_row = numprocs
          WNSD_ky%p_col = 1
        else
          read(*,*)
          read(*,*)
          read(*,*)
          read(*,*)
          read(*,*)
          read(*,*)
        endif ! icalWNSD_ky.gt.0

        read(*,*)
        read(*,*) icalWNSD_kx, WNSD_kx%iwindow
        if(icalWNSD_kx.gt.0) then
          read(*,*)
          read(*,*) WNSD_kx%num_kref
          allocate(WNSD_kx%kref(WNSD_kx%num_kref))
          read(*,*) (WNSD_kx%kref(i),i=1,WNSD_kx%num_kref)
          read(*,*)
          read(*,*) WNSD_kx%ibe, WNSD_kx%iend
          WNSD_kx%jbe  = 1
          WNSD_kx%jend = jmax

          WNSD_kx%kbe = 1
          WNSD_kx%kend = kmax
          WNSD_kx%p_row = numprocs
          WNSD_kx%p_col = 1
        else
          read(*,*)
          read(*,*)
          read(*,*)
          read(*,*)
          read(*,*)
          read(*,*)
        endif ! icalWNSD_ky.gt.0

        if(num_output.gt.14) then
          print *, '# of output variables cannot be greater than 14. STOP ... '
          errcode = 1012
          call  MPI_Abort(MPI_COMM_WORLD,errcode,ierr)
        endif


      endif ! myid.eq.0

      call MPI_Bcast(gridfilename, 400, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(datapath,     400, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(groupname,    400, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(imax,      1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(jmax,      1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(kmax,      1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(iil,       1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(iih,       1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(stride,    1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(itimeave,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(Rm,        1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr) ! molecular
      call MPI_Bcast(num_output,1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      if(.not.allocated(varindex)) allocate(varindex(num_output))
      call MPI_Bcast(varindex,  num_output, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(icalcorrxz,      1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(icalcorrxy,      1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(icalcorryz,      1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(icalcorr3D,      1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(icalWNSD_ky,     1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(icalWNSD_kx,     1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

      if(icalcorrxz.gt.0) then
         call MPI_Bcast(corrxz%num_kref,          1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         if(.not.allocated(corrxz%kref)) allocate(corrxz%kref(corrxz%num_kref))
         call MPI_Bcast(corrxz%kref, corrxz%num_kref, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         call MPI_Bcast(corrxz%ibe,               1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         call MPI_Bcast(corrxz%iend,              1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         call MPI_Bcast(corrxz%jbe,               1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         call MPI_Bcast(corrxz%jend,              1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         call MPI_Bcast(corrxz%kbe,               1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         call MPI_Bcast(corrxz%kend,              1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         call MPI_Bcast(corrxz%iwinl,             1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         call MPI_Bcast(corrxz%iwinr,             1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         call MPI_Bcast(corrxz%p_row,             1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         call MPI_Bcast(corrxz%p_col,             1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

         if( (corrxz%ibe-corrxz%iwinl.lt.1).or.(corrxz%iend+corrxz%iwinr.gt.imax).or. &
            (corrxz%ibe.gt.corrxz%iend) ) then
            if(myid.eq.0) then
               print *, 'i index out of range for corr_xz!'
               print *, 'ibe =', corrxz%ibe, 'iend =', corrxz%iend
               print *, 'iwinl =', corrxz%iwinl, 'iwinr =', corrxz%iwinr
               print *, 'ibe-iwinl =', corrxz%ibe-corrxz%iwinl
               print *, 'iend+iwinr =', corrxz%iend+corrxz%iwinr
            endif
            errcode = 235
            call  MPI_Abort(MPI_COMM_WORLD,errcode,ierr)
         endif   

         if( any(corrxz%kref.lt.corrxz%kbe).or.any(corrxz%kref.gt.corrxz%kend).or. &
            (corrxz%kbe.lt.1).or.(corrxz%kend.gt.kmax) ) then
            if(myid.eq.0) then
               print *, 'k index out of range (kbe<=kref<=kend required)!'
               print *, 'kref =', corrxz%kref
               print *, 'kbe =', corrxz%kbe, 'kend =', corrxz%kend
            endif
            errcode = 246
            call  MPI_Abort(MPI_COMM_WORLD,errcode,ierr)
         endif
      endif

      if(icalcorryz.gt.0) then
         call MPI_Bcast(corryz%num_kref,          1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         if(.not.allocated(corryz%kref)) allocate(corryz%kref(corryz%num_kref))
         call MPI_Bcast(corryz%kref,corryz%num_kref, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         call MPI_Bcast(corryz%ibe,               1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         call MPI_Bcast(corryz%iend,              1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         call MPI_Bcast(corryz%jbe,               1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         call MPI_Bcast(corryz%jend,              1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         call MPI_Bcast(corryz%kbe,               1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         call MPI_Bcast(corryz%kend,              1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         !call MPI_Bcast(corryz%iwinl,             1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         !call MPI_Bcast(corryz%iwinr,             1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         call MPI_Bcast(corryz%p_row,             1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         call MPI_Bcast(corryz%p_col,             1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

         if( any(corryz%kref.lt.corryz%kbe).or.any(corryz%kref.gt.corryz%kend).or. &
            (corryz%kbe.lt.1).or.(corryz%kend.gt.kmax) ) then
            if(myid.eq.0) then
               print *, 'k index out of range (kbe<=kref<=kend required)!'
               print *, 'kref =', corryz%kref
               print *, 'kbe =', corryz%kbe, 'kend =', corryz%kend
            endif
            errcode = 246
            call  MPI_Abort(MPI_COMM_WORLD,errcode,ierr)
         endif

      endif


      if(icalcorrxy.gt.0) then
         call MPI_Bcast(corrxy%num_kref,          1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         if(.not.allocated(corrxy%kref)) allocate(corrxy%kref(corrxy%num_kref))
         call MPI_Bcast(corrxy%kref, corrxy%num_kref, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         call MPI_Bcast(corrxy%ibe,               1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         call MPI_Bcast(corrxy%iend,              1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         call MPI_Bcast(corrxy%jbe,               1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         call MPI_Bcast(corrxy%jend,              1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         call MPI_Bcast(corrxy%kbe,               1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         call MPI_Bcast(corrxy%kend,              1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         call MPI_Bcast(corrxy%iwinl,             1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         call MPI_Bcast(corrxy%iwinr,             1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         call MPI_Bcast(corrxy%p_row,             1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         call MPI_Bcast(corrxy%p_col,             1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

         if( (mod(numprocs,corrxy%num_kref).ne.0).and.(mod(corrxy%num_kref,numprocs).ne.0) ) then
             if(myid.eq.0) print *, 'number of processors/cores must be a multpiles/fraction of number of kplanes for corr_xy'
             errcode = 359
             call MPI_Abort(MPI_COMM_WORLD,errcode,ierr)
         endif


         if( (corrxy%ibe-corrxy%iwinl.lt.1).or.(corrxy%iend+corrxy%iwinr.gt.imax).or. &
            (corrxy%ibe.gt.corrxy%iend) ) then
            if(myid.eq.0) then
               print *, 'i index out of range for corr_xy!'
               print *, 'ibe =', corrxy%ibe, 'iend =', corrxy%iend
               print *, 'iwinl =', corrxy%iwinl, 'iwinr =', corrxy%iwinr
               print *, 'ibe-iwinl =', corrxy%ibe-corrxy%iwinl
               print *, 'iend+iwinr =', corrxy%iend+corrxy%iwinr
            endif
            errcode = 373
            call  MPI_Abort(MPI_COMM_WORLD,errcode,ierr)
         endif

         if( any(corrxy%kref.lt.corrxy%kbe).or.any(corrxy%kref.gt.corrxy%kend).or. &
            (corrxy%kbe.lt.1).or.(corrxy%kend.gt.kmax) ) then
            if(myid.eq.0) then
               print *, 'k index out of range for corr_xy (kbe<=kref<=kend required)!'
               print *, 'kref =', corrxy%kref
               print *, 'kbe =', corrxy%kbe, 'kend =', corrxy%kend
            endif
            errcode = 384
            call  MPI_Abort(MPI_COMM_WORLD,errcode,ierr)
         endif
      endif ! end icalcorrxy.gt.0

      if(icalcorr3D.gt.0) then
         call MPI_Bcast(corr3D%num_kref,          1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         if(.not.allocated(corr3D%kref)) allocate(corr3D%kref(corr3D%num_kref))
         call MPI_Bcast(corr3D%kref,corr3D%num_kref, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         call MPI_Bcast(corr3D%ibe,               1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         call MPI_Bcast(corr3D%iend,              1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         call MPI_Bcast(corr3D%jbe,               1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         call MPI_Bcast(corr3D%jend,              1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         call MPI_Bcast(corr3D%kbe,               1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         call MPI_Bcast(corr3D%kend,              1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         call MPI_Bcast(corr3D%iwinl,             1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         call MPI_Bcast(corr3D%iwinr,             1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         call MPI_Bcast(corr3D%p_row,             1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         call MPI_Bcast(corr3D%p_col,             1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

!         if( (mod(numprocs,(corr3D%iend-corr3D%ibe+1) ).ne.0).and.(mod((corr3D%iend-corr3D%ibe+1),numprocs).ne.0) ) then
!             if(myid.eq.0) print *, 'number of processors/cores must be a multpiles/fraction of (corr3D%iend - corr3D%ibe + 1) for corr3D'
!             errcode = 1220
!             call MPI_Abort(MPI_COMM_WORLD,errcode,ierr)
!         endif


         if( (corr3D%ibe-corr3D%iwinl.lt.1).or.(corr3D%iend+corr3D%iwinr.gt.imax).or. &
            (corr3D%ibe.gt.corr3D%iend) ) then
            if(myid.eq.0) then
               print *, 'i index out of range for corr_xy!'
               print *, 'ibe =', corr3D%ibe, 'iend =', corr3D%iend
               print *, 'iwinl =', corr3D%iwinl, 'iwinr =', corr3D%iwinr
               print *, 'ibe-iwinl =', corr3D%ibe-corr3D%iwinl
               print *, 'iend+iwinr =', corr3D%iend+corr3D%iwinr
            endif
            errcode = 373
            call  MPI_Abort(MPI_COMM_WORLD,errcode,ierr)
         endif

         if( any(corr3D%kref.lt.corr3D%kbe).or.any(corr3D%kref.gt.corr3D%kend).or. &
            (corr3D%kbe.lt.1).or.(corr3D%kend.gt.kmax) ) then
            if(myid.eq.0) then
               print *, 'k index out of range for corr_xy (kbe<=kref<=kend required)!'
               print *, 'kref =', corr3D%kref
               print *, 'kbe =', corr3D%kbe, 'kend =', corr3D%kend
            endif
            errcode = 384
            call  MPI_Abort(MPI_COMM_WORLD,errcode,ierr)
         endif
      endif ! end icalcorr3D.gt.0

      if(icalWNSD_ky.gt.0) then
         call MPI_Bcast(WNSD_ky%num_kref,          1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         if(.not.allocated(WNSD_ky%kref)) allocate(WNSD_ky%kref(WNSD_ky%num_kref))
         call MPI_Bcast(WNSD_ky%kref, WNSD_ky%num_kref, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         call MPI_Bcast(WNSD_ky%ibe,               1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         call MPI_Bcast(WNSD_ky%iend,              1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         call MPI_Bcast(WNSD_ky%jbe,               1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         call MPI_Bcast(WNSD_ky%jend,              1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         call MPI_Bcast(WNSD_ky%kbe,               1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         call MPI_Bcast(WNSD_ky%kend,              1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         call MPI_Bcast(WNSD_ky%p_row,             1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         call MPI_Bcast(WNSD_ky%p_col,             1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      endif

      if(icalWNSD_kx.gt.0) then
         call MPI_Bcast(WNSD_kx%num_kref,          1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         if(.not.allocated(WNSD_kx%kref)) allocate(WNSD_kx%kref(WNSD_kx%num_kref))
         call MPI_Bcast(WNSD_kx%kref, WNSD_kx%num_kref, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         call MPI_Bcast(WNSD_kx%ibe,               1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         call MPI_Bcast(WNSD_kx%iend,              1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         call MPI_Bcast(WNSD_kx%jbe,               1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         call MPI_Bcast(WNSD_kx%jend,              1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         call MPI_Bcast(WNSD_kx%kbe,               1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         call MPI_Bcast(WNSD_kx%kend,              1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         call MPI_Bcast(WNSD_kx%p_row,             1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         call MPI_Bcast(WNSD_kx%p_col,             1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         call MPI_Bcast(WNSD_kx%iinterp_x,         1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      endif

      rbar = R/Rm
      pi = 4.d0*atan(1.d0)

      return
    end subroutine Input


    subroutine SaveCorryz_Ascii(fn)
      character(*), intent(in) :: fn
      character(4) :: znum
      integer :: m, k, i, j
      real(8) :: ydiff(jmax), dy
      real(8), dimension(:,:,:,:), allocatable :: buffer_corrtmp, buffer_Rearrange
      integer :: jhalf, wavenum

      dy = yy(corryz%kbe,corryz%ibe,2) - yy(corryz%kbe,corryz%ibe,1)
      jhalf = jmax/2+1
      allocate(buffer_Rearrange(1:jmax,corryz%kbe:corryz%kend,1:nvarcorr_yz,corryz%num_kref))

      do j = 1, jmax
         wavenum = jhalf-1+j-jmax
         ydiff(j) = dble(wavenum)*dy
         if(j.le.jmax-jhalf) then
            buffer_Rearrange(j,:,:,:) = buffer_corr(jhalf+j,:,:,:)
         else
            buffer_Rearrange(j,:,:,:) = buffer_corr(j-jmax+jhalf,:,:,:)
         endif
      enddo

      open(21,file=trim(fn),status='unknown')
      rewind(21)

      write(21,'(a)') 'variables=ydiff,z,Corr_dim,varrms_ref,varrms'
      do m = 1, corryz%num_kref
         write(unit=znum,fmt='(I04.4)') corryz%kref(m)
         write(21,'("Zone T=kref_",I04.4, ",I=",I4, ",J=",I4 )') corryz%kref(m), jcorryz, kcorryz
         do k = corryz%kbe, corryz%kend
           do j = corryz%jbe, corryz%jend
             write(21,*) ydiff(j), zz(k,corryz%ibe,corryz%jbe), buffer_Rearrange(j,k,1:nvarcorr_yz,m)
           enddo
         enddo
      enddo
      close(21)

    end subroutine SaveCorryz_Ascii


    subroutine SaveCorrxz_Ascii(fn)
      character(*), intent(in) :: fn
      character(4) :: znum
      integer :: m, k, i
      real(8) :: xdiff

      allocate( buffer_rms1(-corrxz%iwinl:corrxz%iwinr, corrxz%kbe:corrxz%kend, nvarcorr_xz, corrxz%num_kref) )
      allocate( buffer_rms2(-corrxz%iwinl:corrxz%iwinr, corrxz%kbe:corrxz%kend, nvarcorr_xz) )

      do m = 1, corrxz%num_kref
        kloc = corrxz%kref(m)
        do k = corrxz%kbe, corrxz%kend
          do ii = -corrxz%iwinl,corrxz%iwinr
            buffer_rms1(ii,k,1:nvarcorr_xz,m) = sqrt( abs(buffer_var1(ii,kloc,1:nvarcorr_xz) - buffer_ave1(ii,kloc,1:nvarcorr_xz)**2) )
          enddo
        enddo
      enddo  ! end m loop

      do k = corrxz%kbe, corrxz%kend
        do ii = -corrxz%iwinl,corrxz%iwinr
          buffer_rms2(ii,k,1:nvarcorr_xz) = sqrt( abs(buffer_var2(ii,k,1:nvarcorr_xz) - buffer_ave2(ii,k,1:nvarcorr_xz)**2) )
        enddo
      enddo

      open(21,file=trim(fn),status='unknown')
      rewind(21)
      write(21,'(a)') 'variables=xdiff,z,Corr_dim,varrms,varrms_ref'
      do m = 1, corrxz%num_kref
         write(unit=znum,fmt='(I04.4)') corrxz%kref(m)
         write(21,'("Zone T=kref_",I04.4, ",I=",I4, ",J=",I4 )') corrxz%kref(m), icorrxz, kcorrxz

         do k = corrxz%kbe, corrxz%kend
         do i = -corrxz%iwinl,corrxz%iwinr
             xdiff = xx(k,iref+i,1)-xx(k,iref,1)
             write(21,*) xdiff, zz(k,iref,1), buffer_corr(i,k,1:nvarcorr_xz,m),&
                         buffer_rms2(i,k,1:nvarave), buffer_rms1(i,k,1:nvarave,m)
         enddo
         enddo
     enddo
     close(21)

     deallocate(buffer_rms1, buffer_rms2)

   end subroutine SaveCorrxz_Ascii

    subroutine SaveCorrxy_Ascii(fn)
      character(*), intent(in) :: fn
      integer :: k, kloc, i, j, m
      real(8) :: xdiff(-corrxy%iwinl:corrxy%iwinr), ydiff(jmax), dy
      integer, dimension(MPI_STATUS_SIZE) :: status
      real(8), dimension(:,:,:,:), allocatable :: buffer_corrtmp, buffer_Rearrange
      real(8), dimension(:,:), allocatable :: buffer_avetmp, buffer_vartmp
      real(8), dimension(:,:,:,:), allocatable :: buffer1_rms
      real(8), dimension(:,:,:), allocatable :: buffer2_rms
      integer :: kstid, kendid, klenid, istid
      integer :: jhalf, wavenum

      do i = -corrxy%iwinl,corrxy%iwinr
         xdiff(i) = xx(corrxy%decomp%zst(1),iref+i,1)-xx(corrxy%decomp%zst(1),iref,1)
      enddo

      dy = yy(corrxy%decomp%zst(1),ist_rd,2) - yy(corrxy%decomp%zst(1),ist_rd,1)
      jhalf = jmax/2+1
      allocate(buffer_Rearrange(-corrxy%iwinl:corrxy%iwinr,1:jmax,1:nvarcorr_xy,corrxy%decomp%zst(1):corrxy%decomp%zen(1)))
      do j = 1, jmax
         wavenum = jhalf-1+j-jmax
         ydiff(j) = dble(wavenum)*dy
         if(j.le.jmax-jhalf) then
            buffer_Rearrange(:,j,:,:) = buffer_corr(:,jhalf+j,:,:)
         else
            buffer_Rearrange(:,j,:,:) = buffer_corr(:,j-jmax+jhalf,:,:)
         endif
      enddo

      if(myid.eq.0) then
         open(22,file=trim(fn),status='unknown')
         rewind(22)

         write(22,'(a)') 'variables=xdiff,ydiff,Corr_dim,varrms_ref,varrms'
         do k = corrxy%decomp%zst(1), corrxy%decomp%zen(1)
            kloc = corrxy%kref(k)
            write(22,'("Zone T=kloc_",I04.4, ",I=",I4, ",J=",I4)') kloc, icorrxy, jmax
            do j = 1, jmax
               do i = -corrxy%iwinl,corrxy%iwinr
                  write(22,*) xdiff(i), ydiff(j), buffer_Rearrange(i,j,1:nvarcorr_xy,k)
               enddo
            enddo
         enddo

         do m = 1, numprocs-1
             call MPI_RECV(istid, 1,MPI_INTEGER,m, m+numprocs+1,MPI_COMM_WORLD,status,ierr)
             call MPI_RECV(kstid, 1,MPI_INTEGER,m, m+numprocs+2,MPI_COMM_WORLD,status,ierr)
             call MPI_RECV(kendid,1,MPI_INTEGER,m, m+numprocs+3,MPI_COMM_WORLD,status,ierr)
             call MPI_RECV(klenid,1,MPI_INTEGER,m, m+numprocs+4,MPI_COMM_WORLD,status,ierr)
             allocate(buffer_corrtmp(-corrxy%iwinl:corrxy%iwinr,1:jmax,1:nvarcorr_xy,kstid:kendid))
             call MPI_RECV(buffer_corrtmp,icorrxy*jmax*nvarcorr_xy*klenid,MPI_DOUBLE_PRECISION,m,m+numprocs,MPI_COMM_WORLD,status, ierr)

             if(istid.eq.1) then
                do k = kstid, kendid
                   kloc = corrxy%kref(k)
                   write(22,'("Zone T=kloc_",I04.4, ",I=",I4, ",J=",I4)') kloc, icorrxy, jmax
                   do j = 1, jmax
                      do i = -corrxy%iwinl,corrxy%iwinr
                         write(22,*) xdiff(i), ydiff(j), buffer_corrtmp(i,j,1:nvarcorr_xy,k)
                      enddo
                  enddo
                enddo
             endif
             deallocate(buffer_corrtmp)
         enddo
         close(22)
       else
            call MPI_SEND(corrxy%decomp%zst(2), 1, MPI_INTEGER, 0, myid+numprocs+1, MPI_COMM_WORLD, ierr)
            call MPI_SEND(corrxy%decomp%zst(1), 1, MPI_INTEGER, 0, myid+numprocs+2, MPI_COMM_WORLD, ierr)
            call MPI_SEND(corrxy%decomp%zen(1), 1, MPI_INTEGER, 0, myid+numprocs+3, MPI_COMM_WORLD, ierr)
            call MPI_SEND(corrxy%decomp%zsz(1), 1, MPI_INTEGER, 0, myid+numprocs+4, MPI_COMM_WORLD, ierr)
            call MPI_SEND(buffer_Rearrange, icorrxy*jmax*nvarcorr_xy*kcorrxy, MPI_DOUBLE_PRECISION, 0, myid+numprocs, MPI_COMM_WORLD, ierr)

       endif
       deallocate(buffer_Rearrange)

   end subroutine SaveCorrxy_Ascii

   subroutine SaveWNSD_ky_Ascii(fn)
     character(*), intent(in) :: fn
      integer :: k, kloc, i, j, m
      integer, dimension(MPI_STATUS_SIZE) :: status
      integer :: kstid, kendid, klenid, istid
      integer :: jhalf, wavenum
      real(8) :: ff, fperiod, dy
      real(8), dimension(:,:,:), allocatable :: spect1d_wt, spect1d_wttmp, buffer_var1_tmp
      real(8), dimension(:), allocatable :: tmpenergy, tmpenergy_tmp

      dy = yy(WNSD_ky%decomp%zst(1),ist_rd,2) - yy(WNSD_ky%decomp%zst(1),ist_rd,1)
      fperiod = dble(jmax)*dy

      allocate( tmpenergy(WNSD_ky%decomp%zst(1):WNSD_ky%decomp%zen(1)) )
      allocate(   spect1d_wt(WNSD_ky%decomp%zst(1):WNSD_ky%decomp%zen(1), 1, 1:jmax/2+1) )
      spect1d_wt(WNSD_ky%decomp%zst(1):WNSD_ky%decomp%zen(1), 1, 1)          =      real( spect1dave(WNSD_ky%decomp%zst(1):WNSD_ky%decomp%zen(1), 1, 1) )*fperiod
      spect1d_wt(WNSD_ky%decomp%zst(1):WNSD_ky%decomp%zen(1), 1, 2:jmax/2+1) = 2.d0*real( spect1dave(WNSD_ky%decomp%zst(1):WNSD_ky%decomp%zen(1), 1, 2:jmax/2+1) )*fperiod
      do k=WNSD_ky%decomp%zst(1), WNSD_ky%decomp%zen(1)
        tmpenergy(k) = sum(spect1d_wt(k,1,1:jmax/2+1))/fperiod
      enddo

      if(myid.eq.0) then
         open(22,file=trim(fn),status='unknown')
         rewind(22)
         write(22,'(a)') 'variables=ky, spectrum, rms2'
         do k = WNSD_ky%decomp%zst(1), WNSD_ky%decomp%zen(1)
            kloc = WNSD_ky%kref(k)
            write(22,'("Zone T=kloc_",I04.4, ",I=",I4)') kloc, jmax/2+1
            do j = 1, jmax/2+1
              ff = dble(j-1)/fperiod*2.d0*pi
              write(22,*) ff, spect1d_wt(k,1,j)*buffer_var1(k,1,1)/tmpenergy(k)/(2.d0*pi), buffer_var1(k,1,1)
            enddo
         enddo

         do m = 1, numprocs-1
             call MPI_RECV(istid, 1,MPI_INTEGER,m, m+numprocs+1,MPI_COMM_WORLD,status,ierr)
             call MPI_RECV(kstid, 1,MPI_INTEGER,m, m+numprocs+2,MPI_COMM_WORLD,status,ierr)
             call MPI_RECV(kendid,1,MPI_INTEGER,m, m+numprocs+3,MPI_COMM_WORLD,status,ierr)
             call MPI_RECV(klenid,1,MPI_INTEGER,m, m+numprocs+4,MPI_COMM_WORLD,status,ierr)
             allocate(buffer_var1_tmp(kstid:kendid,1,1))
             call MPI_RECV(buffer_var1_tmp,kWNSD_ky,MPI_DOUBLE_PRECISION,m, m+numprocs+5,MPI_COMM_WORLD,status,ierr)
             allocate(tmpenergy_tmp(kstid:kendid))
             call MPI_RECV(tmpenergy_tmp,kWNSD_ky,MPI_DOUBLE_PRECISION,m, m+numprocs+6,MPI_COMM_WORLD,status,ierr)
             allocate(spect1d_wttmp(kstid:kendid,1,1:jmax/2+1))
             call MPI_RECV(spect1d_wttmp,(jmax/2+1)*klenid,MPI_DOUBLE_PRECISION,m,m+numprocs,MPI_COMM_WORLD,status, ierr)
             if(istid.eq.1) then
                do k = kstid, kendid
                   kloc = WNSD_ky%kref(k)
                   write(22,'("Zone T=kloc_",I04.4, ",I=",I4)') kloc, jmax/2+1
                   do j = 1, jmax/2+1
                     ff = dble(j-1)/fperiod*2.d0*pi
                      write(22,*) ff, spect1d_wttmp(k,1,j)*buffer_var1_tmp(k,1,1)/tmpenergy_tmp(k)/(2.d0*pi), buffer_var1_tmp(k,1,1)
                   enddo
                enddo
             endif
             deallocate(spect1d_wttmp)
         enddo
         close(22)
       else
            call MPI_SEND(WNSD_ky%decomp%zst(2), 1, MPI_INTEGER, 0, myid+numprocs+1, MPI_COMM_WORLD, ierr)
            call MPI_SEND(WNSD_ky%decomp%zst(1), 1, MPI_INTEGER, 0, myid+numprocs+2, MPI_COMM_WORLD, ierr)
            call MPI_SEND(WNSD_ky%decomp%zen(1), 1, MPI_INTEGER, 0, myid+numprocs+3, MPI_COMM_WORLD, ierr)
            call MPI_SEND(WNSD_ky%decomp%zsz(1), 1, MPI_INTEGER, 0, myid+numprocs+4, MPI_COMM_WORLD, ierr)
            call MPI_SEND(buffer_var1, kWNSD_ky, MPI_DOUBLE_PRECISION, 0, myid+numprocs+5, MPI_COMM_WORLD, ierr)
            call MPI_SEND(  tmpenergy, kWNSD_ky, MPI_DOUBLE_PRECISION, 0, myid+numprocs+6, MPI_COMM_WORLD, ierr)
            call MPI_SEND(spect1d_wt, (jmax/2+1)*kWNSD_ky, MPI_DOUBLE_PRECISION, 0, myid+numprocs, MPI_COMM_WORLD, ierr)
       endif


   end subroutine SaveWNSD_ky_Ascii

   subroutine SaveWNSD_kx_Ascii(fn)
     character(*), intent(in) :: fn
      integer :: k, kloc, i, j, m
      integer, dimension(MPI_STATUS_SIZE) :: status
      integer :: kstid, kendid, klenid, istid
      integer :: jhalf, wavenum
      real(8) :: ff, fperiod, dx
      real(8), dimension(:,:,:), allocatable :: spect1d_wt, spect1d_wttmp, buffer_var1_tmp
      real(8), dimension(:), allocatable :: tmpenergy, tmpenergy_tmp

      dx = xx(WNSD_kx%decomp%yst(1),ist_rd+1,1) - xx(WNSD_kx%decomp%yst(1),ist_rd,1)
      fperiod = dble(ilen_rd)*dx

      allocate( tmpenergy(WNSD_kx%decomp%yst(1):WNSD_kx%decomp%yen(1)) )
      allocate(   spect1d_wt(WNSD_kx%decomp%yst(1):WNSD_kx%decomp%yen(1), 1:ilen_rd/2+1, 1) )
      spect1d_wt(WNSD_kx%decomp%yst(1):WNSD_kx%decomp%yen(1), 1, 1)          =         real( spect1dave(WNSD_kx%decomp%yst(1):WNSD_kx%decomp%yen(1), 1, 1) ) !*fperiod
      spect1d_wt(WNSD_kx%decomp%yst(1):WNSD_kx%decomp%yen(1), 2:ilen_rd/2+1, 1) = 2.d0*real( spect1dave(WNSD_kx%decomp%yst(1):WNSD_kx%decomp%yen(1), 2:ilen_rd/2+1, 1) ) !*fperiod
      do k=WNSD_kx%decomp%yst(1), WNSD_kx%decomp%yen(1)
        tmpenergy(k) = sum(spect1d_wt(k,1:ilen_rd/2+1,1))/fperiod
      enddo

      if(myid.eq.0) then
         open(22,file=trim(fn),status='unknown')
         rewind(22)
         write(22,'(a)') 'variables=kx, spectrum, rms2'
         do k = WNSD_kx%decomp%yst(1), WNSD_kx%decomp%yen(1)
            kloc = WNSD_kx%kref(k)
            write(22,'("Zone T=kloc_",I04.4, ",I=",I4)') kloc, ilen_rd/2+1
            do i = 1, ilen_rd/2+1
              ff = dble(i-1)/fperiod*2.d0*pi
              write(22,*) ff, spect1d_wt(k,i,1)*buffer_var1(k,1,1)/tmpenergy(k)/(2.d0*pi), buffer_var1(k,1,1)
            !  write(22,*) ff, spect1d_wt(k,i,1), buffer_var1(k,1,1)
            enddo
         enddo

         do m = 1, numprocs-1
             call MPI_RECV(istid, 1,MPI_INTEGER,m, m+numprocs+1,MPI_COMM_WORLD,status,ierr)
             call MPI_RECV(kstid, 1,MPI_INTEGER,m, m+numprocs+2,MPI_COMM_WORLD,status,ierr)
             call MPI_RECV(kendid,1,MPI_INTEGER,m, m+numprocs+3,MPI_COMM_WORLD,status,ierr)
             call MPI_RECV(klenid,1,MPI_INTEGER,m, m+numprocs+4,MPI_COMM_WORLD,status,ierr)
             allocate(buffer_var1_tmp(kstid:kendid,1,1))
             call MPI_RECV(buffer_var1_tmp,kWNSD_kx,MPI_DOUBLE_PRECISION,m, m+numprocs+5,MPI_COMM_WORLD,status,ierr)
             allocate(tmpenergy_tmp(kstid:kendid))
             call MPI_RECV(tmpenergy_tmp,kWNSD_kx,MPI_DOUBLE_PRECISION,m, m+numprocs+6,MPI_COMM_WORLD,status,ierr)
             allocate(spect1d_wttmp(kstid:kendid,1:ilen_rd/2+1,1))
             call MPI_RECV(spect1d_wttmp,(ilen_rd/2+1)*klenid,MPI_DOUBLE_PRECISION,m,m+numprocs,MPI_COMM_WORLD,status, ierr)
             if(istid.eq.1) then
                do k = kstid, kendid
                   kloc = WNSD_kx%kref(k)
                   write(22,'("Zone T=kloc_",I04.4, ",I=",I4)') kloc, ilen_rd/2+1
                   do i = 1, ilen_rd/2+1
                     ff = dble(i-1)/fperiod*2.d0*pi
                      write(22,*) ff, spect1d_wttmp(k,i,1)*buffer_var1_tmp(k,1,1)/tmpenergy_tmp(k)/(2.d0*pi), buffer_var1_tmp(k,1,1)
                   enddo
                enddo
             endif
             deallocate(spect1d_wttmp,buffer_var1_tmp)
         enddo
         close(22)
       else
            call MPI_SEND(WNSD_kx%decomp%yst(2), 1, MPI_INTEGER, 0, myid+numprocs+1, MPI_COMM_WORLD, ierr)
            call MPI_SEND(WNSD_kx%decomp%yst(1), 1, MPI_INTEGER, 0, myid+numprocs+2, MPI_COMM_WORLD, ierr)
            call MPI_SEND(WNSD_kx%decomp%yen(1), 1, MPI_INTEGER, 0, myid+numprocs+3, MPI_COMM_WORLD, ierr)
            call MPI_SEND(WNSD_kx%decomp%ysz(1), 1, MPI_INTEGER, 0, myid+numprocs+4, MPI_COMM_WORLD, ierr)
            call MPI_SEND(buffer_var1, kWNSD_kx, MPI_DOUBLE_PRECISION, 0, myid+numprocs+5, MPI_COMM_WORLD, ierr)
            call MPI_SEND(  tmpenergy, kWNSD_kx, MPI_DOUBLE_PRECISION, 0, myid+numprocs+6, MPI_COMM_WORLD, ierr)
            call MPI_SEND(spect1d_wt, (ilen_rd/2+1)*kWNSD_kx, MPI_DOUBLE_PRECISION, 0, myid+numprocs, MPI_COMM_WORLD, ierr)
       endif

   end subroutine SaveWNSD_kx_Ascii

    subroutine SaveCorr3D_Ascii(fn)
      character(*), intent(in) :: fn
      integer :: k, kloc, i, j, m
      real(8) :: xdiff(-corr3D%iwinl:corr3D%iwinr), ydiff(jmax), dy
      integer, dimension(MPI_STATUS_SIZE) :: status
      real(8), dimension(:,:,:,:), allocatable :: buffer_corrtmp, buffer_Rearrange
      real(8), dimension(:), allocatable :: zzz
      integer :: kstid, kendid, klenid, istid
      integer :: jhalf, wavenum

      do i = -corr3D%iwinl,corr3D%iwinr
         xdiff(i) = xx(corr3D%decomp%zst(1),iref+i,1)-xx(corr3D%decomp%zst(1),iref,1)
      enddo

      dy = yy(corr3D%decomp%zst(1),ist_rd,2) - yy(corr3D%decomp%zst(1),ist_rd,1)
      jhalf = jmax/2+1
      allocate(buffer_Rearrange(-corr3D%iwinl:corr3D%iwinr,1:jmax,1:nvarcorr_3D,corr3D%decomp%zst(1):corr3D%decomp%zen(1)))
      do j = 1, jmax
         wavenum = jhalf-1+j-jmax
         ydiff(j) = dble(wavenum)*dy
         if(j.le.jmax-jhalf) then
            buffer_Rearrange(:,j,:,:) = buffer_corr(:,jhalf+j,:,:)
         else
            buffer_Rearrange(:,j,:,:) = buffer_corr(:,j-jmax+jhalf,:,:)
         endif
      enddo

      if(myid.eq.0) then
         open(22,file=trim(fn),status='unknown')
         rewind(22)
         write(22,'(a)') 'variables=xdiff,ydiff,z,Corr_dim,varrms_ref,varrms'
         kloc = kstart ! ref k loc
         write(22,'("Zone T=kloc_",I04.4, ",I=",I4, ",J=",I4, ",K=", I4)') kloc, icorr3D, jmax_rd, kmax_rd
         do k = corr3D%decomp%zst(1), corr3D%decomp%zen(1)
           !kloc = kstart ! ref k loc
           !write(22,'("Zone T=kloc_",I04.4, ",I=",I4, ",J=",I4, ",K=", I4)') kloc, icorr3D, jmax_rd, kmax_rd
           do j = 1, jmax_rd
             do i = -corr3D%iwinl,corr3D%iwinr
               write(22,*) xdiff(i), ydiff(j), zz(k,ist_rd,1), buffer_Rearrange(i,j,1:nvarcorr_3D,k)
             enddo
           enddo
         enddo

         do m = 1, numprocs-1
             call MPI_RECV(istid, 1,MPI_INTEGER,m, m+numprocs+1,MPI_COMM_WORLD,status,ierr)
             call MPI_RECV(kstid, 1,MPI_INTEGER,m, m+numprocs+2,MPI_COMM_WORLD,status,ierr)
             call MPI_RECV(kendid,1,MPI_INTEGER,m, m+numprocs+3,MPI_COMM_WORLD,status,ierr)
             call MPI_RECV(klenid,1,MPI_INTEGER,m, m+numprocs+4,MPI_COMM_WORLD,status,ierr)
             allocate(buffer_corrtmp(-corr3D%iwinl:corr3D%iwinr,1:jmax,1:nvarcorr_3D,kstid:kendid))
             allocate(zzz(kstid:kendid))
             call MPI_RECV(buffer_corrtmp,icorr3D*jmax*nvarcorr_3D*klenid,MPI_DOUBLE_PRECISION,m,m+numprocs,MPI_COMM_WORLD,status, ierr)
             call MPI_RECV(zzz,klenid,MPI_DOUBLE_PRECISION,m,m+numprocs+5,MPI_COMM_WORLD,status, ierr)
             if(istid.eq.1) then
                do k = kstid, kendid
                   do j = 1, jmax_rd
                      do i = -corr3D%iwinl,corr3D%iwinr
                         write(22,*) xdiff(i), ydiff(j), zzz(k), buffer_corrtmp(i,j,1:nvarcorr_3D,k)
                      enddo
                  enddo
                enddo
             endif
             deallocate(buffer_corrtmp)
             deallocate(zzz)
         enddo
         close(22)
       else
            call MPI_SEND(corr3D%decomp%zst(2), 1, MPI_INTEGER, 0, myid+numprocs+1, MPI_COMM_WORLD, ierr)
            call MPI_SEND(corr3D%decomp%zst(1), 1, MPI_INTEGER, 0, myid+numprocs+2, MPI_COMM_WORLD, ierr)
            call MPI_SEND(corr3D%decomp%zen(1), 1, MPI_INTEGER, 0, myid+numprocs+3, MPI_COMM_WORLD, ierr)
            call MPI_SEND(corr3D%decomp%zsz(1), 1, MPI_INTEGER, 0, myid+numprocs+4, MPI_COMM_WORLD, ierr)
            call MPI_SEND(buffer_Rearrange, icorr3D*jmax*nvarcorr_3D*kcorr3D, MPI_DOUBLE_PRECISION, 0, myid+numprocs, MPI_COMM_WORLD, ierr)
            call MPI_SEND(zz(:,ist_rd,1), corr3D%decomp%zsz(1), MPI_DOUBLE_PRECISION, 0, myid+numprocs+5, MPI_COMM_WORLD, ierr)

       endif
       deallocate(buffer_Rearrange)

   end subroutine SaveCorr3D_Ascii


  !subroutine to rearrange correlation output computed by Wiener-Khinchin Theorem (that makes use of FFT)
  !so that the zero-delay (ydiff = 0)  is at the center of the series.
  !np (input): number of points in the correlation along the direction computed using FFT
  !corrin (input): the correlation output from 1D FFT
  !corryout (output): the rearranged correlation
  !wavenum (output): the corresponding wave number of the output spectrum
  subroutine RearrangeFFT1D(np,corrin,corrout,wavenum)
     integer, intent(in) :: np
     real(8), dimension(np), intent(inout) :: corrin
     real(8), dimension(np), intent(out) :: corrout
     integer, dimension(np), intent(out) :: wavenum
     integer :: ihalf, i
     ihalf = np/2+1
     do i=1,np
        wavenum(i) = ihalf-1+i-np
        if (i.le.np-ihalf) then
           corrout(i) = corrin(ihalf+i)
        else
           corrout(i) = corrin(i-np+ihalf)
        end if
     end do
  end subroutine RearrangeFFT1D

end program correlation
