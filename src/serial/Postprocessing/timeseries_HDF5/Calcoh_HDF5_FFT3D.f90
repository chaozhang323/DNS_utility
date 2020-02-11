program Calcoh_FFT3D
use MTSKplane
use modSpect
use modTecbin
use MSYSTIME
use MFFT2D
implicit none
real(8),parameter :: pi = 4.d0*atan(1.)
integer,parameter :: npath = 1
integer,parameter :: nvar = 1, index_var(1) = 4  ! read pressure, p
! input
type(fprop), dimension(:), allocatable :: fileprop
type(tp_DNSIndex) :: DNSIndex_k
integer :: nsection, nsection_kx, ioverlap, iwindow, ntpoint_total, &
           ibe_ave, iend_ave, jbe_ave, jend_ave, ixwindow_half, num_kplane, icylin
integer,dimension(:),allocatable :: kplane
real(8) :: dx_kp, dy_kp, dt_sample, radius = 1.0
real(8),dimension(:,:,:,:),allocatable :: tsdata

call InitHDF5()
call Input()

call main()

contains

subroutine Input()
    logical :: PlaneNameMatch = .false.
    integer :: n, k, kk

    allocate(fileprop(npath))
    ! input parameters
    read(*,*)
    read(*,*) iFileType, icylin
    read(*,*)
    do n=1, npath
        read(*,*)
        read(*,*) fileprop(n)%file_be, fileprop(n)%file_end, fileprop(n)%file_skip
        read(*,*)
        read(*,'(a)') fileprop(n)%filepath
    enddo
    read(*,*)
    read(*,*)
    read(*,*) nsection, nsection_kx, ioverlap, iwindow
    read(*,*)
    read(*,*) ibe_ave, iend_ave, jbe_ave, jend_ave, ixwindow_half
    read(*,*)
    read(*,*) num_kplane
    allocate(kplane(num_kplane))
    read(*,*) (kplane(n), n = 1, num_kplane)
    if(iFileType.eq.0) then
        print *, 'Reading DNS_index information from input file ... '
        read(*,*)
        read(*,*)
        read(*,*) fileprop(1)%ntpoint, fileprop(1)%nskip, dx_kp, dy_kp, dt_sample
        read(*,*)
        read(*,*) DNSIndex_k%ibe, DNSIndex_k%iend, DNSIndex_k%iskip, &
                  DNSIndex_k%jbe, DNSIndex_k%jend, DNSIndex_k%jskip, DNSIndex_k%ibuffer
    else
        call ReadDNS_index_kplane(fileprop(1)%ntpoint, dx_kp, dy_kp, dt_sample, DNSIndex_k,fileprop(1)%filepath)
        fileprop(1)%nskip = 0
    endif
    ! end of reading input file

    if(num_kplane.ne.1) then
        print*,'num_kplane MUST be 1. STOP!'
        stop
    endif

    ! check kplane location index
    if(iFileType.ne.0) then
      if(num_kplane.le.DNSIndex_k%nkplane) then
        do k=1, num_kplane
          do kk=1, DNSIndex_k%nkplane
            if(kplane(k).eq.DNSIndex_k%kplane(kk)) then
              PlaneNameMatch = .true.
            endif
          enddo ! end kk loop
          if(.not.PlaneNameMatch) then
            print *, '#############################################################'
            print *, 'The input kplane name does not match the name in DNSIndex.h5 '
            print *, 'Input kplane name ', kplane(k) ,' does not exist in file DNSIndex.h5 '
            print *, 'Stop ... '
            stop
          else
            PlaneNameMatch = .false.
          endif
        enddo ! end k loop
      else ! num_kplane.gt.DNSIndex_k%nkplane
        print *, 'The input num_kplane does not match the num_kplane in DNSIndex.h5. '
        print *, 'The input num_kplane = ', num_kplane, 'num_kplane in DNSIndex.h5 = ', DNSIndex_k%nkplane
        stop
      endif
    endif

    ! calculate the total number of time point
    ntpoint_total = 0
    if(iFileType.eq.0) then
      do n=1, npath
        fileprop(n)%ntpoint = fileprop(n)%ntpoint - fileprop(n)%nskip
        ntpoint_total = ntpoint_total + fileprop(n)%ntpoint
      enddo
    elseif(iFileType.ge.1) then
      do n=1, npath
        ntpoint_total = ntpoint_total + ( (fileprop(n)%file_end-fileprop(n)%file_be)/fileprop(n)%file_skip + 1 )*fileprop(n)%ntpoint
      enddo
    endif

    if(ioverlap.ne.0.and.ioverlap.ne.1) then
      print*,'ioverlap can ONLY be 0 or 1'
      stop
    end if
    print *, 'number of wall-normal locations =', num_kplane
    print *, 'Wall-normal plane indexes k =', kplane
    print *, 'ntpoint_total =', ntpoint_total, 'dt_sample = ', dt_sample
    print *, 'half ixwindow size =', ixwindow_half

end subroutine Input

subroutine init_data3d()
    integer :: n1_data, n2_data, n3_data, k, kk, kplane_be
    real(8),dimension(:,:,:,:,:),allocatable :: tsdata_tmp
    if(iFileType.eq.0) then
        call InitTSkplane(ntpoint_total, npath, DNSIndex_k, ibe_ave, iend_ave, jbe_ave, jend_ave, &
                        ixwindow_half, ixwindow_half, nvar, index_var, fileprop)
    elseif(iFileType.eq.1 .or. iFileType.eq.2) then
        call InitVOLTSkplane(ntpoint_total, npath, DNSIndex_k, ibe_ave, iend_ave, jbe_ave, jend_ave, &
                           ixwindow_half, ixwindow_half, nvar, index_var, fileprop)
    else
        print *, "iFileType should be 0 or 1, 2"
        stop
    endif

    ! read data
    n1_data=ntpoint_total; n2_data=nypoint_k; n3_data=imaxts_k-imints_k+1
    allocate(tsdata_tmp(ntpoint_total,nypoint_k,imints_k:imaxts_k,1,1),tsdata(n1_data,n2_data,n3_data,num_kplane))
    do k = 1,num_kplane
        do kk = 1,DNSIndex_k%nkplane
            if(kplane(k).eq.DNSIndex_k%kplane(kk)) then
                kplane_be = kk
                call ReadVOLTSkplane(kplane_be,kplane(k),tsdata_tmp)
                tsdata(:,:,:,k) = tsdata_tmp(:,:,:,1,1)
            endif
        enddo
    enddo
    deallocate(tsdata_tmp)

    if(icylin.ge.1) radius = DNSIndex_k%zloc(kplane_be)
end subroutine init_data3d

subroutine main()
    integer :: i,j,k
    integer :: n1_fft0,n3_fft0,n1_fft,n2_fft,n3_fft, nsec1,nsec3, &
               ibe_fft,iend_fft,idim_fft
    real(8) :: dx1,dx2,dx3, lx1,lx2,lx3
!    real(8),dimension(:,:,:),allocatable :: wavefreq_spectra
    complex,dimension(:,:,:),allocatable :: wavefreq_spectra
    complex,dimension(:,:,:),allocatable :: cross_spectra
    real(8),dimension(:),allocatable :: auto_spectra
    integer,dimension(:),allocatable :: wave1,wave2,wave3
    character(4) :: fnum_kloc
    character(500) :: tecname, file_out
    integer :: c1_t, c2_t, TID, NTHREADS
    real(8) :: rate

    call InitSysTime(rate)

    print*,' '
    print*,'read 3D data ...'
    call system_clock(c1_t)
    ! prepare the 3D data: tsdata(ntpoint_total,nypoint_k,imaxts_k-imints_k+1)
    call init_data3d()
!    print*,tsdata(1:10,1,1,1)
    call system_clock(c2_t)
    print*,'    read data time [min]: ',dble(c2_t-c1_t)/rate/60.d0

    print*,'initcoherence ...'
    write(unit=fnum_kloc,fmt='(I04.4)') kplane(1)
    ! initialize FFTpack subroutines
    n1_fft0 = ntpoint_total; n2_fft = nypoint_k; n3_fft0 = 2*ixwindow_half
    nsec1 = nsection; nsec3 = nsection_kx
    dx1 = dt_sample; dx2 = dy_kp; dx3 = dx_kp
    call InitCoherence(n1_fft0,n2_fft,n3_fft0, nsec1,nsec3, iwindow,ioverlap, dx1,dx2,dx3, n1_fft,n3_fft)
    print*,'n1_fft,n2_fft,n3_fft: ',n1_fft,n2_fft,n3_fft
    if(icylin.ge.1) then
        print*,'cylindrical coordinate, radius at kref_'//fnum_kloc//': ',radius
        dx2 = dx2*radius
    endif
    print*,'dt, dy, dx: ', dx1,dx2,dx3

    lx1 = dble(n1_fft-1)*dx1
    lx2 = dble(n2_fft-1)*dx2
    lx3 = dble(n3_fft-1)*dx3

    print*,'calculate wavenumber frequency spectra ...'
    ! wavenumber-frequency spectra by FFT3DF
    ibe_fft = 1; iend_fft = n3_fft0; idim_fft = 0
    allocate(wavefreq_spectra(n1_fft,n2_fft,n3_fft))
    wavefreq_spectra = cmplx(0.d0,0.d0)
!    wavefreq_spectra = 0.d0
    call system_clock(c1_t)
    do i=ibe_ave,iend_ave
        idim_fft = idim_fft+1
        wavefreq_spectra = wavefreq_spectra+spectra3d_ave(tsdata(:,:,ibe_fft:iend_fft,1))
!        wavefreq_spectra = wavefreq_spectra+real(spectra3d_ave(tsdata(:,:,ibe_fft:iend_fft,1)))
        ibe_fft = ibe_fft+1; iend_fft = iend_fft+1
        if(iend_fft.gt.nxpoint_k) exit
    enddo
    wavefreq_spectra = wavefreq_spectra/dble(idim_fft)/(2.d0*pi)**3
    deallocate(tsdata)
    call system_clock(c2_t)
    print*,'    wavenumber frequency spectra time [min]: ',dble(c2_t-c1_t)/rate/60.d0

    print*,'calculate cross-spectrum ...'
    ! cross-spectrum by FFT2DB
    allocate(cross_spectra(n1_fft,n2_fft,n3_fft))
    call system_clock(c1_t)
    call InitFFT2D(n2_fft,n3_fft)
    do i=1,n1_fft
        cross_spectra(i,:,:) = cross_spectra_fft2db(wavefreq_spectra(i,:,:))
    end do
    call system_clock(c2_t)
    print*,'    cross_spectra time [min]: ',dble(c2_t-c1_t)/rate/60.d0

    ! write wavenumber_frequency spectrum
    allocate(wave1(n1_fft),wave2(n2_fft),wave3(n3_fft))
!    call ShiftFFT3D_real(n1_fft,n2_fft,n3_fft,wavefreq_spectra,wavefreq_spectra,wave1,wave2,wave3)
    call ShiftFFT3D(n1_fft,n2_fft,n3_fft,wavefreq_spectra,wavefreq_spectra,wave1,wave2,wave3)

    ! write to tecplot binary
    file_out = 'wavefreq2d_kx_omeg_kref'//fnum_kloc//'.plt'
    print*,'write ',trim(file_out)
    tecname = 'kx, omega, wavefreq_spectra'
    call InitTec(1,n3_fft,n1_fft,1,3,tecname,0)
    forall(i=1:n3_fft) vartmp(i,:,1,1) = 2.*pi*dble(wave3(i))/lx3
    forall(j=1:n1_fft) vartmp(:,j,1,2) = 2.*pi*dble(wave1(j))/lx1
!    forall(i=1:n3_fft,j=1:n1_fft) vartmp(i,j,1,3) = wavefreq_spectra(j,n2_fft/2+1,i)
    forall(i=1:n3_fft,j=1:n1_fft) vartmp(i,j,1,3) = real(wavefreq_spectra(j,n2_fft/2+1,i))
    call WriteTec(trim(file_out))

    file_out = 'wavefreq3d_kref'//fnum_kloc//'.plt'
    print*,'write ',trim(file_out)
    tecname = 'kx, ky, omega, wavefreq_spectra'
    call InitTec(1,n3_fft,n2_fft,n1_fft,4,tecname,0)
    forall(i=1:n3_fft) vartmp(i,:,:,1) = 2.*pi*dble(wave3(i))/lx3
    forall(j=1:n2_fft) vartmp(:,j,:,2) = 2.*pi*dble(wave2(j))/lx2
    forall(k=1:n1_fft) vartmp(:,:,k,3) = 2.*pi*dble(wave1(k))/lx1
!    forall(i=1:n3_fft,j=1:n2_fft,k=1:n1_fft) vartmp(i,j,k,4) = wavefreq_spectra(k,j,i)
    forall(i=1:n3_fft,j=1:n2_fft,k=1:n1_fft) vartmp(i,j,k,4) = real(wavefreq_spectra(k,j,i))
    call WriteTec(trim(file_out))
    if(allocated(vartmp)) deallocate(vartmp)
    deallocate(wavefreq_spectra,wave1,wave2,wave3)

    print*,'shift 3D spectra ...'
    allocate(wave1(n1_fft),wave2(n2_fft),wave3(n3_fft))
    call ShiftFFT3D(n1_fft,n2_fft,n3_fft,cross_spectra,cross_spectra,wave1,wave2,wave3)

    allocate(auto_spectra(n1_fft))
    auto_spectra = sqrt( real(cross_spectra(:,n2_fft/2+1,n3_fft/2+1))**2+imag(cross_spectra(:,n2_fft/2+1,n3_fft/2+1))**2 )

    ! write to ascii
!    open(33,file='coherence_test.dat')
!    write(33,'(a)') 'Variables = f, y, x, spect12_r, spect12_i, spect_ref'
!    write(33,*) 'Zone T = k'//fnum_kloc//', i = ',n1_fft, ', j = ', n2_fft, ', k = ', n3_fft
!    do k=1,n3_fft
!    do j=1,n2_fft
!    do i=1,n1_fft
!        write(33,*) wave1(i)/lx1,wave2(j)*dx2,wave3(k)*dx3,real(cross_spectra(i,j,k)),imag(cross_spectra(i,j,k)), auto_spectra(i)
!    enddo
!    enddo
!    enddo
!    close(33)
    ! write to tecplot binary
    file_out = 'coherence_dy_kref'//fnum_kloc//'.plt'
    print*,'write ',trim(file_out)
    tecname = 'f, y, spect12_r, spect12_i, spect_ref'
    call InitTec(1,n1_fft,n2_fft,1,5,tecname,0)
    forall(i=1:n1_fft) vartmp(i,:,1,1) = dble(wave1(i))/lx1
    forall(j=1:n2_fft) vartmp(:,j,1,2) = dble(wave2(j))*dx2
    forall(i=1:n1_fft,j=1:n2_fft) vartmp(i,j,1,3) = real(cross_spectra(i,j,n3_fft/2+1))
    forall(i=1:n1_fft,j=1:n2_fft) vartmp(i,j,1,4) = imag(cross_spectra(i,j,n3_fft/2+1))
    forall(i=1:n1_fft) vartmp(i,:,1,5) = auto_spectra(i)
    call WriteTec(trim(file_out))

    file_out = 'coherence_dx_kref'//fnum_kloc//'.plt'
    print*,'write ',trim(file_out)
    tecname = 'f, x, spect12_r, spect12_i, spect_ref'
    call InitTec(1,n1_fft,1,n3_fft,5,tecname,0)
    forall(i=1:n1_fft) vartmp(i,1,:,1) = dble(wave1(i))/lx1
    forall(k=1:n3_fft) vartmp(:,1,k,2) = dble(wave3(k))*dx3
    forall(i=1:n1_fft,k=1:n3_fft) vartmp(i,1,k,3) = real(cross_spectra(i,n2_fft/2+1,k))
    forall(i=1:n1_fft,k=1:n3_fft) vartmp(i,1,k,4) = imag(cross_spectra(i,n2_fft/2+1,k))
    forall(i=1:n1_fft) vartmp(i,1,:,5) = auto_spectra(i)
    call WriteTec(trim(file_out))

    file_out = 'coherence_3d_kref'//fnum_kloc//'.plt'
    print*,'write ',trim(file_out)
    tecname = 'f, y, x, spect12_r, spect12_i, spect_ref'
    call InitTec(1,n1_fft,n2_fft,n3_fft,6,tecname,0)
    forall(i=1:n1_fft) vartmp(i,:,:,1) = dble(wave1(i))/lx1
    forall(j=1:n2_fft) vartmp(:,j,:,2) = dble(wave2(j))*dx2
    forall(k=1:n3_fft) vartmp(:,:,k,3) = dble(wave3(k))*dx3
    forall(i=1:n1_fft,j=1:n2_fft,k=1:n3_fft) vartmp(i,j,k,4) = real(cross_spectra(i,j,k))
    forall(i=1:n1_fft,j=1:n2_fft,k=1:n3_fft) vartmp(i,j,k,5) = imag(cross_spectra(i,j,k))
    forall(i=1:n1_fft) vartmp(i,:,:,6) = auto_spectra(i)
    call WriteTec(trim(file_out))
    deallocate(wave1,wave2,wave3)
end subroutine main

end program Calcoh_FFT3D
