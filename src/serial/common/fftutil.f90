
! Signal filtering using fl and fu
! ns: number of sample
! dts: time interval/length interval
! fl: lower bound frequency/wavenumber
! fu: upper bound frequencywavenumber
subroutine FFTFilter(ns,varin,dts,fl,fu,varout)
    implicit none
    integer, intent(in) :: ns
    real(8), intent(in) :: varin(ns)
    real(8), intent(in) :: dts, fl, fu
    real(8), intent(out) :: varout(ns)
    complex :: ctmp(ns)!, ctmp2(ns)
    real(8) :: ff, fperiod
    integer :: n, i, npoint
    integer :: lensave, lenwork, ierr
    real, dimension(:), allocatable :: wsave, wwork

    npoint = ns
    lensave = 2*npoint+int(log(real(npoint))/log(2.))+4
    lenwork = 2*npoint
    fperiod = dble(npoint)*dts

    if(allocated(wsave)) deallocate(wsave)
    if(allocated(wwork)) deallocate(wwork)
    allocate(wsave(lensave),wwork(lenwork))
    ! initalize fft
    call cfft1i(npoint,wsave,lensave,ierr)
    if(ierr.ne.0) then
      print*, 'FFT1D initialization error!'
      stop
    endif

    ctmp(1:npoint) = varin(1:npoint)
    call cfft1f(npoint,1,ctmp,npoint,wsave,lensave,wwork,lenwork,ierr)
    if(ierr.ne.0) then
      print*, 'Error occured in FFT1DF!'
      stop
    endif
    
    do i=1, npoint/2
      ff = dble(i)/fperiod
      if(ff.le.fl.or.ff.ge.fu) then
        ctmp(i) = 0.
      endif
    enddo
    do i=npoint/2+1, npoint
      ff =  abs(dble(i-npoint-1)/fperiod)
      if(ff.le.fl.or.ff.ge.fu) then
        ctmp(i) = 0.
      endif
    enddo

    call cfft1b(npoint,1,ctmp,npoint,wsave,lensave,wwork,lenwork,ierr)
    if(ierr.ne.0) then
      print*, 'Error occured in FFT1DB!'
      stop
    endif
    varout = real(ctmp)

end subroutine FFTFilter


!subroutine to rearrange spectra output by 1D FFT so that
!the zero-frequency component is at the center of the series.
!
!np (input): number of points in the spectrum
!spectin (input): the spectrum output by 1D FFT
!spectout (output): the rearranged spectrum
!wavenum (output): the corresponding wave number of the output spectrum
subroutine RearrangeFFT1D(np,spectin,spectout,wavenum)
  integer, intent(in) :: np
  complex, dimension(np), intent(inout) :: spectin
  complex, dimension(np), intent(out) :: spectout
  integer, dimension(np), intent(out) :: wavenum
  integer :: ihalf, i
  ihalf = np/2+1
  do i=1,np
    wavenum(i) = ihalf-1+i-np
    if (i.le.np-ihalf) then
      spectout(i) = spectin(ihalf+i)
    else
      spectout(i) = spectin(i-np+ihalf)
    end if
  end do
end subroutine RearrangeFFT1D

!subroutine to rearrange spectra output by 2D FFT so that
!the zero-frequency component is at the center of the series.
!
!np1 (input): number of points in 1st dimension of the spectrum
!np2 (input): number of points in 2nd dimension of the spectrum
!spectin (input): the spectrum output by 2D FFT
!spectout (output): the rearranged spectrum
!wavenum1 (output): the corresponding wave number in the 1st dimension of the output spectrum
!wavenum2 (output): the corresponding wave number in the 2nd dimension of the output spectrum
subroutine RearrangeFFT2D(np1,np2,spectin,spectout,wavenum1,wavenum2)
  integer, intent(in) :: np1, np2
  complex, dimension(np1,np2), intent(inout) :: spectin
  complex, dimension(np1,np2), intent(out) :: spectout
  complex, dimension(np2) :: tmp
  integer, intent(out) :: wavenum1(np1),wavenum2(np2)
  integer :: ihalf, jhalf, i, j
  ihalf = np1/2+1; jhalf = np2/2+1
  do i=1,np1
    wavenum1(i) = ihalf-1+i-np1
    if (i.le.np1-ihalf) then
      spectout(i,:) = spectin(ihalf+i,:)
    else
      spectout(i,:) = spectin(i-np1+ihalf,:)
    end if
  end do
  do j=1,np2
    wavenum2(j) = jhalf-1+j-np2
  end do
  do i=1,np1
    tmp = spectout(i,:)
    do j=1,np2
    if (j.le.np2-jhalf) then
      spectout(i,j) = tmp(jhalf+j)
    else
      spectout(i,j) = tmp(j-np2+jhalf)
    end if
    end do
  end do
end subroutine RearrangeFFT2D


subroutine RearrangeFFT2D_new(np1,np2,spectin,spectout,wavenum1,wavenum2)
  integer, intent(in) :: np1, np2
  complex, dimension(np1,np2), intent(inout) :: spectin
  complex, dimension(np1,np2), intent(out) :: spectout
  complex, dimension(np2) :: tmp
  integer, intent(out) :: wavenum1(np1),wavenum2(np2)
  integer :: ihalf, jhalf, i, j
  ihalf = np1/2+1; jhalf = np2/2+1
  do i=1,np1
    wavenum1(i) = ihalf-1+i-np1
   ! if (i.le.np1-ihalf) then
      spectout(i,:) = spectin(i,:)
   ! else
   !   spectout(i,:) = spectin(i-np1+ihalf,:)
   ! end if
  end do
  do j=1,np2
    wavenum2(j) = jhalf-1+j-np2
  end do
  do i=1,np1
    tmp = spectout(i,:)
    do j=1,np2
    if (j.le.np2-jhalf) then
      spectout(i,j) = tmp(jhalf+j)
    else
      spectout(i,j) = tmp(j-np2+jhalf)
    end if
    end do
  end do
end subroutine RearrangeFFT2D_new

!subroutine to rearrange spectra output by 3D FFT so that
!the zero-frequency component is at the center of the series.
!
!np1 (input): number of points in 1st dimension of the spectrum
!np2 (input): number of points in 2nd dimension of the spectrum
!np3 (input): number of points in 2nd dimension of the spectrum
!spectin (input): the spectrum output by 3D FFT
!spectout (output): the rearranged spectrum
!wavenum1 (output): the corresponding wave number in the 1st dimension of the output spectrum
!wavenum2 (output): the corresponding wave number in the 2nd dimension of the output spectrum
!wavenum3 (output): the corresponding wave number in the 3nd dimension of the output spectrum
subroutine RearrangeFFT3D(np1,np2,np3,spectin,spectout,wavenum1,wavenum2,wavenum3)
  integer, intent(in) :: np1, np2, np3
  complex, dimension(np1,np2,np3), intent(inout) :: spectin
  complex, dimension(np1,np2,np3), intent(out) :: spectout
  complex, dimension(np2,np3) :: tmp
  complex, dimension(np3) :: tmp2
  integer, intent(out) :: wavenum1(np1),wavenum2(np2),wavenum3(np3)
  integer :: ihalf, jhalf, khalf, i, j, k
  ihalf = np1/2+1; jhalf = np2/2+1; khalf = np3/2+1
  do i=1,np1
    wavenum1(i) = ihalf-1+i-np1
    if (i.le.np1-ihalf) then
  !print *, 'i = ', i
      spectout(i,:,:) = spectin(ihalf+i,:,:)
    else
      spectout(i,:,:) = spectin(i-np1+ihalf,:,:)
    endif
  enddo
  do j=1,np2
    wavenum2(j) = jhalf-1+j-np2
  enddo
  do i=1,np1
    tmp = spectout(i,:,:)
    do j=1,np2
    if (j.le.np2-jhalf) then
      spectout(i,j,:) = tmp(jhalf+j,:)
    else
      spectout(i,j,:) = tmp(j-np2+jhalf,:)
    endif
    enddo
  enddo

  do k=1,np3
    wavenum3(k) = khalf-1+k-np3
  enddo
  do i=1,np1
    do j=1, np2
      tmp2 = spectout(i,j,:)
      do k=1,np3
        if (k.le.np3-khalf) then
          spectout(i,j,k) = tmp2(khalf+k)
        else
          spectout(i,j,k) = tmp2(k-np3+khalf)
        endif
      enddo ! end k loop
    enddo ! end j loop
  enddo ! end i loop

end subroutine RearrangeFFT3D
