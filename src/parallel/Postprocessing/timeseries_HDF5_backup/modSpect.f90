module modSpect
   use MFFT1D
   use MFFTWindow
   use MTSHDF5
   implicit none
   integer, private :: nsection
   integer, private :: iwindow,  ioverlap, nperiod
   integer, private :: ntpoint,  ntrans
   real(8), private :: dts
   real(8), private :: pi
   logical :: IsSpectInit = .FALSE.

   type(tp_hyperslab) :: fsol_corrxt, fsol_corrxtcross
   real(8), dimension(:,:,:), allocatable :: buffer_corrxt, buffer_corrxtcross

!                  nsection (# of bins)
!                  iwindow ( windowing type)
!                  ioverlap (whether or not 1/2 overlap between neighboring bins)
!                  ntpoint (total number of points), ntrans (# of points per bin)
!                  ntrans (size of bins used for Fourier transformation)
!                  nsection (number of bins; nsection = 8 for standard Welch method)
contains
     subroutine InitSpect1d(nt, npr, iwin, ihflap, dt_sample, npbin, nsect)
        implicit none
        integer, intent(in) :: nt, npr, iwin, ihflap
        real(8), intent(in) :: dt_sample  ! sampling time interval
        integer, intent(out) :: npbin, nsect

        ntpoint = nt
        nperiod = npr
        iwindow = iwin
        ioverlap = ihflap
        dts = dt_sample

        pi = 4.*atan(1.)
        ntrans = nperiod*(ioverlap+1)
        npbin = ntrans
        if (ntrans.gt.ntpoint) then
             print*, 'specified nperiod with ioverlap flag greater than total points available'
             print*, 'ntpoint = ', ntpoint, ', nperiod = ', nperiod, 'ioverlap=',ioverlap
             stop
        else if (nperiod.le.0) then
             nperiod = ntpoint
             ntrans = nperiod
             ioverlap = 0
        end if
        nsection = ntpoint/nperiod - ioverlap
        nsect = nsection
        print *, '# of data points per segments =', ntrans
        print *, '# of segments =', nsection
        print *, 'Sampling frequency (Herz) =', 1./dts
        print *, 'Interval of sampling (sec) =', dts
        print *, 'Length per segments (sec) =', (ntrans-1)*dts
        if(ioverlap.eq.0) print *, 'Segments with no overlap'
        if(ioverlap.eq.1) print *, 'Segments with half overlap'
        call InitFFT1D(ntrans)
        call InitFFTWindow(ntrans,iwindow)
        IsSpectInit = .TRUE.
    end subroutine InitSpect1d

    function spect1d(tfun)
      implicit none
      real(8), intent(in) :: tfun(ntpoint)
      complex :: spect1d(ntrans,nsection)
      real(8) :: tmean
      integer :: nn
      complex, dimension(ntrans) :: ctmp

      ! Calculate the mean and energy of time signal
      tmean = sum(tfun(1:ntpoint))/dble(ntpoint)

      ! Compute spectrum for time data length of ntans
      do nn = 1, nsection
         ctmp(1:ntrans) = wcoeff(1:ntrans)*(tfun((nn-1)*nperiod+1:(nn-1)*nperiod+ntrans)-tmean)
         call FFT1DF(ctmp, ctmp)
         spect1d(:,nn) = ctmp
      enddo
    end function spect1d

    function autospect1d(tfun)
      implicit none
      real(8), intent(in) :: tfun(ntpoint)
      complex :: autospect1d(ntrans)
      complex :: ctmp(ntrans,nsection)

      ctmp = spect1d(tfun)
      autospect1d = sum(ctmp*conjg(ctmp),dim=2)/dble(nsection)
    end function autospect1d

    function crossspect1d(tfun1,tfun2)
      implicit none
      real(8), intent(in) :: tfun1(ntpoint), tfun2(ntpoint)
      complex :: crossspect1d(ntrans)
      complex, dimension(ntrans,nsection) :: ctmp1, ctmp2

      ctmp1 = spect1d(tfun1)
      ctmp2 = spect1d(tfun2)
      crossspect1d = sum(ctmp1*conjg(ctmp2),dim=2)/dble(nsection)

    end function crossspect1d

    subroutine CalPowerspect(spect1d,tmeanenergy,fn)
      implicit none
      complex, intent(in) :: spect1d(ntrans)
      real(8), intent(in) :: tmeanenergy
      real(8) :: powerspect(ntrans/2), tmpenergy
      character(*), intent(in) :: fn
      real(8) :: ff, fperiod
      integer :: i, n
      character(100) :: vars

      ! Ref: Numerical Recipe, Chapter 13.4
      ! tmeanenergy: mean squre amplitude (1/N)*sum_1^N |c_i|^2
      ! PSD is normalized to preserve time-integral squred amplitude int_0^T |c(t)|^2 dt

      fperiod = dble(ntrans-1)*dts  ! period per segments
      powerspect(1) = spect1d(1)
      powerspect(2:ntrans/2) = 2.*spect1d(2:ntrans/2)
      tmpenergy = sum(powerspect(1:ntrans/2))/fperiod  ! sum (fi_i * df), with df = 1/fperiod
      open(21, file='powerspect'//trim(fn)//'.dat')
      rewind(21)
      write(21,'(a)')'variables=f,spectrum'
      write(21,'("Zone T=PSD, I=", I8.8, ", AUXDATA rms2=""", E15.9,"""")') ntrans/2, tmeanenergy
      write(21,'(a,I8.8,a,I3.3)') '# Data points per seg = ', ntrans, '   Num of seg = ', nsection
      write(21,'(a,I1.1,a,I1.1,a,E15.8)') '# Window type = ', iwindow, '   ioverlap = ', ioverlap, '   dt_sample = ', dts
      write(21,'(a,E15.8)') '# meanenergy (mean square amplitude) = ', tmeanenergy
      do i=1,ntrans/2
        ff = real(i-1)/fperiod
        write(21,'(2E20.12)') ff, powerspect(i)*tmeanenergy/tmpenergy
      end do
      close(21)
    end subroutine CalPowerspect

    subroutine CalCoherence(autosp1,autosp2,crosssp,fn)
      implicit none
      complex, intent(in) :: autosp1(ntrans), autosp2(ntrans), crosssp(ntrans)
      character(*), intent(in) :: fn
      complex :: ctmp(ntrans)
      real(8) :: Coh(ntrans)
      real(8) :: ff, fperiod
      integer :: i, n
      character(100) :: vars

      fperiod = dble(ntrans-1)*dts  ! period per segments
      ctmp = abs(crosssp)**2/(abs(autosp1)*abs(autosp2))
      Coh(1) = ctmp(1)
      Coh(2:ntrans/2) = 2.*ctmp(2:ntrans/2)
      open(21, file='powerspect'//trim(fn)//'.dat')
      write(21,'(a)')'variables=f,coherence'
      write(21,'(a,I8.8,a,I3.3)') '# Data points per seg = ', ntrans, '   Num of seg = ', nsection
      write(21,'(a,I1.1,a,I1.1,a,E15.8)') '# Window type = ', iwindow, '   ioverlap = ', ioverlap, '   dt_sample = ', dts
      do i=1,ntrans/2
        ff = real(i-1)/fperiod
        write(21,'(2E20.12)')ff, Coh(i)
      end do
      close(21)
    end subroutine CalCoherence

    subroutine CalCorrt(spect1d,fn)
      implicit none
      complex, intent(in) :: spect1d(ntrans)
      character(*), intent(in) :: fn
      complex(8), dimension(ntrans) :: corrt, corrout
      integer :: twave(ntrans)
      integer ::  n

!      call InitFFT1D(ntrans)
      call FFT1DB(spect1d, corrt)
      call RearrangeFFT1D(ntrans,corrt,corrout,twave)
      open(21, file='corrt'//trim(fn)//'.dat')
      write(21,'(a)')'variables=t,corr'
      write(21,'(a,I8.8,a,I3.3)') '# Data points per seg = ', ntrans, '   Num of seg = ', nsection
      write(21,'(a,I1.1,a,I1.1,a,E15.8)') '# Window type = ', iwindow, '   ioverlap = ', ioverlap, '   dt_sample = ', dts
      do n=1,ntrans
        write(21,'(2E20.12)')real(twave(n))*dts, real(corrout(n))/maxval(real(corrout))
      end do
      close(21)
    end subroutine CalCorrt


    subroutine CalCorrxt_HDF5(ixhalflenl,ixhalflenr,iiloc,dxs,spect1d,fn,Current_n,Current_i)
      implicit none
      integer, intent(in) :: ixhalflenl,ixhalflenr, iiloc, Current_n, Current_i
      real(8), intent(in) :: dxs
      complex, intent(in) :: spect1d(ntrans,-ixhalflenl:ixhalflenr)
      character(*), intent(in) :: fn
      complex(8), dimension(ntrans) :: ctmp
      real(8), dimension(ntrans,-ixhalflenl:ixhalflenr) :: corrxt
      complex(8), dimension(ntrans,-ixhalflenl:ixhalflenr) :: corrout
      integer :: twave(ntrans)
      integer ::  n,i, iloc(1)
      real(8), dimension(ntrans) :: uc
      character(4) :: fnum
      integer :: buffer_attr_integer(4)
      real(8) :: buffer_attr_real(2)

      if(Current_i.eq.0) then
        fsol_corrxt%IsMultiGroup = .false.
      else
        fsol_corrxt%IsMultiGroup = .true.
      endif
      write(unit=fnum,fmt='(I04.4)') iiloc
      fsol_corrxt%gname = '/i'//fnum

      if(Current_n.eq.1) then
!        write(unit=fnum,fmt='(I04.4)') iiloc
!        fsol_corrxt%gname = '/i'//fnum
        fsol_corrxt%dnum = 3
        fsol_corrxt%rank = 2
        if(allocated(fsol_corrxt%dname)) deallocate(fsol_corrxt%dname,fsol_corrxt%dimsf,fsol_corrxt%dimsm, &
                    fsol_corrxt%count,fsol_corrxt%offset,fsol_corrxt%block,fsol_corrxt%stride)
        allocate(fsol_corrxt%dname(fsol_corrxt%dnum),fsol_corrxt%dimsf(fsol_corrxt%rank))
        allocate(fsol_corrxt%dimsm(fsol_corrxt%rank),fsol_corrxt%count(fsol_corrxt%rank))
        allocate(fsol_corrxt%offset(fsol_corrxt%rank),fsol_corrxt%block(fsol_corrxt%rank))
        allocate(fsol_corrxt%stride(fsol_corrxt%rank))
        fsol_corrxt%dname(1) = 't'
        fsol_corrxt%dname(2) = 'x'
        fsol_corrxt%dname(3) = 'corr'

        fsol_corrxt%dimsf(1) = ixhalflenr+ixhalflenl+1
        fsol_corrxt%dimsf(2) = ntrans
        fsol_corrxt%dimsm(1) = fsol_corrxt%dimsf(1)
        fsol_corrxt%dimsm(2) = fsol_corrxt%dimsf(2)
        fsol_corrxt%block  = fsol_corrxt%dimsm
        fsol_corrxt%count  = 1; fsol_corrxt%offset = 0
        fsol_corrxt%stride = 1

        fsol_corrxt%IsHSInitialized = .true.
        if(allocated(buffer_corrxt)) deallocate(buffer_corrxt)
        allocate(buffer_corrxt(fsol_corrxt%dimsf(1),fsol_corrxt%dimsf(2),3))
        if(allocated(fsol_corrxt%attr_name)) deallocate(fsol_corrxt%attr_name)
        allocate(fsol_corrxt%attr_name(6))
        fsol_corrxt%attr_name(1) = '# Data points per seg '
        fsol_corrxt%attr_name(2) = 'Num of seg '
        fsol_corrxt%attr_name(3) = '# Window type '
        fsol_corrxt%attr_name(4) = 'ioverlap '
        fsol_corrxt%attr_name(5) = 'maxval '
        fsol_corrxt%attr_name(6) = 'dt_sample '

        fsol_corrxt%dims_integer = 4
        fsol_corrxt%dims_real = 2

      endif

      if(ixhalflenl.lt.0.or.ixhalflenr.lt.0) then
         print *, 'ixhalflenl, ixhalflenr need to be >= 0. STOP'
      endif
      do i=-ixhalflenl,ixhalflenr
         call FFT1DB(spect1d(:,i), ctmp)
         call RearrangeFFT1D(ntrans,ctmp,corrout(:,i),twave)
     enddo

!      open(21, file='corrxt'//trim(fn)//'.dat')
!        write(21,'(a)')'variables=t,x,corr'
!        write(21,'("Zone T=corrxt, I=", I8, ",J=",I8)') ntrans, ixhalflenl+ixhalflenr+1
!        write(21,'(a,I8.8,a,I3.3,a,E15.8)') '# Data points per seg = ', ntrans, '   Num of seg = ', nsection, '   maxval = ', maxval(abs(real(corrout)))
!        write(21,'(a,I1.1,a,I1.1,a,E15.8)') '# Window type = ', iwindow, '   ioverlap = ', ioverlap, '   dt_sample = ', dts
!      close(21)

!     corrxt = real(corrout)/maxval(abs(real(corrout)))
     print *, 'ixhalflenl =', ixhalflenl, 'ixhalflenr =', ixhalflenr
     print *, 'real(corrout(ntrans/2,0)) =', real(corrout(ntrans/2,0))
     corrxt = real(corrout)/real(corrout(ntrans/2,0))

!print *, 'ntrans = ', ntrans

     do i=-ixhalflenl,ixhalflenr
       do n=1,ntrans
       !     write(21,'(3E20.12)') real(twave(n))*dts, real(i)*dxs, corrxt(n,i)
         buffer_corrxt(i+ixhalflenl+1,n,1) = real(twave(n))*dts
         buffer_corrxt(i+ixhalflenl+1,n,2) = real(i)*dxs
         buffer_corrxt(i+ixhalflenl+1,n,3) = corrxt(n,i)
       enddo
     enddo

     buffer_attr_integer(1) = ntrans;  buffer_attr_integer(2) = nsection
     buffer_attr_integer(3) = iwindow; buffer_attr_integer(4) = ioverlap
     buffer_attr_real(1) = maxval(abs(real(corrout)))
     buffer_attr_real(2) = dts

     fsol_corrxt%fname = 'corrxt'//trim(fn)//'.h5'
     print *, 'Writing file: ', trim(fsol_corrxt%fname)
     call WriteTSHDF5_2D(fsol_corrxt,buffer_corrxt)
     call WriteTSHDF5_Attribute(fsol_corrxt,buffer_attr_integer,buffer_attr_real)



     ! Calculate convection velocity
     ! Ref: Bernardini & Pirozzoli, PhysFluids, 23, 085102 (2011)
      open(21, file='uconv_Pirozzoli_'//trim(fn)//'.dat')
      write(21,'(a)')'variables=t,ucorr'
      write(21,'(a,I8.8,a,I3.3)') '# Data points per seg = ', ntrans, '   Num of seg = ', nsection
      write(21,'(a,I1.1,a,I1.1,a,E15.8)') '# Window type = ', iwindow, '   ioverlap = ', ioverlap, '   dt_sample = ', dts
      do n = ntrans/2+1, ntrans
        iloc = maxloc(corrxt(n,:))
        uc(n) = real(iloc(1))*dxs/( real(twave(n))*dts )
        write(21,'(2E20.12)')real(twave(n))*dts, uc(n)
      enddo
      close(21)

    end subroutine CalCorrxt_HDF5




    subroutine CalCorrxt(ixhalflenl,ixhalflenr,dxs,spect1d,fn)
      implicit none
      integer, intent(in) :: ixhalflenl,ixhalflenr
      real(8), intent(in) :: dxs
      complex, intent(in) :: spect1d(ntrans,-ixhalflenl:ixhalflenr)
      character(*), intent(in) :: fn
      complex(8), dimension(ntrans) :: ctmp
      real(8), dimension(ntrans,-ixhalflenl:ixhalflenr) :: corrxt
      complex(8), dimension(ntrans,-ixhalflenl:ixhalflenr) :: corrout
      integer :: twave(ntrans)
      integer ::  n,i, iloc(1)
      real(8), dimension(ntrans) :: uc

      if(ixhalflenl.lt.0.or.ixhalflenr.lt.0) then
         print *, 'ixhalflenl, ixhalflenr need to be >= 0. STOP'
      endif
      do i=-ixhalflenl,ixhalflenr
         call FFT1DB(spect1d(:,i), ctmp)
         call RearrangeFFT1D(ntrans,ctmp,corrout(:,i),twave)
     enddo

      open(21, file='corrxt'//trim(fn)//'.dat')
      write(21,'(a)')'variables=t,x,corr'
      write(21,'("Zone T=corrxt, I=", I8, ",J=",I8)') ntrans, ixhalflenl+ixhalflenr+1
      write(21,'(a,I8.8,a,I3.3,a,E15.8)') '# Data points per seg = ', ntrans, '   Num of seg = ', nsection, '   maxval = ', maxval(abs(real(corrout)))
      write(21,'(a,I1.1,a,I1.1,a,E15.8)') '# Window type = ', iwindow, '   ioverlap = ', ioverlap, '   dt_sample = ', dts

!     corrxt = real(corrout)/maxval(abs(real(corrout)))
     print *, 'ixhalflenl =', ixhalflenl, 'ixhalflenr =', ixhalflenr
     print *, 'real(corrout(ntrans/2,0)) =', real(corrout(ntrans/2,0))
     corrxt = real(corrout)/real(corrout(ntrans/2,0))
     do i=-ixhalflenl,ixhalflenr
         do n=1,ntrans
            write(21,'(3E20.12)') real(twave(n))*dts, real(i)*dxs, corrxt(n,i)
         enddo
      enddo
      close(21)

     ! Calculate convection velocity
     ! Ref: Bernardini & Pirozzoli, PhysFluids, 23, 085102 (2011)
      open(21, file='uconv_Pirozzoli_'//trim(fn)//'.dat')
      write(21,'(a)')'variables=t,ucorr'
      write(21,'(a,I8.8,a,I3.3)') '# Data points per seg = ', ntrans, '   Num of seg = ', nsection
      write(21,'(a,I1.1,a,I1.1,a,E15.8)') '# Window type = ', iwindow, '   ioverlap = ', ioverlap, '   dt_sample = ', dts
      do n = ntrans/2+1, ntrans
        iloc = maxloc(corrxt(n,:))
        uc(n) = real(iloc(1))*dxs/( real(twave(n))*dts )
        write(21,'(2E20.12)')real(twave(n))*dts, uc(n)
      enddo
      close(21)

    end subroutine CalCorrxt

    subroutine CalCorrxtcross(ixhalflenl,ixhalflenr,dxs,spect12,ratio,fn)
      implicit none
      integer, intent(in) :: ixhalflenl,ixhalflenr
      real(8), intent(in) :: dxs
      complex, intent(in) :: spect12(ntrans,-ixhalflenl:ixhalflenr)
      real(8), intent(in) :: ratio(-ixhalflenl:ixhalflenr)
      character(*), intent(in) :: fn
      complex(8), dimension(ntrans) :: ctmp
      real(8), dimension(ntrans,-ixhalflenl:ixhalflenr) :: corrxt
      complex(8), dimension(ntrans,-ixhalflenl:ixhalflenr) :: corrout
      integer :: twave(ntrans)
      integer ::  n,i, iloc(1)
      real(8), dimension(ntrans) :: uc

      open(21, file='corrxtcross'//trim(fn)//'.dat')
      write(21,'(a)')'variables=t,x,corr'
      write(21,'("Zone T=corrxt, I=", I8, ",J=",I8)') ntrans, ixhalflenl+ixhalflenr+1
      write(21,'(a,I8.8,a,I3.3)') '# Data points per seg = ', ntrans, '   Num of seg = ', nsection
      write(21,'(a,I1.1,a,I1.1,a,E15.8)') '# Window type = ', iwindow, '   ioverlap = ', ioverlap, '   dt_sample = ', dts
      do i=-ixhalflenl,ixhalflenr
         call FFT1DB(spect12(:,i), ctmp)
         call RearrangeFFT1D(ntrans,ctmp,corrout(:,i),twave)
     enddo

     corrxt = 0.d0
     do i=-ixhalflenl,ixhalflenr
         if(abs(corrout(ntrans/2,i)).ne.0.d0) corrxt(:,i) = abs(corrout(:,i))/abs(corrout(ntrans/2,i))*ratio(i)
         do n=1,ntrans
            write(21,'(3E20.12)') real(twave(n))*dts, real(i)*dxs, corrxt(n,i)
         enddo
      enddo
      close(21)
!      print *, 'corrxt=',corrxt(ntrans/2,0)

    end subroutine CalCorrxtcross







    subroutine CalCorrxtcross_HDF5(ixhalflenl,ixhalflenr,dxs,spect12,ratio,fn)
      implicit none
      integer, intent(in) :: ixhalflenl,ixhalflenr
      real(8), intent(in) :: dxs
      complex, intent(in) :: spect12(ntrans,-ixhalflenl:ixhalflenr)
      real(8), intent(in) :: ratio(-ixhalflenl:ixhalflenr)
      character(*), intent(in) :: fn
      complex(8), dimension(ntrans) :: ctmp
      real(8), dimension(ntrans,-ixhalflenl:ixhalflenr) :: corrxt
      complex(8), dimension(ntrans,-ixhalflenl:ixhalflenr) :: corrout
      integer :: twave(ntrans)
      integer ::  n,i, iloc(1)
      real(8), dimension(ntrans) :: uc


      fsol_corrxtcross%gname = '/'
      fsol_corrxtcross%dnum = 3
      fsol_corrxtcross%rank = 2
      allocate(fsol_corrxtcross%dname(fsol_corrxtcross%dnum),fsol_corrxtcross%dimsf(fsol_corrxtcross%rank))
      allocate(fsol_corrxtcross%dimsm(fsol_corrxtcross%rank),fsol_corrxtcross%count(fsol_corrxtcross%rank))
      allocate(fsol_corrxtcross%offset(fsol_corrxtcross%rank),fsol_corrxtcross%block(fsol_corrxtcross%rank))
      allocate(fsol_corrxtcross%stride(fsol_corrxtcross%rank))

      fsol_corrxtcross%dname(1) = 't'
      fsol_corrxtcross%dname(2) = 'x'
      fsol_corrxtcross%dname(3) = 'corr'

      fsol_corrxtcross%dimsf(1) = ixhalflenr+ixhalflenl+1
      fsol_corrxtcross%dimsf(2) = ntrans
      fsol_corrxtcross%dimsm(1) = fsol_corrxtcross%dimsf(1)
      fsol_corrxtcross%dimsm(2) = fsol_corrxtcross%dimsf(2)
      fsol_corrxtcross%block  = fsol_corrxtcross%dimsm
      fsol_corrxtcross%count  = 1; fsol_corrxtcross%offset = 0
      fsol_corrxtcross%stride = 1

      fsol_corrxtcross%IsHSInitialized = .true.

      allocate(buffer_corrxtcross(fsol_corrxtcross%dimsf(1),fsol_corrxtcross%dimsf(2),3))


      do i=-ixhalflenl,ixhalflenr
        if(abs(corrout(ntrans/2,i)).ne.0.d0) corrxt(:,i) = abs(corrout(:,i))/abs(corrout(ntrans/2,i))*ratio(i)
        do n=1,ntrans
          buffer_corrxtcross(i+ixhalflenl+1,n,1) = real(twave(n))*dts
          buffer_corrxtcross(i+ixhalflenl+1,n,2) = real(i)*dxs
          buffer_corrxtcross(i+ixhalflenl+1,n,3) = corrxt(n,i)
!       write(21,'(3E20.12)') real(twave(n))*dts, real(i)*dxs, corrxt(n,i)
        enddo
      enddo

      fsol_corrxtcross%fname = 'corrxtcross'//trim(fn)//'.h5'
      print *, 'Writing file: ', trim(fsol_corrxtcross%fname)

      call WriteTSHDF5_2D(fsol_corrxtcross,buffer_corrxtcross)
!      open(21, file='corrxtcross'//trim(fn)//'.dat')
!      write(21,'(a)')'variables=t,x,corr'
!      write(21,'("Zone T=corrxt, I=", I8, ",J=",I8)') ntrans, ixhalflenl+ixhalflenr+1
!      write(21,'(a,I8.8,a,I3.3)') '# Data points per seg = ', ntrans, '   Num of seg = ', nsection
!      write(21,'(a,I1.1,a,I1.1,a,E15.8)') '# Window type = ', iwindow, '   ioverlap = ', ioverlap, '   dt_sample = ', dts
!      do i=-ixhalflenl,ixhalflenr
!         call FFT1DB(spect12(:,i), ctmp)
!         call RearrangeFFT1D(ntrans,ctmp,corrout(:,i),twave)
!     enddo
!
!     corrxt = 0.d0
!     do i=-ixhalflenl,ixhalflenr
!         if(abs(corrout(ntrans/2,i)).ne.0.d0) corrxt(:,i) = abs(corrout(:,i))/abs(corrout(ntrans/2,i))*ratio(i)
!         do n=1,ntrans
!            write(21,'(3E20.12)') real(twave(n))*dts, real(i)*dxs, corrxt(n,i)
!         enddo
!      enddo
!      close(21)


    end subroutine CalCorrxtcross_HDF5





    subroutine CalPhaseSpeed_VanAtta(spect1d,dxs,fn)
      ! Ref: Stegen & Van Atta, JFM, 42, pp 689-699, 1970
      implicit none
      complex, intent(in) :: spect1d(ntrans)
      real(8), intent(in) :: dxs
      character(*), intent(in) :: fn
      real(8), dimension(ntrans/2) :: phasespeed, phase
      integer :: i
      real(8) :: fperiod, ff, ang

      fperiod = dble(ntrans-1)*dts  ! period per segments
      print *, 'fperiod=',fperiod
      phasespeed=0.; phase=0.
      ang = 0.
      open(21, file='phasespeed_VanAtta_'//trim(fn)//'.dat')
      write(21,'(a)')'variables=f,phasespeed,phase'
      write(21,'(a,I8.8,a,I3.3)') '# Data points per seg = ', ntrans, '   Num of seg = ', nsection
      write(21,'(a,I1.1,a,I1.1,a,E15.8)') '# Window type = ', iwindow, '   ioverlap = ', ioverlap, '   dt_sample = ', dts
      do i=2,ntrans/2
        ff  = (i-1)/fperiod
        phase(i) = atan2(imag(spect1d(i)),real(spect1d(i))) - 2.d0*pi*min(sign(1.,imag(spect1d(i))),0.)
        if (abs(phase(i-1)-phase(i)).ge.1.0) ang = ang+2.d0*pi
!        if(imag(spect1d(i)).lt.0.d0) ang = 2.*pi
        if (abs(phase(i)+ang).gt.1.d-14) phasespeed(i)=2.*pi*ff*dxs/abs((phase(i)+ang))
!        if (abs(real(spect1d(i))).gt.1.d-14) phase(i) = atan(imag(spect1d(i))/real(spect1d(i)))

!        !phase can be larger than pi/2 due to large distance
!        if (phase(i-1).ge.0.and.phase(i).lt.0) ang = ang+pi
!        if (abs(phase(i)+ang).gt.1.d-14) phasespeed(i)=2.*pi*ff*dxs/abs((phase(i)+ang))
!        if (abs(phase(i)+ang).gt.1.d-14) phasespeed(i)=2.*pi*ff*dxs/abs((phase(i)+ang))
        write(21,'(3E20.12)')ff, phasespeed(i), phase(i)
      end do
      close(21)
    end subroutine CalPhaseSpeed_VanAtta

    subroutine CalPhaseSpeed_Jimenez(spect1d,spect1dx,fn)
      ! Ref: Juan C. Alamo & Javier Jimenez, JFM, 640, pp 5-26, 2009
      implicit none
      complex, intent(in) :: spect1d(ntrans), spect1dx(ntrans)
      character(*), intent(in) :: fn
      real(8) :: phasespeed
      integer :: i
      real(8) :: fperiod, ff

      fperiod = dble(ntrans-1)*dts  ! period per segments
      print *, 'fperiod=',fperiod
      open(21, file='phasespeed_Jimenez_'//trim(fn)//'.dat')
      write(21,'(a)')'variables=f,phasespeed'
      write(21,'(a,I8.8,a,I3.3)') '# Data points per seg = ', ntrans, '   Num of seg = ', nsection
      write(21,'(a,I1.1,a,I1.1,a,E15.8)') '# Window type = ', iwindow, '   ioverlap = ', ioverlap, '   dt_sample = ', dts
      do i=2,ntrans/2
        ff  = (i-1)/fperiod
        phasespeed = -2.d0*pi*ff*(spect1d(i)*conjg(spect1d(i)))/imag(spect1dx(i)*conjg(spect1d(i)))
        write(21,'(3E20.12)')ff, phasespeed
      end do
      close(21)

    end subroutine CalPhaseSpeed_Jimenez

end module modSpect

