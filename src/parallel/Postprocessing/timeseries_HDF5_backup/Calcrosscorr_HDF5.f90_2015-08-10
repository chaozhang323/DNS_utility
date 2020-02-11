

    ! program to compute cross-correlation using time series data
    program Calcrosscorr
      use MTSKplane
      use modSpect
      implicit none
      integer :: icrossx, icrossxt
      integer :: nperiod, ioverlap, iwindow, ntrans, nsection
      integer :: nskip, kloc_ref, ixwindowl, ixwindowr, varindex_ref(1)
      integer :: num_kplane
      integer :: kplane_be, kplane_end
      integer, dimension(:), allocatable :: kplane
      integer :: n, ibe_k, iend_k, jbe_k, jend_k
      ! integer :: imints, imaxts, nxpoint, nypoint
      type(fprop), dimension(:), allocatable :: fileprop
      character(4) :: fnum_ref, fnum
      character(8) :: fnum_8
      integer :: nvar_output
      integer, dimension(:), allocatable :: varindex_output
!      character(10), dimension(:), allocatable :: varname_output
      character(10) :: varname_ref(1)
      integer :: nxave, nyave
      integer :: ii, i, j, k, m
      real(8), dimension(:,:,:,:,:), allocatable :: tsdata_ref, tsdata
      real(8), dimension(:), allocatable :: corr12, ave1, ave2, var1, var2, rms1, rms2
      real(8), dimension(:), allocatable :: corrx
      integer, dimension(:), allocatable :: numpt
      real(8) :: dt_sample, dx_kp, dy_kp
      character(400) :: filename
      complex, dimension(:,:), allocatable :: spect12, spect11, spect22
      real(8) :: zloc_kp, dzdk_kp, zloc_ref, dzdk_ref

      type(tp_DNSIndex) :: DNSIndex_k

      integer :: npath, ireadVOL, num_file, ntpoint_total
      integer :: kloc_ref_index

      call Input()
      call InitHDF5()

      print *, 'ntpoint_total =', ntpoint_total, 'dt_sample =', dt_sample

      call InitSpect1d(ntpoint_total, nperiod, iwindow, ioverlap, dt_sample, ntrans, nsection)
      call InitFFT1D(ntrans)
      call InitFFTWindow(ntrans,iwindow)

      call InitVOLTSkplane(ntpoint_total, npath, DNSIndex_k, ibe_k, iend_k, jbe_k, jend_k, &
                           ixwindowl, ixwindowr, nvar_output, varindex_output, fileprop)

      print *, 'Reference wall-normal location k = ', kloc_ref
      write(unit=fnum_ref,fmt='(I04.4)') kloc_ref

      varname_ref(1) = varname_kplane(varindex_ref(1))
      print *, 'Reference variable name is ', trim(varname_ref(1))
      print *, 'Reading Reference kplane ... '

      allocate(corr12(-ixwindowl:ixwindowr),ave1(-ixwindowl:ixwindowr),ave2(-ixwindowl:ixwindowr))
      allocate(var1(-ixwindowl:ixwindowr),var2(-ixwindowl:ixwindowr),numpt(-ixwindowl:ixwindowr))
      allocate(rms1(-ixwindowl:ixwindowr),rms2(-ixwindowl:ixwindowr))
      allocate(corrx(-ixwindowl:ixwindowr))
      allocate(spect12(ntrans,-ixwindowl:ixwindowr),spect11(ntrans,-ixwindowl:ixwindowr),spect22(ntrans,-ixwindowl:ixwindowr))


      num_file = (fileprop(1)%file_end-fileprop(1)%file_be)/fileprop(1)%file_skip + 1

!      call ReadVOLTSkplane_perFile(filename,kloc_ref_index,tsdata_ref) !!!!

      ! Reading k plane data
      allocate(tsdata_ref(ntpoint_total,nypoint_k,imints_k:imaxts_k,1,1))
      allocate(tsdata(ntpoint_total,nypoint_k,imints_k:imaxts_k,1,nvar_output) )  !!!

      do k=1, num_kplane

        do n=1, num_file
          write(unit=fnum_8,fmt='(I08.8)') fileprop(1)%file_be + (n-1)*fileprop(1)%file_skip
          filename = trim(fileprop(1)%filepath)//'timeseries_'//fnum_8//'.h5'
          print *, 'Reading file: ', trim(filename)
          call ReadVOLTSkplane_perFile(filename,kloc_ref_index,tsdata_ref)

          print *, 'wall-normal location k = ', kplane(k)
          call ReadVOLTSkplane_perFile(filename,k+kplane_be-1,tsdata) !!

          write(unit=fnum,fmt='(I04.4)') kplane(k)
          do m=1, nvar_output
            print *, 'variable name is ', trim(varout_k(m))
            corr12 = 0.d0
            ave1 = 0.d0
            ave2 = 0.d0
            var1 = 0.d0
            var2 = 0.d0
            numpt = 0
            do i = 1, nxpoint_k
              do ii = max(-ixwindowl,imints_k-i), min(ixwindowr,imaxts_k-i)

                corr12(ii) = corr12(ii) + (sum(tsdata_ref(:,:,i,1,1)*tsdata(:,:,i+ii,1,m)))
                ave1(ii) = ave1(ii) + sum(tsdata_ref(:,:,i,1,1))
                ave2(ii) = ave2(ii) + sum(tsdata(:,:,i+ii,1,m))
                var1(ii) = var1(ii) + sum(tsdata_ref(:,:,i,1,1)**2)
                var2(ii) = var2(ii) + sum(tsdata(:,:,i+ii,1,m)**2)
!                corr12(ii) = corr12(ii) + (sum(tsdata_ref(1,i,:,:)*tsdata(m,i+ii,:,:)))
!                ave1(ii) = ave1(ii) + sum(tsdata_ref(1,i,:,:))
!                ave2(ii) = ave2(ii) + sum(tsdata(m,i+ii,:,:))
!                var1(ii) = var1(ii) + sum(tsdata_ref(1,i,:,:)**2)
!                var2(ii) = var2(ii) + sum(tsdata(m,i+ii,:,:)**2)
                numpt(ii) = numpt(ii) + 1
              enddo
            enddo
            forall(ii=-ixwindowl:ixwindowr,numpt(ii).gt.0)
               corr12(ii) = corr12(ii)/(ntpoint_total*nypoint_k*numpt(ii))
               ave1(ii) = ave1(ii)/(ntpoint_total*nypoint_k*numpt(ii))
               ave2(ii) = ave2(ii)/(ntpoint_total*nypoint_k*numpt(ii))
               var1(ii) = var1(ii)/(ntpoint_total*nypoint_k*numpt(ii))
               var2(ii) = var2(ii)/(ntpoint_total*nypoint_k*numpt(ii))
               rms1(ii) = sqrt( abs(var1(ii)-ave1(ii)**2) )
               rms2(ii) = sqrt( abs(var2(ii)-ave2(ii)**2) )
               corrx(ii) = (corr12(ii) - ave1(ii)*ave2(ii))/(rms1(ii)*rms2(ii))
            end forall

            if(icrossx.eq.1) then
              filename = 'corrxcross_kref'//fnum_ref//'_k'//fnum//'_'//trim(varname_ref(1))//trim(varout_k(m))//'.h5'
              call outputcrossx(corrx,numpt,trim(filename))
            endif

            if(icrossxt.eq.1) then
              spect12 = (0.d0,0.d0)
              numpt = 0
              do j = 1, nypoint_k
                do i = 1, nxpoint_k
                  do ii = max(-ixwindowl,imints_k-i), min(ixwindowr,imaxts_k-i)
!                    spect12(:,ii) = spect12(:,ii) + crossspect1d(tsdata_ref(1,i,j,:),tsdata(m,i+ii,j,:))
                    spect12(:,ii) = spect12(:,ii) + crossspect1d(tsdata_ref(:,j,i,1,1),tsdata(:,j,i+ii,1,m))
                    numpt(ii) = numpt(ii) + 1
                  enddo
                enddo
              enddo
              filename = '_kref'//fnum_ref//'_k'//fnum//'_'//trim(varname_ref(1))//trim(varout_k(m))
              call Calcorrxtcross_HDF5(ixwindowl,ixwindowr,dx_kp,spect12,corrx,filename)
           endif  ! end calculate x-t cross correlation

          enddo ! end m loop
        enddo ! end n loop



!        print *, 'wall-normal location k = ', kplane(k)
!        call ReadVOLTSkplane_perFile(    ) !!
!        write(unit=fnum,fmt='(I04.4)') kplane(k)
!        do m=1, nvar_output
!          print *, 'variable name is ', trim(varname_output(m))
!          corr12 = 0.d0
!          ave1 = 0.d0
!          ave2 = 0.d0
!          var1 = 0.d0
!          var2 = 0.d0
!          numpt = 0
!          do i = 1, nxpoint
!             do ii = max(-ixwindowl,imints-i), min(ixwindowr,imaxts-i)
!                corr12(ii) = corr12(ii) + (sum(tsdata_ref(1,i,:,:)*tsdata(m,i+ii,:,:)))
!                ave1(ii) = ave1(ii) + sum(tsdata_ref(1,i,:,:))
!                ave2(ii) = ave2(ii) + sum(tsdata(m,i+ii,:,:))
!                var1(ii) = var1(ii) + sum(tsdata_ref(1,i,:,:)**2)
!                var2(ii) = var2(ii) + sum(tsdata(m,i+ii,:,:)**2)
!                numpt(ii) = numpt(ii) + 1
!             enddo
!          enddo
!          forall(ii=-ixwindowl:ixwindowr,numpt(ii).gt.0)
!             corr12(ii) = corr12(ii)/(ntpoint*nypoint*numpt(ii))
!             ave1(ii) = ave1(ii)/(ntpoint*nypoint*numpt(ii))
!             ave2(ii) = ave2(ii)/(ntpoint*nypoint*numpt(ii))
!             var1(ii) = var1(ii)/(ntpoint*nypoint*numpt(ii))
!             var2(ii) = var2(ii)/(ntpoint*nypoint*numpt(ii))
!             rms1(ii) = sqrt( abs(var1(ii)-ave1(ii)**2) )
!             rms2(ii) = sqrt( abs(var2(ii)-ave2(ii)**2) )
!             corrx(ii) = (corr12(ii) - ave1(ii)*ave2(ii))/(rms1(ii)*rms2(ii))
!          end forall
!
!        enddo ! end m loop
!

      enddo ! end k loop



contains

      subroutine Input()
        integer :: i, n

        ! read input file
        read(*,*)
        read(*,*) icrossx, icrossxt
        read(*,*)
        read(*,*) npath, ireadVOL
        allocate(fileprop(npath))
        read(*,*)
        do n=1, npath
          read(*,*)
          read(*,*) fileprop(n)%ntpoint, fileprop(n)%nskip
          read(*,*)
          read(*,*) fileprop(n)%file_be, fileprop(n)%file_end, fileprop(n)%file_skip
          read(*,*)
          read(*,'(a)') fileprop(n)%filepath
        enddo
        read(*,*)
        read(*,*)
        read(*,*) nperiod, ioverlap, iwindow, dt_sample
        read(*,*)
        read(*,*)
        read(*,*) nvar_output
        allocate(varindex_output(nvar_output))
        read(*,*) (varindex_output(i), i = 1, nvar_output)
        read(*,*)
        read(*,*) kloc_ref, kloc_ref_index, varindex_ref(1)
        read(*,*)
        read(*,*) ibe_k, iend_k, jbe_k, jend_k, ixwindowl, ixwindowr, dx_kp, dy_kp
        read(*,*)
        read(*,*) DNSIndex_k%ibe, DNSIndex_k%iend, DNSIndex_k%iskip, &
                  DNSIndex_k%jbe, DNSIndex_k%jend, DNSIndex_k%jskip, DNSIndex_k%ibuffer
        read(*,*)
        read(*,*) num_kplane, kplane_be, kplane_end
        allocate(kplane(num_kplane))
        read(*,*) (kplane(i), i = 1, num_kplane)
        ! end read input file

        print *, 'ixwindowl =', ixwindowl,  'ixwindowr =', ixwindowr
        print *, 'number of wall-normal locations =', num_kplane
        print *, 'Wall-normal plane indexes k =', kplane
        write(*,*)

        ntpoint_total = 0
        if(ireadVOL.eq.0) then
          do n=1, npath
            fileprop(n)%ntpoint = fileprop(n)%ntpoint - fileprop(n)%nskip
            ntpoint_total = ntpoint_total + fileprop(n)%ntpoint
          enddo
        else
          do n=1, npath
            ntpoint_total = ntpoint_total + ( (fileprop(n)%file_end-fileprop(n)%file_be)/fileprop(n)%file_skip + 1 )*fileprop(n)%ntpoint
          enddo
        endif


      end subroutine Input


      subroutine outputcrossx(corr,npt,fn)
        implicit none
        real(8), intent(in) :: corr(-ixwindowl:ixwindowr)
        integer, intent(in) :: npt(-ixwindowl:ixwindowr)
        character(*), intent(in) :: fn

        type(tp_hyperslab) :: fsol_crossx
        real(8), dimension(:,:), allocatable :: buffer_crossx

        fsol_crossx%gname = '/'
        fsol_crossx%dnum = 2
        fsol_crossx%rank = 1
        allocate(fsol_crossx%dname(fsol_crossx%dnum),fsol_crossx%dimsf(fsol_crossx%rank))
        allocate(fsol_crossx%dimsm(fsol_crossx%rank),fsol_crossx%count(fsol_crossx%rank))
        allocate(fsol_crossx%offset(fsol_crossx%rank),fsol_crossx%block(fsol_crossx%rank))
        allocate(fsol_crossx%stride(fsol_crossx%rank))
        fsol_crossx%dname(1) = 'x'
        fsol_crossx%dname(2) = 'corr'

        fsol_crossx%dimsf(1) = ixwindowr + ixwindowl + 1
        fsol_crossx%dimsm(1) = fsol_crossx%dimsf(1)
        fsol_crossx%block = fsol_crossx%dimsm
        fsol_crossx%count = 1; fsol_crossx%offset = 0
        fsol_crossx%stride = 1
        fsol_crossx%IsHSInitialized = .true.

        allocate(buffer_crossx(fsol_crossx%dimsf(1),2))

        do ii=-ixwindowl, ixwindowr
          if(npt(ii).gt.0) then
            buffer_crossx(ii+ixwindowl+1,1) = dble(ii)*dx_kp
            buffer_crossx(ii+ixwindowl+1,2) = corr(ii)
          endif
        enddo
        fsol_crossx%fname = trim(fn)
        print *, 'Writing file: ', trim(fsol_crossx%fname)
        call WriteTSHDF5_1D(fsol_crossx,buffer_crossx)

      end subroutine outputcrossx






    end program Calcrosscorr
