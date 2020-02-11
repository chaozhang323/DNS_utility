

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
      ! character(10), dimension(:), allocatable :: varname_output
      character(10) :: varname_ref(1)
      integer :: nxave, nyave
      integer :: ii, i, j, k, kk, m
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

      call InitHDF5()
      call Input()

      print *, 'ntpoint_total =', ntpoint_total, 'dt_sample =', dt_sample

      call InitSpect1d(ntpoint_total, nsection, iwindow, ioverlap, dt_sample, ntrans, nperiod)
      !call InitFFT1D(ntrans)
      !call InitFFTWindow(ntrans,iwindow)

      if(ireadVOL.eq.0) then
        call InitTSkplane(ntpoint_total, npath, DNSIndex_k, ibe_k, iend_k, jbe_k, jend_k, &
                          ixwindowl, ixwindowr, nvar_output, varindex_output, fileprop)
      elseif(ireadVOL.eq.1) then
        call InitVOLTSkplane(ntpoint_total, npath, DNSIndex_k, ibe_k, iend_k, jbe_k, jend_k, &
                             ixwindowl, ixwindowr, nvar_output, varindex_output, fileprop)
      else
        print *, 'ireadVOL should be 0 or 1. STOP'
        stop
      endif

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

      ! Reading k plane data
      allocate(tsdata_ref(ntpoint_total,nypoint_k,imints_k:imaxts_k,1,1))
      allocate(tsdata(ntpoint_total,nypoint_k,imints_k:imaxts_k,1,nvar_output) )  !!!

    if(ireadVOL.eq.0) then
      call ReadTSkplane(kloc_ref,tsdata_ref)

      do k=1, num_kplane
        print *, 'wall-normal location k = ', kplane(k)
        print *, 'streamwise average range i = ', ibe_k, iend_k
        write(unit=fnum,fmt='(I04.4)') kplane(k)


        call ReadTSkplane(kplane(k),tsdata)

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
                numpt(ii) = numpt(ii) + 1
              enddo
            enddo
            forall(ii=-ixwindowl:ixwindowr,numpt(ii).gt.0)
               corr12(ii) = corr12(ii)/dble(ntpoint_total*nypoint_k*numpt(ii))
               ave1(ii) = ave1(ii)/dble(ntpoint_total*nypoint_k*numpt(ii))
               ave2(ii) = ave2(ii)/dble(ntpoint_total*nypoint_k*numpt(ii))
               var1(ii) = var1(ii)/dble(ntpoint_total*nypoint_k*numpt(ii))
               var2(ii) = var2(ii)/dble(ntpoint_total*nypoint_k*numpt(ii))
               rms1(ii) = sqrt( abs(var1(ii)-ave1(ii)**2) )
               rms2(ii) = sqrt( abs(var2(ii)-ave2(ii)**2) )
               corrx(ii) = (corr12(ii) - ave1(ii)*ave2(ii))/(rms1(ii)*rms2(ii))
            end forall

            if(icrossx.eq.1) then
              filename = 'corrxcross_kref'//fnum_ref//'_k'//fnum//'_'//trim(varname_ref(1))//trim(varout_k(m))//'.dat'
              call outputcrossx(corrx,numpt,trim(filename))
            endif

            !allocate(tsdata(ntpoint_total,nypoint_k,imints_k:imaxts_k,1,nvar_output) )  !!!
            !allocate(tsdata_ref(ntpoint_total,nypoint_k,imints_k:imaxts_k,1,1))
            if(icrossxt.eq.1) then
                spect12 = (0.d0,0.d0)
                numpt = 0
                do j = 1, nypoint_k
                do i = 1, nxpoint_k
                   do ii = max(-ixwindowl,imints_k-i), min(ixwindowr,imaxts_k-i)
                      !spect12(:,ii) = spect12(:,ii) + crossspect1d(tsdata_ref(1,i,j,:),tsdata(m,i+ii,j,:))
                      spect12(:,ii) = spect12(:,ii) + crossspect1d(tsdata_ref(:,j,i,1,1),tsdata(:,j,i+ii,1,m))
                      numpt(ii) = numpt(ii) + 1
                   enddo
                enddo
                enddo
                filename = '_kref'//fnum_ref//'_k'//fnum//'_'//trim(varname_ref(1))//trim(varout_k(m))
                call Calcorrxtcross(ixwindowl,ixwindowr,dx_kp,spect12,corrx,filename)
             endif  ! end calculate x-t cross correlation


          enddo ! end m loop
      enddo ! end k loop

    else ! ireadVOL.eq.1
      ! Read reference plane data
      do kk=1, DNSIndex_k%nkplane
        if(kloc_ref.eq.DNSIndex_k%kplane(kk)) then
          kplane_be = kk
        endif
      enddo ! end kk loop
      call ReadVOLTSkplane(kplane_be,kloc_ref,tsdata_ref)

      do k=1, num_kplane
        print *, 'wall-normal location k = ', kplane(k)
        print *, 'streamwise average range i = ', ibe_k, iend_k
        write(unit=fnum,fmt='(I04.4)') kplane(k)

        ! Read k-plane data
        do kk=1, DNSIndex_k%nkplane
          if(kplane(k).eq.DNSIndex_k%kplane(kk)) then
            kplane_be = kk
          endif
        enddo ! end kk loop
        call ReadVOLTSkplane(kplane_be,kplane(k),tsdata)

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
              filename = 'corrxcross_kref'//fnum_ref//'_k'//fnum//'_'//trim(varname_ref(1))//trim(varout_k(m))//'.dat'
              call outputcrossx(corrx,numpt,trim(filename))
            endif
          enddo ! end m loop

      enddo ! end k loop

    endif ! end ireadVOL.eq.0

contains

      subroutine Input()
        integer :: i, n

        npath = 1
        read(*,*)
        read(*,*) ireadVOL
        read(*,*)
        read(*,*) icrossx, icrossxt
        allocate(fileprop(npath))
        read(*,*)
        do n=1, npath
          read(*,*)
          read(*,*) fileprop(n)%file_be, fileprop(n)%file_end, fileprop(n)%file_skip
          read(*,*)
          read(*,'(a)') fileprop(n)%filepath
        enddo
        read(*,*)
        read(*,*)
        read(*,*) nsection, ioverlap, iwindow
        read(*,*)
        read(*,*)
        read(*,*) nvar_output
        allocate(varindex_output(nvar_output))
        read(*,*) (varindex_output(i), i = 1, nvar_output)
        read(*,*)
        read(*,*) kloc_ref, varindex_ref(1)
        read(*,*)
        read(*,*) ibe_k, iend_k, jbe_k, jend_k, ixwindowl, ixwindowr
        read(*,*)
        read(*,*) num_kplane
        allocate(kplane(num_kplane))
        read(*,*) (kplane(i), i = 1, num_kplane)
      if(ireadVOL.eq.0) then
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

        open(unit=13,file=fn)
        rewind(13)
        write(13,'(a)') 'variables= x corr'
        do ii = -ixwindowl, ixwindowr
          if(npt(ii).gt.0) write(13,*)   dble(ii)*dx_kp, corr(ii)
        enddo
        close(13)
      end subroutine outputcrossx

!      subroutine outputcrossx_HDF5(corr,npt,fn)
!        implicit none
!        real(8), intent(in) :: corr(-ixwindowl:ixwindowr)
!        integer, intent(in) :: npt(-ixwindowl:ixwindowr)
!        character(*), intent(in) :: fn
!
!        type(tp_hyperslab) :: fsol_crossx
!        real(8), dimension(:,:), allocatable :: buffer_crossx
!
!        fsol_crossx%gname = '/'
!        fsol_crossx%dnum = 2
!        fsol_crossx%rank = 1
!        allocate(fsol_crossx%dname(fsol_crossx%dnum),fsol_crossx%dimsf(fsol_crossx%rank))
!        allocate(fsol_crossx%dimsm(fsol_crossx%rank),fsol_crossx%count(fsol_crossx%rank))
!        allocate(fsol_crossx%offset(fsol_crossx%rank),fsol_crossx%block(fsol_crossx%rank))
!        allocate(fsol_crossx%stride(fsol_crossx%rank))
!        fsol_crossx%dname(1) = 'x'
!        fsol_crossx%dname(2) = 'corr'
!
!        fsol_crossx%dimsf(1) = ixwindowr + ixwindowl + 1
!        fsol_crossx%dimsm(1) = fsol_crossx%dimsf(1)
!        fsol_crossx%block = fsol_crossx%dimsm
!        fsol_crossx%count = 1; fsol_crossx%offset = 0
!        fsol_crossx%stride = 1
!        fsol_crossx%IsHSInitialized = .true.
!
!        allocate(buffer_crossx(fsol_crossx%dimsf(1),2))
!
!        do ii=-ixwindowl, ixwindowr
!          if(npt(ii).gt.0) then
!            buffer_crossx(ii+ixwindowl+1,1) = dble(ii)*dx_kp
!            buffer_crossx(ii+ixwindowl+1,2) = corr(ii)
!          endif
!        enddo
!        fsol_crossx%fname = trim(fn)
!        print *, 'Writing file: ', trim(fsol_crossx%fname)
!        call WriteTSHDF5_1D(fsol_crossx,buffer_crossx)
!
!      end subroutine outputcrossx_HDF5






    end program Calcrosscorr
