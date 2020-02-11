

    ! program to compute cross-correlation using time series data
    program Calcrosscorr
      use MTSKplane
      use modSpect
      use modSpecty
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
      character(10) :: varname_ref(1)
      integer :: nxave, nyave
      integer :: iii, ii, i, j, k, m, nn
      real(8), dimension(:,:,:,:,:), allocatable :: tsdata_ref, tsdata, tsdatatmp
!      real(8), dimension(:), allocatable :: corr12, ave1, ave2, var1, var2, rms1, rms2
!      real(8), dimension(:), allocatable :: corrx
      integer, dimension(:), allocatable :: numpt
      real(8) :: dt_sample, dx_kp, dy_kp
      character(400) :: filename, fileoutput
!      complex, dimension(:,:), allocatable :: spect12, spect11, spect22
      real(8) :: zloc_kp, dzdk_kp, zloc_ref, dzdk_ref

      type(tp_DNSIndex) :: DNSIndex_k

      integer :: npath, ireadVOL, num_file, ntpoint_total, ntpoint_tmp, buffer_MulVol
      complex, dimension(:,:), allocatable :: spect1dtmp3, spect1dave3
      integer :: kloc_ref_index

      call Input()
      call InitHDF5()

      print *, 'ntpoint_total =', ntpoint_total, 'dt_sample =', dt_sample

!      call InitSpect1d(ntpoint_total, nperiod, iwindow, ioverlap, dt_sample, ntrans, nsection)
!      call InitFFT1D(ntrans)
!      call InitFFTWindow(ntrans,iwindow)

      call InitVOLTSkplane(ntpoint_total, npath, DNSIndex_k, ibe_k, iend_k, jbe_k, jend_k, &
                           ixwindowl, ixwindowr, nvar_output, varindex_output, fileprop)
      call InitFFT1D(nypoint_k)
      call InitCorry(nypoint_k, dy_kp)

      print *, 'Reference wall-normal location k = ', kloc_ref
      write(unit=fnum_ref,fmt='(I04.4)') kloc_ref

      varname_ref(1) = varname_kplane(varindex_ref(1))
      print *, 'Reference variable name is ', trim(varname_ref(1))
      print *, 'Reading Reference kplane ... '

      num_file = (fileprop(1)%file_end-fileprop(1)%file_be)/fileprop(1)%file_skip + 1

      ! Reading k plane data
      allocate(tsdata_ref(ntpoint_total,nypoint_k,imints_k:imaxts_k,1,1))
      allocate(tsdata(ntpoint_total,nypoint_k,imints_k:imaxts_k,1,nvar_output) )
      allocate(tsdatatmp(fileprop(1)%ntpoint,nypoint_k,imints_k:imaxts_k,1,nvar_output) )


      do k=1, num_kplane
        print *, 'wall-normal location k = ', kplane(k)
        write(unit=fnum,fmt='(I04.4)') kplane(k)
        do n=1, num_file

          do iii=ibe_k, iend_k
            if(iii-ibe_k.eq.0) then
              call InitVOLTSkplane_NoPrintInfo(ntpoint_total, npath, DNSIndex_k, iii, iii, jbe_k, jend_k, &
                                   ixwindowl, ixwindowr, nvar_output, varindex_output, fileprop)
            else
              call InitVOLTSkplane_NoPrintInfo(ntpoint_total, npath, DNSIndex_k, iii, iii, jbe_k, jend_k, &
                                               ixwindowl, ixwindowr, nvar_output, varindex_output, fileprop)
            endif

            ntpoint_tmp = 0
            do nn=-buffer_MulVol, buffer_MulVol
              write(unit=fnum_8,fmt='(I08.8)') fileprop(1)%file_be + (n+nn-1)*fileprop(1)%file_skip
              filename = trim(fileprop(1)%filepath)//'timeseries_'//fnum_8//'.h5'
              print *, 'Reading file: ', trim(filename)
              call ReadVOLTSkplane_perFile(filename,kloc_ref_index,tsdata_ref)

              call ReadVOLTSkplane_perFile(filename,k+kplane_be-1,tsdatatmp) !!!!!!!!!
              tsdata((ntpoint_tmp+1):(ntpoint_tmp+fileprop(1)%ntpoint),1:nypoint_k,imints_k:imaxts_k,1,1:nvar_output) = &
              tsdatatmp(1:fileprop(1)%ntpoint,1:nypoint_k,imints_k:imaxts_k,1,1:nvar_output)
              ntpoint_tmp = ntpoint_tmp + fileprop(1)%ntpoint
            enddo ! end nn loop

            allocate(spect1dtmp3(-ixwindowl:ixwindowr,nypoint_k), spect1dave3(-ixwindowl:ixwindowr,nypoint_k))
            if(.not.allocated(numpt)) allocate(numpt(-ixwindowl:ixwindowr))

            do m=1, nvar_output
              fileoutput = '_corrxcross_kref'//fnum_ref//'_k'//fnum//'_'//trim(varname_ref(1))//trim(varout_k(m))
              numpt = 0
              spect1dtmp3 = (0.,0.); spect1dave3 = (0.,0.)
              do i = 1, nxpoint_k
                do j = 1, ntpoint_total
                  do ii = -min(ixwindowl,i-imints_k), min(ixwindowr,imaxts_k-i)
                    ! print *, 'i = ', i, 'j = ', j, 'ii = ', ii
                    spect1dtmp3(ii,:) = crossspect1dy(tsdata(j,1:nypoint_k,i+ii,1,m),tsdata_ref(j,1:nypoint_k,i,1,m))
                    spect1dave3(ii,:) = spect1dave3(ii,:) + spect1dtmp3(ii,:)
                    numpt(ii) = numpt(ii) + 1
                  enddo
                enddo
              enddo
              do i = 1, nxpoint_k
                do j = 1, ntpoint_total
                  do ii = -min(ixwindowl,i-imints_k), min(ixwindowr,imaxts_k-i)
                    ! print *, 'i = ', i, 'j = ', j, 'ii = ', ii
                    spect1dtmp3(ii,:) = crossspect1dy(tsdata_ref(j,1:nypoint_k,i+ii,1,m),tsdata(j,1:nypoint_k,i,1,m))
                    spect1dave3(ii,:) = spect1dave3(ii,:) + spect1dtmp3(ii,:)
                    numpt(ii) = numpt(ii) + 1
                  enddo
                enddo
              enddo
              forall(ii=-ixwindowl:ixwindowr,numpt(ii).gt.0)
                spect1dave3(ii,:) = spect1dave3(ii,:)/dble(numpt(ii))
              end forall
              call CalCorrxy_HDF5(ixwindowl,ixwindowr,iii,dx_kp,spect1dave3,fileoutput,n,iii-ibe_k)
            enddo ! end m loop

          enddo ! end iii loop

        enddo ! end n loop

      enddo ! end k loop

contains

      subroutine Input()
        integer :: i, n

        ! read input file
        read(*,*)
        read(*,*) icrossx, icrossxt, buffer_MulVol
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


!      subroutine outputcrossx(corr,npt,fn)
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
!      end subroutine outputcrossx






    end program Calcrosscorr
