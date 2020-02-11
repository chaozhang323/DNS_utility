module modSpecty
   use MFFT1D
   use MTSHDF5
   implicit none
   integer, private ::   nypoint
   real(8), private ::   dys
   real(8), private :: pi
   logical, save :: IsCorryInit =  .false.
!                  iwindow ( windowing type)
   type(tp_hyperslab) :: fsol_corrxy
   real(8), dimension(:,:,:), allocatable :: buffer_corrxy


contains
     subroutine InitCorry(ny, dy_sample)
        implicit none
        integer, intent(in) :: ny
        real(8), intent(in) :: dy_sample

        nypoint = ny
        dys = dy_sample
        pi = 4.*atan(1.)
!        call InitFFT1D(nypoint)
        IsCorryInit = .true.

     end subroutine InitCorry

     function spect1dy(tfun)
      implicit none
      real(8), intent(in) :: tfun(nypoint)
      complex :: spect1dy(nypoint)
      real(8) :: tmean
      integer :: i, j, n, nn, nt
      complex, dimension(nypoint) :: ctmp
      real(8) :: tbin(nypoint)


      ! Calculate the mean and energy of time signal
      tmean = sum(tfun(1:nypoint))/dble(nypoint)

      ! Compute spectrum for time data length of ntans
      tbin(1:nypoint) = tfun(1:nypoint) - tmean
      ctmp(1:nypoint) = tbin(1:nypoint)
      call FFT1DF(ctmp, ctmp)
      spect1dy =  ctmp

    end function spect1dy

    function autospect1dy(tfun)
      implicit none
      real(8), intent(in) :: tfun(nypoint)
      complex :: autospect1dy(nypoint)
      complex :: ctmp(nypoint)

      ctmp = spect1dy(tfun)
      autospect1dy = ctmp*conjg(ctmp)
    end function autospect1dy

    function crossspect1dy(tfun1,tfun2)
      implicit none
      real(8), intent(in) :: tfun1(nypoint), tfun2(nypoint)
      complex :: crossspect1dy(nypoint)
      complex, dimension(nypoint) :: ctmp1, ctmp2

      ctmp1 = spect1dy(tfun1)
      ctmp2 = spect1dy(tfun2)
      crossspect1dy = ctmp1*conjg(ctmp2)

    end function crossspect1dy

    subroutine CalCorrxy(ixhalflenl,ixhalflenr,dxs,spect1d,fn)
      implicit none
      integer, intent(in) :: ixhalflenl, ixhalflenr
      real(8), intent(in) :: dxs
      complex, intent(in) :: spect1d(-ixhalflenl:ixhalflenr,nypoint)
      character(*), intent(in) :: fn
      complex(8), dimension(nypoint) :: ctmp
      real(8), dimension(-ixhalflenl:ixhalflenr,nypoint) :: corrxy
      complex(8), dimension(-ixhalflenl:ixhalflenr,nypoint) :: corrout
      integer :: ywave(nypoint)
      integer ::  i, j

        if(.not.IsCorryInit) then
         print *, 'call InitCorry first'
         stop
      endif
      if(ixhalflenl.lt.0.or.ixhalflenr.lt.0) then
         print *, 'ixhalflenl, ixhalflenr need to be > 0. STOP'
      endif
      open(21, file='corrxy'//trim(fn)//'.dat')
      write(21,'(a)')'variables=x,y,corr'
      write(21,'("Zone T=corrxy, I=", I8, ",J=",I8)') ixhalflenl+ixhalflenr+1,nypoint
      do i=-ixhalflenl,ixhalflenr
         call FFT1DB(spect1d(i,:), ctmp)
         call RearrangeFFT1D(nypoint,ctmp,corrout(i,:),ywave)
     enddo
!     corrxy = real(corrout)/maxval(real(corrout))
     corrxy = real(corrout)/abs(real(corrout(0,nypoint/2))+1.e-30)

     do j=1,nypoint
        do i=-ixhalflenl,ixhalflenr
           write(21,'(3E20.12)') real(i)*dxs, real(ywave(j))*dys, corrxy(i,j)
        enddo
     enddo
     close(21)

    end subroutine CalCorrxy


    subroutine CalCorrxy_HDF5(ixhalflenl,ixhalflenr,iloc,dxs,spect1d,fn,Current_n,Current_i)
      implicit none
      integer, intent(in) :: ixhalflenl, ixhalflenr, iloc, Current_n, Current_i
      real(8), intent(in) :: dxs
      complex, intent(in) :: spect1d(-ixhalflenl:ixhalflenr,nypoint)
      character(*), intent(in) :: fn
      complex(8), dimension(nypoint) :: ctmp
      real(8), dimension(-ixhalflenl:ixhalflenr,nypoint) :: corrxy
      complex(8), dimension(-ixhalflenl:ixhalflenr,nypoint) :: corrout
      integer :: ywave(nypoint)
      integer ::  i, j
      character(4) :: fnum

      if(Current_i.eq.0) then
        fsol_corrxy%IsMultiGroup = .false.
      else
        fsol_corrxy%IsMultiGroup = .true.
      endif

      write(unit=fnum,fmt='(I04.4)') iloc
      print *, 'iloc = ', iloc
      fsol_corrxy%gname = '/i'//fnum

      if(Current_n.eq.1) then
        fsol_corrxy%dnum = 3
        fsol_corrxy%rank = 2
        if(allocated(fsol_corrxy%dname)) deallocate(fsol_corrxy%dname,fsol_corrxy%dimsf,fsol_corrxy%dimsm, &
                    fsol_corrxy%count,fsol_corrxy%offset,fsol_corrxy%block,fsol_corrxy%stride)
        allocate(fsol_corrxy%dname(fsol_corrxy%dnum),fsol_corrxy%dimsf(fsol_corrxy%rank))
        allocate(fsol_corrxy%dimsm(fsol_corrxy%rank),fsol_corrxy%count(fsol_corrxy%rank))
        allocate(fsol_corrxy%offset(fsol_corrxy%rank),fsol_corrxy%block(fsol_corrxy%rank))
        allocate(fsol_corrxy%stride(fsol_corrxy%rank))
        fsol_corrxy%dname(1) = 'x'
        fsol_corrxy%dname(2) = 'y'
        fsol_corrxy%dname(3) = 'corr'

        fsol_corrxy%dimsf(1) = ixhalflenr+ixhalflenl+1
        fsol_corrxy%dimsf(2) = nypoint
        fsol_corrxy%dimsm(1) = fsol_corrxy%dimsf(1)
        fsol_corrxy%dimsm(2) = fsol_corrxy%dimsf(2)
        fsol_corrxy%block  = fsol_corrxy%dimsm
        fsol_corrxy%count  = 1; fsol_corrxy%offset = 0
        fsol_corrxy%stride = 1

        fsol_corrxy%IsHSInitialized = .true.
        if(allocated(buffer_corrxy)) deallocate(buffer_corrxy)
        allocate(buffer_corrxy(fsol_corrxy%dimsf(1),fsol_corrxy%dimsf(2),3))
       ! open(21, file='corrxy'//trim(fn)//'.dat')
       ! write(21,'(a)')'variables=x,y,corr'
       ! write(21,'("Zone T=corrxy, I=", I8, ",J=",I8)') ixhalflenl+ixhalflenr+1,nypoint
       ! close(21)
      endif

      if(.not.IsCorryInit) then
         print *, 'call InitCorry first'
         stop
      endif
      if(ixhalflenl.lt.0.or.ixhalflenr.lt.0) then
         print *, 'ixhalflenl, ixhalflenr need to be > 0. STOP'
      endif

      do i=-ixhalflenl,ixhalflenr
         call FFT1DB(spect1d(i,:), ctmp)
         call RearrangeFFT1D(nypoint,ctmp,corrout(i,:),ywave)
     enddo
     ! corrxy = real(corrout)/maxval(real(corrout))
     corrxy = real(corrout)/abs(real(corrout(0,nypoint/2))+1.e-30)

     do j=1,nypoint
        do i=-ixhalflenl,ixhalflenr
          ! write(21,'(3E20.12)') real(i)*dxs, real(ywave(j))*dys, corrxy(i,j)
           buffer_corrxy(i+ixhalflenl+1,j,1) = real(i)*dxs
           buffer_corrxy(i+ixhalflenl+1,j,2) = real(ywave(j))*dys
           buffer_corrxy(i+ixhalflenl+1,j,3) = corrxy(i,j)
        enddo
     enddo
     fsol_corrxy%fname = 'corrxy'//trim(fn)//'.h5'
     print *, 'Writing file: ', trim(fsol_corrxy%fname)
     call WriteTSHDF5_2D(fsol_corrxy,buffer_corrxy)

    end subroutine CalCorrxy_HDF5



end module modSpecty
