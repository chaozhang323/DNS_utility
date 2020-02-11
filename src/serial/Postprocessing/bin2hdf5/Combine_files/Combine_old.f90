program Combine
  use modParallel
  use modPHDF5
  implicit none

!  integer, parameter :: nvar_output = 10
  integer :: nvar_output, i26Variable
  integer :: nskip

  integer :: i, j, k, n
  integer :: ibe_k, iend_k, jbe_k, jend_k
  integer :: jbe_i, jend_i, kbe_i, kend_i
  integer :: ibe_j, iend_j, kbe_j, kend_j
  integer :: ibe_k_inp, iend_k_inp, iskip_k_inp, jbe_k_inp, jend_k_inp, jskip_k_inp
  integer :: jbe_i_inp, jend_i_inp, jskip_i_inp, kbe_i_inp, kend_i_inp, kskip_i_inp
  integer :: ibe_j_inp, iend_j_inp, iskip_j_inp, kbe_j_inp, kend_j_inp, kskip_j_inp
  integer :: ibuffer_k, jbuffer_i, ibuffer_j
  integer :: num_kplane, num_iplane, num_jplane
  integer :: num_file_kp, num_file_ip, num_file_jp

  integer, dimension(:), allocatable :: kplane, iplane, jplane
  real(8), dimension(:,:,:,:,:), allocatable :: buffer_kp, buffer_ip, buffer_jp
  real(8), dimension(:), allocatable :: buffer_time
  character(4) :: fnum, fnum1, fnum2
  character(8) :: fnum_8

  integer :: num_file_k, num_file_i, num_file_j, num_file
  integer :: iConvert_kp, iConvert_ip, iConvert_jp
  integer :: inew_kp, inew_ip, inew_jp
  integer :: iwrite_time_kp, iwrite_time_ip, iwrite_time_jp

  type(tp_rdwt_hdf5) :: TSData_kp, TSData_ip, TSData_jp, TSTime

  integer :: ntpoint, nt_buffer !!!
  integer :: nxpoint_kp, nypoint_kp, nypoint_ip, nzpoint_ip, nxpoint_jp, nzpoint_jp !!!


  character(400) :: filepath
  character(400) :: fileinput, fileoutput

!  (/'u','v','w','p','t','ux','uy','uz','vx','vy','vz','wx','wy','wz',&
!    'uxx','uxy','uxz','vxy','vyy','vyz','wxz','wyz','wzz'/))
  call Init2decomp()
  call InitHDF5()

  call Input()

  call InitTSDataTime(TSTime)

  num_file = ntpoint/nt_buffer

  if(iConvert_kp.eq.1) then
    if(myid.lt.num_file_kp) then
      if(myid.eq.num_file_kp-1) then
        allocate(buffer_kp(nt_buffer,nypoint_kp,nxpoint_kp - ibuffer_k*(num_file_kp-1),1,nvar_output))
      else
        allocate(buffer_kp(nt_buffer,nypoint_kp,ibuffer_k,1,nvar_output))
      endif

      call InitTSDatakplane(TSData_kp)
      do n=nskip/nt_buffer+1, num_file                      ! num_file
        write(unit=fnum_8,fmt='(I08.8)') n*nt_buffer
        TSData_kp%IsSameGroup = .false.
        TSData_kp%IsSameData  = .false.
        do k=1, num_kplane
          ibe_k = (myid*ibuffer_k)*iskip_k_inp + ibe_k_inp
          iend_k = ibe_k + (ibuffer_k - 1)*iskip_k_inp
          if(iend_k.gt.iend_k_inp) then
            iend_k = iend_k_inp
          endif

!          TSData_kp%offset(1) = (n-1)*nt_buffer + nskip
          TSData_kp%offset(1) = (n-1)*nt_buffer
          TSData_kp%offset(2) = 0
          TSData_kp%offset(3) = 0
          TSData_kp%offset(4) = 0

          write(unit=fnum,fmt='(I04.4)')  kplane(k)
          write(unit=fnum1,fmt='(I04.4)') ibe_k
          write(unit=fnum2,fmt='(I04.4)') iend_k
          TSData_kp%fname = trim(filepath)//'timeseries_kplane'//fnum//'_i'//fnum1//'-'//fnum2//'.h5'

          call serial_ReadHDF5_4D(TSData_kp, buffer_kp)

!print *, 'p = ', buffer_kp(1,1,1,1,4)

          if(k.eq.1.and.iwrite_time_kp.eq.1) then
            TSTime%fname = trim(TSData_kp%fname)
!            TSTime%offset(1) = (n-1)*nt_buffer + nskip
            TSTime%offset(1) = (n-1)*nt_buffer
            call ReadHDF5_1D_1v(TSTime, buffer_time)
          endif

          TSData_kp%offset(1) = 0
          TSData_kp%offset(2) = 0
          TSData_kp%offset(3) = ibuffer_k*(myid)
          TSData_kp%offset(4) = k-1

          TSData_kp%fname = trim(fileoutput)//'timeseries_'//fnum_8//'.h5'
          if(myid.eq.0.and.k.eq.1) print *, 'writing file: ', trim(TSData_kp%fname)
          if(inew_kp.eq.0) TSData_kp%IsMultiGroup = .true.
          if(k.gt.1) then
            TSData_kp%IsSameGroup = .true.
            TSData_kp%IsSameData  = .true.
          endif

          call WriteHDF5_4D(TSData_kp, buffer_kp)
          if(k.eq.1) then
            call WriteAttribute(TSData_kp)
            TSTime%offset(1) = 0
            TSTime%fname = trim( TSData_kp%fname)
            if(iwrite_time_kp.eq.1) call WriteHDF5_1D_1v(TSTime, buffer_time)
          endif ! end if
        enddo ! end k loop
      enddo ! end n loop
    endif ! end myid.lt.mum_file_kplane

  endif ! end iConvert_kplane.eq.1


  if(iConvert_ip.eq.1) then
    if(myid.lt.num_file_ip) then
      if(myid.eq.num_file_ip-1) then
        allocate(buffer_ip(nt_buffer,nypoint_ip - jbuffer_i*(num_file_ip-1),nzpoint_ip,1,nvar_output))
      else
        allocate(buffer_ip(nt_buffer,jbuffer_i,nzpoint_ip,1,nvar_output))
      endif

      if(i26Variable.eq.0) then
        call InitTSDataiplane(TSData_ip)
      else
        call InitTSDataiplane_26Variables(TSData_ip)
      endif
      do n=nskip/nt_buffer+1, num_file                      ! num_file
        write(unit=fnum_8,fmt='(I08.8)') n*nt_buffer
        TSData_ip%IsSameGroup = .false.
        TSData_ip%IsSameData  = .false.
        do i=1, num_iplane

          jbe_i = (myid*jbuffer_i)*jskip_i_inp + jbe_i_inp
          jend_i = jbe_i + (jbuffer_i - 1)*jskip_i_inp
          if(jend_i.gt.jend_i_inp) then
            jend_i = jend_i_inp
          endif

!          TSData_ip%offset(1) = (n-1)*nt_buffer + nskip
          TSData_ip%offset(1) = (n-1)*nt_buffer
          TSData_ip%offset(2) = 0
          TSData_ip%offset(3) = 0
          TSData_ip%offset(4) = 0

          write(unit=fnum,fmt='(I04.4)')  iplane(i)
          write(unit=fnum1,fmt='(I04.4)') jbe_i
          write(unit=fnum2,fmt='(I04.4)') jend_i
          TSData_ip%fname = trim(filepath)//'timeseries_iplane'//fnum//'_j'//fnum1//'-'//fnum2//'.h5'

          call serial_ReadHDF5_4D(TSData_ip, buffer_ip)

          if(i.eq.1.and.iwrite_time_ip.eq.1) then
            TSTime%fname = trim(TSData_ip%fname)
!            TSTime%offset(1) = (n-1)*nt_buffer + nskip
            TSTime%offset(1) = (n-1)*nt_buffer
            call ReadHDF5_1D_1v(TSTime, buffer_time)
          endif

          TSData_ip%offset(1) = 0
          TSData_ip%offset(2) = jbuffer_i*myid
          TSData_ip%offset(3) = 0
          TSData_ip%offset(4) = i-1

          TSData_ip%fname = trim(fileoutput)//'timeseries_'//fnum_8//'.h5'
          if(myid.eq.0.and.i.eq.1) print *, 'writing file: ', trim(TSData_ip%fname)
          if(inew_ip.eq.0) TSData_ip%IsMultiGroup = .true.
          if(i.gt.1) then
            TSData_ip%IsSameGroup = .true.
            TSData_ip%IsSameData  = .true.
          endif
          call WriteHDF5_4D(TSData_ip, buffer_ip)
          if(i.eq.1) then
            call WriteAttribute(TSData_ip)
            TSTime%offset(1) = 0
            TSTime%fname = trim( TSData_ip%fname)
            if(iwrite_time_ip.eq.1) call WriteHDF5_1D_1v(TSTime, buffer_time)
          endif
        enddo ! end i loop
      enddo ! end n loop
    endif
  endif ! end iConvert_ip.eq.1

  if(iConvert_jp.eq.1) then
    if(myid.lt.num_file_jp) then
      if(myid.eq.num_file_jp-1) then
        allocate(buffer_jp(nt_buffer,nxpoint_jp - ibuffer_j*(num_file_jp-1),nzpoint_jp,1,nvar_output))
      else
        allocate(buffer_jp(nt_buffer,ibuffer_j,nzpoint_jp,1,nvar_output))
      endif

      call InitTSDatajplane(TSData_jp)
      do n=nskip/nt_buffer+1, num_file                      ! num_file
        write(unit=fnum_8,fmt='(I08.8)') n*nt_buffer
        TSData_jp%IsSameGroup = .false.
        TSData_jp%IsSameData  = .false.
        do j=1, num_jplane

          ibe_j = (myid*ibuffer_j)*iskip_j_inp + ibe_j_inp
          iend_j = ibe_j + (ibuffer_j - 1)*iskip_j_inp
          if(iend_j.gt.iend_j_inp) then
            iend_j = iend_j_inp
          endif

!          TSData_jp%offset(1) = (n-1)*nt_buffer + nskip
          TSData_jp%offset(1) = (n-1)*nt_buffer
          TSData_jp%offset(2) = 0
          TSData_jp%offset(3) = 0
          TSData_jp%offset(4) = 0


          write(unit=fnum,fmt='(I04.4)')  jplane(j)
          write(unit=fnum1,fmt='(I04.4)') ibe_j
          write(unit=fnum2,fmt='(I04.4)') iend_j
          TSData_jp%fname = trim(filepath)//'timeseries_jplane'//fnum//'_i'//fnum1//'-'//fnum2//'.h5'
!          print *, 'Reading file: ', trim(TSData_jp%fname)

          call serial_ReadHDF5_4D(TSData_jp, buffer_jp)

          if(j.eq.1.and.iwrite_time_jp.eq.1) then
            TSTime%fname = trim(TSData_jp%fname)
!            TSTime%offset(1) = (n-1)*nt_buffer + nskip
            TSTime%offset(1) = (n-1)*nt_buffer
            call ReadHDF5_1D_1v(TSTime, buffer_time)
          endif

          TSData_jp%offset(1) = 0
          TSData_jp%offset(2) = ibuffer_j*myid
          TSData_jp%offset(3) = 0
          TSData_jp%offset(4) = j-1

          TSData_jp%fname = trim(fileoutput)//'timeseries_'//fnum_8//'.h5'
          if(myid.eq.0.and.j.eq.1) print *, 'writing file: ', trim(TSData_jp%fname)
          if(inew_jp.eq.0) TSData_jp%IsMultiGroup = .true.
          if(j.gt.1) then
            TSData_jp%IsSameGroup = .true.
            TSData_jp%IsSameData  = .true.
          endif

          call WriteHDF5_4D(TSData_jp, buffer_jp)
          if(j.eq.1) then
            call WriteAttribute(TSData_jp)
            TSTime%offset(1) = 0
            TSTime%fname = trim( TSData_jp%fname)
            if(iwrite_time_jp.eq.1) call WriteHDF5_1D_1v(TSTime, buffer_time)
          endif
        enddo ! end j loop
      enddo ! end n loop
    endif

  endif ! end iConvert_jp.eq.1



  call FinalizeHDF5()
  call Finalize2decomp()



  contains


    subroutine Input()
      integer, parameter :: nid=5
      integer :: i,j

      if(myid.eq.0) then

        read(*,*)
        read(*,*) ntpoint, nskip, nt_buffer, i26Variable
        read(*,*)
        read(*,'(a)') filepath
        read(*,*)
        read(*,'(a)') fileoutput
        read(*,*)
        read(*,*)
        read(*,*) iConvert_kp, inew_kp, iwrite_time_kp
        read(*,*)
        read(*,*) ibe_k_inp, iend_k_inp, iskip_k_inp, jbe_k_inp, jend_k_inp, jskip_k_inp, ibuffer_k
        read(*,*)
        read(*,*) num_kplane
        allocate(kplane(num_kplane))
        read(*,*) (kplane(i), i = 1, num_kplane)
        read(*,*)
        read(*,*)
        read(*,*) iConvert_ip, inew_ip, iwrite_time_ip
        read(*,*)
        read(*,*) jbe_i_inp, jend_i_inp, jskip_i_inp, kbe_i_inp, kend_i_inp, kskip_i_inp, jbuffer_i
        read(*,*)
        read(*,*) num_iplane
        allocate(iplane(num_iplane))
        read(*,*) (iplane(i), i = 1, num_iplane)
        read(*,*)
        read(*,*)
        read(*,*) iConvert_jp, inew_jp, iwrite_time_jp
        read(*,*)
        read(*,*) ibe_j_inp, iend_j_inp, iskip_j_inp, kbe_j_inp, kend_j_inp, kskip_j_inp, ibuffer_j
        read(*,*)
        read(*,*) num_jplane
        allocate(jplane(num_jplane))
        read(*,*) (jplane(j), j = 1, num_jplane)

      endif ! end myid.eq.0

      call MPI_Bcast(ntpoint,      1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(nskip,        1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(nt_buffer,    1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(i26Variable,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(filepath,   400, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(fileoutput, 400, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
      ! kplane
      call MPI_Bcast(iConvert_kp,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(inew_kp,      1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(iwrite_time_kp,1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(ibe_k_inp,    1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(iend_k_inp,   1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(iskip_k_inp,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(jbe_k_inp,    1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(jend_k_inp,   1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(jskip_k_inp,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(ibuffer_k,    1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(num_kplane,   1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      if(.not.allocated(kplane)) allocate(kplane(num_kplane))
      call MPI_Bcast(kplane, num_kplane, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      ! iplane
      call MPI_Bcast(iConvert_ip,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(inew_ip,      1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(iwrite_time_ip,1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(jbe_i_inp,    1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(jend_i_inp,   1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(jskip_i_inp,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(kbe_i_inp,    1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(kend_i_inp,   1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(kskip_i_inp,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(jbuffer_i,    1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(num_iplane,   1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      if(.not.allocated(iplane)) allocate(iplane(num_iplane))
      call MPI_Bcast(iplane, num_iplane, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      ! jplane
      call MPI_Bcast(iConvert_jp,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(inew_jp,      1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(iwrite_time_jp,1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(ibe_j_inp,    1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(iend_j_inp,   1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(iskip_j_inp,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(kbe_j_inp,    1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(kend_j_inp,   1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(kskip_j_inp,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(ibuffer_j,    1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(num_jplane,   1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      if(.not.allocated(jplane)) allocate(jplane(num_jplane))
      call MPI_Bcast(jplane, num_jplane, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

      if(ntpoint.lt.nt_buffer) then
        print *, 'ntpoint should be greater than or equal to nu_buffer'
        print *, 'ntpoint = ', ntpoint, 'nt_buffer = ', nt_buffer
        errcode = 318
        call  MPI_Abort(MPI_COMM_WORLD,errcode,ierr)
      endif

      if(i26Variable.eq.0) then
        nvar_output = 10
      else
        nvar_output = 26
      endif

      if(iConvert_kp.eq.1) then
         if(myid.eq.0) then
         print *, 'number of wall-normal locations =', num_kplane
         print *, 'Wall-normal plane indexes k =', kplane
         endif
         if(mod((iend_k_inp-ibe_k_inp)/iskip_k_inp + 1, ibuffer_k).eq.0) then
           num_file_kp = ((iend_k_inp-ibe_k_inp)/iskip_k_inp +1)/ibuffer_k
         else
           num_file_kp = ((iend_k_inp-ibe_k_inp)/iskip_k_inp +1)/ibuffer_k + 1
         endif
         if(nodes.ne.num_file_kp) then
           print *, 'nodes should be  equal to num_file_kp'
           print *, 'nodes = ', nodes, 'num_file_kp = ', num_file_kp
           errcode = 336
           call  MPI_Abort(MPI_COMM_WORLD,errcode,ierr)
         endif
         nxpoint_kp = (iend_k_inp - ibe_k_inp)/iskip_k_inp + 1
         nypoint_kp = (jend_k_inp - jbe_k_inp)/jskip_k_inp + 1
      endif

      if(iConvert_ip.eq.1) then
         if(myid.eq.0) then
           print *, 'number of const-i locations =', num_iplane
           print *, 'const-i plane indexes i =', iplane
         endif
         if(mod((jend_i_inp-jbe_i_inp)/jskip_i_inp + 1, jbuffer_i).eq.0) then
           num_file_ip = ((jend_i_inp-jbe_i_inp)/jskip_i_inp + 1)/jbuffer_i
         else
           num_file_ip = ((jend_i_inp-jbe_i_inp)/jskip_i_inp + 1)/jbuffer_i + 1
         endif
         if(nodes.ne.num_file_ip) then
           print *, 'nodes should be  equal to num_file_ip'
           print *, 'nodes = ', nodes, 'num_file_ip = ', num_file_ip
           errcode = 356
           call  MPI_Abort(MPI_COMM_WORLD,errcode,ierr)
         endif
         nypoint_ip = (jend_i_inp - jbe_i_inp)/jskip_i_inp + 1
         nzpoint_ip = (kend_i_inp - kbe_i_inp)/kskip_i_inp + 1
      endif

      if(iConvert_jp.eq.1) then
         if(myid.eq.0) then
           print *, 'number of const-j locations =', num_jplane
           print *, 'const-j plane indexes j =', jplane
         endif
         if(mod((iend_j_inp-ibe_j_inp)/iskip_j_inp + 1, ibuffer_j).eq.0) then
           num_file_jp = ((iend_j_inp-ibe_j_inp)/iskip_j_inp + 1)/ibuffer_j
         else
           num_file_jp = ((iend_j_inp-ibe_j_inp)/iskip_j_inp + 1)/ibuffer_j + 1
         endif
         if(nodes.ne.num_file_jp) then
           print *, 'nodes should be  equal to num_file_jp'
           print *, 'nodes = ', nodes, 'num_file_jp = ', num_file_jp
           errcode = 376
           call  MPI_Abort(MPI_COMM_WORLD,errcode,ierr)
         endif
         nxpoint_jp = (iend_j_inp - ibe_j_inp)/iskip_j_inp + 1
         nzpoint_jp = (kend_j_inp - kbe_j_inp)/kskip_j_inp + 1
      endif

    end subroutine input

    subroutine InitTSDataTime(hslab)
      implicit none
      type(tp_rdwt_hdf5), intent(out) :: hslab

      hslab%comm = MPI_COMM_WORLD
      hslab%info = MPI_INFO_NULL
      hslab%gname = '/'
      hslab%dnum = 1
      hslab%rank = 1
      allocate(hslab%dname(hslab%dnum),hslab%dimsf(hslab%rank),hslab%dimsm(hslab%rank))
      allocate(hslab%count(hslab%rank),hslab%offset(hslab%rank),hslab%block(hslab%rank))
      allocate(hslab%stride(hslab%rank))
      hslab%dname(1) = 'time'

      hslab%dimsf = (/nt_buffer/) !!

      hslab%dimsm(1) = nt_buffer
      hslab%block(1) = hslab%dimsm(1)
      hslab%stride = 1
      hslab%count  = 1
      hslab%IsHSInitialized = .true.

      allocate(buffer_time(nt_buffer))

    end subroutine InitTSDataTime

    !
    subroutine InitTSDatakplane_26Variables(hslab)
      implicit none
      type(tp_rdwt_hdf5), intent(out) :: hslab
      integer :: nv = 26

      hslab%comm = MPI_COMM_WORLD
      hslab%info = MPI_INFO_NULL
      hslab%gname = '/kplane'
      hslab%dnum = 26
      hslab%rank = 4
      allocate(hslab%dname(hslab%dnum),hslab%dimsf(hslab%rank),hslab%dimsm(hslab%rank))
      allocate(hslab%count(hslab%rank),hslab%offset(hslab%rank),hslab%block(hslab%rank))
      allocate(hslab%stride(hslab%rank))
      hslab%dname(1)  = "u"
      hslab%dname(2)  = "v"
      hslab%dname(3)  = "w"
      hslab%dname(4)  = "p"
      hslab%dname(5)  = "T"
      hslab%dname(6)  = "ux"
      hslab%dname(7)  = "uy"
      hslab%dname(8)  = "uz"
      hslab%dname(9)  = "vx"
      hslab%dname(10) = "vy"
      hslab%dname(11) = "vz"
      hslab%dname(12) = "wx"
      hslab%dname(13) = "wy"
      hslab%dname(14) = "wz"
      hslab%dname(15) = "uxx"
      hslab%dname(16) = "uxy"
      hslab%dname(17) = "uxz"
      hslab%dname(18) = "vxy"
      hslab%dname(19) = "vyy"
      hslab%dname(20) = "vyz"
      hslab%dname(21) = "wxz"
      hslab%dname(22) = "wyz"
      hslab%dname(23) = "wzz"
      hslab%dname(24) = "uk"
      hslab%dname(25) = "vk"
      hslab%dname(26) = "wk"



      hslab%dimsf = (/nt_buffer,nypoint_kp,nxpoint_kp,num_kplane/)

      hslab%dimsm(1) = nt_buffer
      hslab%dimsm(2) = nypoint_kp
      hslab%dimsm(3) = ibuffer_k
      hslab%dimsm(4) = 1

      if(myid.eq.num_file_kp-1) hslab%dimsm(3) = nxpoint_kp - ibuffer_k*(num_file_kp-1)

      hslab%block(1) = hslab%dimsm(1)
      hslab%block(2) = hslab%dimsm(2)
      hslab%block(3) = hslab%dimsm(3)
      hslab%block(4) = hslab%dimsm(4)

      hslab%stride = 1
      hslab%count  = 1

      hslab%IsHSInitialized = .true.

      ! for attribute
      hslab%attr_name(1) = "number of time steps"
      hslab%attr_name(2) = "istart, iend, iskip"
      hslab%attr_name(3) = "jstart, jend, jskip"
      hslab%attr_name(4) = "K locations"

      hslab%attr_time = nt_buffer
      hslab%info_1st_array(1) = ibe_k_inp
      hslab%info_1st_array(2) = iend_k_inp
      hslab%info_1st_array(3) = iskip_k_inp
      hslab%info_2nd_array(1) = jbe_k_inp
      hslab%info_2nd_array(2) = jend_k_inp
      hslab%info_2nd_array(3) = jskip_k_inp
      hslab%attr_dims_info = 3
      hslab%attr_dims_loc = num_kplane
      allocate(hslab%info_loc(num_kplane))
      hslab%info_loc = kplane

    end subroutine InitTSDatakplane_26Variables



    !
    subroutine InitTSDatakplane(hslab)
      implicit none
      type(tp_rdwt_hdf5), intent(out) :: hslab
      integer :: nv = 10
 
      hslab%comm = MPI_COMM_WORLD
      hslab%info = MPI_INFO_NULL
      hslab%gname = '/kplane'
      hslab%dnum = 10
      hslab%rank = 4
      allocate(hslab%dname(hslab%dnum),hslab%dimsf(hslab%rank),hslab%dimsm(hslab%rank))
      allocate(hslab%count(hslab%rank),hslab%offset(hslab%rank),hslab%block(hslab%rank))
      allocate(hslab%stride(hslab%rank))
      hslab%dname(1) = "u"
      hslab%dname(2) = "v"
      hslab%dname(3) = "w"
      hslab%dname(4) = "p"
      hslab%dname(5) = "T"
      hslab%dname(6) = "uk"
      hslab%dname(7) = "vk"
      hslab%dname(8) = "wk"
      hslab%dname(9) = "pk"
      hslab%dname(10) = "Tk"

      hslab%dimsf = (/nt_buffer,nypoint_kp,nxpoint_kp,num_kplane/)

      hslab%dimsm(1) = nt_buffer
      hslab%dimsm(2) = nypoint_kp
      hslab%dimsm(3) = ibuffer_k
      hslab%dimsm(4) = 1

      if(myid.eq.num_file_kp-1) hslab%dimsm(3) = nxpoint_kp - ibuffer_k*(num_file_kp-1)

      hslab%block(1) = hslab%dimsm(1)
      hslab%block(2) = hslab%dimsm(2)
      hslab%block(3) = hslab%dimsm(3)
      hslab%block(4) = hslab%dimsm(4)

      hslab%stride = 1
      hslab%count  = 1

      hslab%IsHSInitialized = .true.

      ! for attribute
      hslab%attr_name(1) = "number of time steps"
      hslab%attr_name(2) = "istart, iend, iskip"
      hslab%attr_name(3) = "jstart, jend, jskip"
      hslab%attr_name(4) = "K locations"

      hslab%attr_time = nt_buffer
      hslab%info_1st_array(1) = ibe_k_inp
      hslab%info_1st_array(2) = iend_k_inp
      hslab%info_1st_array(3) = iskip_k_inp
      hslab%info_2nd_array(1) = jbe_k_inp
      hslab%info_2nd_array(2) = jend_k_inp
      hslab%info_2nd_array(3) = jskip_k_inp
      hslab%attr_dims_info = 3
      hslab%attr_dims_loc = num_kplane
      allocate(hslab%info_loc(num_kplane))
      hslab%info_loc = kplane

    end subroutine InitTSDatakplane


    subroutine InitTSDataiplane_26Variables(hslab)
      implicit none
      type(tp_rdwt_hdf5), intent(out) :: hslab
      integer :: nv = 26

      hslab%comm = MPI_COMM_WORLD
      hslab%info = MPI_INFO_NULL
      hslab%gname = '/iplane'
      hslab%dnum = 26
      hslab%rank = 4
      allocate(hslab%dname(hslab%dnum),hslab%dimsf(hslab%rank),hslab%dimsm(hslab%rank))
      allocate(hslab%count(hslab%rank),hslab%offset(hslab%rank),hslab%block(hslab%rank))
      allocate(hslab%stride(hslab%rank))
      hslab%dname(1)  = "u"
      hslab%dname(2)  = "v"
      hslab%dname(3)  = "w"
      hslab%dname(4)  = "p"
      hslab%dname(5)  = "T"
      hslab%dname(6)  = "ux"
      hslab%dname(7)  = "uy"
      hslab%dname(8)  = "uz"
      hslab%dname(9)  = "vx"
      hslab%dname(10) = "vy"
      hslab%dname(11) = "vz"
      hslab%dname(12) = "wx"
      hslab%dname(13) = "wy"
      hslab%dname(14) = "wz"
      hslab%dname(15) = "uxx"
      hslab%dname(16) = "uxy"
      hslab%dname(17) = "uxz"
      hslab%dname(18) = "vxy"
      hslab%dname(19) = "vyy"
      hslab%dname(20) = "vyz"
      hslab%dname(21) = "wxz"
      hslab%dname(22) = "wyz"
      hslab%dname(23) = "wzz"
      hslab%dname(24) = "ui"
      hslab%dname(25) = "vi"
      hslab%dname(26) = "wi"

      hslab%dimsf = (/nt_buffer,nypoint_ip,nzpoint_ip,num_iplane/)

      hslab%dimsm(1) = nt_buffer
      hslab%dimsm(2) = jbuffer_i
      hslab%dimsm(3) = nzpoint_ip
      hslab%dimsm(4) = 1

      if(myid.eq.num_file_ip-1) hslab%dimsm(2) = nypoint_ip - jbuffer_i*(num_file_ip-1)

      hslab%block(1) = hslab%dimsm(1)
      hslab%block(2) = hslab%dimsm(2)
      hslab%block(3) = hslab%dimsm(3)
      hslab%block(4) = hslab%dimsm(4)

      hslab%stride = 1
      hslab%count  = 1

      hslab%IsHSInitialized = .true.

      ! for attribute
      hslab%attr_name(1) = "number of time steps"
      hslab%attr_name(2) = "jstart, jend, jskip"
      hslab%attr_name(3) = "kstart, kend, kskip"
      hslab%attr_name(4) = "I locations"

      hslab%attr_time = nt_buffer
      hslab%info_1st_array(1) = jbe_i_inp
      hslab%info_1st_array(2) = jend_i_inp
      hslab%info_1st_array(3) = jskip_i_inp
      hslab%info_2nd_array(1) = kbe_i_inp
      hslab%info_2nd_array(2) = kend_i_inp
      hslab%info_2nd_array(3) = kskip_i_inp
      hslab%attr_dims_info = 3
      hslab%attr_dims_loc = num_iplane
      allocate(hslab%info_loc(num_iplane))
      hslab%info_loc = iplane


    end subroutine InitTSDataiplane_26Variables



    subroutine InitTSDataiplane(hslab)
      implicit none
      type(tp_rdwt_hdf5), intent(out) :: hslab
      integer :: nv = 10

      hslab%comm = MPI_COMM_WORLD
      hslab%info = MPI_INFO_NULL
      hslab%gname = '/iplane'
      hslab%dnum = 10
      hslab%rank = 4
      allocate(hslab%dname(hslab%dnum),hslab%dimsf(hslab%rank),hslab%dimsm(hslab%rank))
      allocate(hslab%count(hslab%rank),hslab%offset(hslab%rank),hslab%block(hslab%rank))
      allocate(hslab%stride(hslab%rank))
      hslab%dname(1) = "u"
      hslab%dname(2) = "v"
      hslab%dname(3) = "w"
      hslab%dname(4) = "p"
      hslab%dname(5) = "T"
      hslab%dname(6) = "ui"
      hslab%dname(7) = "vi"
      hslab%dname(8) = "wi"
      hslab%dname(9) = "pi"
      hslab%dname(10) = "Ti"

      hslab%dimsf = (/nt_buffer,nypoint_ip,nzpoint_ip,num_iplane/)

      hslab%dimsm(1) = nt_buffer
      hslab%dimsm(2) = jbuffer_i
      hslab%dimsm(3) = nzpoint_ip
      hslab%dimsm(4) = 1

      if(myid.eq.num_file_ip-1) hslab%dimsm(2) = nypoint_ip - jbuffer_i*(num_file_ip-1)

      hslab%block(1) = hslab%dimsm(1)
      hslab%block(2) = hslab%dimsm(2)
      hslab%block(3) = hslab%dimsm(3)
      hslab%block(4) = hslab%dimsm(4)

      hslab%stride = 1
      hslab%count  = 1

      hslab%IsHSInitialized = .true.

      ! for attribute
      hslab%attr_name(1) = "number of time steps"
      hslab%attr_name(2) = "jstart, jend, jskip"
      hslab%attr_name(3) = "kstart, kend, kskip"
      hslab%attr_name(4) = "I locations"

      hslab%attr_time = nt_buffer
      hslab%info_1st_array(1) = jbe_i_inp
      hslab%info_1st_array(2) = jend_i_inp
      hslab%info_1st_array(3) = jskip_i_inp
      hslab%info_2nd_array(1) = kbe_i_inp
      hslab%info_2nd_array(2) = kend_i_inp
      hslab%info_2nd_array(3) = kskip_i_inp
      hslab%attr_dims_info = 3
      hslab%attr_dims_loc = num_iplane
      allocate(hslab%info_loc(num_iplane))
      hslab%info_loc = iplane


    end subroutine InitTSDataiplane

    subroutine InitTSDatajplane_26Variables(hslab)
      implicit none
      type(tp_rdwt_hdf5), intent(out) :: hslab
      integer :: nv = 26

      hslab%comm = MPI_COMM_WORLD
      hslab%info = MPI_INFO_NULL
      hslab%gname = '/jplane'
      hslab%dnum = 26
      hslab%rank = 4
      allocate(hslab%dname(hslab%dnum),hslab%dimsf(hslab%rank),hslab%dimsm(hslab%rank))
      allocate(hslab%count(hslab%rank),hslab%offset(hslab%rank),hslab%block(hslab%rank))
      allocate(hslab%stride(hslab%rank))
      hslab%dname(1)  = "u"
      hslab%dname(2)  = "v"
      hslab%dname(3)  = "w"
      hslab%dname(4)  = "p"
      hslab%dname(5)  = "T"
      hslab%dname(6)  = "ux"
      hslab%dname(7)  = "uy"
      hslab%dname(8)  = "uz"
      hslab%dname(9)  = "vx"
      hslab%dname(10) = "vy"
      hslab%dname(11) = "vz"
      hslab%dname(12) = "wx"
      hslab%dname(13) = "wy"
      hslab%dname(14) = "wz"
      hslab%dname(15) = "uxx"
      hslab%dname(16) = "uxy"
      hslab%dname(17) = "uxz"
      hslab%dname(18) = "vxy"
      hslab%dname(19) = "vyy"
      hslab%dname(20) = "vyz"
      hslab%dname(21) = "wxz"
      hslab%dname(22) = "wyz"
      hslab%dname(23) = "wzz"
      hslab%dname(24) = "uj"
      hslab%dname(25) = "vj"
      hslab%dname(26) = "wj"


      hslab%dimsf = (/nt_buffer,nxpoint_jp,nzpoint_jp,num_jplane/)

      hslab%dimsm(1) = nt_buffer
      hslab%dimsm(2) = ibuffer_j
      hslab%dimsm(3) = nzpoint_jp
      hslab%dimsm(4) = 1

      if(myid.eq.num_file_jp-1) hslab%dimsm(2) = nxpoint_jp - ibuffer_j*(num_file_jp-1)

      hslab%block(1) = hslab%dimsm(1)
      hslab%block(2) = hslab%dimsm(2)
      hslab%block(3) = hslab%dimsm(3)
      hslab%block(4) = hslab%dimsm(4)

      hslab%stride = 1
      hslab%count  = 1

      hslab%IsHSInitialized = .true.

      ! for attribute
      hslab%attr_name(1) = "number of time steps"
      hslab%attr_name(2) = "istart, iend, iskip"
      hslab%attr_name(3) = "kstart, kend, kskip"
      hslab%attr_name(4) = "J locations"

      hslab%attr_time = nt_buffer
      hslab%info_1st_array(1) = ibe_j_inp
      hslab%info_1st_array(2) = iend_j_inp
      hslab%info_1st_array(3) = iskip_j_inp
      hslab%info_2nd_array(1) = kbe_j_inp
      hslab%info_2nd_array(2) = kend_j_inp
      hslab%info_2nd_array(3) = kskip_j_inp
      hslab%attr_dims_info = 3
      hslab%attr_dims_loc = num_jplane
      allocate(hslab%info_loc(num_jplane))
      hslab%info_loc = jplane


    end subroutine InitTSDatajplane_26Variables


    subroutine InitTSDatajplane(hslab)
      implicit none
      type(tp_rdwt_hdf5), intent(out) :: hslab
      integer :: nv = 10

      hslab%comm = MPI_COMM_WORLD
      hslab%info = MPI_INFO_NULL
      hslab%gname = '/jplane'
      hslab%dnum = 10
      hslab%rank = 4
      allocate(hslab%dname(hslab%dnum),hslab%dimsf(hslab%rank),hslab%dimsm(hslab%rank))
      allocate(hslab%count(hslab%rank),hslab%offset(hslab%rank),hslab%block(hslab%rank))
      allocate(hslab%stride(hslab%rank))
      hslab%dname(1) = "u"
      hslab%dname(2) = "v"
      hslab%dname(3) = "w"
      hslab%dname(4) = "p"
      hslab%dname(5) = "T"
      hslab%dname(6) = "uj"
      hslab%dname(7) = "vj"
      hslab%dname(8) = "wj"
      hslab%dname(9) = "pj"
      hslab%dname(10) = "Tj"

      hslab%dimsf = (/nt_buffer,nxpoint_jp,nzpoint_jp,num_jplane/)

      hslab%dimsm(1) = nt_buffer
      hslab%dimsm(2) = ibuffer_j
      hslab%dimsm(3) = nzpoint_jp
      hslab%dimsm(4) = 1

      if(myid.eq.num_file_jp-1) hslab%dimsm(2) = nxpoint_jp - ibuffer_j*(num_file_jp-1)

      hslab%block(1) = hslab%dimsm(1)
      hslab%block(2) = hslab%dimsm(2)
      hslab%block(3) = hslab%dimsm(3)
      hslab%block(4) = hslab%dimsm(4)

      hslab%stride = 1
      hslab%count  = 1

      hslab%IsHSInitialized = .true.

      ! for attribute
      hslab%attr_name(1) = "number of time steps"
      hslab%attr_name(2) = "istart, iend, iskip"
      hslab%attr_name(3) = "kstart, kend, kskip"
      hslab%attr_name(4) = "J locations"

      hslab%attr_time = nt_buffer
      hslab%info_1st_array(1) = ibe_j_inp
      hslab%info_1st_array(2) = iend_j_inp
      hslab%info_1st_array(3) = iskip_j_inp
      hslab%info_2nd_array(1) = kbe_j_inp
      hslab%info_2nd_array(2) = kend_j_inp
      hslab%info_2nd_array(3) = kskip_j_inp
      hslab%attr_dims_info = 3
      hslab%attr_dims_loc = num_jplane
      allocate(hslab%info_loc(num_jplane))
      hslab%info_loc = jplane


    end subroutine InitTSDatajplane




end program Combine
