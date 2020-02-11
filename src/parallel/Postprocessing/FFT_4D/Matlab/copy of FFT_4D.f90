

program corr4D
  use decomp_2d
  use MFileIO
  use MPRWHDF5
  use decomp_2d_fft
  use decomp2d_fftw
  implicit none

  real(8), parameter :: R=8314.3D0
  real(8) :: rbar, Rm, pi
  character(8) :: fnum
  character(400) :: fname, datapath, gridfilename, groupname
  integer :: imax, jmax, kmax, kmax_rd, imax_rd, jmax_rd, n, i, ii, j, k, kk, m, nn
  integer :: file_be, file_end, file_skip, num_file

  integer :: num_output, varidx
  integer, dimension(:), allocatable :: varindex
  real(8), dimension(:), allocatable :: ttime
  real(8), dimension(:,:,:), allocatable :: xx, yy, zz, kx, ky, kz, buffer_wt
  real(8), dimension(:,:,:,:,:), allocatable :: buffer_flow
  complex, dimension(:,:,:,:,:), allocatable :: buffer11_comp, buffer22_comp

  integer :: ierr, errcode, myid, numprocs
  real(8) :: dx, dy, dz, dt, Lx, Ly, Lz, Ttotal, omega
  real(8) :: mean_tmp, mean
  integer :: iFFT4D, iFFT3D, itimeave

  type tp_corr
     integer :: ibe, iend, jbe, jend, kbe, kend
     integer :: iwinl, iwinr
     integer :: p_row, p_col ! processor grid
     type(DECOMP_INFO) :: decomp, ph, sp
  end type tp_corr
  type(tp_corr) :: corr3D
  type(tp_rdwt_hdf5) :: fsol, fgrd, fsol_wt, fsol_wt2

  ! initialize MPI
  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
  if(myid.eq.0) print *, 'Started reading parameters'
  call Input()
  if(myid.eq.0) print *, 'Finished reading parameters'

  ! initialize HDF5
  call InitHDF5()
  call InitGridHDF5(fgrd)
  call InitFlowHDF5(fsol)
  call InitHDF5_wt(fsol_wt)
  call InitHDF5_wt2(fsol_wt2)
 
  num_file = (file_end-file_be)/file_skip + 1
  if(myid.eq.0) then
    print *, '##################################################'
    print *, 'Total number of files: ', num_file
    print *, '##################################################'
  endif

  ! do automatic domain decomposition using 2DECOMP&FFT
  kmax_rd = corr3D%kend - corr3D%kbe + 1
  imax_rd = corr3D%iend - corr3D%ibe + 1
  jmax_rd = corr3D%jend - corr3D%jbe + 1
  call decomp_2d_init(kmax_rd, imax_rd, jmax_rd, corr3D%p_row, corr3D%p_col)
  call decomp_info_init(kmax_rd, imax_rd, jmax_rd, corr3D%decomp)

  ! z-pencil is used to read flowdata
  fsol%dimsf  = (/kmax, imax, jmax/)
  fsol%dimsm  = (/corr3D%decomp%zsz(1), corr3D%decomp%zsz(2), corr3D%decomp%zsz(3)/)
  fsol%offset(1) = corr3D%decomp%zst(1)-1 + corr3D%kbe - 1
  fsol%offset(2) = corr3D%decomp%zst(2)-1 + corr3D%ibe - 1
  fsol%offset(3) = corr3D%decomp%zst(3)-1 + corr3D%jbe - 1
  fsol%block = fsol%dimsm

!!########################################################################################################
!!########################### This part is used for grid information #####################################
  fgrd%fname = trim(gridfilename)
  call DetectHDF5(fgrd)
  if(fgrd%dimsf(1).ne.kmax.or.fgrd%dimsf(2).ne.imax.or.fgrd%dimsf(3).ne.jmax) then
    if(myid.eq.0) then
     print *, 'error in file dimension: '
     print *, '###########################'
     print *, 'input dimension: '
     print *, 'imax = ', imax, 'jmax = ', jmax, 'kmax = ', kmax
     print *, '###########################'
     print *, 'file dimension: '
     print *, 'imax = ', fgrd%dimsf(2), 'jmax = ', fgrd%dimsf(3), 'kmax = ', fgrd%dimsf(1)
    endif
    errcode = 94
    call  MPI_Abort(MPI_COMM_WORLD,errcode,ierr)
  endif

  if(myid.eq.0) then
    print *, 'file dimension: imax = ', fgrd%dimsf(2), 'jmax = ', fgrd%dimsf(3), 'kmax = ', fgrd%dimsf(1)
    print *, 'Using x-pencil to read grid file ... '
  endif
  ! x-pencil is used to read grid file
  fgrd%dimsm  = (/corr3D%decomp%xsz(1), corr3D%decomp%xsz(2), corr3D%decomp%xsz(3)/)
  fgrd%offset(1) =  corr3D%decomp%xst(1)-1 + corr3D%kbe - 1
  fgrd%offset(2) =  corr3D%decomp%xst(2)-1 + corr3D%ibe - 1
  fgrd%offset(3) =  corr3D%decomp%xst(3)-1 + corr3D%jbe - 1
  fgrd%block  = fgrd%dimsm

  allocate( xx(corr3D%decomp%xsz(1), corr3D%decomp%xsz(2), corr3D%decomp%xsz(3)), &
            yy(corr3D%decomp%xsz(1), corr3D%decomp%xsz(2), corr3D%decomp%xsz(3)), &
            zz(corr3D%decomp%xsz(1), corr3D%decomp%xsz(2), corr3D%decomp%xsz(3)) )

  if(myid.eq.0) print *, 'Reading grid: ', trim(fgrd%fname)
  call ReadHDF5_3D_1V(fgrd,xx,1)
  call ReadHDF5_3D_1V(fgrd,yy,2)
  call ReadHDF5_3D_1V(fgrd,zz,3)

  if(myid.eq.0) then
    dx = xx(1,2,1) - xx(1,1,1)
    dy = yy(1,1,2) - yy(1,1,1)
    dz = zz(2,1,1) - zz(1,1,1)
    Lx = dx*dble(imax)
    Ly = dy*dble(jmax)
    Lz = dz*dble(kmax)
  endif
  call MPI_Bcast(Lx,        1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(Ly,        1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(Lz,        1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  deallocate(xx,yy,zz)

  allocate( kx( corr3D%decomp%xst(1):corr3D%decomp%xen(1), corr3D%decomp%xst(2):corr3D%decomp%xen(2), corr3D%decomp%xst(3):corr3D%decomp%xen(3) ), &
            ky( corr3D%decomp%xst(1):corr3D%decomp%xen(1), corr3D%decomp%xst(2):corr3D%decomp%xen(2), corr3D%decomp%xst(3):corr3D%decomp%xen(3) ), &
            kz( corr3D%decomp%xst(1):corr3D%decomp%xen(1), corr3D%decomp%xst(2):corr3D%decomp%xen(2), corr3D%decomp%xst(3):corr3D%decomp%xen(3) )  )

  do i=corr3D%decomp%xst(2), corr3D%decomp%xen(2)
    kx(:,i,:)= dble(i-1)*2.d0*pi/Lx
  enddo
  do j=corr3D%decomp%xst(3), corr3D%decomp%xen(3)
    ky(:,:,j)= dble(j-1)*2.d0*pi/Ly
  enddo
  do k=corr3D%decomp%xst(1), corr3D%decomp%xen(1)
    kz(k,:,:)= dble(k-1)*2.d0*pi/Lz
  enddo

  ! using x-pencil to write the dataset
  fsol_wt2%fname = 'WaveNumber_3D.h5'
  if(myid.eq.0)  print *, 'Writing file: ', trim(fsol_wt2%fname)
  fsol_wt2%dimsf = (/kmax_rd,imax_rd,jmax_rd/)
  fsol_wt2%dimsm = (/corr3D%decomp%xsz(1), corr3D%decomp%xsz(2), corr3D%decomp%xsz(3)/)
  fsol_wt2%block = fsol_wt2%dimsm
  fsol_wt2%offset = (/corr3D%decomp%xst(1)-1, corr3D%decomp%xst(2)-1, corr3D%decomp%xst(3)-1/)
  call WriteHDF5_3D_1V(fsol_wt2,kx,1)
  call WriteHDF5_3D_1V(fsol_wt2,ky,2)
  call WriteHDF5_3D_1V(fsol_wt2,kz,3)
  deallocate(kx,ky,kz)
!!########################################################################################################

  fsol%gname = trim(groupname)
  fsol_wt%dimsf = (/kmax_rd,imax_rd,jmax_rd/)
  fsol_wt%dimsm = (/corr3D%decomp%xsz(1), corr3D%decomp%xsz(2), corr3D%decomp%xsz(3)/)
  fsol_wt%block = fsol_wt%dimsm
  fsol_wt%offset = (/corr3D%decomp%xst(1)-1, corr3D%decomp%xst(2)-1, corr3D%decomp%xst(3)-1/)


  allocate( buffer_flow(corr3D%decomp%zsz(1),corr3D%decomp%zsz(2),corr3D%decomp%zsz(3), num_file,1) )
  allocate( buffer11_comp(1:corr3D%decomp%zsz(1),1:corr3D%decomp%zsz(2),corr3D%decomp%zsz(3),1:num_file/2+1,1) )
  allocate( buffer_wt(corr3D%decomp%xsz(1), corr3D%decomp%xsz(2), corr3D%decomp%xsz(3)) )
  allocate( buffer22_comp(corr3D%decomp%xsz(1), corr3D%decomp%xsz(2), corr3D%decomp%xsz(3),num_file/2+1,1) )
  allocate( ttime(num_file) )

  do nn=1, num_output ! u, v, w, p, T
    varidx = varindex(nn)
    if(myid.eq.0) then
      print *, '###########################'
      print *, 'Reading variable: ', trim(fsol%dname(varidx))
      print *, '###########################'
    endif
    ! using z-pencil to read the dataset
    do n=1, num_file
      write(unit=fnum,fmt='(I08.8)') file_be + (n-1)*file_skip
      fsol%fname = trim(datapath)//fnum//'.h5'
      if(myid.eq.0)  print *, 'Reading file: ', trim(fsol%fname)
      call ReadHDF5_3D_1V(fsol, buffer_flow(:,:,:,n,1),varidx) ! for test, only read pressure
      if(myid.eq.0) call ReadHDF5_scalar(fsol,ttime(n))
    enddo ! end n loop

    if(myid.eq.0) then
      dt = ttime(2) - ttime(1)
      Ttotal = dt*dble(num_file)
    endif
    call MPI_Bcast(Ttotal,        1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    ! calculate and reduce the mean value
    mean_tmp = sum(buffer_flow)
    call MPI_ALLREDUCE(mean_tmp,mean,1,MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)
    mean = mean/dble(imax_rd*jmax_rd*kmax_rd*num_file)
    buffer_flow = buffer_flow - mean
    if(myid.eq.0) then
      print *, 'mean value for variable (', trim(fsol%dname(varidx)), ') = ', mean
    endif

!!########################################################################################################
!!########################### This part is the 4D forward FFT ############################################
    if(myid.eq.0) print *, 'Doing 4D fft '

    ! compute 1D real-to-complex transformation in time direction
    call fftw_init_zpencil(kmax_rd, imax_rd, num_file)
    do j=1, corr3D%decomp%zsz(3)
      call fftw_r2c_1m_z( buffer_flow(:,:,j,:,1), buffer11_comp(:,:,j,:,1)  )
    enddo
    call fftw_finalize
    ! normalization
    buffer11_comp = buffer11_comp/dble(num_file)

    call decomp_2d_fft_init(PHYSICAL_IN_Z)
    ! compute 3D complex-to-complex transformation in i, j, k directions
    do n=1, num_file/2+1
      call decomp_2d_fft_3d(buffer11_comp(:,:,:,n,1),buffer22_comp(:,:,:,n,1),DECOMP_2D_FFT_FORWARD)
    enddo
    call decomp_2d_fft_finalize
    ! normalization
    buffer22_comp = buffer22_comp/dble(imax_rd*jmax_rd*kmax_rd)
!!########################################################################################################

!!########################################################################################################
!!########################### using x-pencil to write the dataset ########################################
    do n=1, num_file/2+1
      write(unit=fnum,fmt='(I08.8)') n-1
      fsol_wt%fname = trim(fsol%dname(varidx))//'_mode'//fnum//'.h5'
      if(myid.eq.0)  print *, 'Writing file: ', trim(fsol_wt%fname)

      buffer_wt(:,:,:) = real(buffer22_comp(:,:,:,n,1))
      call WriteHDF5_3D_1V(fsol_wt,buffer_wt,1)
      buffer_wt(:,:,:) = imag(buffer22_comp(:,:,:,n,1))
      call WriteHDF5_3D_1V(fsol_wt,buffer_wt,2)
      buffer_wt(:,:,:) = atan2(imag(buffer22_comp(:,:,:,n,1)),real(buffer22_comp(:,:,:,n,1)))
      call WriteHDF5_3D_1V(fsol_wt,buffer_wt,3)
      buffer_wt(:,:,:) = 2.d0*sqrt( real(buffer22_comp(:,:,:,n,1))**2 + imag(buffer22_comp(:,:,:,n,1))**2 )
      call WriteHDF5_3D_1V(fsol_wt,buffer_wt,4)

      if(myid.eq.0) then
        omega = dble(n-1)/Ttotal
        fsol_wt%sname = 'omega'
        call WriteHDF5_scalar(fsol_wt,omega)
        fsol_wt%sname = 'mean'
        call WriteHDF5_scalar(fsol_wt,mean)
      endif
    enddo ! end n loop
!!########################################################################################################

  enddo ! end nn loop u, v, w, p, T


  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! compute real-to-complex transform
  !  input is z-pencil real    data whose global size is (kmax_rd*imax_rd*jmax_rd)
  ! output is x-pencil complex data whose global size is (kmax_rd*imax_rd*(jmax_rd/2+1))
  ! call decomp_2d_fft_3d(buffer_flow(:,:,:,1,1),spect11tmp(:,:,:,1,1))
  ! compute complex-to-real transform
  ! call decomp_2d_fft_3d(spect11tmp,buffer_flow(:,:,:,1,2))
  ! normalization
  ! buffer_flow = buffer_flow/real(kmax_rd)/real(imax_rd)/(jmax_rd)
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  deallocate( buffer_flow )
  deallocate( buffer11_comp, buffer22_comp )

  call decomp_info_finalize(corr3D%decomp)
  call decomp_2d_finalize

  call FinalizeHDF5()
  call MPI_FINALIZE(ierr)

contains

    subroutine CalWaveNumber()



    end subroutine CalWaveNumber

    subroutine Input()
      integer :: i

      if(myid.eq.0) then
        read(*,*)
        read(*,*) iFFT4D, iFFT3D
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
        read(*,*) file_be, file_end, file_skip, itimeave
        read(*,*)
        read(*,*) Rm
        read(*,*)
        read(*,*)
        read(*,*) corr3D%p_row, corr3D%p_col
        read(*,*)
        read(*,*) corr3D%ibe, corr3D%iend, corr3D%jbe, corr3D%jend, corr3D%kbe, corr3D%kend
        read(*,*)
        read(*,*) num_output
        allocate(varindex(num_output))
        read(*,*) (varindex(i),i=1,num_output)

      endif ! myid.eq.0
      call MPI_Bcast(iFFT4D,    1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(iFFT3D,    1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(gridfilename, 400, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(datapath,     400, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(groupname,    400, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(imax,      1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(jmax,      1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(kmax,      1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(file_be,   1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(file_end,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(file_skip, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(itimeave,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(Rm,        1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr) ! molecular

      call MPI_Bcast(corr3D%p_row,             1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(corr3D%p_col,             1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(corr3D%ibe,               1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(corr3D%iend,              1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(corr3D%jbe,               1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(corr3D%jend,              1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(corr3D%kbe,               1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(corr3D%kend,              1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(num_output,               1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      if(.not.allocated(varindex)) allocate(varindex(num_output))
      call MPI_Bcast(varindex,        num_output, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

      if(numprocs.ne.(corr3D%p_row*corr3D%p_col)) then
        if(myid.eq.0) then
          print *, 'Number of processors allocated: ', numprocs
          print *, 'Please try again and give the processor number as: ', corr3D%p_row*corr3D%p_col
        endif
        call MPI_FINALIZE(ierr)
        stop
      endif

      if(iFFT3D+iFFT4D.gt.1) then
        if(myid.eq.0) then
          print *, 'Only one of the following can be set to 1 ... STOP!!! '
          print *, 'iFFT4D =', iFFT4D, 'iFFT3D =', iFFT3D
        endif
        errcode = 344
        call  MPI_Abort(MPI_COMM_WORLD,errcode,ierr)
      endif

      if((corr3D%ibe.gt.corr3D%iend).or.(corr3D%jbe.gt.corr3D%jend).or.(corr3D%kbe.gt.corr3D%kend)) then
        if(myid.eq.0) then
          print *, 'The index out of range ... '
          print *, 'ibe =', corr3D%ibe, 'iend =', corr3D%iend
          print *, 'jbe =', corr3D%jbe, 'jend =', corr3D%jend
          print *, 'kbe =', corr3D%kbe, 'kend =', corr3D%kend
        endif
        errcode = 336
        call  MPI_Abort(MPI_COMM_WORLD,errcode,ierr)
      endif
      rbar = R/Rm
      pi = 4.d0*atan(1.d0)

      return
    end subroutine Input

    subroutine InitHDF5_wt(hslab)
       type(tp_rdwt_hdf5), intent(out) :: hslab

       hslab%comm = MPI_COMM_WORLD
       hslab%info = MPI_INFO_NULL
       hslab%gname = '/'
       hslab%dnum = 4
       hslab%rank = 3

       allocate(hslab%dname(hslab%dnum), hslab%dimsf(hslab%rank), hslab%dimsm(hslab%rank))
       allocate(hslab%count(hslab%rank), hslab%offset(hslab%rank), hslab%block(hslab%rank), hslab%stride(hslab%rank))

       hslab%dname(1) = 'Re'
       hslab%dname(2) = 'Im'
       hslab%dname(3) = 'phase'
       hslab%dname(4) = 'amp'
       hslab%sname = 'omega'
       hslab%IsHSInitialized = .true.
       ! default value
       hslab%count = 1
       hslab%stride = 1
       hslab%offset = 0
     end subroutine InitHDF5_wt

     subroutine InitHDF5_wt2(hslab)
       type(tp_rdwt_hdf5), intent(out) :: hslab

       hslab%comm = MPI_COMM_WORLD
       hslab%info = MPI_INFO_NULL
       hslab%gname = '/'
       hslab%dnum = 3
       hslab%rank = 3

       allocate(hslab%dname(hslab%dnum), hslab%dimsf(hslab%rank), hslab%dimsm(hslab%rank))
       allocate(hslab%count(hslab%rank), hslab%offset(hslab%rank), hslab%block(hslab%rank), hslab%stride(hslab%rank))

       hslab%dname(1) = 'kx'
       hslab%dname(2) = 'ky'
       hslab%dname(3) = 'kz'
       hslab%IsHSInitialized = .true.
       ! default value
       hslab%count = 1
       hslab%stride = 1
       hslab%offset = 0
     end subroutine InitHDF5_wt2

end program corr4D
