program correlation
  use decomp_2d
  use MPRWHDF5
  implicit none

  real(8), parameter :: R=8314.3D0
  character(400) :: fname, datapath, gridfilename, groupname
  character(8) :: fnum, fnum3, fnum4
  character(4) :: fnum1, fnum2, fnum5, fnum6
  integer :: imax, jmax, kmax
  integer :: file_be, file_end, file_skip, num_file
  real(8) :: rbar, Rm, pi

  type tp_corr
     integer :: ibe, iend, jbe, jend, kbe, kend
     integer :: num_kref
     integer, dimension(:), allocatable :: kref
     integer :: iwinl, iwinr
     integer :: p_row, p_col ! processor grid
     type(DECOMP_INFO) :: decomp, ph, sp
     integer :: iinterp_x, imax_new, iinterp_y, jmax_new, iwindow
  end type tp_corr
  type(tp_corr) :: corr3D

  integer, parameter :: nvarcorr_3D = 3
  integer, dimension(:), allocatable :: varindex
  integer :: num_output, varidx
  character(10) :: dname(14)
  !>                   1   2   3   4   5    6     7   8     9   10   11   12   13   14
  parameter(dname = (/'u','v','w','p','T','rho','P0','T0','ru','rv','rw','uv','uw','vw'/))
  real(8), dimension(:,:,:,:,:), allocatable :: buffer_flow
  real(8), dimension(:,:,:), allocatable :: xx, yy, zz
  complex, dimension(:,:,:,:,:), allocatable :: buffer11_comp, buffer22_comp

  integer :: ierr, errcode, myid, numprocs, kloc
  integer :: kmax_rd, imax_rd, jmax_rd, n, i, ii, j, k, kk, m, nn
  logical, dimension(3) :: periodic_bc
  integer :: dim1_line, dim2_line
  integer :: coords1(2), coords2(2)
  integer :: sumtmp
  type(tp_rdwt_hdf5) :: fsol, fgrd
  integer, dimension(3) :: fft_start, fft_end, fft_size

  ! initialize MPI
  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
  if(myid.eq.0) print *, 'Started reading parameters'
  call Input()
  if(myid.eq.0) print *, 'Finished reading parameters'
  call InitHDF5()

  call InitGridHDF5(fgrd)
  call InitFlowHDF5(fsol)
 
  num_file = (file_end-file_be)/file_skip + 1

  if(myid.eq.0) then
    print *, '##################################################'
    print *, 'Total number of files: ', num_file
    print *, '##################################################'
  endif

  ! do automatic domain decomposition using 2DECOMP&FFT
  call decomp_2d_init(kmax, imax, jmax, corr3D%p_row, corr3D%p_col)
  call decomp_info_init(kmax, imax, jmax, corr3D%decomp)

  fsol%dimsf  = (/kmax, imax, jmax/)
  fsol%dimsm  = (/corr3D%decomp%xsz(1), corr3D%decomp%xsz(2), corr3D%decomp%xsz(3)/)
  fsol%offset(1) = corr3D%decomp%xst(1)-1
  fsol%offset(2) = corr3D%decomp%xst(2)-1
  fsol%offset(3) = corr3D%decomp%xst(3)-1
  fsol%block = fsol%dimsm

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
  fgrd%dimsm  = fsol%dimsm
  fgrd%offset = fsol%offset
  fgrd%block  = fgrd%dimsm

  allocate( xx(corr3D%decomp%xsz(1), corr3D%decomp%xsz(2), corr3D%decomp%xsz(3)), &
            yy(corr3D%decomp%xsz(1), corr3D%decomp%xsz(2), corr3D%decomp%xsz(3)), &
            zz(corr3D%decomp%xsz(1), corr3D%decomp%xsz(2), corr3D%decomp%xsz(3)) )

  if(myid.eq.0) print *, 'Reading grid: ', trim(fgrd%fname)
  print *, 'fgrd%dimsf = ', fgrd%dimsf

!  call ReadHDF5_3D_1V(fgrd,xx,1)
!  call ReadHDF5_3D_1V(fgrd,yy,2)
!  call ReadHDF5_3D_1V(fgrd,zz,3)

!  allocate( buffer_flow(corr3D%decomp%zsz(1),corr3D%decomp%zsz(2),corr3D%decomp%zsz(3), num_file,5) )
!
!  fsol%gname = trim(groupname)
!  do n=1, num_file
!    write(unit=fnum,fmt='(I08.8)') file_be + (n-1)*file_skip
!    fsol%fname = trim(datapath)//fnum//'.h5'
!    if(myid.eq.0)  print *, 'Reading file: ', trim(fsol%fname)
!    !call ReadHDF5_3D(fsol,buffer_flow(:,:,:,n,:))
!    call ReadHDF5_3D_1V(fsol, buffer_flow(:,:,:,n,1),4) ! for test, only read pressure
!  enddo



!  deallocate( xx, yy, zz, buffer_flow )


  call decomp_info_finalize(corr3D%decomp)
  call decomp_2d_finalize
 
  call FinalizeHDF5()
  call MPI_FINALIZE(ierr)

  contains

    subroutine Input()
      integer :: i

      if(myid.eq.0) then
        read(*,*)
        read(*,*) corr3D%p_row, corr3D%p_col
        read(*,*)
        read(*,*)
        read(*,*) imax, jmax, kmax
        read(*,*)
        read(*,*) file_be, file_end, file_skip
        read(*,*)
        read(*,*)
        read(*,'(A)')gridfilename
        read(*,*)
        read(*,'(A)')datapath  !path where data files are stored
        read(*,*)
        read(*,'(A)')groupname
      endif ! myid.eq.0

      call MPI_Bcast(corr3D%p_row,             1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(corr3D%p_col,             1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

      call MPI_Bcast(gridfilename, 400, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(datapath,     400, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(groupname,    400, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(imax,      1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(jmax,      1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(kmax,      1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(file_be,   1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(file_end,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(file_skip, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

      if(numprocs.ne.(corr3D%p_row*corr3D%p_col)) then
        if(myid.eq.0) then
          print *, 'Number of processors allocated: ', numprocs
          print *, 'Please try again and give the processor number as: ', corr3D%p_row*corr3D%p_col
        endif
        call MPI_FINALIZE(ierr)
        stop
      endif

      return
    end subroutine Input


end program correlation
