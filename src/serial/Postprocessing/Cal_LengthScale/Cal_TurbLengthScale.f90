
!> calculate turbulent Mach number Mt, Reynolds number based on Taylor length scale, Taylor length scale
!>  dissipation length scale

program Cal_TurbLengthScale
  use modRWHDF5
  implicit none

  integer :: i, j, k, ilen, jlen, klen
  integer :: nxpoint, nzpoint, nxpoint_grid, nzpoint_grid
  integer :: ibe, iend, kbe, kend, ibe_grid, iend_grid, kbe_grid, kend_grid

  character(400) :: filepath, gridpath, fname
  real(8), dimension(:,:,:), allocatable :: buffer
  real(8), dimension(:,:,:), allocatable :: buffer_grid
  real(8), dimension(:,:,:), allocatable :: buffer_output

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       integer :: TotZone
       type(Pointer3D), dimension(:), allocatable :: vars
       real(8), dimension(:,:,:,:), allocatable, target :: vartmp
       integer, dimension(:), allocatable :: varidx
       logical :: IsFileNew, IsFileClose
       integer :: ilen, jlen, klen, nvarout, numshare
       integer :: nvars
       integer, dimension(:), allocatable :: ShareVar
       real(8) :: time
       character(600) :: varoutname
       !real(8) :: aa
       character(4) :: fnum

       type(tp_rdwt_hdf5) :: grd, fsol

  call InitHDF5()
  call Input()

  call InitReadInfo()


  call FinalizeHDF5()

contains

       subroutine InitReadInfo()
         implicit none

         fsol%rank = 2
         fsol%dnum = 8
         fsol%gname = 'Stat2d'
         allocate(fsol%dname(fsol%dnum))
         allocate(fsol%dimsf(fsol%rank) )
         fsol%dname(1) = 'uave'
         fsol%dname(2) = 'vave'
         fsol%dname(3) = 'wave'
         fsol%dname(4) = 'u2'
         fsol%dname(5) = 'v2'
         fsol%dname(6) = 'w2'
         fsol%dname(7) = 'rhoave'
         fsol%dname(8) = 'tave'

         grd%rank = 2
         grd%dnum = 3
         grd%gname = '/'
         allocate(grd%dname(grd%dnum))
         allocate(grd%dimsf(grd%rank))
         grd%dname(1) = 'x'
         grd%dname(2) = 'z'
         grd%dname(3) = 'dkdz'
         ! grd%dname(3) = 'dzdz'

         fsol%dimsf = (/nzpoint_grid,nxpoint_grid/)
         fsol%IsHSInitialized = .true.

         grd%dimsf = (/nzpoint_grid,nxpoint_grid/)
         grd%IsHSInitialized = .true.

       end subroutine InitReadInfo


  subroutine Input()
    implicit none

    read(*,*)
    read(*,'(a)') gridpath
    read(*,*)
    read(*,'(a)') filepath
    read(*,*)
    read(*,*) ibe, iend, kbe, kend

    ibe_grid = ibe; iend_grid = iend
    kbe_grid = kbe; kend_grid = kend

    nxpoint = iend - ibe + 1
    nzpoint = kend - kbe + 1
    nxpoint_grid = iend_grid - ibe_grid + 1
    nzpoint_grid = kend_grid - kbe_grid + 1

    print *, '*************************************************************'
    print *, 'average file info'
    print *, 'nxpoint =       ', nxpoint,      '  nzpoint =      ', nzpoint
    print *, '*************************************************************'

  end subroutine Input




end program Cal_TurbLengthScale
