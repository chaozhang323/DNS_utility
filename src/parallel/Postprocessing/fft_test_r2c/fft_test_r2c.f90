!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Main test program for the FFT r2c/c2r interface
!  using z-pencil
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program fft_test_r2c

  use decomp_2d
  use decomp_2d_fft
  use glassman
  use decomp_2d_io
  
  use MPI
  
  implicit none
  !include "fftw3.f"
  
  integer, parameter :: nx=1, ny=1, nz=11
  integer, parameter :: p_row=1, p_col=1
  

  real(8), parameter :: pi=4.d0*dble(atan(1.0))
  real(8), parameter :: Lx = 2.d0*pi
  real(8) :: dx = Lx/nz !!!!!!!!!!!!!!!!!!!!

  real(mytype), allocatable, dimension(:,:,:) :: in, in2
  
  real(mytype), dimension(nx,ny,nz) :: in_global, in_g2, in_g3
  complex(mytype), dimension(nx/2+1,ny,nz) :: out_global
  

  real(mytype) :: err
  integer :: fh, ierror, i,j,k, n,iol
  complex(8) :: tmp
  
  call MPI_INIT(ierror)
  call decomp_2d_init(nx,ny,nz,p_row,p_col)

  ! generage the initial value for in_global
!  do k=1, nz
!    do j=1, ny
!      do i=1, nx
!        in_global(i,j,k) = dble(sin((i-1)*dx))  !!
!      enddo
!    enddo
!  enddo


10 format(1x,6(:,'(',F5.2,',',F5.2,')'))
20 format(1x,6F5.2)


!  allocate (in(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)))
!  allocate (in2(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)))
!  ! each processor gets its local portion of global data
!  do k=xstart(3),xend(3)
!     do j=xstart(2),xend(2)
!        do i=xstart(1),xend(1)
!           in(i,j,k) = in_global(i,j,k)
!        end do
!     end do
!  end do

!  call ddy_fft1st(nx,ny,nz,in,Lx,in2)

!  if (nrank==0) then
!     write(*,*) ' - after backward transform and normalisation'
!     write(*,*) ' real output held by rank 0:'
!     write(*,20) in2

!     write(*,*) in2
!  end if

!  deallocate(in,in2)
  

  do k=1, nz
    do j=1, ny
      do i=1, nx
        in_global(i,j,k) = dble(sin((k-1)*dx))  !!
      enddo
    enddo
  enddo

  allocate(in(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3)))
  allocate(in2(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3)))
  ! each processor gets its local portion of global data
  do k=zstart(3),zend(3)
     do j=zstart(2),zend(2)
        do i=zstart(1),zend(1)
           in(i,j,k) = in_global(i,j,k)
        end do
     end do
  end do

  if (nrank==0) then
     write(*,*) ' '
     write(*,*) '*** Distributed computation (X-pencil input)'
     write(*,*) ' real input held by rank 0:'
     write(*,20) in
  end if


  call ddy_fft3rd(nx,ny,nz,in,Lx,in2)

  if (nrank==0) then
     write(*,*) ' - after backward transform and normalisation'
     write(*,*) ' real output held by rank 0:'
     write(*,20) in2
     write(*,*) in2

  end if


  call decomp_2d_finalize
  call MPI_FINALIZE(ierror)
  
contains

  subroutine ddy_fft1st(nx,ny,nz,in,Lx,out)
    implicit none
    integer, intent(in) :: nx, ny, nz
    real(8), intent(in) :: Lx
    real(8), dimension(:,:,:), intent(in) :: in
    real(8), dimension(:,:,:), intent(out) :: out
    complex(8), dimension(:,:,:), allocatable :: out_fft
    real(8), parameter :: pi=4.d0*dble(atan(1.0))
    complex(8), parameter :: ii = (0.0, 1.0)
    integer, dimension(3) :: fft_start, fft_end, fft_size
    integer :: i, j, k
    real(8) :: kk

    call decomp_2d_fft_init
    call decomp_2d_fft_get_size(fft_start,fft_end,fft_size)
    allocate (out_fft(fft_start(1):fft_end(1), &
                      fft_start(2):fft_end(2), &
                      fft_start(3):fft_end(3)))
    call decomp_2d_fft_3d(in,out_fft)

    do i=fft_start(1), fft_end(1)
      kk = 2.d0*pi*(i-1)/Lx
      out_fft(i,:,:) = out_fft(i,:,:)*ii*kk
    enddo

    call decomp_2d_fft_3d(out_fft,out)

    ! normalization
    out = out/real(nx)/real(ny)/real(nz)

    call decomp_2d_fft_finalize
  end subroutine ddy_fft1st


  subroutine ddy_fft3rd(nx,ny,nz,in,Lz,out)
    implicit none
    integer, intent(in) :: nx, ny, nz
    real(8), intent(in) :: Lz
    real(8), dimension(:,:,:), intent(in) :: in
    real(8), dimension(:,:,:), intent(out) :: out
    complex(8), dimension(:,:,:), allocatable :: out_fft
    real(8), parameter :: pi=4.d0*dble(atan(1.0))
    complex(8), parameter :: ii = (0.0, 1.0)
    integer, dimension(3) :: fft_start, fft_end, fft_size
    integer :: i, j, k
    real(8) :: kk

    call decomp_2d_fft_init(PHYSICAL_IN_Z)
    call decomp_2d_fft_get_size(fft_start,fft_end,fft_size)
    allocate (out_fft(fft_start(1):fft_end(1), &
                      fft_start(2):fft_end(2), &
                      fft_start(3):fft_end(3)))
    call decomp_2d_fft_3d(in,out_fft)

    do i=fft_start(3), fft_end(3)
      kk = 2.d0*pi*(i-1)/Lz
      out_fft(:,:,i) = out_fft(:,:,i)*ii*kk
    enddo

    call decomp_2d_fft_3d(out_fft,out)

    ! normalization
    out = out/real(nx)/real(ny)/real(nz)

    call decomp_2d_fft_finalize
  end subroutine ddy_fft3rd







end program fft_test_r2c
