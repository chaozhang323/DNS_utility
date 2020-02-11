!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Main test program for the FFT r2c/c2r interface
!  also demonstrate the use of the IO library 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program fft_test_r2c

  use decomp_2d
  use decomp_2d_fft
  use glassman
  use decomp_2d_io
  
  use MPI
  
  implicit none
  !include "fftw3.f"
  
  integer, parameter :: nx=11, ny=4, nz=1
  integer, parameter :: p_row=4, p_col=1
  

  real(8), parameter :: pi=4.d0*dble(atan(1.0))
  complex(8), parameter :: ii = (0.0, 1.0)
  real(8), parameter :: Lx = 2.d0*pi
  real(8) :: dx = Lx/nx
  real(8), dimension(nx/2+1,1,1) :: kx


  real(mytype), allocatable, dimension(:,:,:) :: in, in2
  complex(mytype), allocatable, dimension(:,:,:) :: out
  
  integer, dimension(3) :: fft_start, fft_end, fft_size
  
  real(mytype), dimension(nx,ny,nz) :: in_global, in_g2, in_g3
  complex(mytype), dimension(nx/2+1,ny,nz) :: out_global
  
  integer (kind=MPI_OFFSET_KIND) :: filesize, disp
  
  real(mytype) :: err
  !integer*8 :: plan
  integer :: fh, ierror, i,j,k, n,iol
  complex(8) :: tmp
  
  call MPI_INIT(ierror)
  call decomp_2d_init(nx,ny,nz,p_row,p_col)
  call decomp_2d_fft_init
  
!  call random_number(in_global)


  ! generage the initial value for in_global
  do k=1, nz
    do j=1, ny
      do i=1, nx
        in_global(i,j,k) = dble(sin((i-1)*dx))
      enddo
    enddo
  enddo
!do i=1, nx
!  print *, 'sin = ', sin((i-1)*dx)
!enddo
!do i=1, nx
!  print *, 'cos = ', cos((i-1)*dx)
!enddo

10 format(1x,6(:,'(',F5.2,',',F5.2,')'))
20 format(1x,6F5.2)


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Test the real-to-complex interface (r2c) 
  
  !  input is X-pencil real    data whose global size is nx*ny*nz
  ! output is Z-pencil complex data whose global size is (nx/2+1)*ny*nz
  allocate (in(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)))
  
  call decomp_2d_fft_get_size(fft_start,fft_end,fft_size)
  allocate (out(fft_start(1):fft_end(1), &
       fft_start(2):fft_end(2), &
       fft_start(3):fft_end(3)))
  
  ! each processor gets its local portion of global data
  do k=xstart(3),xend(3)
     do j=xstart(2),xend(2)
        do i=xstart(1),xend(1)
           in(i,j,k) = in_global(i,j,k)
        end do
     end do
  end do

!print *, 'nrank = ', nrank
!print *, 'xstart(1) = ', xstart(1)
!print *, 'xend(1) = ', xend(1)
!print *, 'xstart(2) = ', xstart(2)
!print *, 'xend(2) = ', xend(2)
!print *, 'xstart(3) = ', xstart(3)
!print *, 'xend(3) = ', xend(3)
!print *, 'xsize(1) = ', xsize(1)
!print *, 'xsize(2) = ', xsize(2)
!print *, 'xsize(3) = ', xsize(3)
!print *, 'fft_start(1) = ', fft_start(1)
!print *, 'fft_end(1) = ', fft_end(1)
!print *, 'fft_start(2) = ', fft_start(2)
!print *, 'fft_end(2) = ', fft_end(2)
!print *, 'fft_start(3) = ', fft_start(3)
!print *, 'fft_end(3) = ', fft_end(3)
!print *, 'fft_size(1) = ', fft_size(1)
!print *, 'fft_size(2) = ', fft_size(2)
!print *, 'fft_size(3) = ', fft_size(3)

   print *, 'nrank = ', nrank
!  if (nrank==0) then
     write(*,*) ' '
     write(*,*) '*** Distributed computation (X-pencil input)'
     write(*,*) ' real input held by rank 0:'
     write(*,20) in
!  end if
  
  ! compute r2c transform 
  call decomp_2d_fft_3d(in,out)
  

   print *, 'nrank = ', nrank
  if (nrank==0) then
     write(*,*) ' - after forward transform'
     write(*,*) ' complex output held by rank 0:'
     write(*,10) out
!     write(*,*) 'size(out,1) =', size(out,1)
!     write(*,*) 'size(out,2) =', size(out,2)
!     write(*,*) 'size(out,3) =', size(out,3)
!     out(1,1,1) = out(1,1,1)*0
!     out(1,1,2) = out(1,1,2)*(0,1.)
!     out(1,1,3) = out(1,1,3)*(0,2.)

!    do i=1, nx/2
!      kx(i,1,1) = 2.d0*pi*(i-1)/Lx
!    enddo
!    kx(nx/2+1,1,1) = 0.0

  endif ! end nrank==0




    do i=fft_start(1), fft_end(1) !!!!!
 !     tmp = ((i-1)*20*3.14)*(0,1)
 !     tmp = (i-1)*(0,1)
      kx(i,1,1) = 2.d0*pi*(i-1)/Lx

      out(i,:,:) = out(i,:,:)*ii*kx(i,1,1)
      print*, 'out(i,1,1) = ', out(i,1,1)

    enddo

!  end if

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Test the complex-to-real interface (c2r)

  allocate (in2(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)))

  ! compute c2r transform
  call decomp_2d_fft_3d(out,in2)

  ! normalisation
  in2 = in2 / real(nx) / real(ny) / real(nz)

  ! write the data recovered by inverse FFT to file
!  call decomp_2d_write_var(fh,disp,1,in2)

  if (nrank==0) then
     write(*,*) ' - after backward transform and normalisation'
     write(*,*) ' real output held by rank 0:'
     write(*,20) in2

     write(*,*) in2
  end if

   deallocate(in,in2,out)













!  call decomp_2d_fft_finalize
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Repeat the above but using Z-pencil input
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  call decomp_2d_fft_init(PHYSICAL_IN_Z)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Test the real-to-complex interface (r2c)

!  allocate (in(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3)))
!  call decomp_2d_fft_get_size(fft_start,fft_end,fft_size)
!  allocate (out(fft_start(1):fft_end(1), &
!       fft_start(2):fft_end(2), &
!       fft_start(3):fft_end(3)))


!  print *, 'nrank = ', nrank
!  print *, 'zstart(1) = ', zstart(1)
!  print *, 'zend(1) = ', zend(1)
!  print *, 'zstart(2) = ', zstart(2)
!  print *, 'zend(2) = ', zend(2)
!  print *, 'zstart(3) = ', zstart(3)
!  print *, 'zend(3) = ', zend(3)
!  print *, 'fft_start(1) = ', fft_start(1)
!  print *, 'fft_end(1) = ', fft_end(1)
!  print *, 'fft_start(2) = ', fft_start(2)
!  print *, 'fft_end(2) = ', fft_end(2)
!  print *, 'fft_start(3) = ', fft_start(3)
!  print *, 'fft_end(3) = ', fft_end(3)


  ! each processor gets its local portion of global data
!  do k=zstart(3),zend(3)
!     do j=zstart(2),zend(2)
!        do i=zstart(1),zend(1)
!           in(i,j,k) = in_global(i,j,k)
!        end do
!     end do
!  end do
!  if (nrank==0) then
!     write(*,*) ' '
!     write(*,*) '*** Distributed computation (Z-pencil input)'
!     write(*,*) ' real input held by rank 0:'
!     write(*,20) in
!  end if

  ! compute r2c transform
!  call decomp_2d_fft_3d(in,out)

!  if (nrank==0) then
!     write(*,*) ' - after forward transform'
!     write(*,*) ' complex output held by rank 0:'
!     write(*,10) out
!  end if
  
  call decomp_2d_fft_finalize
  call decomp_2d_finalize
  call MPI_FINALIZE(ierror)
  
end program fft_test_r2c
