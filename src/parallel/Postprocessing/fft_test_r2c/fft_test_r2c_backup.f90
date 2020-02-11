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
  
  integer, parameter :: nx=4, ny=2, nz=3
  integer, parameter :: p_row=2, p_col=2
  
  real(mytype), allocatable, dimension(:,:,:) :: in, in2
  complex(mytype), allocatable, dimension(:,:,:) :: out
  
  integer, dimension(3) :: fft_start, fft_end, fft_size
  
  real(mytype), dimension(nx,ny,nz) :: in_global, in_g2, in_g3
  complex(mytype), dimension(nx/2+1,ny,nz) :: out_global
  
  integer (kind=MPI_OFFSET_KIND) :: filesize, disp
  
  real(mytype) :: err
  !integer*8 :: plan
  integer :: fh, ierror, i,j,k, n,iol
  
  call MPI_INIT(ierror)
  call decomp_2d_init(nx,ny,nz,p_row,p_col)
  call decomp_2d_fft_init
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Compute a small problem all on rank 0 as reference
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call random_number(in_global)


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
  
  ! write input to file
  call decomp_2d_write_var(fh,disp,1,in)
  
  if (nrank==0) then
     write(*,*) ' '
     write(*,*) '*** Distributed computation (X-pencil input)'
     write(*,*) ' real input held by rank 0:'
     write(*,20) in
  end if
  
  ! compute r2c transform 
  call decomp_2d_fft_3d(in,out)
  
  if (nrank==0) then
     write(*,*) ' - after forward transform'
     write(*,*) ' complex output held by rank 0:'
     write(*,10) out
  end if

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
  end if



  
  call decomp_2d_fft_finalize
  call decomp_2d_finalize
  call MPI_FINALIZE(ierror)
  
end program fft_test_r2c
