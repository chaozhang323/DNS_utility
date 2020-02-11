!=======================================================================
! This is part of the 2DECOMP&FFT library
! 
! 2DECOMP&FFT is a software framework for general-purpose 2D (pencil) 
! decomposition. It also implements a highly scalable distributed
! one-dimensional Fast Fourier Transform (FFT).
!
! Copyright (C) 2009-2012 Ning Li, the Numerical Algorithms Group (NAG)
!
!=======================================================================

! This is the FFTW (version 3.x) implementation of the FFT library

module decomp2d_fftw
  
  use decomp_2d  ! 2D decomposition module
  use MFFTWindow
  
  implicit none

  include "fftw3.f"
  
  private        ! Make everything private unless declared public

  ! engine-specific global variables
  integer, save :: plan_type = FFTW_MEASURE

  ! FFTW plans
  ! j=1,2,3 corresponds to the 1D FFTs in X,Y,Z direction, respectively
  ! For c2c transforms: 
  !     use plan(-1,j) for  forward transform; 
  !     use plan( 1,j) for backward transform;
  ! For r2c/c2r transforms:
  !     use plan(0,j) for r2c transforms;
  !     use plan(2,j) for c2r transforms;
  integer*8, save :: plan(-1:2,3)

  integer, parameter :: DECOMP_2D_FFT_FORWARD = -1
  integer, parameter :: DECOMP_2D_FFT_BACKWARD = 1

  ! Physical space data can be stored in either X-pencil or Z-pencil
  integer, parameter :: PHYSICAL_IN_X = 1
  integer, parameter :: PHYSICAL_IN_Z = 3
  integer, parameter :: PHYSICAL_IN_Y = 2

  integer, save :: format                 ! input X-pencil or Z-pencil

  ! The libary can only be initialised once
  logical, save :: initialised = .false.

  ! Global size of the FFT
  integer, save :: nx_fft, ny_fft, nz_fft

  ! 2D processor grid
  integer, save, dimension(2) :: dims

  ! Decomposition objects
  TYPE(DECOMP_INFO), save :: ph  ! physical space
  TYPE(DECOMP_INFO), save :: sp  ! spectral space

  public :: fftw_init_xpencil, fftw_init_zpencil, fftw_init_ypencil, &
            fftw_r2c_1m_x, fftw_r2c_1m_z, fftw_r2c_1m_y, &
            fftw_c2r_1m_x, fftw_c2r_1m_z, fftw_c2r_1m_y, &
            fftw_spect_1m_x, fftw_spect_1m_z, fftw_spect_1m_y,&
            fftw_crossspect_1m_x, fftw_crossspect_1m_z, fftw_crossspect_1m_y, &
            fftw_autospect_1m_x, fftw_autospect_1m_z, fftw_autospect_1m_y,&
            fftw_finalize


contains

  ! Initialise the FFT library to perform arbitrary size transforms for x-pencil 3D data
  subroutine fftw_init_xpencil(nx, ny, nz)

    implicit none
    integer, intent(IN) :: nx, ny, nz
    logical, dimension(2) :: dummy_periods
    integer, dimension(2) :: dummy_coords
    integer :: status, errorcode, ierror

    if (initialised) then
       errorcode = 4
       call decomp_2d_abort(errorcode, &
            'FFT library should only be initialised once')
    end if

    format = 1
    nx_fft = nx
    ny_fft = ny
    nz_fft = nz

    ! determine the processor grid in use
    call MPI_CART_GET(DECOMP_2D_COMM_CART_X, 2, &
         dims, dummy_periods, dummy_coords, ierror)

    ! for c2r/r2c interface:
    ! if in physical space, a real array is of size: nx*ny*nz
    ! in spectral space, the complex array is of size:
    !         (nx/2+1)*ny*nz, if PHYSICAL_IN_X
    !      or nx*ny*(nz/2+1), if PHYSICAL_IN_Z

    call decomp_info_init(nx, ny, nz, ph)
    call decomp_info_init(nx/2+1, ny, nz, sp)

    call init_fftw_engine
    initialised = .true.

    return
  end subroutine fftw_init_xpencil

  ! Initialise the FFT library to perform arbitrary size transforms for z-pencil 3D data
  subroutine fftw_init_zpencil(nx, ny, nz)

    implicit none
    integer, intent(IN) :: nx, ny, nz

    logical, dimension(2) :: dummy_periods
    integer, dimension(2) :: dummy_coords
    integer :: status, errorcode, ierror

    if (initialised) then
       errorcode = 4
       call decomp_2d_abort(errorcode, &
            'FFT library should only be initialised once')
    end if

    format = 3
    nx_fft = nx
    ny_fft = ny
    nz_fft = nz

    ! determine the processor grid in use
    call MPI_CART_GET(DECOMP_2D_COMM_CART_X, 2, &
         dims, dummy_periods, dummy_coords, ierror)

    ! for c2r/r2c interface:
    ! if in physical space, a real array is of size: nx*ny*nz
    ! in spectral space, the complex array is of size:
    !         (nx/2+1)*ny*nz, if PHYSICAL_IN_X
    !      or nx*ny*(nz/2+1), if PHYSICAL_IN_Z

    call decomp_info_init(nx, ny, nz, ph)
    call decomp_info_init(nx, ny, nz/2+1, sp)

    call init_fftw_engine
    initialised = .true.

    return
  end subroutine fftw_init_zpencil

  subroutine fftw_init_ypencil(nx, ny, nz)

    implicit none
    integer, intent(IN) :: nx, ny, nz
    logical, dimension(2) :: dummy_periods
    integer, dimension(2) :: dummy_coords
    integer :: status, errorcode, ierror

    if (initialised) then
       errorcode = 4
       call decomp_2d_abort(errorcode, &
            'FFT library should only be initialised once')
    end if

    format = 2
    nx_fft = nx
    ny_fft = ny
    nz_fft = nz

    ! determine the processor grid in use
    call MPI_CART_GET(DECOMP_2D_COMM_CART_X, 2, &
         dims, dummy_periods, dummy_coords, ierror)

    ! for c2r/r2c interface:
    ! if in physical space, a real array is of size: nx*ny*nz
    ! in spectral space, the complex array is of size:
    !         (nx/2+1)*ny*nz, if PHYSICAL_IN_X
    !      or nx*ny*(nz/2+1), if PHYSICAL_IN_Z

    call decomp_info_init(nx, ny, nz, ph)
    call decomp_info_init(nx, ny/2+1, nz, sp)

    call init_fftw_engine
    initialised = .true.

    return
  end subroutine fftw_init_ypencil

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Return the size, starting/ending index of the distributed array 
  !  whose global size is (nx/2+1)*ny*nz, for defining data structures
  !  in r2c and c2r interfaces
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine fftw_get_size(istart, iend, isize)

    implicit none
    integer, dimension(3), intent(OUT) :: istart, iend, isize

    if (format==PHYSICAL_IN_X) then
       istart = sp%zst
       iend   = sp%zen
       isize  = sp%zsz
    else if (format==PHYSICAL_IN_Z) then
       istart = sp%xst
       iend   = sp%xen
       isize  = sp%xsz
    end if

    return
  end subroutine fftw_get_size


 ! Return a FFTW3 plan for multiple 1D r2c FFTs in X direction
  subroutine r2c_1m_x_plan(plan1, decomp_ph, decomp_sp)

    implicit none

    integer*8, intent(OUT) :: plan1
    TYPE(DECOMP_INFO), intent(IN) :: decomp_ph
    TYPE(DECOMP_INFO), intent(IN) :: decomp_sp

    real(mytype), allocatable, dimension(:,:,:) :: a1
    complex(mytype), allocatable, dimension(:,:,:) :: a2

    allocate(a1(decomp_ph%xsz(1),decomp_ph%xsz(2),decomp_ph%xsz(3)))
    allocate(a2(decomp_sp%xsz(1),decomp_sp%xsz(2),decomp_sp%xsz(3)))
    call dfftw_plan_many_dft_r2c(plan1, 1, decomp_ph%xsz(1), &
         decomp_ph%xsz(2)*decomp_ph%xsz(3), a1, decomp_ph%xsz(1), 1, &
         decomp_ph%xsz(1), a2, decomp_sp%xsz(1), 1, decomp_sp%xsz(1), &
         plan_type)
    deallocate(a1,a2)    

    return
  end subroutine r2c_1m_x_plan


  ! Return a FFTW3 plan for multiple 1D c2r FFTs in X direction
  subroutine c2r_1m_x_plan(plan1, decomp_sp, decomp_ph)

    implicit none

    integer*8, intent(OUT) :: plan1
    TYPE(DECOMP_INFO), intent(IN) :: decomp_sp
    TYPE(DECOMP_INFO), intent(IN) :: decomp_ph

    complex(mytype), allocatable, dimension(:,:,:) :: a1
    real(mytype), allocatable, dimension(:,:,:) :: a2

    allocate(a1(decomp_sp%xsz(1),decomp_sp%xsz(2),decomp_sp%xsz(3)))
    allocate(a2(decomp_ph%xsz(1),decomp_ph%xsz(2),decomp_ph%xsz(3)))
    call dfftw_plan_many_dft_c2r(plan1, 1, decomp_ph%xsz(1), &
         decomp_ph%xsz(2)*decomp_ph%xsz(3), a1, decomp_sp%xsz(1), 1, &
         decomp_sp%xsz(1), a2, decomp_ph%xsz(1), 1, decomp_ph%xsz(1), &
         plan_type)
    deallocate(a1,a2)

    return
  end subroutine c2r_1m_x_plan


  ! Return a FFTW3 plan for multiple 1D r2c FFTs in Z direction
  subroutine r2c_1m_z_plan(plan1, decomp_ph, decomp_sp)

    implicit none

    integer*8, intent(OUT) :: plan1
    TYPE(DECOMP_INFO), intent(IN) :: decomp_ph
    TYPE(DECOMP_INFO), intent(IN) :: decomp_sp

    real(mytype), allocatable, dimension(:,:,:) :: a1
    complex(mytype), allocatable, dimension(:,:,:) :: a2

    allocate(a1(decomp_ph%zsz(1),decomp_ph%zsz(2),decomp_ph%zsz(3)))
    allocate(a2(decomp_sp%zsz(1),decomp_sp%zsz(2),decomp_sp%zsz(3)))
    call dfftw_plan_many_dft_r2c(plan1, 1, decomp_ph%zsz(3), &
         decomp_ph%zsz(1)*decomp_ph%zsz(2), a1, decomp_ph%zsz(3), &
         decomp_ph%zsz(1)*decomp_ph%zsz(2), 1, a2, decomp_sp%zsz(3), &
         decomp_sp%zsz(1)*decomp_sp%zsz(2), 1, plan_type)
    deallocate(a1,a2)

    return
  end subroutine r2c_1m_z_plan


  ! Return a FFTW3 plan for multiple 1D c2r FFTs in Z direction
  subroutine c2r_1m_z_plan(plan1, decomp_sp, decomp_ph)

    implicit none

    integer*8, intent(OUT) :: plan1
    TYPE(DECOMP_INFO), intent(IN) :: decomp_sp
    TYPE(DECOMP_INFO), intent(IN) :: decomp_ph

    complex(mytype), allocatable, dimension(:,:,:) :: a1
    real(mytype), allocatable, dimension(:,:,:) :: a2

    allocate(a1(decomp_sp%zsz(1),decomp_sp%zsz(2),decomp_sp%zsz(3)))
    allocate(a2(decomp_ph%zsz(1),decomp_ph%zsz(2),decomp_ph%zsz(3)))

    call dfftw_plan_many_dft_c2r(plan1, 1, decomp_ph%zsz(3), &
         decomp_ph%zsz(1)*decomp_ph%zsz(2), a1, decomp_sp%zsz(3), &
         decomp_sp%zsz(1)*decomp_sp%zsz(2), 1, a2, decomp_ph%zsz(3), &
         decomp_ph%zsz(1)*decomp_ph%zsz(2), 1, plan_type)
    deallocate(a1,a2)     

    return
  end subroutine c2r_1m_z_plan

 ! Return a FFTW3 plan for multiple 1D r2c FFTs in Y direction
  subroutine r2c_1m_y_plan(plan1, decomp_ph, decomp_sp)

    implicit none

    integer*8, intent(OUT) :: plan1
    TYPE(DECOMP_INFO), intent(IN) :: decomp_ph
    TYPE(DECOMP_INFO), intent(IN) :: decomp_sp

    real(mytype), allocatable, dimension(:,:,:) :: a1
    complex(mytype), allocatable, dimension(:,:,:) :: a2

    ! Due to memory pattern of 3D arrays, 1D FFTs along Y have to be
    ! done one Z-plane at a time. So plan for 2D data sets here.
    allocate(a1(decomp_ph%ysz(1),decomp_ph%ysz(2),decomp_ph%ysz(3)))
    allocate(a2(decomp_sp%ysz(1),decomp_sp%ysz(2),decomp_sp%ysz(3)))
    call dfftw_plan_many_dft_r2c(plan1, 1, decomp_ph%ysz(2), decomp_ph%ysz(1), &
         a1, decomp_ph%ysz(2), decomp_ph%ysz(1), 1, a2, decomp_sp%ysz(2), &
         decomp_sp%ysz(1), 1, plan_type)
    deallocate(a1,a2)    

    return
  end subroutine r2c_1m_y_plan


  ! Return a FFTW3 plan for multiple 1D c2r FFTs in X direction
  subroutine c2r_1m_y_plan(plan1, decomp_sp, decomp_ph)

    implicit none

    integer*8, intent(OUT) :: plan1
    TYPE(DECOMP_INFO), intent(IN) :: decomp_sp
    TYPE(DECOMP_INFO), intent(IN) :: decomp_ph

    complex(mytype), allocatable, dimension(:,:,:) :: a1
    real(mytype), allocatable, dimension(:,:,:) :: a2

    allocate(a1(decomp_sp%ysz(1),decomp_sp%ysz(2),decomp_sp%ysz(3)))
    allocate(a2(decomp_ph%ysz(1),decomp_ph%ysz(2),decomp_ph%ysz(3)))
    call dfftw_plan_many_dft_c2r(plan1, 1, decomp_ph%ysz(2), decomp_ph%ysz(1), &
         a1, decomp_sp%ysz(2), decomp_sp%ysz(1), 1, a2, decomp_ph%ysz(2), &
         decomp_ph%ysz(1), 1, plan_type)
    deallocate(a1,a2)   

    return
  end subroutine c2r_1m_y_plan

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  This routine performs one-time initialisations for the FFT engine
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_fftw_engine

    implicit none

    if (nrank==0) then
       write(*,*) ' '
       write(*,*) '***** Using the FFTW (version 3.x) engine *****'
       write(*,*) ' '
    end if

    if (format == PHYSICAL_IN_X) then

       ! For R2C/C2R tranforms
       call r2c_1m_x_plan(plan(0,1), ph, sp)
       call c2r_1m_x_plan(plan(2,1), sp, ph)

    else if (format == PHYSICAL_IN_Z) then

       ! For R2C/C2R tranforms
       call r2c_1m_z_plan(plan(0,3), ph, sp)
       call c2r_1m_z_plan(plan(2,3), sp, ph)

    else if (format == PHYSICAL_IN_Y) then

       ! For R2C/C2R tranforms
       call r2c_1m_y_plan(plan(0,2), ph, sp)
       call c2r_1m_y_plan(plan(2,2), sp, ph)

    end if

    return
  end subroutine init_fftw_engine

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  This routine performs one-time finalisations for the FFT engine
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine fftw_finalize

    implicit none

    call decomp_info_finalize(ph)
    call decomp_info_finalize(sp)

    call dfftw_destroy_plan(plan(0,1))
    call dfftw_destroy_plan(plan(2,1))
    call dfftw_destroy_plan(plan(0,3))
    call dfftw_destroy_plan(plan(2,3))

    initialised = .false.

    return
  end subroutine fftw_finalize


  ! Following routines calculate multiple one-dimensional FFTs to form 
  ! the basis of three-dimensional FFTs.

  ! r2c transform, multiple 1D FFTs in x direction
  subroutine fftw_r2c_1m_x(input, output)

    implicit none

    real(mytype), dimension(:,:,:), intent(IN)  ::  input
    complex(mytype), dimension(:,:,:), intent(OUT) :: output

    call dfftw_execute_dft_r2c(plan(0,1), input, output)

    return

  end subroutine fftw_r2c_1m_x

  ! r2c transform, multiple 1D FFTs in z direction
  subroutine fftw_r2c_1m_z(input, output)

    implicit none

    real(mytype), dimension(:,:,:), intent(IN)  ::  input
    complex(mytype), dimension(:,:,:), intent(OUT) :: output

    call dfftw_execute_dft_r2c(plan(0,3), input, output)

    return

  end subroutine fftw_r2c_1m_z

  ! r2c transform, multiple 1D FFTs in y direction
  subroutine fftw_r2c_1m_y(input, output)

    implicit none

    real(mytype), dimension(:,:,:), intent(IN)  ::  input
    complex(mytype), dimension(:,:,:), intent(OUT) :: output

    call dfftw_execute_dft_r2c(plan(0,2), input, output)

    return

  end subroutine fftw_r2c_1m_y


  ! c2r transform, multiple 1D FFTs in x direction
  subroutine fftw_c2r_1m_x(input, output)

    implicit none

    complex(mytype), dimension(:,:,:), intent(IN)  ::  input
    real(mytype), dimension(:,:,:), intent(OUT) :: output

    call dfftw_execute_dft_c2r(plan(2,1), input, output)

    return

  end subroutine fftw_c2r_1m_x

  ! c2r transform, multiple 1D FFTs in z direction
  subroutine fftw_c2r_1m_z(input, output)

    implicit none

    complex(mytype), dimension(:,:,:), intent(IN) :: input
    real(mytype), dimension(:,:,:), intent(OUT) :: output

    call dfftw_execute_dft_c2r(plan(2,3), input, output)

    return

  end subroutine fftw_c2r_1m_z

  ! c2r transform, multiple 1D FFTs in y direction
  subroutine fftw_c2r_1m_y(input, output)

    implicit none

    complex(mytype), dimension(:,:,:), intent(IN)  ::  input
    real(mytype), dimension(:,:,:), intent(OUT) :: output

    call dfftw_execute_dft_c2r(plan(2,2), input, output)

    return

  end subroutine fftw_c2r_1m_y

  ! Calculate wave-number spectrum along x direction
  subroutine fftw_spect_1m_x(input, output)
    implicit none

    real(mytype), dimension(:,:,:), intent(IN)  ::  input
    complex(mytype), dimension(:,:,:), intent(OUT) :: output
    real(mytype) :: input_mean(size(input,dim=2),size(input,dim=3))
    real(mytype) ::  input_fluct(size(input,dim=1),size(input,dim=2),size(input,dim=3))
    integer :: nx, i

    nx = size(input,dim=1)
    input_mean = sum(input,dim=1)/dble(nx)
    forall(i=1:nx)
       input_fluct(i,:,:) = input(i,:,:) - input_mean(:,:)
    end forall
    call fftw_r2c_1m_x(input_fluct, output)

    return
  end subroutine fftw_spect_1m_x

  subroutine fftw_autospect_1m_x(input, output)
      implicit none
      real(8), dimension(:,:,:), intent(in) :: input
      complex, dimension(:,:,:), intent(out) :: output
      complex :: ctmp(size(output,dim=1),size(output,dim=2),size(output,dim=3))

      call fftw_spect_1m_x(input, ctmp)
      output = ctmp*conjg(ctmp)

      return
  end subroutine fftw_autospect_1m_x

  subroutine fftw_crossspect_1m_x(input1, input2, output)
      implicit none
      real(8), dimension(:,:,:), intent(in) :: input1
      real(8), dimension(:,:,:), intent(in) :: input2
      complex, dimension(:,:,:), intent(out) :: output
      complex :: ctmp1(size(output,dim=1),size(output,dim=2),size(output,dim=3))
      complex :: ctmp2(size(output,dim=1),size(output,dim=2),size(output,dim=3))

      call fftw_spect_1m_x(input1, ctmp1)
      call fftw_spect_1m_x(input2, ctmp2)
      output = ctmp1*conjg(ctmp2)
  end subroutine fftw_crossspect_1m_x

  ! Calculate wave-number spectrum along z direction
  subroutine fftw_spect_1m_z(input, output)
    implicit none

    real(mytype), dimension(:,:,:), intent(IN)  ::  input
    complex(mytype), dimension(:,:,:), intent(OUT) :: output
    real(mytype) :: input_mean(size(input,dim=1),size(input,dim=2))
    real(mytype)  ::  input_fluct(size(input,dim=1),size(input,dim=2),size(input,dim=3))
    integer :: nz, k

    nz = size(input,dim=3)
    input_mean = sum(input,dim=3)/dble(nz)
    forall(k=1:nz)
       input_fluct(:,:,k) = input(:,:,k) - input_mean(:,:)
    end forall
    call fftw_r2c_1m_z(input_fluct, output)

    return
  end subroutine fftw_spect_1m_z

  subroutine fftw_autospect_1m_z(input, output)
      implicit none
      real(8), dimension(:,:,:), intent(in) :: input
      complex, dimension(:,:,:), intent(out) :: output
      complex :: ctmp(size(output,dim=1),size(output,dim=2),size(output,dim=3))

      call fftw_spect_1m_z(input, ctmp)
      output = ctmp*conjg(ctmp)

      return
  end subroutine fftw_autospect_1m_z

  subroutine fftw_crossspect_1m_z(input1, input2, output)
      implicit none
      real(8), dimension(:,:,:), intent(in) :: input1
      real(8), dimension(:,:,:), intent(in) :: input2
      complex, dimension(:,:,:), intent(out) :: output
      complex :: ctmp1(size(output,dim=1),size(output,dim=2),size(output,dim=3))
      complex :: ctmp2(size(output,dim=1),size(output,dim=2),size(output,dim=3))


     !print *, size(input1,dim=1), size(input1,dim=2), size(input1,dim=3)
     !print *, size(input2,dim=1), size(input2,dim=2), size(input2,dim=3)
      call fftw_spect_1m_z(input1, ctmp1)

      call fftw_spect_1m_z(input2, ctmp2)

      output = ctmp1*conjg(ctmp2)

  end subroutine fftw_crossspect_1m_z

  ! Calculate wave-number spectrum along y direction
  subroutine fftw_spect_1m_y(input, output)
    implicit none

    real(mytype), dimension(:,:,:), intent(IN)  ::  input
    complex(mytype), dimension(:,:,:), intent(OUT) :: output
    real(mytype) :: input_mean(size(input,dim=1),size(input,dim=3))
    real(mytype) ::  input_fluct(size(input,dim=1),size(input,dim=2),size(input,dim=3))
    integer :: ny, j

    ny = size(input,dim=2)
    input_mean = sum(input,dim=2)/dble(ny)
    forall(j=1:ny)
       input_fluct(:,j,:) = input(:,j,:) - input_mean(:,:)
    end forall

    call fftw_r2c_1m_y(input_fluct, output)

    return
  end subroutine fftw_spect_1m_y

  subroutine fftw_autospect_1m_y(input, output)
      implicit none
      real(8), dimension(:,:,:), intent(in) :: input
      complex, dimension(:,:,:), intent(out) :: output
      complex :: ctmp(size(output,dim=1),size(output,dim=2),size(output,dim=3))

      call fftw_spect_1m_y(input, ctmp)
      output = ctmp*conjg(ctmp)

      return
  end subroutine fftw_autospect_1m_y

  subroutine fftw_crossspect_1m_y(input1, input2, output)
      implicit none
      real(8), dimension(:,:,:), intent(in) :: input1
      real(8), dimension(:,:,:), intent(in) :: input2
      complex, dimension(:,:,:), intent(out) :: output
      complex :: ctmp1(size(output,dim=1),size(output,dim=2),size(output,dim=3))
      complex :: ctmp2(size(output,dim=1),size(output,dim=2),size(output,dim=3))

      call fftw_spect_1m_y(input1, ctmp1)
      call fftw_spect_1m_y(input2, ctmp2)
      output = ctmp1*conjg(ctmp2)
  end subroutine fftw_crossspect_1m_y

!  !subroutine to rearrange spectra output by 1D FFT so that
!  !the zero-frequency component is at the center of the series.
!  !np (input): number of points in the spectrum
!  !spectin (input): the spectrum output by 1D FFT
!  !spectout (output): the rearranged spectrum
!  !wavenum (output): the corresponding wave number of the output spectrum
!  subroutine RearrangeFFT1D_real(input,output,wavenum)
!     complex, dimension(:,:,:), intent(inout) :: input
!     complex, dimension(np), intent(out) :: spectout
!     integer, dimension(np), intent(out) :: wavenum
!     integer :: ihalf, i
!     ihalf = np/2+1
!     do i=1,np
!        wavenum(i) = ihalf-1+i-np
!        if (i.le.np-ihalf) then
!           spectout(i) = spectin(ihalf+i)
!        else
!           spectout(i) = spectin(i-np+ihalf)
!        end if
!     end do
!  end subroutine RearrangeFFT1D

end module decomp2d_fftw


!module for computing window functions
module MFFTWindow
  use decomp_2d

  integer, private :: npoint, ntype
  real, dimension(:), allocatable :: wcoeff !window function values
  real :: energy_ratio  !=1/window energy
  contains
    !initialize window function computing
    !np (input): number of points used
    !nt (input): type of the window to be computed
    !            0: flat top window
    !            1: Hanning  window
    !            2: modified Hanning  window
    !            3: Hamming window
    subroutine InitFFTWindow(np, nt)
      implicit none
      integer, intent(in) :: np, nt
      npoint = np
      ntype = nt
      if (allocated(wcoeff)) deallocate(wcoeff)
      allocate(wcoeff(npoint))
      call CalWindow()
      if(nrank.eq.0) call WriteWindow('window.dat')
    end subroutine InitFFTWindow

    !subroutine to compute the window functions
    subroutine CalWindow()
      integer :: i
      real :: pi
      pi  = 4.*atan(1.)
      select case(ntype)
      case(0)  !flat top window
        if(nrank.eq.0) print*,'calculating flat top window'
        wcoeff = 1.
      case(1)  !Hanning
        if(nrank.eq.0) print*,'calculating hanning window'
        do i=1, npoint
          wcoeff(i) = 0.5*(1.-cos(2.*pi*real(i-1)/real(npoint)))
        end do
      case(2)  !modified Hanning with flat top
        if(nrank.eq.0) print*,'calculating modified hanning window'
        wcoeff = 1.
        do i=1, npoint/8-1
          wcoeff(i) = 0.5*(1.-cos(8.*pi*real(i-1)/real(npoint)))
        end do
        do i=7*npoint/8+1, npoint
          wcoeff(i) = 0.5*(1.-cos(8.*pi*real(i-1)/real(npoint)))
        end do
      case(3) ! Hamming
        if(nrank.eq.0) print*,'calculating hamming window'
        do i=1, npoint
          wcoeff(i) = 0.54 - 0.46*cos(2.*pi*real(i-1)/real(npoint))
        end do
      case default
        if(nrank.eq.0) print*,'unknow window type with ntype = ', ntype
        stop
      end select
      energy_ratio = 1./(sum(wcoeff*wcoeff)/npoint)
      wcoeff = wcoeff*sqrt(energy_ratio)
    end subroutine CalWindow

    !output window function to a file
    subroutine WriteWindow(fn)
      character(*),intent(in) :: fn
      integer :: i
      open(11,file=fn)
      write(11,'(a)')'variables=i,w'
      do i=1,npoint
        write(11,'(I8, E20.11)')i,wcoeff(i)
      end do
      close(11)
    end subroutine WriteWindow
end module MFFTWindow




