program main
! C binding
use, intrinsic :: iso_c_binding
implicit none

double precision, parameter :: pi = 4*ATAN(1.0)
complex, parameter :: ii =(0.0,1.0)

integer(C_INT), parameter :: Nx = 32
integer(C_INT), parameter :: Ny = Nx
double precision, parameter :: Lx = 2*pi, Ly = 2*pi
! Derived paramenter
double precision, parameter :: dx = Lx/Nx, dy = Ly/Ny

real(C_DOUBLE), dimension(Nx,Ny) :: x,y, u0,in,dudx,dudxE, errdU
real(C_DOUBLE), dimension(Nx/2+1,Ny) :: kx, ky

! Fourier space variables
complex(C_DOUBLE_COMPLEX), dimension(Nx/2+1,Ny) :: u_hat_x, out
! indices
integer :: i, j
!---FFTW plans
type(C_PTR) :: pf, pb

! FFTW include
include 'fftw3.f03'

write(*,'(A)') 'Starting...'

! Grid
forall(i=1:Nx,j=1:Ny)
    x(i,j) = (i-1)*dx
    y(i,j) = (j-1)*dy
end forall

! Compute the wavenumbers
forall(i=1:Nx/2,j=1:Ny) kx(i,j)=2*pi*(i-1)/Lx
kx(Nx/2+1,:) = 0.0
forall(i=1:Nx/2+1,j=1:Ny/2) ky(i,j)=2*pi*(j-1)/Ly
forall(i=1:Nx/2+1,j=Ny/2+1:Ny) ky(i,j)=2*pi*(j-Ny-1)/Ly

! Initial Condition
u0 = sin(2*x)
dudxE = 2*cos(2*x)

! Go to Fourier Space
in = u0
pf = fftw_plan_dft_r2c_2d(Ny, Nx, in,out ,FFTW_ESTIMATE)
call fftw_execute_dft_r2c(pf,in,out)
u_hat_x = out

! Derivative
out = ii*kx*out

! Back to physical space
pb = fftw_plan_dft_c2r_2d(Ny, Nx, out,in,FFTW_ESTIMATE)
call fftw_execute_dft_c2r(pb,out,in)

! rescale
dudx = in/Nx/Ny

! check the derivative
errdU = dudx - dudxE

! Write file
write(*,*) 'Writing to files...'

OPEN(UNIT=1, FILE="out_for.dat", ACTION="write", STATUS="replace", &
   FORM="unformatted")
WRITE(1) kx,u0,dudx,errdU,abs(u_hat_x)
CLOSE(UNIT=1)

call fftw_destroy_plan(pf)
call fftw_destroy_plan(pb)

write(*,'(A)') 'Done!'

end program main
