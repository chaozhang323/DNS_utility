program main
use modRWHDF5
use coord_converter
implicit none

!!----------declare varirables----------!!
type(tp_rdwt_hdf5) :: hgrid, hflow
integer :: imax, jmax, kmax
integer :: icyl2cart
character(len=256) :: filename_grid_out, filename_flow_out, filename_grid_in, filename_flow_in
real(kind=8), dimension(:,:,:,:), allocatable :: grid, flow
real(kind=8), dimension(:,:,:), allocatable :: y,z, v, w 
real(kind=8), dimension(2) :: origin
real(kind=8) :: time
!!----------declare namelist-----------------!!
namelist /inputs/ icyl2cart, &
                  imax, jmax, kmax, &
                  origin


!! ---------- initialization -----------!!
call InitHDF5()
call InitGridHDF5(hgrid)
call InitFlowHDF5(hflow)

open(unit=10, file='coord_converter.inp')
read(10, nml=inputs)
read(10, *)
read(10, *) filename_grid_in
read(10, *) filename_flow_in
read(10, *)
read(10, *) filename_grid_out
read(10, *) filename_flow_out

!filename_grid_out = 'gir'
!write(10, nml=inputs)
close(10)
print*, "kmax, imax, jmax = ", kmax, imax, jmax
print*, "Origin is at ", origin

hgrid%fname = trim(filename_grid_in)
hflow%fname = trim(filename_flow_in)
hgrid%dimsf = (/kmax, imax, jmax/)
hflow%dimsf = (/kmax, imax, jmax/)
hflow%sname = 'time'

allocate(grid(kmax,imax,jmax,3))
allocate(flow(kmax,imax,jmax,5))
allocate(y(kmax,imax,jmax))
allocate(z(kmax,imax,jmax))
allocate(v(kmax,imax,jmax))
allocate(w(kmax,imax,jmax))

!!------------read----------------------!!
print*, "Reading ", filename_grid_in
call ReadHDF5_3D(hgrid, grid)
print*, "Reading ", filename_flow_in
call ReadHDF5_3D(hflow, flow)
call ReadHDF5_scalar(hflow, time)
!!-----------convert--------------------!!
!origin = 0.d0
if(icyl2cart==1) then
    print*, "Starting converting from cylindrical system to Cartesian..."
    !print*, grid(:,:,:,1)
    call coord_cyl_to_cart(grid(:,:,:,2), grid(:,:,:,3), y,z, origin)
    !                      theta          r              y z
    call vel_cyl_to_cart(flow(:,:,:,2), flow(:,:,:,3), v,w, grid(:,:,:,2))
    !                      ut           ur             v w  theta
elseif(icyl2cart==0) then
    print*, "Starting converting from Cartesian system to cylindrical..."
    call coord_cart_to_cyl(grid(:,:,:,2), grid(:,:,:,3), y,    z, origin)
    !                      y              z              theta r
    call vel_cart_to_cyl(flow(:,:,:,2), flow(:,:,:,3), v, w, y)
    !                    v              w              ut ur theta
else
    print*, 'Wrong flag icyl2cart. It is either 0 or 1.'
    stop
endif

!!---------write----------------------!!
grid(:,:,:,2) = y
grid(:,:,:,3) = z 
flow(:,:,:,2) = v
flow(:,:,:,3) = w

hgrid%fname = trim(filename_grid_out)
hflow%fname = trim(filename_flow_out)

print*, "Writing ", filename_grid_out
call WriteHDF5_3D(hgrid, grid)
print*, "Writing ", filename_flow_out
call WriteHDF5_3D(hflow, flow)
call WriteHDF5_scalar(hflow, time)
print*, 'Converted successfully.'

end program
