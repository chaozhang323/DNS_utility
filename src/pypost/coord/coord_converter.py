from util import util_cyl, IO_util

## params
dir_in = '../'                                      # directory in
filename_in_grid = 'grid.h5'
filename_in_data = 'flowdata_00000000.h5'
dir_out = './'                                      # directory out
filename_out_grid = 'grid.h5'
filename_out_data = 'flowdata_00000000.h5'
if_cart_to_cyl = True                               # if converting cartesian to cylidrical? set True/False
## end params

## read files
grid = IO_util.read_hdf5(dir_in+filename_in_grid)
data = IO_util.read_hdf5(dir_in+filename_in_data)

## convert
if if_cart_to_cyl:
    x, y, z = util_cyl.coord_cart_to_cyl(grid['x'], grid['y'], grid['z'])
    u, v, w = util_cyl.vel_cart_to_cyl(data['u'], data['v'], data['w'], x, y, z)
else:
    x, y, z = util_cyl.coord_cyl_to_cart(grid['x'], grid['y'], grid['z'])
    u, v, w = util_cyl.vel_cyl_to_cart(data['u'], data['v'], data['w'], grid['x'], grid['y'], grid['z'])

##out
grid['x'] = x
grid['y'] = y
grid['z'] = z
data['u'] = u
data['v'] = v
data['w'] = w

IO_util.write_hdf5(dir_out+filename_out_grid, grid)
IO_util.write_hdf5(dir_out+filename_out_data, data)
