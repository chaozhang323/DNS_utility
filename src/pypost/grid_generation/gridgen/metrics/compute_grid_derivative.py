import numpy as np
from util import IO_util
from util import util_f

def main():
    dir_in = '/usr/local/home/yl8bc/duannas/yl8bc/data/cone/run3_polar/REST/'
    filename_in_grid = 'grid.h5'
    stencilsize = 6
    bc = [0,0, 1,1, 0,1]
    ## end params

    ## read grid
    grid_in = IO_util.read_hdf5(dir_in+filename_in_grid)
    x = grid_in['x']
    y = grid_in['y']
    z = grid_in['z']
    ny,nx,nz = x.shape

    ## extend grid
    x,y,z = util_f.mod_metrics.extend_grid(x,y,z, stencilsize)

    ## compute grid der
    grid_derivative = util_f.mod_metrics.compute_grid_derivative(x,y,z, ny,nx,nz,stencilsize, bc)
    grid_derivative = grid_derivative[:,(1,0,2),:,:,:]
    
    ## condition gridder
    grid_derivative[(0,2), 1, :,:,:] = 0.
    grid_derivative[1, (0,2), :,:,:] = 0.
    grid_derivative = np.repeat(grid_derivative[:,:,stencilsize+ny/2:stencilsize+ny/2+1,:,:], grid_derivative.shape[2], axis=2)
    grid_derivative[:,1, :2,:,:] = 0.
    grid_derivative[:,1,-2:,:,:] = 0.

    ## out
    grid_out = {'dxdi':grid_derivative[0,0,:,:,:], \
                'dxdj':grid_derivative[0,1,:,:,:], \
                'dxdk':grid_derivative[0,2,:,:,:], \
                'dydi':grid_derivative[1,0,:,:,:], \
                'dydj':grid_derivative[1,1,:,:,:], \
                'dydk':grid_derivative[1,2,:,:,:], \
                'dzdi':grid_derivative[2,0,:,:,:], \
                'dzdj':grid_derivative[2,1,:,:,:], \
                'dzdk':grid_derivative[2,2,:,:,:]}
    IO_util.write_hdf5(dir_in+'gridDerivative_f.h5', grid_out)


##
def examine_equal(data, axis):
    tmp = np.swapaxes(data, 0, axis)
    equals = [tmp[(0,Ellipsis)]==tmp[(n,Ellipsis)] for n in range(tmp.shape[0])]
    if np.all(equals):
        print(' All equal along axis=%d'%(axis))
    else:
        print(' Not all equal along axis=%d'%(axis))
        
    return


if __name__=='__main__':
    main()
