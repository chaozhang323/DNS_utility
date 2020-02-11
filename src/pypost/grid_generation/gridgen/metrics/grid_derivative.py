import numpy as np
from util import util_f

def native(grid_in, stencilsize=6, bc_g=[0,0,0,0,0,0], bc_gd=[0,0,0,0,0,0]):
    x = grid_in['x']
    y = grid_in['y']
    z = grid_in['z']
    nx,ny,nz = x.shape
    x_work,y_work,z_work = util_f.mod_metrics.extend_grid(x,y,z, stencilsize, bc_g)
    grid_out = {'x':x_work, 'y':y_work, 'z':z_work}
    gd_out = util_f.mod_metrics.compute_grid_derivative(x_work,y_work,z_work, nx,ny,nz, 6, bc_gd)
    
    if ny==1:
        grid_derivative[1,1,Ellipsis] = 1.

    return grid_out, gd_out

def analytically_2d_(grid_in, stencilsize=6, bc_g=[0,0,0,0,0,0], bc_gd=[0,0,0,0,0,0]):
    '''
        x,y,z(i,j,k) -- assume j is independent
    '''

    grid_tmp, gd = native(grid_in, stencilsize, bc_g, bc_gd)
    
    # condition
    gd[(0,2), 1, :,:,:] = 0.
    gd[1, (0,2), :,:,:] = 0.
    ny = gd.shape[3]
    gd = np.repeat(gd[:,:,:,ny//2:ny//2+1,:], ny, axis=3)
    if ny==1+2*stencilsize:
        gd[1,1, :,:,:] = 1.

    return gd


def analytically_2d(grid_in, stencilsize=6, bc_g=[0,0,0,0,0,0], bc_gd=[0,0,0,0,0,0]):
    '''
        Assume y is the spanwise/axisymmetric direction.
        In:
            grid_in (dict)
            stencilsize=6 /* for grid extension
            bc=[0,0,0,0,0,0] /* 0 for periodic, 1 for one sided. Order: j,i,k 
        Out:
            x,y,z /* extended grid
            grid_derivative
    '''        
    x = grid_in['x']
    y = grid_in['y']
    z = grid_in['z']
    ny,nx,nz = x.shape

    ## extend grid
    y,x,z = util_f.mod_metrics.extend_grid(y,x,z, stencilsize, bc_g)

    ## compute grid der
    grid_derivative = util_f.mod_metrics.compute_grid_derivative(x,y,z, ny,nx,nz,stencilsize, bc_gd)
    grid_derivative = grid_derivative[:,(1,0,2),:,:,:] # ddy,ddx,ddz => ddx,ddy,ddz
    
    ## condition gridder
    grid_derivative[(0,2), 1, :,:,:] = 0.
    grid_derivative[1, (0,2), :,:,:] = 0.
    grid_derivative = np.repeat(grid_derivative[:,:,stencilsize+ny//2:stencilsize+ny//2+1,:,:], grid_derivative.shape[2], axis=2)
    grid_derivative[:,1, :2,:,:] = 0.
    grid_derivative[:,1,-2:,:,:] = 0.
    if ny==1:
        grid_derivative[1,1,Ellipsis] = 1.

    ## out
    grid_out = {'x':x, 'y':y, 'z':z}
       

    return grid_out, grid_derivative

def pack_grid_derivative(gdd_in):
    gdd_out = {'dxdi':gdd_in[0,0,:,:,:], \
              'dxdj':gdd_in[0,1,:,:,:], \
              'dxdk':gdd_in[0,2,:,:,:], \
              'dydi':gdd_in[1,0,:,:,:], \
              'dydj':gdd_in[1,1,:,:,:], \
              'dydk':gdd_in[1,2,:,:,:], \
              'dzdi':gdd_in[2,0,:,:,:], \
              'dzdj':gdd_in[2,1,:,:,:], \
              'dzdk':gdd_in[2,2,:,:,:]}
    return gdd_out


def examine_equal(data, axis):
    tmp = np.swapaxes(data, 0, axis)
    equals = [tmp[(0,Ellipsis)]==tmp[(n,Ellipsis)] for n in range(tmp.shape[0])]
    if np.all(equals):
        print(' All equal along axis=%d'%(axis))
    else:
        print(' Not all equal along axis=%d'%(axis))
        
    return


