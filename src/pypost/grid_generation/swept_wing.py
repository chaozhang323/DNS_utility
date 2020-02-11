import numpy as np
from util import IO_util

def main():
    ## params
    dims = (32, 100, 140)
    dir_out = '/usr/local/home/yl8bc/yl8bc/data/swept_wing/REST/'
    filename_out_grid = 'grid'
    ## end params
    
    ## define x and z
    xi = np.expand_dims(np.linspace(0., 1., dims[1]), 1)
    zeta = np.expand_dims(np.linspace(0., 10., dims[2]), 0)
    x = xi**2 + 0.5*(1.-zeta)
    z = np.sqrt(2.)*xi*zeta

    ## out
    grid = {'x':x, 'z':z}
    IO_util.write_hdf5(dir_out+filename_out_grid, grid)
    
if __name__=='__main__':
    main()
