from util import IO_util, util_interp
import numpy as np

if __name__ == '__main__':
    ##---params --##
    shape = (512, 128)
    dir_in = './InitBody/'
    filename_delta = 'BL.dat'
    ##
    dir_out = '../../data/grids/'
    filename_in_grid = 'cone_grid_parabolic_part2.h5'
    filename_out = 'cone_para_part2'
    ##---params end----##
    
    ## Input
    grid_in, names = IO_util.read_hdf5(dir_out+filename_in_grid)
    delta = np.loadtxt(dir_in+filename_delta)
    
    ## assess delta
    dist = np.linspace( grid_in['z'][0,0], grid_in['z'][0,-1], shape[1] )
    for i in range(shape[0]):
        grid_in['z'][i,:] = np.interp( 
    

    ## output
    IO_util.write_hdf5(dir_out+filename_out, grid_new)

    
    
