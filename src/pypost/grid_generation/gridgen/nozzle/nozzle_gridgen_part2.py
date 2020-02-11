from util import IO_util, util_tp
import numpy as np
from util import gridgen, interp2d
from nozzle_gridgen_part1 import normalized_flipped
import time

def main():
    ## parameters
    dir_in = '/usr/local/home/yl8bc/duannas/'
    filename_in_grid = 'nozzle_grid_test.h5'
    filename_outlet = 'Test_Incompact3d/HLB_Nozzle/AcousticZgrid_k500.dat'
    ## output
    dir_out = '../../data/nozzle/'
    filename_out = 'nozzle_grid_test_try1'
    ## BELOW NO PARAMS ARE DEFINED
    tstart = time.time()

    ## Read grid(uniform)
    grid, vars_name = IO_util.read_hdf5(dir_out+filename_in_grid)
    u = grid['x']
    v = grid['z']
    dset, zones_name, vars_name = util_tp.read_tp(dir_in+filename_outlet)
    data_outlet, nodemap = util_tp.read_tp_zone(dset, 'SubZone', vars_name)

    ## intersections
    nx,ny = u.shape
    v[-1,:] = v[-1,-1]*normalized_flipped(data_outlet['z'], ny)
    u,v = gridgen.cm_gridgen.gridgen_checker(u,v)
    
    print 'Checker done. Elapsed time: %f'%(time.time()-tstart)

    ## Write grid
    grid['x'] = u
    grid['z'] = v
    IO_util.write_hdf5(dir_out+filename_out, grid)



def draw_lines(u,v, grid_outlet):
    nx,ny = u.shape
    v[-1,:] = v[-1,-1]*normalized_flipped(grid_outlet['SubZone']['z'], ny)
    u1,v1 = gridgen.cm_gridgen.init_tfi(u,v)
    for i in range(1,nx-1):
        for j in range(1,ny-1):
            print i,j
            for m in range(0,nx-1):
                for n in range(0,ny-1):
                    p00 = np.array([u[i,n],v[i,n]])
                    p01 = np.array([u[i,n+1],v[i,n+1]])
                    p10 = np.array([u1[m,j],v1[m,j]])
                    p11 = np.array([u1[m+1,j],v1[m+1,j]])
                    inter, point = interp2d.interp2d.get_intersect_point(p00,p01,p10,p11)
                    if inter==1:
                        u[i,j] = point[0]
                        v[i,j] = point[1]
                        break
                if inter==1:
                    break
    return u,v


if __name__ == '__main__':
    main() 
