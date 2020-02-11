from util import IO_util, util_tp, util_interp, util_data, util_flow
import numpy as np
from scipy.interpolate import griddata
import time


def main():
    ###params###
    
    mrange = range(1,99,1)
    mmax = len(mrange)
    nmax = 75
    

    dir_in = '/usr/local/home/yl8bc/data_local/cone/run3_Col_native/REST/'
    filename_in_grid = 'grid.h5'
    filename_in_data = 'flowdata_00024000.h5'
    ## files output: 
    dir_out = '/usr/local/home/yl8bc/yl8bc/data/cone/FLUENT/'
    filename_out_profiles= 'profiles_dns_GCL'
    filename_out_integrals= 'integrals_dns_GCL'
    ##--below no param is defined--##
    tstart = time.time()
    keys_change = {'X':'x', 'Y':'z', 'X Velocity':'uave', 'Y Velocity':'wave', 'Pressure':'pave', 'Density':'rhoave', 'Temperature':'tave'}
    grid_keys = ['x','z'] 


    ## read grid and data
    grid = IO_util.read_hdf5(dir_in+filename_in_grid)
    data = IO_util.read_hdf5(dir_in+filename_in_data)
    print('Data read. Time elapsed: %f secs'%(time.time()-tstart))

    ## slice 
    data.update(grid); data.pop('time')
    sdata = util_data.structured_data(data, ['y','x','z'])
    sdata = sdata.slice_data(0, 0)

    ## wall normal profiles
    profiles_out = sdata.get_wall_normal_profiles('right', mrange, nmax)

    ## integral
    delta = np.zeros(mmax)
    delta_star = np.zeros(mmax)
    theta = np.zeros(mmax)
    for n in range(mmax):
        wd = profiles_out['wd'][n,:]
        x = profiles_out['x'][n,:]
        z = profiles_out['z'][n,:]
        u = profiles_out['u'][n,:]
        w = profiles_out['w'][n,:]
        rho = util_flow.get_rho(profiles_out['p'][n,:], profiles_out['T'][n,:])
        up = util_flow.get_up_2d(u,w, x,z)
        delta[n] = util_flow.get_delta(wd, up)
        delta_star[n] = util_flow.get_delta_star(wd, up, rho)
        theta[n] = util_flow.get_theta(wd, up, rho)
     
    int_out = {}
    int_out['x'] = profiles_out['x'][:,0]
    int_out['z'] = profiles_out['z'][:,0]
    int_out['delta'] = delta
    int_out['delta*'] = delta_star
    int_out['theta'] = theta
     

    ## out
    IO_util.write_hdf5(dir_out+filename_out_profiles, profiles_out)
    IO_util.write_hdf5(dir_out+filename_out_integrals, int_out)



if __name__ == '__main__':
    main()
