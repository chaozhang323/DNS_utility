from util import IO_util, util_tp, util_interp, util_data, util_flow
import numpy as np
from scipy.interpolate import griddata
import time


def main():
    ###params###
    mrange = range(1,99,1)
    mmax = len(mrange)
    nmax = 90 
    
    dir_in = '/usr/local/home/yl8bc/duannas/jhyt7/Acoustics/HIFiRE_Cone/FluentData/'#'/usr/local/home/yl8bc/duannas/jhyt7/Acoustics/HLB/Interp_Fluent_DNS/'
    filename_in_data = 'CircularCone2D.dat'# 'nozzle2d_fluent.dat'
    filename_in_case = 'CircularCone2D.cas'# 'nozzle2d_fluent.dat'
    ## names of the zones below
    interior = 'unspecified'
    wall_right = 'wall'    
    ## files output: 
    dir_out = '/usr/local/home/yl8bc/yl8bc/data/cone/FLUENT/'
    filename_out_data = 'data_fluent'
    filename_out_profiles= 'profiles_fluent'
    filename_out_integrals= 'integrals_fluent'
    ##--below no param is defined--##
    tstart = time.time()
    keys_change = {'X':'x', 'Y':'z', 'X Velocity':'u', 'Y Velocity':'w', 'Pressure':'p', 'Density':'rho', 'Temperature':'T'}
    grid_keys = ['x','z'] 


    ## read grid and data
    data = util_tp.read_tp(dir_in+filename_in_data, file_type='fluent', case_filenames=dir_in+filename_in_case)
    print('Data read. Time elapsed: %f secs'%(time.time()-tstart))

    ## recover structured
    print('Rebuilding...')
    data_int = util_data.change_dict(data[interior], keys_change)
    data_int = util_data.recover_structured_data(data_int, grid_keys=grid_keys)
    #data_int = recovery.recover_interior(data_int, grid_keys)

    R = util_flow.get_R(data_int['p'], data_int['rho'], data_int['T'])
    print(R)

    ## boundary
    data_right = util_data.change_dict(data[wall_right], keys_change)
    data_right = util_data.recover_structured_data(data_right, grid_keys=grid_keys, standard=False)
    data_right['rho'] = data_right['p'] / (R*data_right['T'])
    data = {name:np.concatenate((data_int[name], data_right[name][np.newaxis,:]), axis=0) for name in data_int.keys()}
    
    print('Rebuilt. Time elapsed: %f secs'%(time.time()-tstart))

    ## wall normal profiles

    sdata = util_data.structured_data(data, grid_keys) 
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
        rho = profiles_out['rho'][n,:]
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
    IO_util.write_hdf5(dir_out+filename_out_data, data)
    IO_util.write_hdf5(dir_out+filename_out_profiles, profiles_out)
    IO_util.write_hdf5(dir_out+filename_out_integrals, int_out)



if __name__ == '__main__':
    main()
