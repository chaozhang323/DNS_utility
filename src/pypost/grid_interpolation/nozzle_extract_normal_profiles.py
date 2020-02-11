from util import IO_util, util_tp, util_interp, util_data, util_flow
import numpy as np
import time


def main():
    '''
    Note:
        The boundary thickness calculation relies on a clear profile  that converges. current code cannot intelligently tell which wall normal profile is "perfect". Thus boundary layer thickness might yield to poor quality if the profile is badly extracted.
    '''
    ###params## 
    mrange = range(100,300,100) ## along the wall at which index to extract
    nmax = 75 ## along wall normal direction how many nodes to extract
    fluid_type = 'nitrogen'

    dir_in = './'
    filename_in_data = 'structured_tec.plt'
    is_in_tecplot =  True # Input is of tecplot format ??? False to be HDF5, True to be Tecplot format

    ## files output---------------------------------! 
    dir_out = './'
    filename_out_profiles = 'profiles'
    filename_out_integrals= 'integrals'
    ## names of variables---------------------------!
    name_x = 'X'           # name of these variables!
    name_z = 'Y'
    name_u = 'u'
    name_w = 'w'
    name_p = 'p'
    name_T = 'T'
    ## --------------------------------------
    keys_change = {name_x:'x', name_z:'z'}
    grid_keys = ['x','z'] 
    ## end params

    tstart = time.time()
    ## read grid and data
    if is_in_tecplot:
        data = util_tp.read_tp(dir_in+filename_in_data, order='F')
        data = next(iter(data.values()))
    else:
        data = IO_util.read_hdf5(dir_in+filename_in_data)
    print('Data read. Time elapsed: %f secs'%(time.time()-tstart))

    ## slice 
    data = util_data.change_dict(data, keys_change=keys_change)
    sdata = util_data.structured_data(data, grid_keys)

    ## wall normal profiles
    profiles_out = sdata.get_wall_normal_profiles('top', mrange, nmax)

    ## integral
    mmax = len(mrange)
    delta = np.zeros(mmax)
    delta_star = np.zeros(mmax)
    theta = np.zeros(mmax)
    tauw = np.zeros(mmax)
    utau = np.zeros(mmax)
    ztau = np.zeros(mmax)
    for n in range(mmax):
        wd = profiles_out['wd'][n,:]
        x = profiles_out['x'][n,:]
        z = profiles_out['z'][n,:]
        u = profiles_out[name_u][n,:]
        w = profiles_out[name_w][n,:]
        p = profiles_out[name_p][n,:]
        T = profiles_out[name_T][n,:]
        
        rho = util_flow.get_rho(profiles_out[name_p][n,:], profiles_out[name_T][n,:], fluid_type=fluid_type)
        up = util_flow.get_up_2d(u,w, x,z)
        delta[n] = util_flow.get_delta(wd, up)
        delta_star[n] = util_flow.get_delta_star(wd, up, rho)
        theta[n] = util_flow.get_theta(wd, up, rho)
        mu = util_flow.get_mu_Sutherland(T[0],fluid_type=fluid_type)
        tauw[n] = np.abs( mu*up[1]/wd[1] )
        utau[n] = np.sqrt(tauw[n]/rho[0])
        ztau[n] = mu / rho[0] / utau[n]

    int_out = {}
    int_out['x'] = profiles_out['x'][:,0]
    int_out['z'] = profiles_out['z'][:,0]
    int_out['delta'] = delta
    int_out['delta*'] = delta_star
    int_out['theta'] = theta
    int_out['tauw'] = tauw
    int_out['utau'] = utau
    int_out['ztau'] = ztau
     

    ## out
    IO_util.write_hdf5(dir_out+filename_out_profiles, profiles_out)
    IO_util.write_hdf5(dir_out+filename_out_integrals, int_out)
    IO_util.write_ascii_point(dir_out+filename_out_profiles, profiles_out)
    IO_util.write_ascii_point(dir_out+filename_out_integrals, int_out)



if __name__ == '__main__':
    main()
