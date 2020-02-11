from util import IO_util, util_interp
import numpy as np
from scipy import interpolate

if __name__ == '__main__':
    ##---params --##
    shape = (2, 100, 140)
    dir_in_old = '/usr/local/home/yl8bc/duannas/yl8bc/data/cone/run6_release/REST/'#'./InitBody/'
    filename_grid_old = 'grid.h5'
    filename_data_old = 'flowdata_00200000.h5'
    dir_in_new = '/usr/local/home/yl8bc/duannas/yl8bc/data/cone/'
    filename_grid_new = 'cone_grid_flat_part2.h5'
    ## out
    dir_out = '/usr/local/home/yl8bc/duannas/yl8bc/data/cone/run7_2d/REST/'
    #filename_out =
    ##---params end----##
    ## read data/grid 
    grid_old = IO_util.read_hdf5(dir_in_old+filename_grid_old)
    data_old = IO_util.read_hdf5(dir_in_old+filename_data_old)
    grid_new = IO_util.read_hdf5(dir_in_new+filename_grid_new)
   
    grid_old = {'x':grid_old['x'], 'z':grid_old['z']}
    
    ## get IC
    data_new = util_interp.structured_interp(grid_old, data_old, grid_new, robust=True, method='linear')
    data_new.update(grid_new)
    
    ## distribution 
    x0 = data_new['x'][0,0]
    xm = data_new['x'][-1,0]
    dist_old = data_new['x'][:,0].copy()
    dist_new = np.linspace(x0, xm, data_new['x'].shape[0])

    ## interp to new dist
    for key,var in data_new.iteritems():
        for n in range(0, data_new['x'].shape[1]):
            data_new[key][:,n] = np.interp( dist_new, dist_old, var[:,n])

    ## intep to new shape
    data1 = {}
    for key,var in data_new.iteritems():
        data1[key] = np.zeros((data_new['x'].shape[0],shape[2]))
        for n in range(0, data_new['x'].shape[0]):
            data1[key][n,:] = util_interp.interp_by_index(var[n,:], shape[2])
    data2 = {}
    for key,var in data1.iteritems():
        data2[key] = np.zeros(shape[1:])
        for n in range(0, shape[2]):
            data2[key][:,n] = util_interp.interp_by_index(var[:,n], shape[1])

    ## fix axis boundary
    if True:
        x = data2['x']
        z = data2['z']
        L = z[:,-1]
        d = z[:,1]
        r = L / (L+0.5*d)
        r = r[:,np.newaxis]
        z = L[:,np.newaxis]*(1.-r)   + r*z
        #x[:,0] = x[:,1]
        data2['x'] = x
        data2['z'] = z
    
    data2['x'][:,0] = data2['x'][:,1]


    ## expand to 3d
    for key, var in data2.iteritems():
        data2[key] =  np.repeat(var[np.newaxis,:,:], shape[0], axis=0)

    y = 2.*np.pi/shape[0] * np.arange(0, shape[0], 1.) 
    #np.linspace(0., 2.*np.pi, shape[0]+1)[:-1]
    y = np.repeat(y[:,np.newaxis], shape[1], axis=1)
    y = np.repeat(y[:,:,np.newaxis], shape[2], axis=2)
    data2['y'] = y
    #data2['w'][:] = 0.
    data2['v'][:] = 0.
    data2['u'][:,-1,:] = 0.
    data2['w'][:,-1,:] = 0.
    data2['T'][:,-1,:] = 300. 
    
    
    ## output
    grid = {key:data2[key] for key in ['x','y','z']}
    data = {key:data2[key] for key in ['T','p','u', 'v','w']}
    data['time'] = np.array([0.])
    data_inlet = {key:data2[key][:,0:1,:] for key in ['T','p','u', 'v','w']}


    IO_util.write_hdf5(dir_out+'grid', grid)
    IO_util.write_hdf5(dir_out+'flowdata_00000000', data)
    IO_util.write_hdf5(dir_out+'inlet', data_inlet)
    IO_util.write_hdf5(dir_out+'vis', data2)


    
    
