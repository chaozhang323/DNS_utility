from util import util_interp
import numpy as np
from . import sensor, dist

def reshape_by_index(shape, data_old):
    
    data_new = {}
    for name,var in data_old.items():
        try:
            shape_old = var.shape
            rank = var.ndim
            data_work = var.copy()
            for dim in range(rank):
                data_work = util_interp.interp_by_index(shape[dim], data_work, dim)
            data_new[name] = data_work
        except:
            print('%s not reshaped.'%name)

    return data_new

def expand_spanwise(shape_new, data_old, dim_new=0, span_name=None, span_len=2.*np.pi):
    '''
        In:
            shape_new (int)
            data_old (dict)
    '''

    data_new = {}
    for name,var in data_old.items():
        try:
            data_new[name] = np.repeat( np.expand_dims(var,dim_new), shape_new, axis=dim_new )
        except:
            print('%s not expanded.'%name)

    ## spanwise 
    if span_name is not None:
        shape_full_new = next(iter(data_new.values())).shape
        y = span_len/shape_new*np.arange(0., shape_new, 1.)
        for npts in shape_full_new[dim_new+1:]:
            y = np.repeat( np.expand_dims(y, -1), npts, axis=-1 )
        for npts in shape_full_new[0:dim_new][::-1]:
            y = np.repeat( np.expand_dims(y, 0), npts, axis=0)
        data_new[span_name] = y
    return data_new

def redistribute(dist_new, dist_old, data_old, dim=0):
    data_new = {}
    for name,var in data_old.items():
        var_work = np.swapaxes(var, -1, dim)
        for idx in np.ndindex(var_work.shape[:-1]):
            var_work[idx] = np.interp(dist_new, dist_old, var_work[idx])

        data_new[name] = np.swapaxes(var_work, -1, dim)

    return data_new

def redistribute_by_ducros(grid_in, data_in, dist_key, ducros, thres=0.15, dim=0, a=2.):
    ducros_work = np.swapaxes(ducros, dim, -1)
    grid_work = {key:np.swapaxes(value, dim, -1) for key,value in grid_in.items()}
    data_work = {key:np.swapaxes(value, dim, -1) for key,value in data_in.items()}
    
    # find max
    ducros_max = np.max(ducros_work)
    if ducros_max > thres:
        idx = np.unravel_index(np.argmax(ducros_work), ducros_work.shape)
        print("Find max-value index = ",idx)
        dist_old = grid_work[dist_key][idx[:-1]].copy()
        dist_new = dist.ptrans(dist_old).poly_3rd(a, idx[-1])
        grid_work = redistribute(dist_new, dist_old, grid_work, dim=-1)
        data_work = redistribute(dist_new, dist_old, data_work, dim=-1)
    
    grid_out = {key:np.swapaxes(value, dim, -1) for key,value in grid_work.items()}
    data_out = {key:np.swapaxes(value, dim, -1) for key,value in data_work.items()}
    return grid_out, data_out





