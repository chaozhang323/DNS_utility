import numpy as np
from scipy import interpolate

def structured_interp(grid_in, data_in, grid_out, robust=False, **kwargs):
    ## pack grid_in/out
    grid_in_pack = []
    grid_out_pack = []
    for name in grid_out.keys():
        grid_in_pack.append(grid_in[name].ravel())
        grid_out_pack.append(grid_out[name].ravel())
    
    grid_in_pack = np.stack(grid_in_pack, axis=1)
    grid_out_pack = np.stack(grid_out_pack, axis=1)
    
    ## remember shape
    shape = grid_out[next(iter(grid_out))].shape
    num_var = len(data_in)

    ## interpolation
    n = 0
    data_out = {}
    for name, var in data_in.items():
        try:
            data_out[name] = interpolate.griddata(grid_in_pack, var.ravel(), grid_out_pack, **kwargs)
            n+=1
            print('Interpolated %d out of %d'%(n, num_var))
        except:
            print('%s not interpolated!'%(name))
            pass

    ## robust control
    key1 = next(iter(data_out))
    I_nan = np.isnan(data_out[key1])
    num_nan = np.sum(I_nan)
    if num_nan>0 and robust==True:
        print('There are points found out of convex! Nearest interpolation enabled!')
        for name in list(data_out.keys()):
            try:
                data_out[name][I_nan] = interpolate.griddata(grid_in_pack, data_in[name].ravel(), grid_out_pack[I_nan], method='nearest')
            except:
                print('Robust control failed for %s.'%(name))
                pass

    ## reshape
    for name in list(data_out.keys()):
        data_out[name] = data_out[name].reshape(shape)

    return data_out


def interp_by_index(shape_new, data_old, dim, kind='linear'):
    shape_old = data_old.shape[dim]
    
    if shape_old==shape_new:
        return data_old.copy()

    index_old = np.linspace(0., shape_old-1., shape_old)
    index_new = np.linspace(0., shape_old-1., shape_new)

    #data_new = np.interp( index_new, index_old, u, axis=dim)
    interp_linear = interpolate.interp1d(index_old, data_old, kind=kind, axis=dim)
    data_new = interp_linear(index_new)
    
    return data_new
    
