import numpy as np
from util import util_f, IO_util
from . import metrics

def ducros(u,v,w, grid_in, stencilsize=6, u_inf=1.,delta_in=1.):
    '''
        x,y,z(j,i,k)
    '''
    u_work = u
    v_work = v
    w_work = w
    nx,ny,nz = u_work.shape
    
    gd = metrics.grid_derivative.analytically_2d_(grid_in, bc_gd=[1,1,1,1,1,1])
    
    mm, Ji = util_f.mod_metrics.compute_mesh_metrics(gd)
    mm = mm[:,:, stencilsize:-stencilsize, stencilsize:-stencilsize, stencilsize:-stencilsize]
    
    div = util_f.flow.div_3d(mm, u_work,v_work,w_work)
    curl = util_f.flow.curl_3d(mm, u_work,v_work,w_work)
    
    div_sq = div**2
    curl_sq = np.sum(curl**2, axis=0)
    
    theta = div_sq / (div_sq + curl_sq + (u_inf/delta_in)**2)
    return theta


def extract_large_ducros(ducros, thres=0.15, dim=0):
    for idx in np.ndindex(ducros.shape[:-1]):
        ducros_line = ducros[idx]
        ducros_max = np.max(ducros_line)
        if ducros_max > thres:
            yield idx, np.argmax(ducros_line)
            

    

