import numpy as np
from util import IO_util
from util import util_f

def compute_mesh_metrics(gdd):
    mm, Ji = util_f.mod_metrics.compute_mesh_metrics(gdd) 
    return mm, Ji

def compute_metric_identities(mm, Ji, stencilsize=6, bc=[0,0,0,0,0,0]):
    '''
        mm(:,:,j,i,k)
    '''
    mm_work = mm
    Ji_work = Ji
    shape = Ji_work.shape
    nx = shape[1] - 2*stencilsize
    ny = shape[2] - 2*stencilsize
    nz = shape[3] - 2*stencilsize
    I = util_f.mod_metrics.compute_metric_identities(mm_work, Ji_work, nx,ny,nz, stencilsize, bc)

    return I
   
def pack_mesh_metrics(mm):
    mm_out = {'didx':mm[0,0,:,:,:], \
              'didy':mm[0,1,:,:,:], \
              'didz':mm[0,2,:,:,:], \
              'djdx':mm[1,0,:,:,:], \
              'djdy':mm[1,1,:,:,:], \
              'djdz':mm[1,2,:,:,:], \
              'dkdx':mm[2,0,:,:,:], \
              'dkdy':mm[2,1,:,:,:], \
              'dkdz':mm[2,2,:,:,:] }
    return mm_out

def pack_metric_identities(mi):
    mi_out = {'Ixi':mi[0,:,:,:], \
              'Ieta':mi[1,:,:,:], \
              'Izeta':mi[2,:,:,:]}
    return mi_out

def examine_equal(data, axis):
    tmp = np.swapaxes(data, 0, axis)
    equals = [tmp[(0,Ellipsis)]==tmp[(n,Ellipsis)] for n in range(tmp.shape[0])]
    if np.all(equals):
        print(' All equal along axis=%d'%(axis))
    else:
        print(' Not all equal along axis=%d'%(axis))
        
    return


