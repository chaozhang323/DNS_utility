import numpy as np
import h5py
from pyevtk.hl import gridToVTK, imageToVTK

def convert_inst_2d(filename_grid, filename_inst, nst, nen, nsk, filename_op):
    fid = h5py.File(filename_grid, 'r')
    x = fid.get('x').__array__().T[0:590,0:190,0]
    y = fid.get('y').__array__().T[0:590,0:190,0]
    ep = fid.get('ep').__array__().T[0:590,0:190,0]
    nx,ny = x.shape

    ## read data from HDF5
    ntot = (nen - nst +1) / nsk + 1
    varname = ['u','v','w','p']
    data = np.zeros([nx,ny,len(varname)])
    n = 0
    for step in range(nst,nen,nsk):
        filename = filename_inst[0] + step.__str__().zfill(8) + filename_inst[1]
        print(filename)
        fid = h5py.File(filename, 'r')
        gid = fid.get(list(fid.keys())[0])
        for idset in range(varname.__len__()):
            data[:,:,idset] = gid.get(varname[idset]).__array__().T[:,:,0]
        fid.close()
        n = n + 1

        filename = filename_op[0] + step.__str__().zfill(8) + filename_op[1]
        fid = h5py.File(filename, 'w')
        for idset in range(4):
            did = fid.create_dataset(varname[idset], [nx, ny], data=data[:, :, idset])
        fid.close()
    return

def check_inst_gkeys(filename):
    fid = h5py.File(filename, 'r')
    gkeys = list(fid.keys())
    print(gkeys)

def check_inst_dkeys(filename, gkey):
    fid = h5py.File(filename, 'r')
    gid = fid.get(gkey)
    dkeys = list(gid.keys())
    print(dkeys)

def check_inst_dset(filename, dkey):
    fid = h5py.File(filename, 'r')
    did = fid.get(dkey)
    return did.shape





