import h5py
import numpy as np
#from pyevtk.hl import gridToVTK, imageToVTK
#import tecplot as tp
#import logging
#from  util.derivative_lele import *


def read_grid(filename_grid, dim=2, slc=None):
    """Read Cartesian grid HDF5 file.
    In:     grid .h5 filename
            dimension of grid, optional
            slc(slice) of specific part, optional
    Out:    grid dictionary
            nx, ny, nz"""
    ## get shape and slice
    fid = h5py.File(filename_grid, 'r')
    if dim==2:
        varnames = ['x', 'y', 'ep']
        if slc is None: slc = np.s_[0,:,:]
    if dim==3:
        varnames = ['x', 'y', 'z', 'ep']
        if slc is None: slc = np.s_[:,:,:]

    dset = fid.get(varnames[0])
    shape = dset[slc].shape
    (nx,ny,nz) = dset.shape
    ## read variables
    grid = {}
    for varname in varnames:
        try:
            dset = fid.get(varname)
            grid[varname] = np.zeros(shape)
            dset.read_direct(grid[varname], source_sel=slc)
            grid[varname] = grid[varname].transpose()
        except:
            pass
    fid.close()
    return grid, nx, ny, nz

def read_inst(filename_inst, gkey, nst, nen, nsk, slc):
    #====get data shape====#
    varname = ['u','v','w','p']
    filename = filename_inst[0] + nst.__str__().zfill(8) + filename_inst[1]
    fid = h5py.File(filename, 'r')
    gid = fid.get(gkey)
    dset = gid.get(varname[0])
    if slc is None:
        shape = dset.shape
    else:
        shape = dset[slc].shape
    fid.close()
    nt = shape[0]
    ntot = ((nen - nst) / nsk + 1) * nt
    ## read variables
    inst = {}
    for name in varname:
        inst[name] = np.zeros((ntot,)+shape[1:])
    #====read inst====#
    n = 0
    for step in range(nst,nen+1,nsk):
        filename = filename_inst[0] + step.__str__().zfill(8) + filename_inst[1]
        print(filename)
        fid = h5py.File(filename, 'r')
        gid = fid.get(gkey)
        for name in varname:
            did = gid.get(name)
            did.read_direct(inst[name], source_sel=slc, dest_sel=np.s_[n*nt:(n+1)*nt,:,:,:])
        fid.close()
        n += 1

    ## process inst
    for name in varname:
        inst[name] = inst[name].transpose()

    return inst

def read_outdated_inst(filename, gkey):
    #====get data shape====#
    varname = ['u','v','w','p']
    fid = h5py.File(filename, 'r')
    gid = fid.get(gkey)
    dset = gid.get(varname[0])
    shape = dset.shape
    fid.close()
    ## read variables
    inst = {}
    for name in varname:
        inst[name] = np.zeros(shape)
    #====read inst====#
    fid = h5py.File(filename, 'r')
    gid = fid.get(gkey)
    for name in varname:
        did = gid.get(name)
        did.read_direct(inst[name])
    fid.close()
    ## process inst
    for name in varname:
        inst[name] = inst[name].transpose()

    return inst

def read_flowdata(filename_flowdata):
    varname = ['u','v','w']
    #====get data shape====#
    fid = h5py.File(filename_flowdata, 'r')
    flowdata = {}

    #====read flowdata====#
    fid = h5py.File(filename_flowdata, 'r')
    for name in varname:
        did = fid.get(name)
        data = did[:,:,:].transpose()
        flowdata[name] = data
    fid.close()

    return flowdata

def write_inst_vtk_structured(filename_vtk,x,y,z,inst):
    for step in range(inst.shape[3]):
        filename = filename_vtk + step.__str__().zfill(8)
        print(filename)
        u = np.array(inst[:,:,:,step,0])
        v = np.array(inst[:,:,:,step,1])
        w = np.array(inst[:,:,:,step,2])
        p = np.array(inst[:,:,:,step,3])
        gridToVTK(filename, x,y,z, pointData={"u": u, "v": v, "w": w, "p": p})
    return

def read_stat_2d_basic(filename_stat, nst, nen, nsk):
    ntot = (nen-nst+1) / nsk
    statname = ['umean', 'vmean', 'wmean', 'pmean', 'uumean', 'vvmean', 'wwmean', 'uvmean', 'uwmean', 'vwmean']
    #==get stat size==#
    filename = filename_stat[0] + nst.__str__().zfill(8) + filename_stat[1] + (nst + nsk - 1).__str__().zfill(8) + filename_stat[2]
    fid = h5py.File(filename, 'r')
    gid = fid.get('stat2d')
    did = gid.get(statname[0])
    ny, nx = did.shape

    stat = np.zeros([nx,ny,len(statname)])
    for step in range(nst,nen,nsk):
        filename = filename_stat[0] + step.__str__().zfill(8) + filename_stat[1] + (step+nsk-1).__str__().zfill(8)+ filename_stat[2]
        print(filename)
        fid = h5py.File(filename, 'r')
        gid = fid.get('stat2d')
        for idset in range(statname.__len__()):
            stat[:,:,idset] += gid.get(statname[idset]).__array__().T
        fid.close()
    ## averaging temporally
    for idset in range(statname.__len__()):
        stat[:,:,idset] /= ntot

    ## get Reynolds stress
    rstressname = ['tau11','tau22','tau33','tau12','tau13','tau23']
    rstress = np.zeros([nx,ny,6])
    rstress[:,:,0] = stat[:,:,4] - stat[:,:,0]*stat[:,:,0]
    rstress[:,:,1] = stat[:,:,5] - stat[:,:,1]*stat[:,:,1]
    rstress[:,:,2] = stat[:,:,6] - stat[:,:,2]*stat[:,:,2]
    rstress[:,:,3] = stat[:,:,7] - stat[:,:,0]*stat[:,:,1]
    rstress[:,:,4] = stat[:,:,8] - stat[:,:,0]*stat[:,:,2]
    rstress[:,:,5] = stat[:,:,9] - stat[:,:,1]*stat[:,:,2]

    stat[:,:,4:10] = rstress
    return stat, statname[0:4]+rstressname

def read_stat_2d(filename_stat, nst, nen, nsk):
    ntot = (nen-nst+1) / nsk
    #==get stat size==#
    filename = filename_stat[0] + nst.__str__().zfill(8) + filename_stat[1] + (nst + nsk - 1).__str__().zfill(8) + filename_stat[2]
    fid = h5py.File(filename, 'r')
    gid = fid.get('stat2d')
    statname = list(gid.keys())
    did = gid.get(statname[0])
    ny, nx = did.shape

    stat = {}
    for name in statname:
        stat[name] = np.zeros([ny,nx])
    tmp = np.zeros([ny,nx])
    for step in range(nst,nen,nsk):
        filename = filename_stat[0] + step.__str__().zfill(8) + filename_stat[1] + (step+nsk-1).__str__().zfill(8)+ filename_stat[2]
        print(filename)
        fid = h5py.File(filename, 'r')
        gid = fid.get('stat2d')
        for name in statname:
            dset = gid.get(name)
            dset.read_direct(tmp)
            stat[name] += tmp
        fid.close()
    ## averaging temporally and transposing
    for name in statname:
        stat[name] /= ntot
        stat[name] = stat[name].transpose()

    return stat

#def read_hdf5(filename, vars_name=None, transposed=False):
#    """Read direct datasets under a H5 file.
#     In:    H5 filename.
#            vars_name(variable name list).
#     Out:   data dictionary.
#            data name in H5 file."""
#    fid = h5py.File(filename, 'r')
#    if vars_name==None: vars_name = fid.keys()
#
#    data = {}
#    for var_name in vars_name:
#        try:
#            dset = fid.get(var_name)
#            shape = dset.shape
#            data[var_name] = np.zeros(shape)
#            dset.read_direct(data[var_name])
#            if transposed: data[var_name] = data[var_name].transpose()
#        except:
#            print var_name, 'not read.'
#            pass
#    fid.close()
#    print 'Read from ', filename
#    print 'Variables names = '
#    print '\n'.join(vars_name)
#
#    return data, vars_name
def read_hdf5(filename, namelist=None, **kwargs):
    """Read direct datasets under a H5 file.
     In:    H5 filename.
            vars_name(variable name list).
     Out:   data dictionary.
            data name in H5 file."""

    print('Reading %s...'%filename)

    fid = h5py.File(filename, mode='r')
    
    data = read_hdf5_tree(fid, namelist, **kwargs)

    fid.close()
    
    print('Finished reading %s.'%filename)
    return data

def read_hdf5_tree(gid, namelist, transposed=False):
    def tinker(gid, nmlst):
        if nmlst is not None:
            nmlst_out = nmlst
        else:
            keys = list(gid.keys())
            nmlst_out = {key:None for key in keys}
        return nmlst_out

    ## prepare keys to access
    nmlst_work = tinker(gid, namelist)

    ## read data
    data = {}
    for key in nmlst_work:
        dset = gid.get(key)
        if type(dset)==h5py.Group:
            print('Accessing group %s...'%key)
            subdata = read_hdf5_tree(dset, nmlst_work[key], transposed=transposed)
            data[key] = subdata
        else:
            print('Accessing dataset %s...'%key)
            shape = dset.shape
            data[key] = np.zeros(shape)
            dset.read_direct(data[key])
    return data


def read_hdf5_group(filename, gname, vars_name=None):
    """Read direct datasets under a H5 group.
     In:    H5 filename
            group name.
     Out:   data dictionary
            data name in H5 file"""
    fid = h5py.File(filename, 'r')
    gid = fid.get(gname)
    if vars_name is None:   vars_name = list(gid.keys())

    data = {}
    for var_name in vars_name:
        try:
            dset = gid.get(var_name)
            shape = dset.shape
            data[var_name] = np.zeros(shape)
            dset.read_direct(data[var_name])
        except:
            pass
    fid.close()
    print('Read from ', ''.join((filename,'/',gname)))
    print('Variables names = ')
    print('\n'.join(vars_name))

    return data, vars_name

def write_ascii_point(filename, data, title='Untitled'):
    """Write ASCII data file by point
    """
    vars_name = list(data.keys())
    num_var = len(vars_name)
    num_point = np.size( data[vars_name[0]] )
    shape = data[vars_name[0]].shape
    data_raveled = np.zeros((num_point,num_var))
    for var_name,n in zip(vars_name,list(range(num_var))):
        data_raveled[:,n] = data[var_name].ravel(order='F')

    f = open(filename+'.dat', mode='w')
    f.write('TITLE = "%s" \n'%title)

    f.write('Variables = ')
    for var_name in vars_name[:-1]:
        f.write('"%s",'%var_name)
    f.write('"%s"\n'%vars_name[-1])

    f.write('Zone I=%d'%shape[0])
    if len(shape)>=2:   f.write(',J=%d'%shape[1])
    f.write('\n')

    for i in range(num_point):
        for n in range(num_var):
            f.write('%.6E\t'%data_raveled[i,n])
        f.write('\n')
    f.close()
    return

def peep_hdf5(filename):
    fid = h5py.File(filename, 'r')
    vars_name = list(fid.keys())
    print(vars_name)
    return vars_name


def write_stat_hdf5_2d(filename,x,y,ep, stat):
    fid = h5py.File(filename+'.h5', 'w')
    nx, ny = x.shape

    fid.create_dataset('x', [nx,ny], data=x)
    fid.create_dataset('y', [nx,ny], data=y)
    fid.create_dataset('ep', [nx,ny], data=ep)
    for name, var in stat.items():
        fid.create_dataset(name, [nx,ny], data=var)

    fid.close()
    return

def write_hdf5(filename, data):
    """Write direct datasets under a H5 file.
    In:    H5 filename
            data dictionary
    Out:    None"""
    
    if '.h5' in filename:
        fid = h5py.File(filename, 'w')
    else:
        filename = filename+'.h5'
        fid = h5py.File(filename, 'w')

    print('Writing %s...'%filename)

    write_hdf5_group(fid, data)

    fid.close()
    print('Finished writting %s.'%filename)
    return

#def write_hdf5_group(filename, data):
#    if '.h5' in filename:
#        fid = h5py.File(filename, 'w')
#    else:
#        fid = h5py.File(filename+'.h5', 'w')
#
#
#    for gname, gvar in data.iteritems():
#        gid = fid.create_group(gname)
#        for dname, dvar in gvar.iteritems():
#            shape = dvar.shape
#            gid.create_dataset(dname, shape, data=dvar)
#
#    fid.close()
#    return

def write_hdf5_group(fid, data):
    for name,var in data.items():
        if type(var) is dict:
            print('Accessing group %s...'%name)
            gid = fid.create_group(name)
            write_hdf5_group(gid, var)
        else:
            print('Writing dataset %s...'%name)
            shape = var.shape
            fid.create_dataset(name, shape, data=var)
    return

def write_stat_hdf5_1d(filename, y, stat, statname):
    fid = h5py.File(filename, 'w')
    ny, = y.shape

    did = fid.create_dataset('y', [ny], data=y)
    for idset in range(len(statname)):
        did = fid.create_dataset(statname[idset], [ny], data=stat[:,idset])
    fid.close()
    return

def write_inst_hdf5(filename, inst):
    fid = h5py.File(filename+'.h5', 'w')
    shape = inst['u'].shape

    for name, var in inst.items():
        fid.create_dataset(name, shape, data=var)

    fid.close()
    return

def read_stat_proc_2d(filename_op):
    fid = h5py.File(filename_op, 'r')
    ## read data from HDF5
    x = fid.get('x').__array__()
    y = fid.get('y').__array__()
    ep = fid.get('ep').__array__()
    nx,ny = x.shape



    statname = list(fid.keys())
    data = {}
    for name in statname:
        data[name] = fid.get(name).__array__()

    fid.close()
    return x,y,ep,data,nx,ny

def read_stat_proc_1d(filename_op):
    fid = h5py.File(filename_op, 'r')
    ## read stat from HDF5
    y = fid.get('y').__array__()
    ny, = y.shape

    statname = ['umean','vmean','wmean','pmean','tau11','tau22','tau33','tau12','tau13','tau23']
    stat = np.zeros([ny,len(statname)])
    for idset in range(statname.__len__()):
        stat[:,idset] += fid.get(statname[idset]).__array__()

    fid.close()
    return y,stat,statname,ny

def write_tp_2d_time(filename, x,y, data, ep=None, dt=1e0, old=False):
    if old:
        filename += '.plt'
    else:
        filename += '.szplt'
    print('Writing '+filename+'...')
    logging.basicConfig(level=logging.INFO)
    tp.new_layout()
    frame = tp.active_frame()
    varname = list(data.keys())
    if ep is None:
        dset = frame.create_dataset('data', var_names = ['x','y']+varname)
    else:
        dset = frame.create_dataset('data', var_names = ['x','y','ep']+varname)

    zone = dset.add_ordered_zone('zone', x.shape)
    zone.values('x')[:] = x.ravel(order='F')
    zone.values('y')[:] = y.ravel(order='F')
    if not ep is None:
        zone.values('ep')[:] = ep.ravel(order='F')
    for name in varname:
        zone.values(name)[:] = data[name].ravel(order='F')

    N_time = data[varname[0]].shape[-1]
    n = 0
    zone = dset.add_ordered_zone(n.__str__().zfill(8), x.shape[:-1])
    zone.solution_time = n*dt
    zone.values('x')[:] = x.ravel(order='F')
    zone.values('y')[:] = y.ravel(order='F')
    for var in data.values():
        zone.values('u')[:] = var[:,:,0,n].ravel(order='F')
        zone.values('v')[:] = var[:,:,0,n].ravel(order='F')
    zone0 = zone
    zonen = []
    for n in range(1, N_time):
        zone = dset.add_ordered_zone(n.__str__().zfill(8), x.shape)
        zone.solution_time = n*dt
        for var in data.values():
            zone.values('u')[:] = var[:,:,0,n].ravel(order='F')
            zone.values('v')[:] = var[:,:,0,n].ravel(order='F')
        zonen.append( zone )
    dset.share_variables(zone0, zonen, [dset.variable('x'), dset.variable('y')])

    if old:
        tp.data.save_tecplot_plt(filename)
    else:
        tp.data.save_tecplot_szl(filename)
    print('Finished writing '+filename+'.')

def write_tp_2d(filename, x,y, data, ep=None, old=False):
    if old:
        filename += '.plt'
    else:
        filename += '.szplt'
    print('Writing '+filename+'...')
    logging.basicConfig(level=logging.INFO)
    tp.new_layout()
    frame = tp.active_frame()
    varname = list(data.keys())
    if ep is None:
        dset = frame.create_dataset('data', var_names = ['x','y']+varname)
    else:
        dset = frame.create_dataset('data', var_names = ['x','y','ep']+varname)

    zone = dset.add_ordered_zone('zone', x.shape)
    zone.values('x')[:] = x.ravel(order='F')
    zone.values('y')[:] = y.ravel(order='F')
    if not ep is None:
        zone.values('ep')[:] = ep.ravel(order='F')
    for name in varname:
        zone.values(name)[:] = data[name].ravel(order='F')

    if old:
        tp.data.save_tecplot_plt(filename)
    else:
        tp.data.save_tecplot_szl(filename)
    print('Finished writing '+filename+'.')

def write_tp_3d(filename, x,y,z, data, varname=['u','v','w'], ep=None, old=False):
    if old:
        filename += '.plt'
    else:
        filename += '.szplt'
    print('Writing '+filename+'...')
    logging.basicConfig(level=logging.INFO)
    tp.new_layout()
    frame = tp.active_frame()
    if ep is None:
        dset = frame.create_dataset('flowdata', var_names = ['x','y','z']+varname)
    else:
        dset = frame.create_dataset('flowdata', var_names = ['x','y','z','ep']+varname)
    n = 0
    zone = dset.add_ordered_zone(n.__str__().zfill(8), x.shape)
    zone.values('x')[:] = x.ravel(order='F')
    zone.values('y')[:] = y.ravel(order='F')
    zone.values('z')[:] = z.ravel(order='F')
    if not ep is None:
        zone.values('ep')[:] = ep.ravel(order='F')
    for name in varname:
        zone.values(name)[:] = data[name].ravel(order='F')

    if old:
        tp.data.save_tecplot_plt(filename)
    else:
        tp.data.save_tecplot_szl(filename)
    print('Finished writing '+filename+'.')
