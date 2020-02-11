import tecplot as tp
import logging
import os
import numpy as np

def read_tp(flnm, zonesname=None, varsname=None, order='C', file_type=None, case_filenames=None):
    tp.new_layout()
    logging.basicConfig(level=logging.INFO)
    if file_type=='plt' or flnm[-4:]=='.plt' or flnm[-4:]=='.dat':
        dset = tp.data.load_tecplot(flnm)
    elif file_type=='szplt' or flnm[-6:]=='.szplt':
        dset = tp.data.load_tecplot_szl(flnm)
    elif file_type=='fluent':
        dset = tp.data.load_fluent(case_filenames=case_filenames, data_filenames=flnm)

    if zonesname==None:
        zonesname = [zone.name for zone in dset.zones()]
    if varsname==None:
        varsname = [var.name for var in dset.variables()]

    ## access data
    data_out = {}
    print('Reading from %s'%flnm)
    for zone_name in zonesname:
        zone = dset.zone(zone_name)
        print('Accessing zone %s of type %s ... \n num of elements = %d, num of points = %d:'%(zone.name, zone.zone_type, zone.num_elements, zone.num_points))
        data_out[zone.name] = {}
        try:
            nodemap = zone.nodemap
            data_out[zone.name]['nodemap'] = nodemap_to_array(nodemap)
        except:
            print('Nodemap not read!')
            pass

        for var_name in varsname:
            val = zone.values(var_name)
            print('Accessing variable %s ...'%var_name)
            data_out[zone_name][var_name] = val.as_numpy_array()
            try:
                data_out[zone_name][var_name] = data_out[zone_name][var_name].reshape(val.shape, order=order)
            except:
                print('%s is not reshaped, and therefore flat.'%var_name)
                pass

    print('Finished reading from %s'%flnm)
            
    
    return data_out



def read_tp_zone(dset, zone_name, vars_name, get_nodemap=False, keep_shape=False, order='C'):
    data = {}
    zone = dset.zone(zone_name)
    print('In zone %s, num of elements = %d, num of points = %d'%(zone_name,zone.num_elements, zone.num_points))
    for var_name in vars_name:
        if keep_shape:
            shape = zone.values(var_name).shape
            data[var_name] = zone.values(var_name).as_numpy_array().reshape(shape, order=order)
        else:
            data[var_name] = zone.values(var_name).as_numpy_array()
    ## nodemap
    if get_nodemap:
        nodemap = nodemap_to_array(zone.nodemap)
    else:
        nodemap = None
    return data, nodemap

def write_tp(filename, data, old=False, order='C'):
    if old:
        filename += '.plt'
    else:
        filename += '.szplt'
    print('Writing '+filename+'...')
    logging.basicConfig(level=logging.INFO)
    tp.new_layout()
    frame = tp.active_frame()

    varname = list(data.keys())
    dset = frame.create_dataset('data', var_names = varname)

    zone = dset.add_ordered_zone('zone', data[varname[0]].shape)
    for name in varname:
        zone.values(name)[:] = data[name].ravel(order=order)

    if old:
        tp.data.save_tecplot_plt(filename)
    else:
        tp.data.save_tecplot_szl(filename)
    print('Finished writing '+filename+'.')

def write_tp_grouped(filename, data, data_shared, old=False, order='C'):
    print('Writing '+filename+'...')
    logging.basicConfig(level=logging.INFO)
    tp.new_layout()
    frame = tp.active_frame()

    ## useful params
    sorted_keys = sorted(data.keys())
    varname_data = list(next(iter(data.values())).keys())
    varname_shared = list(data_shared.keys())
    varname_all = varname_data + varname_shared
    shape = data_shared[varname_shared[0]].shape

    ## create dset
    dset = frame.create_dataset('data', var_names = varname_all)

    ## create data by zone
    zone_list = []
    for zonename in sorted_keys:
        subdata = data[zonename]
        zone = dset.add_ordered_zone(zonename, shape)
        zone_list.append(zone)
        for name,var in subdata.items():
            try:
                zone.values(name)[:] = var.ravel(order=order)
            except:
                print('Failed to write %s ...'%name)
    
    ## share variables
    for name,var in data_shared.items():
        zone_list[0].values(name)[:] = var.ravel(order=order)
    variables_shared = [dset.variable(varname) for varname in varname_shared]
    dset.share_variables(zone_list[0], zone_list[1:], variables_shared)

    ## out
    if old:
        tp.data.save_tecplot_plt(filename+'.plt')
    else:
        tp.data.save_tecplot_szl(filename+'.szplt')
    print('Finished writing '+filename+'.')

def nodemap_to_array(nodemap):
    array = np.zeros(nodemap.shape, dtype=int)
    for i in range(nodemap.shape[0]):
        array[i,:] = nodemap[i]
    return array
