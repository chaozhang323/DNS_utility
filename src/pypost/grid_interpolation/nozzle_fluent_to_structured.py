from util import IO_util, util_tp, util_interp, util_data, util_flow
import numpy as np
from scipy.interpolate import griddata
import time
from collections import OrderedDict


def main():
    '''
    This script converts unstructured cell-centered fluent data to structured data. Tecplot needs to be used to load fluent case and data files and write the data in Tecplot ASCII (*.dat) or binary format (*.plt,*szplt), which is used for variable "filename_in_data".
    run this script: python3 nozzle_fluent.py
    I've added one script that builds structured data from a Fluent/PLT file under /DNSMST_utility/src/pypost/grid_interpolation/nozzle_fluent.py. It has been tested to work with data located at duannas/duanl/Acoustics/TestFluentConversion which is a nozzle with wall on the top.

    If you're trying this on some other dataset, please pay attention to the following potential issues.
    Input data file_type. The data could be of PLT, SZPLT, DAT type, so please specify for the function read_tp()
    Data variable. Those variables could be named differently. By default the x-axis goes by the name 'X' and y-axis 'Y' and all other variables will be reconstructred to structured data according to these too coordinate variables.
    Wall direction. By default the wall is on the top.
    Broken connectivity. If the boundary points are not arranged in a ordered manner(eg. not following descending order with respect to 'X'), those points might need to be reordered. But this is not common anyway.
    The runtime for my test(around 200k datapoints) is ~20mins and the code is not optimized to its best state yet.
    ATTENTION!!! ------------------------------
    POTENTIAL FAILURE:
    1. Unmatached variable names. please go to line 37  keys_change and re-define!
    2. Wrong direction. Please go to line to chooose one of those boundaries!
    3. Input file format. Caution with tecplot files, for the format is unclear(ordered/unordered).    
    -----------------------------------------------'''
    ###params###
    dir_in = './'
    filename_in_data = 'nozzle_cartesian_231x101_kwsst_refine1_2nd_UnstrTecplotCellCentered.plt'
    ## names of the zones below
    interior = 'unspecified'
    wall_top = 'wall'    
    ## files output: 
    dir_out = './'
    filename_out_data = 'nozzle_cartesian_460x200_kwsst_2nd_structured'
    ##--below no param is defined--##
    tstart = time.time()
    keys_change = {'X':'x', 'Y':'z', 'X Velocity':'u', 'Y Velocity':'w', 'Pressure':'p', 'Density':'rho', 'Temperature':'T'}
    grid_keys = ['x','z'] 


    ## read grid and data
    data = util_tp.read_tp(dir_in+filename_in_data)
    print('Data read. Time elapsed: %f secs'%(time.time()-tstart))

    ## recover structured
    print('Rebuilding...')
    data_int = util_data.change_dict(data[interior], keys_change)
    data_int = util_data.recover_structured_data(data_int, grid_keys=grid_keys)
    

    R = util_flow.get_R(data_int['p'], data_int['rho'], data_int['T'])
    print(R)

    ## boundary
    data_top = util_data.change_dict(data[wall_top], keys_change)
    data_top = util_data.recover_structured_data(data_top, grid_keys=grid_keys, standard=False)
    data_top['rho'] = data_top['p'] / (R*data_top['T'])
    data_top = sort_top(data_top, grid_keys[0])
    
    ## one of these boundries below, need to choose!!
    ## top 
    data = {name:np.concatenate((data_int[name], data_top[name][:,np.newaxis]), axis=1) for name in data_int.keys()}
    ## bot
    #data = {name:np.concatenate((data_top[name][:,np.newaxis], data_int[name]), axis=1) for name in data_int.keys()}
    ## right
    #data = {name:np.concatenate((data_int[name], data_top[name][np.newaxis,:]), axis=0) for name in data_int.keys()}
    ## left
    #data = {name:np.concatenate((data_top[name][np.newaxis,:], data_int[name]), axis=0) for name in data_int.keys()}
    
    print('Rebuilt. Time elapsed: %f secs'%(time.time()-tstart))

    ## out
    # Put gridkeys in first two
    def sortkey(x):
        if x[0]==grid_keys[0]:
            return 0
        elif x[0]==grid_keys[1]:
            return 1
        else:
            return 2

    data =OrderedDict( sorted(data.items(), key=sortkey))
    print(    data.keys())
    IO_util.write_hdf5(dir_out+filename_out_data, data)
    util_tp.write_tp(dir_out+filename_out_data, data, old=True, order='F')



def sort_top(data, gridkey):
    I = np.argsort(data[gridkey])
    data_out = {}
    for key,val in data.items():
        data_out[key] = val[I]
    return data_out


if __name__ == '__main__':
    main()


