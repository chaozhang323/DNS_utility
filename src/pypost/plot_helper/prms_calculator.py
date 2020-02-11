import sys
from os import path
sys.path.append( '/usr/local/home/yl8bc/Nic/PH/post/')#path.dirname( path.dirname( path.abspath(__file__) ) ) )
from util import IO_util, util_plot

import numpy as np
import matplotlib.pyplot as plt

def read_hdf5_group_series(filename_list, gname, vars_name):
    data, names =  IO_util.read_hdf5_group(filename_list[0], gname, vars_name)
    for filename in filename_list[1:]:
        data_tmp, names = IO_util.read_hdf5_group(filename, gname, vars_name)
        for name in names:
            data[name] += data_tmp[name]
    for name in names:
        data[name] *= 1./len(filename_list)

    return data, names

if __name__ == '__main__':
    ## Input:
    ##      dir_in: directory for AVERAGE input files
    ##      filename_stat:  a list for filename of the stat files
    ##      filename_grid:  filename of the grid file
    ##      start, end, skip:   control of stat files
    ##      index_k_vs_x: the index k
    ##      ib, ien:    for windowing
    dir_in = '/usr/local/home/yl8bc/duannas/Test_Incompact3d/Sample_AverageAcoustics/AVERAGE/'
    filename_stat_h5 = ['AveAcousticStat_', '.h5']
    filename_stat_plt = 'AveAcousticStat_timeave00032000-00131000_AppendGroup_Stat2d.plt'
    filename_grid = 'AveAcousticStat_GridMetrics.h5'
    start, end, skip = 127000, 131000, 1000
    index_k_vs_x = 10
    ibe, ien = 10, 50
    in_h5 = False
    in_plt = not in_h5
    ##  Out:
    ##      dir_out:    directory for output
    dir_out = './'
    export_lay = True


    ##  read from stats
    if in_h5:
        filename_list = []
        for n in range(start, end+skip, skip):
            filename = dir_in + filename_stat_h5[0] + n.__str__().zfill(8) + filename_stat_h5[1]
            filename_list.append(filename)
        stat_in, varsname = read_hdf5_group_series(filename_list, 'Stat2d', ['pave','p2'])
        stat_bot, varsname = read_hdf5_group_series(filename_list, 'Int_BotWall', ['tauw','delta'])
        stat_top, varsname = read_hdf5_group_series(filename_list, 'Int_TopWall', ['tauw','delta'])
        tauw = 0.5*(stat_bot['tauw'] + stat_top['tauw'])
        delta = 0.5*(stat_bot['delta'] + stat_top['delta'])
    if in_plt:
        from util import util_tp
        dset, zones_name, vars_name = util_tp.read_tp(dir_in+filename_stat_plt)
        stat_in = util_tp.read_tp_zone(dset, zones_name[0], ['pave','p2','tauw','delta'], order='F')
        tauw = 0.5*(stat_in['tauw'][:,0] + stat_in['tauw'][:,-1])
        delta = 0.5*(stat_in['delta'][:,0] + stat_in['delta'][:,-1])


    ##  calculation
    stat_out = {}
    stat_out['prms'] = np.sqrt( np.abs( stat_in['p2']-stat_in['pave']**2 ) )
    prms_tauw_vs_x = stat_out['prms'][:,index_k_vs_x] / tauw
    #
    prms_tauw_vs_z = np.mean(stat_out['prms'][ibe:ien,:]/tauw[ibe:ien,np.newaxis], axis=0)

    ##  read grid atg walls
    grid, gridname = IO_util.read_hdf5(dir_in+filename_grid, vars_name=['x','z'])
    x = grid['x'][:,0]
    z = np.mean(grid['z'][ibe:ien,:]/delta[ibe:ien,np.newaxis], axis=0)

    ##  plots
    fig, axe = plt.subplots(1,1)
    axe.plot(x, prms_tauw_vs_x)
    util_plot.labels(axe, 'x/m', r'$p_{rms}/\tau_{w}$')
    util_plot.title(axe, r'$p_{rms}/\tau_{w}$ along k = %d'%index_k_vs_x)
    fig.savefig(dir_out+'prms_tauw_vs_x.png')
    fig, axe = plt.subplots(1,1)
    axe.plot(z, prms_tauw_vs_z)
    util_plot.labels(axe, r'$z/\delta$', r'$p_{rms}/\tau_{w}$')
    util_plot.title(axe, r'$p_{rms}/\tau_{w}$ in window i = %d~%d'%(ibe,ien))
    fig.savefig(dir_out+'prms_tauw_vs_z.png')

    ## tecplot
    if export_lay:
        import logging
        import tecplot as tp
        logging.basicConfig(level=logging.INFO)
        tp.new_layout()
        frame = tp.active_frame()

        dset = frame.create_dataset('data', var_names = ['x','z','prms_tauw_vs_x', 'prms_tauw_vs_z'])

        zone_x = dset.add_ordered_zone('prms_tauw_vs_x', x.shape)
        zone_x.values('x')[:] = x.ravel(order='F')
        zone_x.values('prms_tauw_vs_x')[:] = prms_tauw_vs_x.ravel(order='F')
        zone_z = dset.add_ordered_zone('prms_tauw_vs_z', z.shape)
        zone_z.values('z')[:] = z.ravel(order='F')
        zone_z.values('prms_tauw_vs_z')[:] = prms_tauw_vs_z.ravel(order='F')

        from tecplot.constant import PlotType, FillMode
        plot = frame.plot(PlotType.XYLine)
        plot.activate()
        plot.add_linemap(zone=zone_x, x=dset.variable('x'), y=dset.variable('prms_tauw_vs_x'))
        plot.add_linemap(zone=zone_z, x=dset.variable('z'), y=dset.variable('prms_tauw_vs_z'))

        tp.data.save_tecplot_plt(dir_out+'prms_tauw')
        tp.save_layout(dir_out+'prms_tauw.lay')

    plt.show()