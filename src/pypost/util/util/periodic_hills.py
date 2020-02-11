import numpy as np
from .stat_util import *
from .IO_util import *
from stat_io import *
import matplotlib.pyplot as plt
import os
from .wall_util import *
import scipy.signal as signal
from util.grid import grid

class periodic_hills:
    """A class of periodic hills for post-processing"""
    dirs_dict = {'AVERAGE':'AVERAGE/', 'TIMESERIES':'TIMESERIES/', 'REST':'REST/'}
    flnms_dict = {'grid':'grid.h5', 'stat':['Avespantime_', '-', '.h5'], 'timeseries':['timeseries_', '.h5'], 'flowdata':['flowdata_', '.h5']}

    def __init__(self, name, path, Re=5600):
        self.name = name
        self.path = path
        self.Re = Re
        self.grid = grid(self.path+self.dirs_dict['REST']+self.flnms_dict['grid'])

    def compare_with_Breuer(self, ax, x_line):
        x_list = [0.05, 0.5, 1., 2., 3., 4., 5., 6., 7., 8.]
        n_line = x_list.index(x_line) + 1
        fid = open('./Breuer/UFR3-30_X_5600_data_CR-'+n_line.__str__().zfill(3)+'.dat', 'r')
        str = fid.readlines()
        fid.close()
        del str[:6]
        ipdata1 = np.loadtxt(str, delimiter=',')
        ax[0, 0].plot(ipdata1[:, 1], ipdata1[:, 0], '.')
        ax[0, 1].plot(ipdata1[:, 2], ipdata1[:, 0], '.')
        ax[1, 0].plot(ipdata1[:, 3], ipdata1[:, 0], '.')
        ax[1, 1].plot(ipdata1[:, 4], ipdata1[:, 0], '.')

        fid = open('./Breuer/UFR3-30_data-NP-Re5600-DNS2-'+n_line.__str__().zfill(2)+'.dat', 'r')
        str = fid.readlines()
        fid.close()
        ipdata2 = np.loadtxt(str)
        ax[0, 0].plot(ipdata2[:, 1], ipdata2[:, 0], '-.', linewidth=3)
        ax[0, 1].plot(ipdata2[:, 2], ipdata2[:, 0], '-.', linewidth=3)
        ax[1, 0].plot(ipdata2[:, 3], ipdata2[:, 0], '-.', linewidth=3)
        ax[1, 1].plot(ipdata2[:, 4], ipdata2[:, 0], '-.', linewidth=3)

        fid = open('./Breuer/UFR3-30_C_5600_data_MB-'+n_line.__str__().zfill(3)+'.dat', 'r')
        str = fid.readlines()
        fid.close()
        del str[:37]
        ipdata3 = np.loadtxt(str)
        ax[0, 0].plot(ipdata3[:, 1], ipdata3[:, 0], '--', linewidth=3)
        ax[0, 1].plot(ipdata3[:, 2], ipdata3[:, 0], '--', linewidth=3)
        ax[1, 0].plot(ipdata3[:, 3], ipdata3[:, 0], '--', linewidth=3)
        ax[1, 1].plot(ipdata3[:, 4], ipdata3[:, 0], '--', linewidth=3)
        ax[1, 1].legend(['present', 'experiments', 'MGLET', 'LESOCC'], fontsize=25)

    def convert_stats(self, nst,nen,nsk, flnm_save):
        x, y, z, ep, nx, ny, nz = self.grid.get_grid(dim=2)
        dx, dy, dz = self.grid.get_dx()

        flnm_stat = [self.path + self.dirs_dict['AVERAGE'] + self.flnms_dict['stat'][0]] + self.flnms_dict['stat'][1:]
        stat = read_stat_2d(flnm_stat, nst,nen,nsk)
        stat_out = compute_stats_2d(stat, dx, dy, self.Re)
        write_stat_hdf5_2d(flnm_save,x,y,ep,stat_out)
        self.stat_save = {'filename':flnm_save, 'interval':(nst,nen,nsk)}

        return x, y, ep, stat_out

    def convert_stats_basic(self, nst,nen,nsk, flnm_save):
        flnm_grid = self.path + self.dirs_dict['REST'] + self.flnms_dict['grid']
        x, y, z, ep, nx, ny, nz = read_grid(flnm_grid)
        flnm_stat = [self.path + self.dirs_dict['AVERAGE'] + self.flnms_dict['stat'][0]] + self.flnms_dict['stat'][1:]
        stat, statname = read_stat_2d_basic(flnm_stat, nst,nen,nsk)
        x = x[:,:,0]
        y = y[:,:,0]
        ep = ep[:,:,0]
        write_stat_hdf5_2d(flnm_save,x,y,ep,stat,statname)
        self.stat_save = {'filename':flnm_save, 'interval':(nst,nen,nsk)}


    def get_stats(self):
        flnm_stat = self.stat_save['filename']
        x,y,ep,stat, nx,ny = read_stat_proc_2d(flnm_stat)

        return x,y,ep, stat, nx,ny

    def draw_line_x(self, x,y,ep, stat, nx,ny, x_line, compare=False):
        stat_line = {}
        for name in stat:
            stat_line[name] = np.zeros([ny])
            for j in range(ny):
                stat_line[name][j] = np.interp(x_line, x[:,0], stat[name][:,j])

        if compare:
            fig, ax = plt.subplots(2, 2)
        else:
            fig, ax = plt.subplots(3, 2)
        fig.suptitle('Averaged variables and Reynolds stress along y-direction at x = ' + x_line.__str__(), fontdict={'size':28})
        ax[0,0].plot(stat_line['U'], y[0,:], linewidth=2)
        ax[0,0].set_xlabel(r'$U/U_{b}$', fontdict={'size':28})
        ax[0,0].set_ylabel(r'$y/h$', fontdict={'size':28})
        ax[0,0].tick_params(axis='both', labelsize=25)
        ax[0,1].plot(stat_line['V'], y[0,:], linewidth=2)
        ax[0,1].set_xlabel(r'$V/U_{b}$', fontdict={'size':28})
        ax[0,1].set_ylabel(r'$y/h$', fontdict={'size':28})
        ax[0,1].tick_params(axis='both', labelsize=25)
        ax[1,0].plot(stat_line['uu'], y[0,:], linewidth=2)
        ax[1,0].set_xlabel(r'$u''u''/U_{b}^{2}$', fontdict={'size':28})
        ax[1,0].set_ylabel(r'$y/h$', fontdict={'size':28})
        ax[1,0].tick_params(axis='both', labelsize=25)
        ax[1,1].plot(stat_line['vv'], y[0,:], linewidth=2)
        ax[1,1].set_xlabel(r'$v''v''/U_{b}^{2}$', fontdict={'size':28})
        ax[1,1].set_ylabel(r'$y/h$', fontdict={'size':28})
        ax[1,1].tick_params(axis='both', labelsize=25)
        if not compare:
            ax[2,0].plot(stat_line['tau33'], y[0,:])
            ax[2,0].set_xlabel('w\'w\'', fontdict={'size':25})
            ax[2,0].set_ylabel('y/h', fontdict={'size':25})
            ax[2,0].tick_params(axis='both', labelsize=25)
            ax[2,1].plot(stat_line['tau12'], y[0,:])
            ax[2,1].set_xlabel('u\'v\'', fontdict={'size':25})
            ax[2,1].set_ylabel('y/h', fontdict={'size':25})
            ax[2,1].tick_params(axis='both', labelsize=25)
        return stat_line, fig, ax

    def draw_law_of_wall(self, x,y,ep, stat, nx,ny, x_line):
        u = np.zeros(ny)
        for j in range(ny):
            u[j] = np.interp(x_line, x[:,0], stat['umean'][:,j])

        tauw = (u[1]-u[0]) / (y[0,1]-y[0,0])
        ut = np.sqrt(tauw/self.Re)
        Ret = np.sqrt(tauw*self.Re)
        yplus = y[0,:] * Ret
        u /= ut

        fig, ax = plt.subplots(1, 1)
        fig.suptitle(r'$u^{+}$ vs. $y^{+}$ at $x/h$ = ' + x_line.__str__(), fontdict={'size':28})
        ax.semilogx(yplus[:(ny-1)], u[:(ny-1)], linewidth=2)
        ax.plot(yplus[yplus<16], yplus[yplus<16], '--', linewidth=2)
        ax.plot(yplus[np.logical_and(yplus<180, yplus>6)], 1./0.41*np.log(yplus[np.logical_and(yplus<180, yplus>6)]) + 5.4, '-.', linewidth=2)
        ax.legend(['present', r'$u^{+}=y^{+}$', r'$u^{+}=\frac{1}{0.41}log(y^{+})+5.4$'])
        ax.set_xlabel(r'$y+$', fontdict={'size':28})
        ax.set_ylabel(r'$u+$', fontdict={'size':28})
        ax.tick_params(axis='both', labelsize=25)
        return u, fig, ax


    def get_inst(self, gkey, nst,nen,nsk, slc, span_average=False):
        filename_inst = [self.path + self.dirs_dict['TIMESERIES'] + self.flnms_dict['timeseries'][0]] + self.flnms_dict['timeseries'][1:]
        inst = read_inst(filename_inst, gkey, nst, nen, nsk, slc)
        if span_average:
            for name in inst.keys():
                inst[name] = np.mean(inst[name], axis=2, keepdims=True)
        return inst

    def get_grid3d(self, slc=None, staggered=False):
        flnm_grid = self.path + self.dirs_dict['REST'] + self.flnms_dict['grid']
        x, y, z, ep, nx, ny, nz = read_grid(flnm_grid, slc=slc, dim=3, staggered=staggered)
        return x, y, z, ep, nx, ny, nz

    def get_flowdata(self, n):
        filename_flowdata = self.path + self.dirs_dict['REST'] + self.flnms_dict['flowdata'][0] + n.__str__().zfill(8) + self.flnms_dict['flowdata'][1]
        flowdata = read_flowdata(filename_flowdata)

        return flowdata

    def draw_inst_point(self, inst, name, fft=True, dt=1.e-1, axs=None, scale=1., nseg=5):
        N_time = inst['u'].shape[-1]
        T = dt*np.arange(0, N_time, 1, dtype=float)
        peak_psd = []
        psd_list = []

        if axs is None:
            fig, ax = plt.subplots(1,1)
            fig2, ax2 = plt.subplots(1,1)
        else:
            ax, ax2 = axs

        data = inst[name]
        point_data = data.reshape([-1,data.shape[-1]])
        for i in range(point_data.shape[0]):
            ax.plot(T/scale, point_data[i,:])
            ax.set_xlabel(r'$tU_{b}/L_{s}$', fontdict={'size':45})
            ax.set_ylabel('dimensionless '+name, fontdict={'size':25})
            ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0), useMathText=True)
            ax.yaxis.offsetText.set_fontsize(25)
            ax.yaxis.offsetText.set_position((-0.08,.5))
            ax.tick_params(axis='both', labelsize=25)

            ax_top = ax.twiny()
            ax_top.plot(T, point_data[i,:], linewidth=0.01, label='_nolegend_')
            ax_top.set_xbound(np.array(ax.get_xbound())*scale)
            ax_top.set_xlabel(r'$tU_{b}/h$', fontdict={'size':45})
            ax_top.tick_params(axis='both', labelsize=25)

            if fft:
                N = N_time / nseg
                M = N/2
                freq, psd = signal.welch(point_data[i,:], fs=1./dt, window='hanning', nperseg=N, noverlap=M)

                idx_maxpsd = np.argmax(psd*freq)
                print(freq[idx_maxpsd], (freq*scale)[idx_maxpsd])
                peak_psd.append(freq[idx_maxpsd])

                ax2.semilogx(freq*scale, psd*freq)
                ax2.set_xlabel(r'$fL_{s}/U_{b}$', fontdict={'size':45})
                ax2.set_ylabel('Premultiplied PSD of dimensionless '+name, fontdict={'size':25})
                ax2.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
                ax2.yaxis.offsetText.set_fontsize(25)
                ax2.tick_params(axis='both', labelsize=25)
                ax_top = ax2.twiny()
                ax_top.semilogx(freq, psd*freq)
                ax_top.set_xlabel(r'$fh/U_{b}$', fontdict={'size':45})
                ax_top.set_xbound(np.array(ax2.get_xbound())/scale)
                ax_top.tick_params(axis='both', labelsize=25)

                psd_list.append(psd)

        # np.save('peak_psd', peak_psd)
        # f = open('/usr/local/home/yl8bc/Nic/PH/validation/psd/u.dat', mode='w')
        # f.write(N_time.__str__()+'\n')
        # for i in range(N_time):
        #     f.write(inst['u'][0,0,0,i].__str__()+'\n')
        # f.close()
        return ax, ax2, freq, psd_list

    def get_wall(self, x,y,ep):
        wu_list = find_wall_2d(x,y,ep)
        x_w = []
        y_w = []
        for wu in wu_list:
            x_w.append( x[wu['index']])
            y_w.append( y[wu['index']])

        return wu_list, x_w, y_w

    def draw_Cf(self, wu_list, x_w, y_w, stat=None, slc=None, inst=None, time_step=[0]):
        Cf_stat = []
        Cf_inst = []
        for wu  in wu_list:
            oy = np.sum(wu['neighbors'][-2:])
            ox = np.sum(wu['neighbors'][:2])
            wx = np.abs(oy) / (np.abs(ox)+np.abs(oy))
            wy = np.abs(ox) / (np.abs(ox)+np.abs(oy))
            if stat:
                u = stat['U'][slc][wu['index']]
                v = stat['V'][slc][wu['index']]
                tauw = u*oy*wx + v*ox*wy
                Cf = tauw*2./self.Re
                Cf_stat.append(Cf)
            if inst:
                u = inst['u'][wu['index']+(0,time_step)]
                v = inst['v'][wu['index']+(0,time_step)]
                tauw = u*oy*wx + v*ox*wy
                Cf = tauw*2./self.Re
                Cf_inst.append(Cf)

        sorted_index = np.argsort(x_w, axis=-1)
        x_w = np.array(x_w)[sorted_index]
        y_w = np.array(y_w)[sorted_index]
        if stat is not None: Cf_stat = np.array(Cf_stat)[sorted_index]
        if inst is not None: Cf_inst = np.array(Cf_inst)[sorted_index]

        fig, ax = plt.subplots()
        # plt.plot(x_w, y_w, '.')
        ax.plot(x_w, np.zeros_like(x_w), '-')
        if stat: ax.plot(x_w, Cf_stat, '-')
        if inst: ax.plot(x_w, Cf_inst, '--')
        ax.set_xlabel('x', fontdict={'size':30})
        ax.set_ylabel('Cf', fontdict={'size':30})
        ax.tick_params(axis='both', labelsize=30)

        return x_w, y_w, Cf_stat, Cf_inst

    def find_reattachment_point_inst(self, x_w, Cf_inst, fft=False, dt=1.e-2, scale=1.):
        N_point, N_time = Cf_inst.shape

        T = dt*np.arange(0, N_time, 1, dtype=float)

        xr_list = []
        for n in range(N_time):
            for i in range(N_point):
                if x_w[i]>2.3 and Cf_inst[i,n] > 0.:
                    if i==0:
                        xr = x_w[i]
                    else:
                        xr = np.interp(0., Cf_inst[i-1:i,n], x_w[i-1:i])
                    xr_list.append(xr)
                    break
                if i==N_point-1: xr_list.append(x_w[i])

        if fft:
            N = N_time
            M = N/2
            freq, psd = signal.welch(xr_list, fs=1./dt, window='hamming', nperseg=N, noverlap=M)
            idx_maxpsd = np.argmax(psd*freq)
            print(freq[idx_maxpsd], (freq*scale)[idx_maxpsd], (psd*freq)[idx_maxpsd])
        if fft:
            fig, ax = plt.subplots(1,1)
            ax.plot(T/scale, xr_list)
            ax.set_xlabel(r'$tU_{b}/L_{s}$', fontdict={'size':45})
            ax.set_ylabel(r'reattachment point $x_{r}/h$', fontdict={'size':25})
            ax.tick_params(axis='both', labelsize=25)
            fig, ax2 = plt.subplots(1,1)
            ax2.semilogx(freq*scale, freq*psd)
            ax2.set_ylabel('Premultiplied PSD of reattachment point', fontdict={'size':25})
            ax2.set_xlabel(r'$fL_{s}/U_{b}$', fontdict={'size':45})
            ax2.tick_params(axis='both', labelsize=25)

            ax_top = ax.twiny()
            ax_top.plot(T, xr_list, linewidth=0.01)
            ax_top.set_xlabel(r'$tU_{b}/h$', fontdict={'size':45})
            ax_top.set_xbound(np.array(ax.get_xbound())*scale)
            ax_top.tick_params(axis='both', labelsize=25)
            ax_top = ax2.twiny()
            ax_top.semilogx(freq, freq*psd, linewidth=0.01)
            ax_top.set_xlabel(r'$fh/U_{b}$', fontdict={'size':45})
            ax_top.set_xbound(np.array(ax2.get_xbound())/scale)
            ax_top.tick_params(axis='both', labelsize=25)
        else:
            fig, ax = plt.subplots(1,1)
            ax.plot(T, xr_list)
            ax.set_xlabel(r'$th/U_{b}$', fontdict={'size':25})
            ax.set_ylabel(r'$x/h$', fontdict={'size':25})
            ax.tick_params(axis='both', labelsize=25)
        if fft:
            pass

    def FIR_filter_inst(self,  inst, dt=5.e-2, ntaps = 50):
        N_time = inst['u'].shape[-1]
        shape = inst['u'].shape
        b = signal.firwin(ntaps, 1./(1./dt))
        inst_filtered = {}
        for name, data in inst.items():
            inst_filtered[name] = np.zeros(shape)
            for i in range(shape[0]):
                for j in range(shape[1]):
                    for k in range(shape[2]):
                        timeseries = data[i,j,k,:]
                        inst_filtered[name][i,j,k,:] = signal.filtfilt(b, [1.],timeseries)#[ntaps:-ntaps]
        return inst_filtered

    def draw_AIM(self, x,y,ep, stat, nx,ny, x_line, ax=None):
        stat_line = {}
        for name in stat:
            stat_line[name] = np.zeros([ny])
            for j in range(ny):
                stat_line[name][j] = np.interp(x_line, x[:,0], stat[name][:,j])
        a = np.zeros([ny,3,3])
        tau = np.zeros([ny,3,3])
        tau[:,0,0] = stat_line['uu']
        tau[:,0,1] = stat_line['uv']
        tau[:,0,2] = stat_line['uw']
        tau[:,1,0] = stat_line['uv']
        tau[:,1,1] = stat_line['vv']
        tau[:,1,2] = stat_line['vw']
        tau[:,2,0] = stat_line['uw']
        tau[:,2,1] = stat_line['vw']
        tau[:,2,2] = stat_line['ww']
        ukuk = tau[:,0,0]+tau[:,1,1]+tau[:,2,2]
        for i in range(3):
            for j in range(3):
                a[:,i,j] = tau[:,i,j] / ukuk
                if i==j: a[:,i,j] -= 1./3.
        # for i in range(3):
        #     for j in range(3-i,3):
        #         a[:,i,j] = a[:,j,i]
        II = np.zeros([ny])
        III = np.zeros([ny])
        for i in range(3):
            for j in range(3):
                II = II + a[:,i,j]*a[:,j,i]
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    III = III + a[:,i,j]*a[:,j,k]*a[:,k,i]
        ## draw stuff
        if ax is None:
            fig, ax = plt.subplots()
            xx = np.arange(-1./36, 2./9., .01)
            ax.plot(xx, 2./9+2.*xx, '--', label='_nolegend_')
            xx = np.arange(0., 2./9., .01)
            ax.plot(xx, 3./2.*(4./3.*xx)**(2./3.), '--', label='_nolegend_')
            fig.suptitle('Anisotropy Invariant Map at x = ' + x_line.__str__(), fontdict={'size':28})

        ax.plot(III, II,'-', linewidth=2)
        ax.set_xlabel('III', fontdict={'size':28})
        ax.set_ylabel('II', fontdict={'size':28})
        ax.tick_params(axis='both', labelsize=25)

        return ax

class hills:
    filename = 'phs.npy'
    def __init__(self):
        if os.path.isfile(self.filename):
            self.hills_list = np.load(self.filename).tolist()
        else:
            self.hills_list = []

    def search_hill(self, name):
        ifind = False
        for ph in self.hills_list:
            if ph.name == name:
                ifind = True
                return ph
        if not ifind:
            print('This case does not exist!')
            return None

    def add_hill(self, ph_new):
        ph = self.search_hill(ph_new.name)
        if ph is None:
            self.hills_list.append(ph_new)
        else:
            print('This case already exists!')

    def delete_hill(self, name):
        ph = self.search_hill(name)
        if ph is None:
            print('This case does not exist!')
        else:
            self.hills_list.remove(ph)

    def save(self):
        np.save(self.filename, self.hills_list)

    def show_hills(self):
        for ph in self.hills_list:
            print('Name: ' + ph.name + ', averaged by ', ph.stat_save['interval'], ', in ', ph.path)