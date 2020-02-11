import h5py
from . import IO_util
import numpy as np
from .derivative_lele import *
from .refine_1d import *

def collect_stat(filename_tuple, namelist=None, group='stat2d'):
    stat_out = IO_util.read_hdf5(filename_tuple[0], namelist=namelist)
    stat_out = stat_out[group]
    for filename in filename_tuple[1:]:
        stat = IO_util.read_hdf5(filename, namelist=namelist)
        for key in stat_out.keys():
            stat_out[key] += stat[group][key]
    
    num_step = float(len(filename_tuple))
    for key in stat_out.keys():
        stat_out[key] /= num_step

    return stat_out
            
    



def convert_stat_2d(filename_grid, filename_stat, nst, nen, nsk, filename_op):
    ## grid info
    fid = h5py.File(filename_grid, 'r')
    x = fid.get('x').__array__().T[:,:,0]
    y = fid.get('y').__array__().T[:,:,0]
    ep = fid.get('ep').__array__().T[:,:,0]
    nx,ny = x.shape

    ## read data from HDF5
    ntot = (nen-nst+1) / nsk
    statname = ['umean','vmean','wmean','pmean','uumean','vvmean','wwmean','uvmean','uwmean','vwmean']
    data = np.zeros([nx,ny,len(statname)])
    for step in range(nst,nen,nsk):
        filename = filename_stat[0] + step.__str__().zfill(8) + filename_stat[1] + (step+nsk-1).__str__().zfill(8)+ filename_stat[2]
        print(filename)
        fid = h5py.File(filename, 'r')
        gid = fid.get('stat2d')
        for idset in range(statname.__len__()):
            data[:,:,idset] += gid.get(statname[idset]).__array__().T[:,:]
        fid.close()

    ## averaging temporally
    for idset in range(statname.__len__()):
        data[:,:,idset] /= ntot

    ## get Reynolds stress
    rstressname = ['tau11','tau22','tau33','tau12','tau13','tau23']
    rstress = np.zeros([nx,ny,6])
    rstress[:,:,0] = data[:,:,4] - data[:,:,0]*data[:,:,0]
    rstress[:,:,1] = data[:,:,5] - data[:,:,1]*data[:,:,1]
    rstress[:,:,2] = data[:,:,6] - data[:,:,2]*data[:,:,2]
    rstress[:,:,3] = data[:,:,7] - data[:,:,0]*data[:,:,1]
    rstress[:,:,4] = data[:,:,8] - data[:,:,0]*data[:,:,2]
    rstress[:,:,5] = data[:,:,9] - data[:,:,1]*data[:,:,2]

    ## output
    fid = h5py.File(filename_op, 'w')
    did = fid.create_dataset('x', [nx,ny], data=x)
    did = fid.create_dataset('y', [nx,ny], data=y)
    did = fid.create_dataset('ep', [nx,ny], data=ep)
    for idset in range(4):
        did = fid.create_dataset(statname[idset], [nx,ny], data=data[:,:,idset])
    for idset in range(6):
        did = fid.create_dataset(rstressname[idset], [nx,ny], data=rstress[:,:,idset])
    fid.close()
    return



def convert_stat_budgets_2d(filename_grid, filename_stat, nst, nen, nsk, filename_op):
    ## grid info
    fid = h5py.File(filename_grid, 'r')
    x = fid.get('x').__array__().T[:,:,0]
    y = fid.get('y').__array__().T[:,:,0]
    ep = fid.get('ep').__array__().T[:,:,0]
    nx,ny = x.shape

    ## read data from HDF5
    ntot = (nen-nst+1) / nsk
    statname = ['umean','vmean','wmean','pmean', \
                'uumean','vvmean','wwmean','uvmean','uwmean','vwmean','pumean','pvmean','pwmean','ppmean', \
                'uuumean','vvvmean','wwwmean','uuvmean','uuwmean','vvumean','vvwmean','wwumean','wwvmean','uvwmean', \
                'pdudxmean','pdudymean','pdudzmean','pdvdxmean','pdvdymean','pdvdzmean','pdwdxmean','pdwdymean','pdwdzmean', \
                'dudumean','dvdvmean','dwdwmean','dudvmean','dudwmean','dvdwmean']#33-38
    data = np.zeros([nx,ny,len(statname)])
    for step in range(nst,nen,nsk):
        filename = filename_stat[0] + step.__str__().zfill(8) + filename_stat[1] + (step+nsk-1).__str__().zfill(8)+ filename_stat[2]
        print(filename)
        fid = h5py.File(filename, 'r')
        gid = fid.get('stat2d')
        for idset in range(statname.__len__()):
            data[:,:,idset] += gid.get(statname[idset]).__array__().T[:,:]
        fid.close()

    ## averaging temporally
    for idset in range(statname.__len__()):
        data[:,:,idset] /= ntot

    ## get Reynolds stress
    uu = data[:,:,4] - data[:,:,0]*data[:,:,0]
    vv = data[:,:,5] - data[:,:,1]*data[:,:,1]
    ww = data[:,:,6] - data[:,:,2]*data[:,:,2]
    uv = data[:,:,7] - data[:,:,0]*data[:,:,1]
    uw = data[:,:,8] - data[:,:,0]*data[:,:,2]
    vw = data[:,:,9] - data[:,:,1]*data[:,:,2]
    uuu = data[:,:,14] - data[:,:,0]^3 - 3*data[:,:,0]*uu
    vvv = data[:,:,15] - data[:,:,1]^3 - 3*data[:,:,1]*vv
    www = data[:,:,16] - data[:,:,2]^3 - 3*data[:,:,2]*ww
    uuv = data[:,:,17] - data[:,:,0]^2*data[:,:,1] - 2*data[:,:,0]*uv - data[:,:,1]*uu
    uuw = data[:,:,18] - data[:,:,0]^2*data[:,:,2] - 2*data[:,:,0]*uw - data[:,:,2]*uu
    vvu = data[:,:,19] - data[:,:,0]*data[:,:,1]^2 - 2*data[:,:,1]*uv - data[:,:,0]*vv
    vvw = data[:,:,20] - data[:,:,2]*data[:,:,1]^2 - 2*data[:,:,1]*vw - data[:,:,2]*vv
    wwu = data[:,:,21] - data[:,:,0]*data[:,:,2]^2 - 2*data[:,:,2]*uw - data[:,:,0]*ww
    wwv = data[:,:,22] - data[:,:,1]*data[:,:,2]^2 - 2*data[:,:,2]*vw - data[:,:,2]*vv
    uvw = data[:,:,23] - data[:,:,0]*data[:,:,1]*data[:,:,2] - data[:,:,0]*vw - data[:,:,1]*uw - data[:,:,2]*uv

    yly = y[0,ny-1] - y[0,0]
    xlx = x[nx-1,0] - x[0,0]
    ncly = 2
    nclx = 0
    istret = 2
    beta = 2
    yp,ypi,ppy,pp2y,pp4y,ppyi,pp2yi,pp4yi = refine_incompact3d_1d(yly, ncly, ny, istret, beta)
    ##
    dy = yly / (ny+1)
    dx = xlx / nx

    dUdx = np.zeros([nx,ny])
    dUdx = np.zeros([nx,ny])
    dUdx = np.zeros([nx,ny])
    dUdx = np.zeros([nx,ny])
    for j in range(ny):
        stat1[:,j,0] = derive_1st_lele(stat[:,j,0], dx, nx, nclx)
        stat1[:,j,2] = derive_1st_lele(stat[:,j,1], dx, nx, nclx)
        stat1[:,j,4] = derive_1st_lele(stat[:,j,2], dx, nx, nclx)
        stat1[:,j,6] = derive_1st_lele(stat[:,j,3], dx, nx, nclx)
        stat2[:,j,0] = derive_2nd_lele(stat[:,j,0], dx, nx, nclx)
        stat2[:,j,2] = derive_2nd_lele(stat[:,j,1], dx, nx, nclx)
        stat2[:,j,4] = derive_2nd_lele(stat[:,j,2], dx, nx, nclx)
    for i in range(nx):
        stat1[i,:,1] = derive_1st_lele(stat[i,:,0], dy, ny, ncly)
        stat1[i,:,3] = derive_1st_lele(stat[i,:,1], dy, ny, ncly)
        stat1[i,:,5] = derive_1st_lele(stat[i,:,2], dy, ny, ncly)
        stat1[i,:,7] = derive_1st_lele(stat[i,:,3], dy, ny, ncly)
        stat1[i,:,1] *= ppy
        stat1[i,:,3] *= ppy
        stat1[i,:,5] *= ppy
        stat1[i,:,7] *= ppy
        stat2[i,:,1] = derive_2nd_lele(stat[i,:,0], dy, ny, ncly)
        stat2[i,:,3] = derive_2nd_lele(stat[i,:,1], dy, ny, ncly)
        stat2[i,:,5] = derive_2nd_lele(stat[i,:,2], dy, ny, ncly)
        stat2[i,:,1] = stat2[i,:,1]*pp2y - pp4y*stat1[i,:,1]
        stat2[i,:,3] = stat2[i,:,3]*pp2y - pp4y*stat1[i,:,3]
        stat2[i,:,5] = stat2[i,:,5]*pp2y - pp4y*stat1[i,:,5]
    ## calculate budgets 2d
    P11 = -(uu*dUdx+uv*dUdy)*2
    P22 = -(uv*dVdx+vv*dVdy)*2
    P33 = -(uw*dWdx+ww*dWdy)*2
    P12 = -(uv*dUdx+vv*dUdy) -(uu*dVdx+uv*dVdy)
    P13 = -(uw*dUdx+vw*dUdy) -(uu*dWdx+uv*dWdy)
    P23 = -(uw*dVdx+vw*dVdy) -(uv*dWdx+vv*dWdy)

    T11 = -(duuudx+duuvdy)
    T22 = -(dvvudx+dvvvdy)
    T33 = -(dwwudx+dwwvdy)
    T12 = -(duuvdx+dvvudy)
    T13 = -(duuwdx+duvwdy)
    T23 = -(duvwdx+dvvwdy)

    D11 = 1/Re * (d2uudx+d2uudy)
    D22 = 1/Re * (d2vvdx+d2vvdy)
    D33 = 1/Re * (d2wwdx+d2wwdy)
    D12 = 1/Re * (d2uvdx+d2uvdy)
    D13 = 1/Re * (d2uwdx+d2uwdy)
    D23 = 1/Re * (d2vwdx+d2vwdy)

    Dp11 = -(dpudx*2)
    Dp22 = -(dpvdy*2)
    Dp33 = np.zeros([nx,ny])
    Dp12 = -(dpvdx + dpudy)
    Dp13 = -(dpwdx)
    Dp23 = -(dpwdy)

    Fi11 = pdudx * 2
    Fi22 = pdvdy * 2
    Fi33 = np.zeros([nx,ny])
    Fi12 = pdudy + pdvdx
    Fi13 = pdwdx
    Fi23 = pdwdy

    ep11 = 2/Re * (dudu)
    ep22 = 2/Re * (dvdv)
    ep33 = 2/Re * (dwdw)
    ep12 = 2/Re * (dudv)
    ep13 = 2/Re * (dudw)
    ep23 = 2/Re * (dvdw)


    ## output
    fid = h5py.File(filename_op, 'w')
    did = fid.create_dataset('x', [nx,ny], data=x)
    did = fid.create_dataset('y', [nx,ny], data=y)
    did = fid.create_dataset('ep', [nx,ny], data=ep)
    for idset in range(4):
        did = fid.create_dataset(statname[idset], [nx,ny], data=data[:,:,idset])
    for idset in range(6):
        did = fid.create_dataset(rstressname[idset], [nx,ny], data=rstress[:,:,idset])
    fid.close()
    return


def compute_stats_simple(stat):
    ## wrap the basic stats
    stat_out = {}
    stat_out['U'] = stat['umean']
    stat_out['V'] = stat['vmean']
    stat_out['W'] = stat['wmean']
    stat_out['P'] = stat['pmean']
    stat_out['uu'] = stat['uumean'] - stat_out['U'] * stat_out['U']
    stat_out['vv'] = stat['vvmean'] - stat_out['V'] * stat_out['V']
    stat_out['ww'] = stat['wwmean'] - stat_out['W'] * stat_out['W']
    stat_out['uv'] = stat['uvmean'] - stat_out['U'] * stat_out['V']
    stat_out['uw'] = stat['uwmean'] - stat_out['U'] * stat_out['W']
    stat_out['vw'] = stat['vwmean'] - stat_out['V'] * stat_out['W']

    ## temperature
    if 'tmean' in stat:
        stat_out['T'] = stat['tmean']
        stat_out['tt'] = stat['ttmean'] - stat_out['T'] * stat_out['T']

    return stat_out

def compute_stats_2d(stat, dx, dy, Re):
    nx, ny = stat['umean'].shape
    ## wrap the basic stats
    stat_out = {}
    stat_out['U'] = stat['umean']
    stat_out['V'] = stat['vmean']
    stat_out['W'] = stat['wmean']
    stat_out['P'] = stat['pmean']
    stat_out['uu'] = stat['uumean'] - stat_out['U']*stat_out['U']
    stat_out['vv'] = stat['vvmean'] - stat_out['V']*stat_out['V']
    stat_out['ww'] = stat['wwmean'] - stat_out['W']*stat_out['W']
    stat_out['uv'] = stat['uvmean'] - stat_out['U']*stat_out['V']
    stat_out['uw'] = stat['uwmean'] - stat_out['U']*stat_out['W']
    stat_out['vw'] = stat['vwmean'] - stat_out['V']*stat_out['W']
    ## get budgets breakdowns
    stat['uuu'] = stat['uuumean'] - 3.*stat_out['U']*stat_out['uu'] - np.power(stat_out['U'],3)
    stat['vvv'] = stat['vvvmean'] - 3.*stat_out['V']*stat_out['vv'] - np.power(stat_out['V'],3)
    stat['www'] = stat['wwwmean'] - 3.*stat_out['W']*stat_out['ww'] - np.power(stat_out['W'],3)
    stat['uuv'] = stat['uuvmean'] - 2.*stat_out['U']*stat_out['uv'] - 1.*stat_out['V']*stat_out['uu'] - np.power(stat_out['U'],2)*stat_out['V']
    stat['uuw'] = stat['uuwmean'] - 2.*stat_out['U']*stat_out['uw'] - 1.*stat_out['W']*stat_out['uu'] - np.power(stat_out['U'],2)*stat_out['W']
    stat['vvu'] = stat['vvumean'] - 2.*stat_out['V']*stat_out['uv'] - 1.*stat_out['U']*stat_out['vv'] - np.power(stat_out['V'],2)*stat_out['U']
    stat['vvw'] = stat['vvwmean'] - 2.*stat_out['V']*stat_out['vw'] - 1.*stat_out['W']*stat_out['vv'] - np.power(stat_out['V'],2)*stat_out['W']
    stat['wwu'] = stat['wwumean'] - 2.*stat_out['W']*stat_out['uw'] - 1.*stat_out['U']*stat_out['ww'] - np.power(stat_out['W'],2)*stat_out['U']
    stat['wwv'] = stat['wwvmean'] - 2.*stat_out['W']*stat_out['vw'] - 1.*stat_out['V']*stat_out['ww'] - np.power(stat_out['W'],2)*stat_out['V']
    stat['uvw'] = stat['uvwmean'] - 1.*stat_out['U']*stat_out['vw'] - 1.*stat_out['V']*stat_out['uw'] - 1.*stat_out['W']*stat_out['uv'] \
                  - stat_out['U']*stat_out['V']*stat_out['W']
    stat['pu'] = stat['pumean'] - stat_out['P']*stat_out['U']
    stat['pv'] = stat['pvmean'] - stat_out['P']*stat_out['V']
    stat['pw'] = stat['pwmean'] - stat_out['P']*stat_out['W']


    stat['dUdx'] = np.zeros_like(stat_out['U'])
    stat['dUdy'] = np.zeros_like(stat_out['U'])
    stat['dVdx'] = np.zeros_like(stat_out['U'])
    stat['dVdy'] = np.zeros_like(stat_out['U'])
    stat['dWdx'] = np.zeros_like(stat_out['U'])
    stat['dWdy'] = np.zeros_like(stat_out['U'])
    stat['duuudx'] = np.zeros_like(stat_out['U'])
    stat['duuvdx'] = np.zeros_like(stat_out['U'])
    stat['duuwdx'] = np.zeros_like(stat_out['U'])
    stat['dvvudx'] = np.zeros_like(stat_out['U'])
    stat['dwwudx'] = np.zeros_like(stat_out['U'])
    stat['duvwdx'] = np.zeros_like(stat_out['U'])
    stat['dvvvdy'] = np.zeros_like(stat_out['U'])
    stat['dvvudy'] = np.zeros_like(stat_out['U'])
    stat['dvvwdy'] = np.zeros_like(stat_out['U'])
    stat['duuvdy'] = np.zeros_like(stat_out['U'])
    stat['dwwvdy'] = np.zeros_like(stat_out['U'])
    stat['duvwdy'] = np.zeros_like(stat_out['U'])
    stat['dduudxx'] = np.zeros_like(stat_out['U'])
    stat['dduudyy'] = np.zeros_like(stat_out['U'])
    stat['ddvvdxx'] = np.zeros_like(stat_out['U'])
    stat['ddvvdyy'] = np.zeros_like(stat_out['U'])
    stat['ddwwdxx'] = np.zeros_like(stat_out['U'])
    stat['ddwwdyy'] = np.zeros_like(stat_out['U'])
    stat['dduvdxx'] = np.zeros_like(stat_out['U'])
    stat['dduvdyy'] = np.zeros_like(stat_out['U'])
    stat['dduwdxx'] = np.zeros_like(stat_out['U'])
    stat['dduwdyy'] = np.zeros_like(stat_out['U'])
    stat['ddvwdxx'] = np.zeros_like(stat_out['U'])
    stat['ddvwdyy'] = np.zeros_like(stat_out['U'])
    stat['duudy'] = np.zeros_like(stat_out['U'])
    stat['dvvdy'] = np.zeros_like(stat_out['U'])
    stat['dwwdy'] = np.zeros_like(stat_out['U'])
    stat['duvdy'] = np.zeros_like(stat_out['U'])
    stat['duwdy'] = np.zeros_like(stat_out['U'])
    stat['dvwdy'] = np.zeros_like(stat_out['U'])
    stat['dpudx'] = np.zeros_like(stat_out['U'])
    stat['dpvdx'] = np.zeros_like(stat_out['U'])
    stat['dpwdx'] = np.zeros_like(stat_out['U'])
    stat['dpudy'] = np.zeros_like(stat_out['U'])
    stat['dpvdy'] = np.zeros_like(stat_out['U'])
    stat['dpwdy'] = np.zeros_like(stat_out['U'])
    yp,ypi,ppy,pp2y,pp4y,ppyi,pp2yi,pp4yi = refine_incompact3d_1d(3.036, 2, ny, 2, 2.)
    for j in range(ny):
        stat['dUdx'][:,j] = derive_1st_lele(stat_out['U'][:,j], dx, nx, 0 )
        stat['dVdx'][:,j] = derive_1st_lele(stat_out['V'][:,j], dx, nx, 0 )
        stat['dWdx'][:,j] = derive_1st_lele(stat_out['W'][:,j], dx, nx, 0 )
        stat['duuudx'][:,j] = derive_1st_lele(stat['uuu'][:,j], dx, nx, 0 )
        stat['dvvudx'][:,j] = derive_1st_lele(stat['vvu'][:,j], dx, nx, 0 )
        stat['dwwudx'][:,j] = derive_1st_lele(stat['wwu'][:,j], dx, nx, 0 )
        stat['duuvdx'][:,j] = derive_1st_lele(stat['uuv'][:,j], dx, nx, 0 )
        stat['duvwdx'][:,j] = derive_1st_lele(stat['uvw'][:,j], dx, nx, 0 )
        stat['dduudxx'][:,j] = derive_2nd_lele(stat_out['uu'][:,j], dx, nx, 0 )
        stat['ddvvdxx'][:,j] = derive_2nd_lele(stat_out['vv'][:,j], dx, nx, 0 )
        stat['ddwwdxx'][:,j] = derive_2nd_lele(stat_out['ww'][:,j], dx, nx, 0 )
        stat['dduvdxx'][:,j] = derive_2nd_lele(stat_out['uv'][:,j], dx, nx, 0 )
        stat['dduwdxx'][:,j] = derive_2nd_lele(stat_out['uw'][:,j], dx, nx, 0 )
        stat['ddvwdxx'][:,j] = derive_2nd_lele(stat_out['vw'][:,j], dx, nx, 0 )
        stat['dpudx'][:,j] = derive_1st_lele(stat['pu'][:,j], dx, nx, 0 )
        stat['dpvdx'][:,j] = derive_1st_lele(stat['pv'][:,j], dx, nx, 0 )
        stat['dpwdx'][:,j] = derive_1st_lele(stat['pw'][:,j], dx, nx, 0 )

    for i in range(nx):
        stat['dUdy'][i,:] = derive_1st_lele(stat_out['U'][i,:], dy, ny, 2 )*ppy
        stat['dVdy'][i,:] = derive_1st_lele(stat_out['V'][i,:], dy, ny, 2 )*ppy
        stat['dWdy'][i,:] = derive_1st_lele(stat_out['W'][i,:], dy, ny, 2 )*ppy
        stat['dvvvdy'][i,:] = derive_1st_lele(stat['vvv'][i,:], dy, ny, 2 )*ppy
        stat['dvvudy'][i,:] = derive_1st_lele(stat['vvu'][i,:], dy, ny, 2 )*ppy
        stat['dvvwdy'][i,:] = derive_1st_lele(stat['vvw'][i,:], dy, ny, 2 )*ppy
        stat['duuvdy'][i,:] = derive_1st_lele(stat['uuv'][i,:], dy, ny, 2 )*ppy
        stat['dwwvdy'][i,:] = derive_1st_lele(stat['wwv'][i,:], dy, ny, 2 )*ppy
        stat['duvwdy'][i,:] = derive_1st_lele(stat['uvw'][i,:], dy, ny, 2 )*ppy
        stat['duudy'][i,:] = derive_1st_lele(stat_out['uu'][i,:], dy, ny, 2 )*ppy
        stat['dduudyy'][i,:] = derive_1st_lele(stat['duudy'][i,:], dy, ny, 2 )*ppy
        stat['dvvdy'][i,:] = derive_1st_lele(stat_out['vv'][i,:], dy, ny, 2 )*ppy
        stat['ddvvdyy'][i,:] = derive_1st_lele(stat['dvvdy'][i,:], dy, ny, 2 )*ppy
        stat['dwwdy'][i,:] = derive_1st_lele(stat_out['ww'][i,:], dy, ny, 2 )*ppy
        stat['ddwwdyy'][i,:] = derive_1st_lele(stat['dwwdy'][i,:], dy, ny, 2 )*ppy
        stat['duvdy'][i,:] = derive_1st_lele(stat_out['uv'][i,:], dy, ny, 2 )*ppy
        stat['dduvdyy'][i,:] = derive_1st_lele(stat['duvdy'][i,:], dy, ny, 2 )*ppy
        stat['duwdy'][i,:] = derive_1st_lele(stat_out['uw'][i,:], dy, ny, 2 )*ppy
        stat['dduwdyy'][i,:] = derive_1st_lele(stat['duwdy'][i,:], dy, ny, 2 )*ppy
        stat['dvwdy'][i,:] = derive_1st_lele(stat_out['vw'][i,:], dy, ny, 2 )*ppy
        stat['ddvwdyy'][i,:] = derive_1st_lele(stat['dvwdy'][i,:], dy, ny, 2 )*ppy
        stat['dpudy'][i,:] = derive_1st_lele(stat['pu'][i,:], dy, ny, 2 )
        stat['dpvdy'][i,:] = derive_1st_lele(stat['pv'][i,:], dy, ny, 2 )
        stat['dpwdy'][i,:] = derive_1st_lele(stat['pw'][i,:], dy, ny, 2 )

    stat['pdudx'] = stat['pdudxmean'] - stat_out['P']*stat['dUdx']
    stat['pdudy'] = stat['pdudymean'] - stat_out['P']*stat['dUdy']
    stat['pdudz'] = stat['pdudzmean']
    stat['pdvdx'] = stat['pdvdxmean'] - stat_out['P']*stat['dVdx']
    stat['pdvdy'] = stat['pdvdymean'] - stat_out['P']*stat['dVdy']
    stat['pdvdz'] = stat['pdvdzmean']
    stat['pdwdx'] = stat['pdwdxmean'] - stat_out['P']*stat['dWdx']
    stat['pdwdy'] = stat['pdwdymean'] - stat_out['P']*stat['dWdy']
    stat['pdwdz'] = stat['pdwdzmean']
    ## get budgets
    stat_out['P11'] = - (stat_out['uu']*stat['dUdx'] + stat_out['uv']*stat['dUdy'] \
                         + stat_out['uu']*stat['dUdx'] + stat_out['uv']*stat['dUdy'])
    stat_out['P22'] = - (stat_out['uv']*stat['dVdx'] + stat_out['vv']*stat['dVdy'] \
                         + stat_out['uv']*stat['dVdx'] + stat_out['vv']*stat['dVdy'])
    stat_out['P33'] = - (stat_out['uw']*stat['dWdx'] + stat_out['vw']*stat['dWdy'] \
                         + stat_out['uw']*stat['dWdx'] + stat_out['vw']*stat['dWdy'])
    stat_out['P12'] = - (stat_out['uu']*stat['dVdx'] + stat_out['uv']*stat['dVdy'] \
                         + stat_out['uv']*stat['dUdx'] + stat_out['vv']*stat['dUdy'])
    stat_out['P13'] = - (stat_out['uu']*stat['dWdx'] + stat_out['uv']*stat['dWdy'] \
                         + stat_out['uw']*stat['dUdx'] + stat_out['vw']*stat['dUdy'])
    stat_out['P23'] = - (stat_out['uv']*stat['dWdx'] + stat_out['vv']*stat['dWdy'] \
                         + stat_out['uw']*stat['dVdx'] + stat_out['vw']*stat['dVdy'])
    stat_out['Pk'] = 0.5 * (stat_out['P11'] + stat_out['P22'] + stat_out['P33'])
    stat_out['T11'] = - (stat['duuudx'] + stat['duuvdy'])
    stat_out['T22'] = - (stat['dvvudx'] + stat['dvvvdy'])
    stat_out['T33'] = - (stat['dwwudx'] + stat['dwwvdy'])
    stat_out['T12'] = - (stat['duuvdx'] + stat['dvvudy'])
    stat_out['T13'] = - (stat['duuwdx'] + stat['duvwdy'])
    stat_out['T23'] = - (stat['duvwdx'] + stat['dvvwdy'])
    stat_out['Tk'] = 0.5 * (stat_out['T11'] + stat_out['T22'] + stat_out['T33'])
    stat_out['D11'] = - (stat['dduudxx'] + stat['dduudyy'])/Re
    stat_out['D22'] = - (stat['ddvvdxx'] + stat['ddvvdyy'])/Re
    stat_out['D33'] = - (stat['ddwwdxx'] + stat['ddwwdyy'])/Re
    stat_out['D12'] = - (stat['dduvdxx'] + stat['dduvdyy'])/Re
    stat_out['D13'] = - (stat['dduwdxx'] + stat['dduwdyy'])/Re
    stat_out['D23'] = - (stat['ddvwdxx'] + stat['ddvwdyy'])/Re
    stat_out['Dk'] = 0.5 * (stat_out['D11'] + stat_out['D22'] + stat_out['D33'])
    stat_out['Dp11'] = - (stat['dpudx'] + stat['dpudx'])
    stat_out['Dp22'] = - (stat['dpvdy'] + stat['dpvdy'])
    stat_out['Dp12'] = - (stat['dpudy'] + stat['dpvdx'])
    stat_out['Dp13'] = - (stat['dpwdx'])
    stat_out['Dp23'] = - (stat['dpwdy'])
    stat_out['Dpk'] = 0.5 * (stat_out['Dp11'] + stat_out['Dp22'])
    stat_out['F11'] = stat['pdudx'] + stat['pdudx']
    stat_out['F22'] = stat['pdvdy'] + stat['pdvdy']
    stat_out['F33'] = stat['pdwdz'] + stat['pdwdz']
    stat_out['F12'] = stat['pdudy'] + stat['pdvdx']
    stat_out['F13'] = stat['pdudz'] + stat['pdwdx']
    stat_out['F23'] = stat['pdvdz'] + stat['pdwdy']
    stat_out['Fk'] = 0.5 * (stat_out['F11'] + stat_out['F22'] + stat_out['F33'])
    stat_out['e11'] = 2./Re * (stat['dudumean'] - stat['dUdx']*stat['dUdx'] - stat['dUdy']*stat['dUdy'])
    stat_out['e22'] = 2./Re * (stat['dvdvmean'] - stat['dVdx']*stat['dVdx'] - stat['dVdy']*stat['dVdy'])
    stat_out['e33'] = 2./Re * (stat['dwdwmean'] - stat['dWdx']*stat['dWdx'] - stat['dWdy']*stat['dWdy'])
    stat_out['e12'] = 2./Re * (stat['dudvmean'] - stat['dUdx']*stat['dVdx'] - stat['dUdy']*stat['dVdy'])
    stat_out['e13'] = 2./Re * (stat['dudwmean'] - stat['dUdx']*stat['dWdx'] - stat['dUdy']*stat['dWdy'])
    stat_out['e23'] = 2./Re * (stat['dvdwmean'] - stat['dVdx']*stat['dWdx'] - stat['dVdy']*stat['dWdy'])
    stat_out['ek'] = 0.5 * (stat_out['e11'] + stat_out['e22'] + stat_out['e33'])

    ## futher
    stat_out['k'] = 0.5*(stat_out['uu'] + stat_out['vv'] + stat_out['ww'])
    stat_out['b'] = np.sqrt( np.power(stat_out['uu']/stat_out['k']/2.-1./3., 2) + np.power(stat_out['vv']/stat_out['k']/2.-1./3., 2) \
                             + np.power(stat_out['ww']/stat_out['k']/2.-1./3., 2) + 2.*np.power(stat_out['uv']/stat_out['k']/2., 2) \
                             + 2.*np.power(stat_out['uw']/stat_out['k']/2., 2) + 2.*np.power(stat_out['vw']/stat_out['k']/2., 2))
    stat_out['d'] = np.sqrt( np.power(stat_out['e11']/stat_out['ek']/2.-1./3., 2) + np.power(stat_out['e22']/stat_out['ek']/2.-1./3., 2) \
                             + np.power(stat_out['e33']/stat_out['ek']/2.-1./3., 2) + 2.*np.power(stat_out['e12']/stat_out['ek']/2., 2) \
                             + 2.*np.power(stat_out['e13']/stat_out['ek']/2., 2) + 2.*np.power(stat_out['e23']/stat_out['ek']/2., 2))
    S = np.sqrt( stat['dUdx']*stat['dUdx'] + stat['dUdy']*stat['dUdy']  \
                 + stat['dVdx']*stat['dVdx'] + stat['dVdy']*stat['dVdy']  \
                 + stat['dWdx']*stat['dWdx'] + stat['dWdy']*stat['dWdy'] )
    stat_out['Sstar'] = S*(2.*stat_out['k'])/stat_out['ek']
    stat_out['Sc'] = S*np.sqrt(1./stat_out['ek']/Re)
    stat_out['eta'] = np.power(1./Re, 3./4.)*np.power(1./stat_out['ek'], 1./4.)

    ## temperature

    return stat_out


def write_stat_2d_vt(x,y,ep,stat,nx,ny, path, filename):
    statname = ['umean','vmean','wmean','pmean','tau11','tau22','tau33','tau12','tau13','tau23']
    #regularization
    ep[:,0] = 1.
    ep[:,-1] = 1.
    wall_coor = []
    for i in range(nx):
        for j in range(np.floor(ny/2).astype(int)):
            if ep[i,j]==1.:
                stat[i,j,:] = 0.
                if (i-1>=0 and ep[i-1,j]==0.) or (i+1<ep.shape[0] and ep[i+1,j]==0.) or \
                    (j-1>=0 and ep[i,j-1]==0.) or (j+1<ep.shape[1] and ep[i,j+1]==0.):
                    wall_coor.append([x[i,j],y[i,j]])
    for i in range(nx):
        for j in range(np.floor((ny+1)/2).astype(int), ny):
            if ep[i,j]==1.:
                stat[i,j,:] = 0.
                if (i-1>=0 and ep[i-1,j]==0.) or (i+1<ep.shape[0] and ep[i+1,j]==0.) or \
                    (j-1>=0 and ep[i,j-1]==0.) or (j+1<ep.shape[1] and ep[i,j+1]==0.):
                    wall_coor.append([x[i,j],y[i,j]])
    wall_coor = np.array(wall_coor)

    fid = open(path+filename, 'w')
    for i in range(nx):
        for j in range(ny):
            fid.write(x[i,j].__str__()+'\t')
            fid.write(y[i,j].__str__()+'\t')

            for idset in range(statname.__len__()):
                fid.write(stat[i,j,idset].__str__()+'\t')
            fid.write('\n')
    fid.close()

    fid_wx = open(path+'wall_x', 'w')
    fid_wy = open(path+'wall_y', 'w')
    for coor in wall_coor:
        fid_wx.write(coor[0].__str__())
        fid_wx.write('\n')
        fid_wy.write(coor[1].__str__())
        fid_wy.write('\n')
    fid_wx.close()
    fid_wy.close()

    return
