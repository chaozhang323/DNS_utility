import numpy as np
from numpy.core.multiarray import normalize_axis_index
from numpy.lib.index_tricks import ndindex
from util import util_f, util_interp

def examine_equal(data, axis):
    tmp = np.swapaxes(data, 0, axis)
    equals = [tmp[(0,Ellipsis)]==tmp[(n,Ellipsis)] for n in range(tmp.shape[0])]
    if np.all(equals):
        print('All equal along axis=%d .'%(axis))
    else:
        print('Warning: Not all equal along axis=%d !'%(axis))
    return

def appy_long_axis_2(a1, a2):
    nd = f.ndim
    axis = normalize_axis_index(axis, nd)
    if axis!=0:
        f = np.swapaxes(f, axis, 0)
        x = np.swapaxes(x, axis, 0)
    
    dfdx = np.empty_like(f)
    inds = ((i,Ellipsis) for i in range(f.shape[0]))
    for ind in inds:
        dfdx[ind] = diff(f[ind], x[ind], **kwargs)

    if axis!=0:
        dfdx = np.swapaxes(dfdx, axis, 0)
    


def diff_1st(f, x, axis=0, **kwargs):
    """ Generally 1st order local finite difference.
    Central for interior points, and one-sided for boundary points."""

    def diff(f, x, bc=0):
        dfdx = np.zeros_like(f)
        dfdx[1:-1] = (f[2:]-f[:-2])/(x[2:]-x[:-2])
        if bc==0:
            dfdx[0] = (f[1]-f[-1])/(2.*(x[1]-x[0]))
            dfdx[-1] = (f[0]-f[-2])/(2.*(x[-1]-x[-2]))
        elif bc==1:
            dfdx[0] = (f[1]-f[0])/(x[1]-x[0])
            dfdx[-1] = (f[-1]-f[-2])/(x[-1]-x[-2])
        return dfdx
    def wrapper_diff(f, **kwargs):
        nf = f.shape[0]/2
        x = f[nf:]
        f = f[:nf]
        dfdx = diff(f, x, **kwargs)
        return dfdx

    f_work = np.concatenate((f,x),axis=axis)
    dfdx = np.apply_along_axis(wrapper_diff, axis, f_work, **kwargs)

    
    
    return dfdx

class structured_data:
    def __init__(self, data, gridkeys):
        self.data = data
        self.gridkeys = gridkeys
        return

    @property
    def data(self):
        return self._data
    @data.setter
    def data(self, data):
        if type(data) is not dict:
            raise TypeError("Input is not a dictionary!")
        if any( [type(x) is not np.ndarray for x in list(data.values())] ):
            raise TypeError("Input is not a dictionary of numpy array!")
        if len(set( [x.shape for x in list(data.values())])) != 1:
            raise ValueError('Variable shapes do not match!')
        
        self._rank = list(data.values())[0].ndim
        self._shape = list(data.values())[0].shape
        
        self._data = data

        return

    @property
    def gridkeys(self):
        return self._gridkeys
    @gridkeys.setter
    def gridkeys(self, keys):
        if type(keys) is not list:
            raise TypeError("Input is not a list!")
        if len(keys)!=self.rank:
            raise ValueError('Num of keys /= rank = %d!'%self.rank)
        if len(set(keys))!=self.rank:
            raise ValueError('Repeated keys!')
        if any( [type(x) is not str for x in keys] ):
            raise TypeError('Input is not a list of strings!')
        if not set(keys)<=set(self.data.keys()):
            raise ValueError('Gridkeys not found in data!')

        self._gridkeys = keys
        self._datakeys = [key for key in list(self.data.keys()) if (key not in self.gridkeys) ]
        return

    @property
    def datakeys(self):
        return self._datakeys

    @property
    def rank(self):
        return self._rank
    
    @property
    def shape(self):
        return self._shape

    def slice_data(self, dim=0, index=0):
        if dim>=self.rank or dim<0:
            raise ValueError('Wrong dim input!')
        if type(index) is not int:
            raise ValueError('Index is not an integer!')
        
        gridkeys = [self.gridkeys[n] for n in range(self.rank) if n!=dim]
        data = {key: np.take(var, index, axis=dim) for key,var in self.data.items()}
 
        sdata = structured_data(data, gridkeys)

        return sdata

    def reposition(self):
        if self.rank != 2:
            print('This method is intended only for 2D data. Nothing will be done.')
            return

        xkey,ykey = self.gridkeys

        x = self.data[xkey]
        y = self.data[ykey]
        allkeys = list(self.data.keys())

        if x[0,0] > x[-1,-1]:
            for key in allkeys:
                self.data[key] = self.data[key][::-1,:]
            print('Flipped %s.'%xkey)
        if y[0,0] > y[-1,-1]:
            for key in allkeys:
                self.data[key] = self.data[key][:,::-1]
            print('Flipped %s.'%ykey)
        if x[-1,0] < x[0,-1]:
            transpose(self)
            #for key in allkeys:
            #    self.data[key] = self.data[key].transpose()
            print('Transposed.')
        return

    def transpose(self):
        for key,var in self.data.items():
            self.data[key] = var.transpose()
        self.gridkeys.reverse()
        return
    
    def get_wall_normal_profiles(self, wall_type, m_range, nmax, method='linear'):
        if self.rank!=2:
            raise ValueError('This only works for rank = 2!')
        f_wall = {'bot':util_f.interp2d.get_wall_normal_bot, \
                  'top':util_f.interp2d.get_wall_normal_top, \
                  'right':util_f.interp2d.get_wall_normal_right}
        x_out,z_out,wd_out = f_wall[wall_type](self.data[self.gridkeys[0]], \
                                self.data[self.gridkeys[1]], \
                                m_range, nmax)
        grid_out = {self.gridkeys[0]:x_out, self.gridkeys[1]:z_out}
        grid_in = {self.gridkeys[0]:self.data[self.gridkeys[0]], \
                   self.gridkeys[1]:self.data[self.gridkeys[1]]}
        data_in = {key:self.data[key] for key in self.datakeys}
        data_out = util_interp.structured_interp(grid_in, data_in, grid_out, robust=True, method=method)
        data_out.update(grid_out)
        data_out.update({'wd':wd_out})
        return data_out
        
class profile(structured_data):
    def __init__(self, data, key):
        super(profile,self).__init__(data, key)
        print(self.data)
        if self.rank!=1:
            raise ValueError('Rank /= 1!')
        return



def recover_structured_data(data_in, nodemap=None, standard=True, \
                            grid_keys=['x','z'] ):
    if nodemap is None:
        nodemap = data_in.pop('nodemap')
    data_out = {}
    x_key = grid_keys[0]
    y_key = grid_keys[1]
    def gen_data_keys(keys):
        for key in keys:
            if key!=x_key and key!=y_key:
                yield key
    data_keys = gen_data_keys(data_in.keys())

    ## calculate element center
    point_coor = np.stack([data_in[x_key], \
                           data_in[y_key]], axis=1)
    grid = util_f.uns2str.get_element_center(nodemap, point_coor)
    
    ## standard restructure
    if standard:
        elemap, num_neighbors = util_f.uns2str.get_rough_elemap(nodemap)
        mask_bdry, num_bdry = util_f.uns2str.get_num_boundary_elements(num_neighbors)
        index_bdry, num_neighbors_bdry = util_f.uns2str.get_boundary_elements(elemap, num_neighbors, mask_bdry, num_bdry)
        shape = util_f.uns2str.get_domain_shape(num_neighbors_bdry)
        grid_map = util_f.uns2str.rebuild_structured_grid(elemap, num_neighbors, index_bdry, num_neighbors_bdry, shape[0], shape[1])
        print('New shape = ', shape)
        ## coordinate rebuilding
        data_out[x_key] = util_f.uns2str.get_structured_data(grid_map, grid[:,0])
        data_out[y_key] = util_f.uns2str.get_structured_data(grid_map, grid[:,1])

        ## decide whether to reposition
        if data_out[x_key][0,0] > data_out[x_key][-1,-1]:
            print('Flipped I.')
            for key in grid_keys:
                data_out[key] = data_out[key][::-1,:]
            grid_map = grid_map[::-1,:]
        if data_out[y_key][0,0] > data_out[y_key][-1,-1]:
            print('Flipped J.')
            for key in grid_keys:
                data_out[key] = data_out[key][:,::-1]
            grid_map = grid_map[:,::-1]
        #if data_out[y_key][-1,0] > data_out[y_key][0,-1]:
        #    print('Transposed.')
        #    for key in grid_keys:
        #        data_out[key] = data_out[key].transpose()
        #    grid_map = grid_map.transpose()
        if data_out[x_key][-1,0] < data_out[x_key][0,-1]:
            print('Transposed.')
            for key in grid_keys:
                data_out[key] = np.transpose(data_out[key]) #data_out[key][::-1,:]
            grid_map = np.transpose(grid_map) #grid_map[::-1,:]
        
        #
        #if x[0,0] > x[-1,-1]:
        #    for key in allkeys:
        #        self.data[key] = self.data[key][::-1,:]
        #    print('Flipped %s.'%xkey)
        #if y[0,0] > y[-1,-1]:
        #    for key in allkeys:
        #        self.data[key] = self.data[key][:,::-1]
        #    print('Flipped %s.'%ykey)
        #if x[-1,0] < x[0,-1]:
        #    transpose(self)
        # 
        ## interpolate variables
        for key in data_keys:
            if key!=x_key and key!=y_key:
                data_out[key] = util_f.uns2str.get_structured_data(grid_map, data_in[key])
    
    ## non-standard: only coordinate
    else:
        data_out = data_in.copy()
        data_out[x_key] = grid[:,0]
        data_out[y_key] = grid[:,1]
        #for key in data_keys:
        #    data_out[key] = data_in[key]

    return data_out


def change_dict(dict_in, keys_change=None):
    dict_out = dict_in
    if keys_change==None:
        pass
    else:
        for key_in, key_out in keys_change.items():
            try:
                dict_out[key_out] = dict_in[key_in]
                dict_out.pop(key_in)
            except:
                print('Failed to convert %s to %s !'%(key_in, key_out))
                pass
    try:
        dict_out['nodemap'] = dict_in['nodemap']
    except:
        print('Nodemap is not passed to new dict!')
        pass
    return dict_out

def reorder(data, ordered_key):
    dim = data[ordered_key].ndim
    assert (dim==1), "Dimension not equal to 1."
    I = np.argsort(data[ordered_key])

    for key,var in data.items():
        data[key] = var[I]
    
    return data
