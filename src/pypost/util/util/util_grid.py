import numpy as np
from util import util_f 
from util import IO_util

class structured_base(object):
    def __init__(self, data, gridkeys):
        self.data = data
        self.gridkeys = gridkeys
        return

    @property
    def data(self):
        return self.__data
    @data.setter
    def data(self, data):
        if type(data) != dict:
            raise ValueError("Input is not a dictionary!")
        if any( [type(x) != np.ndarray for x in list(data.values())] ):
            raise ValueError("Input is not a dictionary of numpy array!")
        if len(set( [x.shape for x in list(data.values())])) != 1:
            raise ValueError('Variable shapes do not match!')
        
        self.__rank = list(data.values())[0].ndim
        self.__data = data
        
        return

    @property
    def gridkeys(self):
        return self.__gridkeys
    @gridkeys.setter
    def gridkeys(self, keys):
        if type(keys) is not list:
            raise ValueError("Input is not a list!")
        if len(keys)!=self.__rank:
            raise ValueError('Num of keys /= %d!'%self.__rank)
        if len(set(keys))!=self.__rank:
            raise ValueError('Repeated keys!')
        if any( [type(x)!=str for x in keys] ):
            raise ValueError('Input is not a list of strings!')
        if not set(keys)<set(self.__data.keys()):
            raise ValueError('Gridkeys not found in data!')

        self.__gridkeys = keys
        self.__datakeys = [key for key in list(self.__data.keys()) if (key not in self.__gridkeys) ]
        return

    @property
    def datakeys(self):
        return self.__datakeys
    
class structured_data(structured_base):
    def slice_data(self, dim, index):
        if dim>=self.__rank or dim<0:
            raise ValueError('Wrong dim input!')
        if type(index) != int:
            raise ValueError('Index is not an integer!')
        
        #if self.__rank==3:
        slice_data = structured_data()
        #else:
        #    slice_data = structured_data_2d()

        slice_data.gridkeys = [self.gridkeys[n] for n in range(self.__rank) if n!=dim]
        slice_data.data = {key: np.take(var, index, axis=dim) for key,var in self.data.items()}
        return slice_data
    
    def reposition(self):
        if self.__rank!=2:
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
        for key,var in self.__data.items():
            self.__data[key] = var.transpose()
        self.__gridkeys.reverse()
        return







        





def gridgen_orth(shape, choice_left, choice_right, choice_bottom, choice_top, \
                filename_out, varname=['x','z'], \
                tol_in=1.e-8, tol_bdry=0, nsave=64, fixed_boundary=[0,0,0,0], bend_boundary=[0,0,0,0]):
    u = np.zeros(shape)
    v = np.zeros(shape)
    
    ## init boundary
    idx_left, u[0,:],v[0,:] = util_f.cm_gridgen.init_index(choice_left, shape[1], 1)
    idx_right, u[-1,:],v[-1,:] = util_f.cm_gridgen.init_index(choice_right, shape[1], 1)
    idx_bottom, u[:,0],v[:,0] = util_f.cm_gridgen.init_index(choice_bottom, shape[0], 1)
    idx_top, u[:,-1],v[:,-1] = util_f.cm_gridgen.init_index(choice_top, shape[0], 1)
    
    ## init via tfi
    u,v = util_f.cm_gridgen.init_tfi(u, v)
    
    grid = {varname[0]:u, varname[1]:v}
    IO_util.write_hdf5(filename_out, grid)
    
    ## computation
    print('Gridgen orth starts computing...')
    iloop = True
    while iloop:
        u,v,idx_left,idx_right,idx_bottom,idx_top,iconverge = \
        util_f.cm_gridgen.compute_grid(u,v,idx_left,idx_right, \
                                            idx_bottom,idx_top, \
                                            choice_left,choice_right, \
                                            choice_bottom,choice_top, \
                                            tol_in, tol_bdry, nsave, \
                                            fixed_boundary,  \
                                            bend_boundary )
        ## out
        grid = {varname[0]:u, varname[1]:v}
        IO_util.write_hdf5(filename_out, grid)
 
        if iconverge==1:
            print('Converged!')
            iloop=False

    return u,v
