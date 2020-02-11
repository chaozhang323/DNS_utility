from util import util_interp
import numpy as np

def add_time(data):
    data['time'] = np.array([0.])
    return

def ICgen(grid_old, data_old, grid_new, robust=True, method='linear'):
    data_new = util_interp.structured_interp(grid_old, data_old, grid_new, robust=True, method='linear')
    #data_new.update(grid_new)

    return data_new


