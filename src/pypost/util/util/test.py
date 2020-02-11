from . import util_data as f
import numpy as np

d = {'x':np.zeros((2))}
k = ['x']

a = f.profile(d,k)
