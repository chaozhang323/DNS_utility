import numpy as np

def PSD_simple(u, *args, **kwargs):
    u_f = np.fft.fftn(u, **kwargs)
    u_f_mag = np.absolute(u_f)
    return u_f_mag
