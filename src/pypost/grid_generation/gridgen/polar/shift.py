import numpy as np

def Colonius_shift(z, dim):
    '''
        Shift z for polar treatment.
    '''
    z_work = np.swapaxes(z, 0, dim)
    
    ## check if z alligns with z = 0.
    if np.any(z_work[0,Ellipsis] != 0.):
        raise ValueError("z does not allign with z = 0.")


    L = z_work[-1:, Ellipsis]
    d = z_work[1:2, Ellipsis]
    r = L / (L+.5*d)

    a0 = L*(1.-r)
    #a0 = np.repeat(a0, z_work.shape[0], axis=0)
    z_work = a0 + r*z_work

    z_work = np.swapaxes(z_work, dim, 0)
    return z_work
