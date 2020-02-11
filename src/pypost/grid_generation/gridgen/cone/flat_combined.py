import numpy as np
from util import util_geom

def gridgen(u0,v0, u1,v1, n0,n1):
    u = np.concatenate((u0[:-1,:],u1), axis=0)
    v = np.concatenate((v0[:-1,:],v1), axis=0)

    smooth_with_poly(u,v, n0,n1)
    
    return u, v

def smooth_with_poly(u,v, n0, n1):
    x = np.arange(n0, n1+1, 1.)
    x_ends = [n0, n1]


    def get_f(u, n0, n1):
        f = np.zeros(6)
        f[0] = u[n0]
        f[1] = u[n1]
        f[2] = u[n0] - u[n0-1]
        f[3] = u[n1+1] - u[n1]
        f[4] = u[n0-2] -2.*u[n0-1] + u[n0]
        f[5] = u[n1+2] -2.*u[n1+1] + u[n1]
        return f

    def smooth_1d(u, v):
        # u
        f_ends = get_f(u, n0, n1)
        p = util_geom.make_poly_from_endpoints(f_ends, x_ends)
        u[n0:n1+1] = np.polyval(p, x)
        # v
        f_ends = get_f(v, n0, n1)
        p = util_geom.make_poly_from_endpoints(f_ends, x_ends)
        v[n0:n1+1] = np.polyval(p, x)
        return


    for n in range(u.shape[1]):
        smooth_1d(u[:,n], v[:,n])

    return 

