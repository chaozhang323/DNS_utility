import numpy as np
from util import util_f

def intersect_1p1v_line(p, v, x, y):
    inter, point = util_f.interp2d.get_intersect_point_1p1v_line(p, v, x, y)
    if inter==1:
        return point
    else:
        raise ValueError('Did not find intersect!')
        return

def unit_vector_from_angle_2d(angle):
    return np.stack( [np.cos(angle), np.sin(angle)], axis=-1 )

def get_rightmost_point_and_angle(x, y):
    p = np.stack( [x[-1],y[-1]], axis=-1 )
    a = np.arctan2(y[-1]-y[-2], x[-1]-x[-2])
    return p, a

def get_leftmost_point_and_angle(x, y):
    p = np.stack( [x[0],y[0]], axis=-1 )
    a = np.arctan2(y[1]-y[0], x[1]-x[0])
    return p, a

def connect_endpoints(p0,v0, p1,v1):
    a = np.zeros((4,4))
    b = np.zeros((4,1))
    
    a[0, 0] = p0[0]**3
    a[0, 1] = p0[0]**2
    a[0, 2] = p0[0]
    a[0, 3] = 1. 
    a[1, 0] = p1[0]**3
    a[1, 1] = p1[0]**2
    a[1, 2] = p1[0]
    a[1, 3] = 1. 
    a[2, 0] = 3.*p0[0]**2*v0[0]
    a[2, 1] = 2.*p0[0]*v0[0]
    a[2, 2] = v0[0]
    a[3, 0] = 3.*p1[0]**2*v1[0]
    a[3, 1] = 2.*p1[0]*v1[0]
    a[3, 2] = v1[0]


    

    #a[0:2,0] = 1.
    #a[0,1] = p0[0]
    #a[1,1] = p1[0]
    #a[0,2] = p0[0]**2
    #a[1,2] = p1[0]**2
    #a[0,3] = p0[0]**3
    #a[1,3] = p1[0]**3
    #
    #a[2,1] = v0[0]
    #a[3,1] = v1[0]
    #a[2,2] = 2.*p0[0]*v0[0]
    #a[3,2] = 2.*p1[0]*v1[0]
    #a[2,3] = 3.*p0[0]**2*v0[0]
    #a[3,3] = 3.*p1[0]**2*v1[0]
    # 
    b[0] = p0[1]
    b[1] = p1[1]
    b[2] = v0[1]
    b[3] = v1[1]
    # print a
    # print b
    x = np.linalg.solve(a,b)
    return x

def make_poly_from_endpoints(f, x):
    dof = f.shape[0]

    ## make poly power/coefs ##
    coef = [ np.ones((dof,)) ]
    power = [ np.arange(dof-1, -1, -1) ]
    for n in range(dof//2-1):
        coef.append( coef[n][:-1]*power[n][:-1] )
        power.append( power[n][:-1] - 1 )
    
    ## compute for polys
    A = np.zeros((dof,dof))
    A[0, :] = coef[0]*(x[0]**power[0])
    A[1, :] = coef[0]*(x[1]**power[0])
    for n in range(1, dof//2):
        A[2*n,:-n] = coef[n]*(x[0]**power[n])
        A[2*n+1,:-n] = coef[n]*(x[1]**power[n])
    ##
    p = np.linalg.solve(A, f)
    return p

def reorder(data, ordered_key):
    dim = data[ordered_key].ndim
    assert (dim==1), "Dimension not equal to 1."
    I = np.argsort(data[ordered_key])

    for key,var in data.items():
        data[key] = var[I]
    
    return data
