from util import IO_util, util_geom, util_grid
import matplotlib.pyplot as plt
import numpy as np
from .. import top


def gridgen(nx, xl,yl, u_in,v_in, filename_out, \
            num_trans1=3,num_trans2=3,num_trans3=8, check=False):
    trans0 = nx
    trans3 = nx - num_trans3 
    trans2 = nx + num_trans2 - 1
    trans1 = trans2 + num_trans1

    ## make transition grid
    u_trans = np.apply_along_axis(extend_left, 0, u_in, num=num_trans1+num_trans2)
    v_trans = np.apply_along_axis(extend_left, 0, v_in, num=num_trans1+num_trans2)

    ## define shock
    y_right = v_trans[0,:]
    x_right = u_trans[0,:]
	
    x_right = u_in[0,:]
    y_right = v_in[0,:]

    ## define inlet
    x_left = xl*np.ones((2**16,))
    y_left = np.linspace(0.,yl, x_left.shape[0])

    ## define outlet
    p_l1,v_l1 = top.get_pv_for_top(x_left, y_left)
    p_r1,v_r1 = top.get_pv_for_top(x_right, y_right)
    x_top, y_top = top.decide_top(p_l1, v_l1, p_r1, v_r1, 2**16)

    ## define central line
    x_bot = np.linspace(x_left[0], x_right[0], 2**16)
    y_bot = np.zeros_like(x_bot) 
    
    if check:
        plt.plot(x_bot,y_bot)
        plt.plot(x_top,y_top)
        plt.plot(x_left,y_left)
        plt.plot(x_right,y_right)
        plt.gca().set_aspect(1.)
        plt.show()

    ## define params for grid gen
    shape = (nx, u_in.shape[1])
    u = np.zeros(shape)
    v = np.zeros(shape)


    choice_left = np.stack([x_left,y_left], axis=1)
    choice_right = np.stack([x_right,y_right], axis=1)
    choice_bottom = np.stack([x_bot, y_bot], axis=1)
    choice_top = np.stack([x_top,y_top], axis=1)

    ## computation
    u,v = util_grid.gridgen_orth(shape, choice_left, choice_right, choice_bottom, choice_top, \
                                    filename_out, varname=['x','z'], \
                                    tol_in=1.e-8, tol_bdry=0, nsave=64, \
                                    fixed_boundary=[0,1,0,0], bend_boundary=[0,0,1,0])
    #u1 = np.concatenate( ( u, u_trans[1:,:], u_in), axis=0 )
    #v1 = np.concatenate( ( v, v_trans[1:,:], v_in), axis=0 )
    #smooth_with_poly(u1, v1, 5, trans1, trans2, trans3)


    return u, v
def extend_left(u, num=4):
    du = u[1] - u[0]
    u_ext = u[0] + np.linspace(-num, -1., num)*du
    return u_ext

def smooth_with_poly(u,v, deg, n1, n2, n3):
    idx13 = np.arange(n3, n1, 1.)#[np.newaxis,:]
    def get_f(u, n1, n3):
        f = np.zeros(6)
        f[0] = u[n1]
        f[1] = u[n3]
        f[2] = u[n1] - u[n1-1]
        f[3] = u[n3+1] - u[n3]
        f[4] = u[n1-2] -2.*u[n1-1] + u[n1]
        f[5] = u[n3+2] -2.*u[n3+1] + u[n3]
        return f

    def smooth_1d(u, v,n1=1,n3=3):
        x = [n3, n1]
        ##
        f = get_f(u, n3, n1)
        p = util_geom.make_poly_from_endpoints(x, f)
        u[n3:n1] = np.polyval(p, idx13)
        f = get_f(v, n3, n1)
        p = util_geom.make_poly_from_endpoints(x, f)
        v[n3:n1] = np.polyval(p, idx13)
        return

    for n in range(u.shape[1]):
        i=1
        smooth_1d(u[:,n], v[:,n], n1=n1,n3=n3)

    return 

