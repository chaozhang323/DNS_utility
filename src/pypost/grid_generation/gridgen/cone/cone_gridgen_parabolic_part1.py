from util import IO_util, util_geom, util_grid
import matplotlib.pyplot as plt
import numpy as np
from util import gridgen, interp2d


def normalized_flipped(u, nx):
	uscale = u[-1] - u[0]
	u -= u[0]
	u /= uscale
	u = np.flip(u, axis=0)
	u = 1. - u
	nx0 = u.shape[0]
	u = np.interp( np.linspace(0.,nx0-1.,nx), np.linspace(0.,nx0-1.,nx0), u)

	return u

def Billig_cone_shock(M, R, beta, y):
    delta = R*0.143*np.exp(3.24/M**2)
    Rc = R*1.143*np.exp(0.54/np.power(M-1.,1.2))
    tan_beta = np.tan(beta)
    x = R + delta - Rc/tan_beta**2*(np.sqrt(1+y**2*tan_beta**2/Rc**2) - 1.)
    return x

def Billig_cone_shock_inverse(M, R, beta, x):
    delta = R*0.143*np.exp(3.24/M**2)
    Rc = R*1.143*np.exp(0.54/np.power(M-1.,1.2))
    tan_beta = np.tan(beta)
    y = np.sqrt( ( ( (x-R-delta)/(-Rc/tan_beta**2)+1. )**2 - 1.)*Rc**2/tan_beta**2 )
    return y

def Billig_cone(R, theta, y):
    Inose = abs(y) < np.cos(theta)*R
    Ibody = np.logical_not( Inose )
    x = np.zeros_like(y)
    x[Inose] = np.sqrt(R**2-y[Inose]**2)
    x[Ibody] = -np.abs( y[Ibody]/np.tan(theta) ) + R/np.sin(theta)
    num_nose = np.count_nonzero(Inose)
    p = np.polyfit(y[num_nose/4*3:num_nose/4*5],x[num_nose/4*3:num_nose/4*5], 10)
    x[num_nose/4*3:num_nose/4*5] = np.polyval(p, y[num_nose/4*3:num_nose/4*5])

    return x

def decide_shock_by_wall(pw,vw, xs,ys):
    print interp2d.interp2d.__doc__
    inter,point = interp2d.interp2d.get_intersect_point_3p1v_line(pw,vw, xs,ys)
    return point

def decide_outlet(p0,v0, p1,v1, num):
    poly = util_geom.orthogonal_connect_spline(p0,v0,p1,v1)
    x = np.linspace(p0[0],p1[0],num)
    y = np.polyval(poly[::-1], x)
    return x, y

def decide_shock(M,R,beta, x0,y0, x_trans=0.):
    x1 = x0
    y1 = Billig_cone_shock_inverse(M,R,beta, x1-x_trans)
    y_top = np.linspace(y1, 0., 2**16)
    x_top  = Billig_cone_shock(M, R, beta, y_top) + x_trans
    p0 = np.array([x0,y0])
    v0 = np.array([np.sin(theta_c), np.cos(theta_c)])
    p1 = decide_shock_by_wall(p0, v0, x_top,y_top)
    x1 = p1[0]
    y1 = p1[1]
    y_top = np.linspace(y1, 0., 2**16)
    x_top  = Billig_cone_shock(M, R, beta, y_top) + x_trans
    return x_top, y_top, p1, x1,y1, p0,v0

def expand_shock(x_bot, y_bot, x_top, y_top, rate = 1.0):
    x_top[:] = (1.-rate)*x_bot + rate*x_top
    y_top[:] = (1.-rate)*y_bot + rate*y_top
    return x_top[0],y_top[0]

def get_vp_for_outlet(x,y):
    p = np.array( [x[-1],y[-1]] )
    angle = np.arctan2( y[-1]-y[-2], x[-1]-x[-2])
    angle -= 0.5*np.pi
    v = util_geom.unit_vector_from_angle_2d(angle)
    return p,v

if __name__ == '__main__':
    ###---parameters----###
    shape = (32,256)
    start_over = True
    # Files required to start interpolation
    dir_in = './InitBody/'
    filename_in_body = 'ConeSurfGeometry.dat'
    filename_in_shock = 'ConeShockGeometry.dat'
    filename_in_dleta = 'BL.dat'
    # Files generated by the program
    dir_out = '../../data/grids/'
    filename_out = 'cone_grid_parabolic_part1'
    ###---no parameters below---#
    # print gridgen.cm_gridgen.__doc__

    ## define wall
    surf = np.loadtxt(dir_in+filename_in_body, skiprows=2)
    x_right = surf[:,0]
    y_right = surf[:,1]

    ## define shock
    shock = np.loadtxt(dir_in+filename_in_shock, skiprows=2)
    x_left = shock[:,0] + shock[0,0]
    y_left = shock[:,1]

    ## define outlet
    p_l1,v_l1 = get_vp_for_outlet(x_left, y_left)
    p_r1,v_r1 = get_vp_for_outlet(x_right, y_right)
    x_top, y_top = decide_outlet(p_l1, v_l1, p_r1, v_r1, 2**16)

    ## define central line
    x_bot = np.linspace(x_left[0], x_right[0], 2**16)
    y_bot = np.zeros_like(x_bot)

    ## for a look
    plt.plot(x_bot,y_bot)
    plt.plot(x_top,y_top)
    plt.plot(x_left,y_left)
    plt.plot(x_right, y_right)
    plt.gca().set_aspect(1.)
    plt.show()

    ## define params for grid gen
    u = np.zeros(shape)
    v = np.zeros(shape)


    choice_left = np.stack([x_left,y_left], axis=1)

    choice_right = np.stack([x_right,y_right], axis=1)

    choice_bottom = np.stack([x_bot, y_bot], axis=1)

    choice_top = np.stack([x_top,y_top], axis=1)

    choice_left = gridgen.cm_gridgen.refine_boundary(choice_left, 2**16/choice_left.shape[0])
    choice_right = gridgen.cm_gridgen.refine_boundary(choice_right, 2**16/choice_right.shape[0])
 
   ## computation
    u,v = util_grid.gridgen_orth(shape, choice_left, choice_right, choice_bottom, choice_top, \
                                    dir_out+filename_out, varname=['x','z'], \
                                    tol_in=1.e-8, tol_bdry=0, nsave=30, \
                                    fixed_boundary=[0,0,0,0], bend_boundary=[0,0,1,0])
    exit()

    ## star over?
    if start_over:
        #IC
        idx_left, u[0,:],v[0,:] = gridgen.cm_gridgen.init_index(choice_left, shape[1], 1)
        idx_right, u[-1,:],v[-1,:] = gridgen.cm_gridgen.init_index(choice_right, shape[1], 1)
        idx_bottom, u[:,0],v[:,0] = gridgen.cm_gridgen.init_index(choice_bottom, shape[0], 1)
        idx_top, u[:,-1],v[:,-1] = gridgen.cm_gridgen.init_index(choice_top, shape[0], 1)
        u,v = gridgen.cm_gridgen.init_tfi(u, v)
        grid = {'x':u, 'y':v}
        IO_util.write_hdf5(dir_out+filename_out, grid)
    else:
        idx_left = np.load(dir_out+'idx_left.npy')
        idx_right = np.load(dir_out+'idx_right.npy')
        idx_bottom = np.load(dir_out+'idx_bottom.npy')
        idx_top = np.load(dir_out+'idx_top.npy')
        grid, vars_name = IO_util.read_hdf5(dir_out+filename_out+'.h5')
        u = grid['x']
        v = grid['y']

    ## computation
    iloop = True
    while iloop:

        u,v,idx_left,idx_right,idx_bottom,idx_top,iconverge = gridgen.cm_gridgen.compute_grid(u,v,idx_left,idx_right,
                                                                                              idx_bottom,idx_top,
                                                                                              choice_left,choice_right,
                                                                                              choice_bottom,choice_top, 1e-8, .0, 32,
                                                                                              np.array([0,0,0,0],dtype=int),
                                                                                              np.array([0,0,0,0],dtype=int))
        #iconverge = 1
        np.save(dir_out+'idx_left', idx_left)
        np.save(dir_out+'idx_right', idx_right)
        np.save(dir_out+'idx_bottom', idx_bottom)
        np.save(dir_out+'idx_top', idx_top)
        grid = {'x':u, 'y':v}
        IO_util.write_hdf5(dir_out+filename_out, grid)
        if iconverge==1:
            print 'Converged!'
            iloop=False

