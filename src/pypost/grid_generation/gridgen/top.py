from util import IO_util, util_geom, util_grid
import matplotlib.pyplot as plt
import numpy as np


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
    print(interp2d.interp2d.__doc__)
    inter,point = interp2d.interp2d.get_intersect_point_3p1v_line(pw,vw, xs,ys)
    return point

def decide_top(p0,v0, p1,v1, num):
    poly = util_geom.connect_endpoints(p0,v0,p1,v1)
    x = np.linspace(p0[0],p1[0],num)
    y = np.polyval(poly, x)
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

def get_pv_for_top(x,y):
    '''
        Input is from left/right boundary points.
    '''

    p, angle = util_geom.get_rightmost_point_and_angle(x, y)
    p = np.array( [x[-1],y[-1]] )
    angle = np.arctan2( y[-1]-y[-2], x[-1]-x[-2])
    angle -= 0.5*np.pi
    v = util_geom.unit_vector_from_angle_2d(angle)
    return p,v


