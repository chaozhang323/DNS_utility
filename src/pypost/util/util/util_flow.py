import numpy as np
from . import util_f, util_geom

def rortex_2d(dudx,dudy, dvdx,dvdy):
    if dudy.shape!=dvdx.shape:
        raise ValueError('Shape of inputs do not match!')
    rortex = util_f.flow.rortex_2d(dudx,dudy,dvdx,dvdy)
    
    return rortex

def get_mu_Sutherland(T, fluid_type='air'):
    if fluid_type=='air':
        c1 = 1.458e-6
        c2 = 110.4
        c3 = 0.
    elif fluid_type=='nitrogen':
        c1 = 1.418e-6
        c2 = 122.1
        c3 = 5.

    omega = .5
    tmp1 = 10.**(-c3/T)
    mu = c1*T**(1.+omega) / (T+c2*tmp1)
    return mu


def get_R(p,rho,T):
    R = p / (rho*T)
    R = np.mean(R)
    return R

def get_rho(p, T, fluid_type='air'):
    if fluid_type=='air':
        R = 287.0
    elif fluid_type=='nitrogen':
        R = 297.0
    rho = p / (R*T)
    return rho

def get_up_2d(u,v, x,y):
    p,a = util_geom.get_rightmost_point_and_angle(x,y)
    a -= np.pi*.5
    tau = util_geom.unit_vector_from_angle_2d(a)

    up = u*tau[0] + v*tau[1]
    return up

def get_delta(y, u):
    ''' compute boundary thickness'''
    percent = u / u[-1]
    idx = np.argwhere(percent>=.99)
    delta = y[idx[0,0,]]
    return delta

def get_delta_star(y, u, rho):
    percent = u / u[-1]
    idx = np.argwhere(percent>=.99)
    limit = idx[0,0]
    h = y[1:limit] - y[:limit-1]
    f = 1. - rho[:limit-1]*u[:limit-1]/(rho[limit]*u[limit])
    delta = np.dot(f, h)
    return delta

def get_theta(y, u, rho):
    percent = u / u[-1]
    idx = np.argwhere(percent>=.99)
    limit = idx[0,0]
    h = y[1:limit] - y[:limit-1]
    f = (1. - u[:limit-1]/u[limit])*rho[:limit-1]*u[limit]/(rho[limit]*u[limit])
    delta = np.dot(f, h)
    return delta

def Helmholtz_decomposition(u, v, w):
    '''
        u,v,w are required to be uniformly distributed.
    '''
    
    ## freq
    nx,ny,nz = u.shape
    kx = np.fft.fftfreq(nx).reshape(nx,1,1)
    ky = np.fft.fftfreq(ny).reshape(1,ny,1)
    kz = np.fft.fftfreq(nz).reshape(1,1,nz)
    k_sq = kx**2 + ky**2 + kz**2
    k_sq[k_sq==0.] = 1.

    ##
    u_f = np.fft.fftn(u)
    v_f = np.fft.fftn(v)
    w_f = np.fft.fftn(w)

    div_f = u_f*kx + v_f*ky + w_f*kz
    div_f_normed = div_f / k_sq
    
    u_dila_f = kx*div_f_normed
    v_dila_f = ky*div_f_normed
    w_dila_f = kz*div_f_normed

    u_sole_f = u_f - kx*div_f_normed
    v_sole_f = v_f - ky*div_f_normed
    w_sole_f = w_f - kz*div_f_normed

    ## div
    divv_f = 1j*np.pi*2.*(kx*u_sole_f + ky*v_sole_f + kz*w_sole_f)
    divv = np.real(np.fft.ifftn(divv_f))
    ## curl
    curl_x_f = 1j*np.pi*2.*ky*w_dila_f - 1j*np.pi*2.*kz*v_dila_f 
    curl_y_f = 1j*np.pi*2.*kz*u_dila_f - 1j*np.pi*2.*kx*w_dila_f 
    curl_z_f = 1j*np.pi*2.*kx*v_dila_f - 1j*np.pi*2.*ky*u_dila_f 
    curl_x = np.real(np.fft.ifftn(curl_x_f))
    curl_y = np.real(np.fft.ifftn(curl_y_f))
    curl_z = np.real(np.fft.ifftn(curl_z_f))
    
    ##
    u_sole = np.real(np.fft.ifftn(u_sole_f))
    v_sole = np.real(np.fft.ifftn(v_sole_f))
    w_sole = np.real(np.fft.ifftn(w_sole_f))

    u_dila = u - u_sole
    v_dila = v - v_sole
    w_dila = w - w_sole

    return u_sole,u_dila, v_sole,v_dila, w_sole,w_dila, divv, curl_x,curl_y,curl_z
