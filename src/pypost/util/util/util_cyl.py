import numpy as np

def coord_cart_to_cyl(x,y,z, pole=np.array([0.,0.])):
    theta = np.arctan2(y-pole[0], z-pole[1])
    r = np.sqrt((y-pole[0])**2 + (z-pole[1])**2)
    x = x.copy()
    return x, theta, r

def vel_cart_to_cyl(u,v,w, x,theta,r):
    ux = u.copy()
    ut = -v*np.sin(theta) + w*np.cos(theta)
    ur = v*np.cos(theta) + w*np.sin(theta)
    return ux, ut, ur

def coord_cyl_to_cart(x, theta, r):
    x = x.copy()
    z = r*np.cos(theta)
    y = r*np.sin(theta)
    return x,y,z

def vel_cyl_to_cart(ux,ut,ur, x,theta,r):
    u = ux.copy()
    v = ur*np.cos(theta) - ut*np.sin(theta)
    w = ur*np.sin(theta) + ut*np.cos(theta)
    return u, v, w


