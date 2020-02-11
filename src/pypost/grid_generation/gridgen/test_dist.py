import numpy as np

def f_tanh(a, t, x):
    f = (np.tanh(a*(x-t)) - np.tanh(a*(-t))) / (np.tanh(a*(1-t)) - np.tanh(a*(-t)))
    return f

def dfdx_tanh(a, t, x):
    dfdx = a/(np.tanh(a*(1.-t))+np.tanh(a*t)) * (1.-np.tanh(a*(x-t))**2)
    return dfdx

def dfdx_3rd(a, t, x):
    dfdx = 3.*a*x**2 - 2.*a*(1.+t)*x + (1+a*t)
    return dfdx

if __name__=='__main__':
    t = np.linspace(0., 1., 32)[:,np.newaxis]
    x = np.linspace(0., 1., 32)[np.newaxis,:]
    a = 1.

    f = f_tanh(a, t, x)
    dfdx = dfdx_tanh(a, t, x)
    dfdx = dfdx_3rd(a, t, x)
