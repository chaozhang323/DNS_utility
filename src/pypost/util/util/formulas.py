import numpy as np

def tbl_scalar_Kader1981(y, yplus, Pr, type='channel'):
    beta = np.power(3.85*np.power(Pr,1./3.) - 1.3, 2) + 2.12*np.log(Pr)
    Gamma = (np.power(10.,-2)*np.power(Pr*yplus,4)) / (1.+5.*np.power(Pr,3)*yplus)
    if type=='channel':
        a1 = Pr*yplus*np.exp(-Gamma)
        b1 = 1.5*(2.-y) / (1.+2.*(1.-y))
        tplus = a1 + (2.12*np.log((1.+yplus)*b1)+beta)*np.exp(-1./Gamma)

    return tplus

def tbl_T_Kim1989():
    data = np.loadtxt('/usr/local/home/yl8bc/Nic/BFS/T.dat')

    return data

def tbl_Trms_Kim1989():
    data = np.loadtxt('/usr/local/home/yl8bc/Nic/BFS/Trms.dat')

    return data

def tbl_u(yplus, c):
    uplus = np.zeros_like(yplus)
    for i in range(yplus.shape[0]):
        if yplus[i] < 8.:
            uplus[i] = yplus[i]
        else:
            uplus[i] = 1./0.41*np.log(yplus[i]) + c
    return uplus