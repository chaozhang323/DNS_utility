import numpy as np
import scipy.ndimage.morphology as spm

def find_wall_vars_2d(x, y, ep, u, v):
    ep = ep.astype(int)
    nx, ny = ep.shape
    e1 = [[1,1], [1,0]]
    e2 = [[1,0], [1,1]]
    e3 = [[1,0], [1,0], [1,0]]
    xw = list()
    yw = list()
    dudnw = list()
    dw = list()

    find_e1 = spm.binary_hit_or_miss(ep, structure1=e1, origin1=[0,0])
    idx_e1 = np.argwhere(find_e1)
    for idx in idx_e1:
        i, j = idx
        dx = x[i,j]-x[i-1,j]
        dy = y[i,j]-y[i,j-1]
        dudn = u[i,j]/dy - v[i,j]/dx
        dudnw.append(dudn)
        dw.append((dx+dy))
        xw.append(x[i,j])
        yw.append(y[i,j])
    find_e2 = spm.binary_hit_or_miss(ep, structure1=e2, origin1=[-1,0])
    idx_e2 = np.argwhere(find_e2)
    for idx in idx_e2:
        i, j = idx
        dx = x[i+1,j]-x[i,j]
        dy = y[i,j]-y[i,j-1]
        dudn = u[i,j]/dy - v[i,j]/dx
        dudnw.append(dudn)
        dw.append((dx+dy))
        xw.append(x[i,j])
        yw.append(y[i,j])
    find_e3 = spm.binary_hit_or_miss(ep, structure1=e3, origin1=[0,0])
    idx_e3 = np.argwhere(find_e3)
    for idx in idx_e3:
        i, j = idx
        dy = y[i,j]-y[i,j-1]
        dudn = u[i,j]/dy
        dudnw.append(dudn)
        dw.append((dy))
        xw.append(x[i,j])
        yw.append(y[i,j-1])

    sorted_index = np.argsort(xw)
    xw = np.array(xw)[sorted_index]
    yw = np.array(yw)[sorted_index]
    dudnw = np.array(dudnw)[sorted_index]
    dw = np.array(dw)[sorted_index]


    return xw, yw, dudnw, dw

def find_wall_2d(x,y, ep_in):
    ep = ep_in.astype(int)
    e1 = [[1,1], [1,0]]
    e2 = [[1,0], [1,1]]
    e3 = [[1,0], [1,0], [1,0]]

    wu_list = []

    find_e1 = spm.binary_hit_or_miss(ep, structure1=e1, origin1=[0,0])
    idx_e1 = np.argwhere(find_e1)
    for idx in idx_e1:
        idx = tuple(idx)
        wu = {'index':idx, 'neighbors':(-1./(x[idx] - x[(idx[0]-1,idx[1])]), 0, 1./(y[idx] - y[(idx[0], idx[1]-1)]), 0)}
        wu_list.append(wu)

    find_e2 = spm.binary_hit_or_miss(ep, structure1=e2, origin1=[-1,0])
    idx_e2 = np.argwhere(find_e2)
    for idx in idx_e2:
        idx = tuple(idx)
        wu = {'index':idx, 'neighbors':(0., 1./(x[(idx[0]+1,idx[1])] - x[idx]), 0, 1./(y[idx] - y[(idx[0], idx[1]-1)]), 0)}
        wu_list.append(wu)

    find_e3 = spm.binary_hit_or_miss(ep, structure1=e3, origin1=[0,0])
    idx_e3 = np.argwhere(find_e3)
    for idx in idx_e3:
        idx = tuple(idx)
        wu = {'index':idx, 'neighbors':(0., 0, 1./(y[idx] - y[(idx[0], idx[1]-1)]), 0)}
        wu_list.append(wu)

    return wu_list

