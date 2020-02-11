import sys
from os import path
sys.path.append( path.dirname( path.dirname( path.abspath(__file__) ) ) )
from util import IO_util, util_plot
import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as signal
import scipy.integrate as intergrate
from numba import jit
import itertools as it
from numpy import f2py
with open("../grid_generation/gridgen.f90") as sourcefile:
    sourcecode = sourcefile.read()
f2py.compile(sourcecode, modulename='gridgen', extension='.f90')
import gridgen


@jit
def gaussian_DxNy(u_in, error, report=64):
    err = 1.
    n = 0
    u = u_in
    u0 = u
    one_over_three = 1./3.
    while err>error:
        u[1:-1,1:-1] = 0.25*(u0[:-2,1:-1] + u0[2:,1:-1] + u0[1:-1,:-2] + u0[1:-1,2:])
        u[1:-1,0] = one_over_three*(u0[:-2,0] + u0[2:,0] + u0[1:-1,1])
        u[1:-1,-1] = one_over_three*(u0[:-2,-1] + u0[2:,-1] + u0[1:-1,-2])
        u0 = u
        err1 = np.mean(np.abs(0.25*(u[:-2,1:-1] + u[2:,1:-1] + u[1:-1,:-2] + u[1:-1,2:]) - u[1:-1,1:-1]))
        err2 = np.mean(np.abs(one_over_three*(u[:-2,0] + u[2:,0] + u[1:-1,1]) - u[1:-1,0]))
        err3 = np.mean(np.abs(one_over_three*(u[:-2,-1] + u[2:,-1] + u[1:-1,-2]) - u[1:-1,-1]))
        err = one_over_three*np.mean(err1+err2+err3)
        n += 1
        if divmod(n,report):
            print err
    return u

@jit
def gaussian_DyNx(u_in, error, report=64):
    err = 1.
    n = 0
    u = u_in
    u0 = u
    one_over_three = 1./3.
    while err>error:
        u[1:-1,1:-1] = 0.25*(u0[:-2,1:-1] + u0[2:,1:-1] + u0[1:-1,:-2] + u0[1:-1,2:])
        u[0,1:-1] = one_over_three*(u0[0, :-2] + u0[0, 2:] + u0[1, 1:-1])
        u[-1, 1:-1] = one_over_three*(u0[-1, :-2] + u0[-1, 2:] + u0[-1, 1:-1])
        u0 = u
        err1 = np.mean(np.abs(0.25*(u[:-2,1:-1] + u[2:,1:-1] + u[1:-1,:-2] + u[1:-1,2:]) - u[1:-1,1:-1]))
        err2 = np.mean(np.abs(one_over_three*(u[0,:-2] + u[0, 2:] + u[1, 1:-1]) - u[0, 1:-1]))
        err3 = np.mean(np.abs(one_over_three*(u[-1, :-2] + u[-1, 2:] + u[-1, 1:-1]) - u[-1, 1:-1]))
        err = one_over_three*(err1+err2+err3)
        n += 1
        if divmod(n,report):
            print err
    return u

@jit
def gaussian_Dxy(u_in, error, report=64):
    err = 1.
    n = 0
    u = u_in
    u0 = u
    while err>error:
        u[1:-1,1:-1] = 0.25*(u0[:-2,1:-1] + u0[2:,1:-1] + u0[1:-1,:-2] + u0[1:-1,2:])

        u0 = u
        err = np.mean(np.abs(0.25*(u[:-2,1:-1] + u[2:,1:-1] + u[1:-1,:-2] + u[1:-1,2:]) - u[1:-1,1:-1]))
        n += 1
        if divmod(n,report):
            print err
    return u

def laplace_Dxy(u, v, error=1.e-6):
    u = gaussian_Dxy(u, error)
    v = gaussian_Dxy(v, error)
    return u,v

def inverse_Laplace(u, v, dx,dy, error = 1.e-6):
    shape = u.shape
    u1 = u
    v1 = v
    err = 1.
    while err>error:
        for i in range(1, shape[0]-1):
            for j in range(1, shape[1]-1):
                dudx = (u[i+1,j]-u[i-1,j]) / dx
                dudy = (u[i,j+1]-u[i,j-1]) / dy
                dvdx = (v[i+1,j]-v[i-1,j]) / dx
                dvdy = (v[i,j+1]-v[i,j-1]) / dy
                alpha = (dudy**2 + dvdy**2) / dx**2
                beta = 0.#(dudx*dudy + dvdx*dvdy) /dx/dy
                gamma = (dudx**2 + dvdx**2) / dy**2
                rhs = alpha*(u[i+1,j]+u[i-1,j])  + gamma*(u[i,j+1]+u[i,j-1]) - 2.*beta*(u[i+1,j+1]+u[i-1,j-1]-u[i+1,j-1]-u[i-1,j+1])
                u1[i,j] = rhs / (alpha*2. + gamma*2.)
                rhs = alpha*(v[i+1,j]+v[i-1,j])  + gamma*(v[i,j+1]+v[i,j-1]) - 2.*beta*(v[i+1,j+1]+v[i-1,j-1]-v[i+1,j-1]-v[i-1,j+1])
                v1[i,j] = rhs / (alpha*2. + gamma*2.)
        u = u1
        v = v1
        err = 0.
        for i in range(1, shape[0]-1):
            for j in range(1, shape[1]-1):
                dudx = (u[i+1,j]-u[i-1,j]) / dx
                dudy = (u[i,j+1]-u[i,j-1]) / dy
                dvdx = (v[i+1,j]-v[i-1,j]) / dx
                dvdy = (v[i,j+1]-v[i,j-1]) / dy
                alpha = (dudy**2 + dvdy**2) / dx**2
                beta = 0.#(dudx*dudy + dvdx*dvdy) /dx/dy
                gamma = (dudx**2 + dvdx**2) / dy**2
                rhs = alpha*(u[i+1,j]+u[i-1,j]) + gamma*(u[i,j+1]+u[i,j-1]) - 2.*beta*(u[i+1,j+1]+u[i-1,j-1]-u[i+1,j-1]-u[i-1,j+1])
                err += abs(rhs - (alpha*2. + gamma*2.)*u[i,j])
                rhs = alpha*(v[i+1,j]+v[i-1,j]) + gamma*(v[i,j+1]+v[i,j-1]) - 2.*beta*(v[i+1,j+1]+v[i-1,j-1]-v[i+1,j-1]-v[i-1,j+1])
                err += abs(rhs - (alpha*2. + gamma*2.)*v[i,j])
        err /= shape[0]*shape[1]
        print err,u[1,1]



if __name__ == '__main__':
    print gridgen.cm_gridgen.__doc__
    shape = (128, 32)
    x = np.linspace(0., 1., shape[0])
    y = np.linspace(0., 1., shape[1])

    u = np.zeros(shape)
    v = np.zeros(shape)

    # num_choice = 1280
    # choice_left = np.zeros((num_choice,2))
    # choice_left[:,0] = -1.
    # choice_left[:,1] = np.linspace(0., 4., num_choice)
    #
    # num_choice = 1280
    # choice_right = np.zeros((num_choice,2))
    # choice_right[:,0] = np.linspace(0., 2., num_choice)
    # choice_right[:,1] = np.linspace(0., 2., num_choice)**2
    #
    # num_choice = 640
    # choice_bottom = np.zeros((num_choice,2))
    # choice_bottom[:,0] = np.linspace(-1., 0., num_choice)
    # choice_bottom[:,1] = np.zeros((num_choice,))
    #
    # num_choice = 640
    # choice_top = np.zeros((num_choice,2))
    # choice_top[:,0] = np.linspace(-1., 2., num_choice)
    # choice_top[:,1] = 4.*np.ones((num_choice,))

    num_choice = 640
    choice_left = np.zeros((num_choice,2))
    choice_left[:,0] = 0.
    choice_left[:,1] = np.linspace(0.2, 2., num_choice)

    num_choice = 640
    choice_right = np.zeros((num_choice,2))
    choice_right[:,0] = 0.
    choice_right[:,1] = np.linspace(-0.2, -2., num_choice)

    num_choice = 640
    choice_bottom = np.zeros((num_choice,2))
    choice_bottom[:,0] = -25.*(np.linspace(0.2, -0.2, num_choice)**2 - 0.04)
    choice_bottom[:,1] = np.linspace(0.2, -0.2, num_choice)

    num_choice = 32000
    choice_top = np.zeros((num_choice,2))
    choice_top[:,0] = -0.5*(np.linspace(2., -2., num_choice)**2 - 4.)
    choice_top[:,1] = np.linspace(2., -2., num_choice)

    idx_left, u[0,:],v[0,:] = gridgen.cm_gridgen.init_index(choice_left, shape[1])
    idx_right, u[-1,:],v[-1,:] = gridgen.cm_gridgen.init_index(choice_right, shape[1])
    idx_bottom, u[:,0],v[:,0] = gridgen.cm_gridgen.init_index(choice_bottom, shape[0])
    idx_top, u[:,-1],v[:,-1] = gridgen.cm_gridgen.init_index(choice_top, shape[0])

    iloop = True
    while iloop:

        # for i in range(2):
        #     u,v = gridgen.cm_gridgen.inverse_laplace(u, v, x[1]-x[0], y[1]-y[0], 1.e-10)
        #
        #     unit_vector, dev = gridgen.cm_gridgen.evaluate_left_boundary(u,v)
        #     idx_left1, u[0,1:-1],v[0,1:-1] = gridgen.cm_gridgen.choose_boundary(choice_left, unit_vector, dev, idx_left)
        #
        #     unit_vector, dev = gridgen.cm_gridgen.evaluate_right_boundary(u, v)
        #     idx_right1, u[-1,1:-1],v[-1,1:-1] = gridgen.cm_gridgen.choose_boundary(choice_right, unit_vector, dev, idx_right)
        #
        #
        #     unit_vector, dev = gridgen.cm_gridgen.evaluate_bottom_boundary(u, v)
        #     idx_bottom1, u[1:-1,0],v[1:-1,0] = gridgen.cm_gridgen.choose_boundary(choice_bottom, unit_vector, dev, idx_bottom)
        #
        #     unit_vector, dev = gridgen.cm_gridgen.evaluate_top_boundary(u, v)
        #     idx_top1, u[1:-1,-1],v[1:-1,-1] = gridgen.cm_gridgen.choose_boundary(choice_top, unit_vector, dev, idx_top)
        #     if np.all(idx_left==idx_left1) \
        #             and np.all(idx_right==idx_right1) \
        #             and np.all(idx_bottom==idx_bottom1)  \
        #             and np.all(idx_top==idx_top1) : iloop=False
        #     idx_right = idx_right1
        #     idx_left = idx_left1
        #     idx_bottom = idx_bottom1
        #     idx_top = idx_top1

        u,v,idx_left,idx_right,idx_bottom,idx_top,iconverge = gridgen.cm_gridgen.compute_grid(u,v,idx_left,idx_right,
                                                                                              idx_bottom,idx_top,
                                                                                              choice_left,choice_right,
                                                                                              choice_bottom,choice_top, 1e-10, 0.5, 100)
        # print dev
        # iconverge=1
        #print left_angle
        grid = {'x':u, 'y':v}
        IO_util.write_hdf5('TEST', grid)
        if iconverge==1:
            print 'Converged!'
            iloop=False

