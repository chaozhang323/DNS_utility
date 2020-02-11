import numpy as np
import h5py
import os
from  .IO_util import *

class grid:
    """Cartesian Grid"""
    def __init__(self, filename):
        self.filename = filename
        self.dim = None

    def read_grid(self, dim=2, slc=None):
        x, y, z, ep, nx, ny, nz = read_grid(self.filename, dim=dim, slc=slc, staggered=False)
        self.dim = dim
        self.x = x
        self.y = y
        self.z = z
        self.ep = ep
        self.nx = nx
        self.ny = ny
        self.nz = nz
        return x, y, z, ep, nx, ny, nz

    def update_grid_info(self, Lx=1., Ly=1., Lz=1.):
        self.Lx = Lx
        self.Ly = Ly
        self.Lz = Lz

    def get_grid(self, dim=2, slc=None):
        if self.dim != dim:
            x, y, z, ep, nx, ny, nz = self.read_grid(dim=dim, slc=slc)
        x = self.x
        y = self.y
        z = self.z
        ep = self.ep
        nx = self.nx
        ny = self.ny
        nz = self.nz
        return x, y, z, ep, nx, ny, nz

    def get_dx(self):
        dx = self.Lx/self.nx
        dy = self.Ly/self.ny
        dz = self.Lz/self.nz
        return dx, dy, dz





