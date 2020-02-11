import numpy as np

class pdist:
    def __init__(self, x0, xm, num):
        self.x0 = x0
        self.xm = xm
        self.num = num
        return

    def tanh(self, a, type='both'):
        x0 = self.x0
        xm = self.xm
        L = self.xm - self.x0
        x = np.linspace(self.x0, self.xm, self.num)
        def tanh_b():
            f = np.tanh(a*(2./L*(x-x0)-1.)) / np.tanh(a) *L*0.5 + (x0+xm)*0.5
            return f
        def tanh_l():
            f = np.tanh(a*(1./L*(x-x0)-1.)) / np.tanh(a) *L + xm
            return f
        def tanh_r():
            f = np.tanh(a*(1./L*(x-x0)   )) / np.tanh(a) * L + x0
            return f
        
        f_dict = {'both':tanh_b, 'left':tanh_l, 'right':tanh_r}
        return f_dict[type]()

class ptrans:
    def __init__(self, x):
        self.x = x
        return
   
    @property
    def x(self):
        return self._x
    @x.setter
    def x(self, x):
        if type(x) is not np.ndarray:
            raise TypeError("Input is not an ndarray.")
        if x.ndim!=1:
            raise ValueError("Input is not of rank 1.")
        if x.shape[0]<=2:
            raise ValueError("Input is less than 3 points.")
        self._x = x
        return



    def tanh(self, a, idx_c):
        x = self.x
        x0 = x[0]
        xm = x[-1]
        xc = x[idx_c]
        L = xm - x0
        t0 = np.tanh(a/L*(x0-xc))
        tm = np.tanh(a/L*(xm-xc))
        
        f = (np.tanh(a*(x-xc)/L) - t0) / (tm-t0) * (xm-x0)  +  x0
        return f

    def poly_3rd(self, r, idx_c):
        x = self.x
        x0 = x[0]
        xm = x[-1]
        xc = x[idx_c]
        L = xm - x0
        beta = (xm-xc) / L
        a = (1.-r) / beta / (1.-beta)

        x = (x-x0) / L
        f = a*x**3 - a*(1+xc)*x**2 + (1+a*xc)*x
        f = f*L + x0
        return f

    def expo(self, r, idx_c):
        x = self.x
        xl = x[0:idx_c+1]
        xr = x[idx_c:]
        dxl = xl[1:] - xl[:-1]
        dxr = xr[1:] - xr[:-1]
        Ll = xl[-1] - xl[0]
        Lr = xr[-1] - xr[0]
        dxl = dxl[::-1]
        
        def get_dx_new(r, L, dx):
            nx = dx.shape[0]
            a1 = 1.
            a2 = 1./r
            a = .5*(a1+a2)
            L_new = 0.
            while np.abs(L_new-L)>1.e-10:
                dx_new = np.array([r*dx[n]*a**n for n in range(nx)])
                L_new = np.sum(dx_new)
                if L_new < L:
                    a1 = a
                else:
                    a2 = a
                a = .5*(a1+a2)
                print(a, L_new)
            return dx_new
        
        dxl_new = get_dx_new(r, Ll, dxl)
        dxr_new = get_dx_new(r, Lr, dxr)
        f = x.copy()
        for n in range(dxl.shape[0]):
            f[idx_c-n-1] = f[idx_c-n] - dxl_new[n]
        for n in range(dxr.shape[0]):
            f[idx_c+n+1] = f[idx_c+n] + dxr_new[n]
        return f
