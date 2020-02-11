import numpy as np
pi = np.pi

def derive_1st_lele(f, dx, nx, nc):
    fbx = np.zeros([nx])
    fcx = np.zeros([nx])
    ffx = np.zeros([nx])
    rhs = np.zeros([nx])
    rhs1 = np.zeros([nx])
    rhs2 = np.zeros([nx])

    al = 1./3.; a = 14./9.; b = 1./9.;
    if nc == 0:
        fbx[0] = 0.
        fcx[0] = 2.
        ffx[0] = al
        fbx[1] = al
        fcx[1] = 1
        ffx[1] = al
        rhs1[0] = a/(2.*dx)*(f[1]-f[nx-1]) + b/(4.*dx)*(f[2]-f[nx-2])
        rhs1[1] = a/(2.*dx)*(f[2]-f[0]) + b/(4.*dx)*(f[3]-f[nx-1])
        rhs2[0] = -1.
        rhs2[1] = 0
        for i in range(2,nx-2):
            fbx[i] = al
            fcx[i] = 1.
            ffx[i] = al
            rhs1[i] = a/(2.*dx)*(f[i+1]-f[i-1]) + b/(4.*dx)*(f[i+2]-f[i-2])
            rhs2[i] = 0
        fbx[nx-2] = al
        fcx[nx-2] = 1.
        ffx[nx-2] = al
        fbx[nx-1] = al
        fcx[nx-1] = 1. + al*al
        ffx[nx-1] = 0.
        rhs1[nx-2] = a/(2.*dx)*(f[nx-1]-f[nx-3]) + b/(4.*dx)*(f[0]-f[nx-4])
        rhs1[nx-1] = a/(2.*dx)*(f[0]-f[nx-2]) + b/(4.*dx)*(f[1]-f[nx-3])
        rhs2[nx-2] = 0.
        rhs2[nx-1] = al

        ans1 = TDMAsolver(fbx,fcx,ffx,rhs1)
        ans2 = TDMAsolver(fbx,fcx,ffx,rhs2)
        ans = ans1 - ans2*(ans1[0] -al*ans1[nx-1])/(1.+ans2[0]-al*ans2[nx-1])
        ans = ans.reshape([nx])
        return ans
    else:
        fbx[0] = 0.
        fcx[0] = 1.
        ffx[0] = 2.
        fbx[1] = 1./4.
        fcx[1] = 1.
        ffx[1] = 1./4.
        rhs[0] = 1./(2.*dx) * (-5.*f[0]+4.*f[1]+f[2])
        rhs[1] = 3./2./(2.*dx) * (f[2] - f[0])

        for i in range(2,nx-2):
            fbx[i] = al
            fcx[i] = 1.
            ffx[i] = al
            rhs[i] = a/(2.*dx)*(f[i+1]-f[i-1]) + b/(4*dx)*(f[i+2]-f[i-2])

        fbx[nx-2] = 1./4.
        fcx[nx-2] = 1.
        ffx[nx-2] = 1./4.
        fbx[nx-1] = 2.
        fcx[nx-1] = 1.
        ffx[nx-1] = 0.
        rhs[nx-2] = 3./2./(2.*dx) * (-f[nx-3] + f[nx-1])
        rhs[nx-1] = 1./(2.*dx) * (5.*f[nx-1]-4.*f[nx-2]-f[nx-3])

        ans = TDMAsolver(fbx,fcx,ffx,rhs)
        ans = ans.reshape([nx])
        return ans

def derive_2nd_lele(f, dx, nx, nc):
    fbx = np.zeros([nx])
    fcx = np.zeros([nx])
    ffx = np.zeros([nx])
    rhs = np.zeros([nx])
    rhs1 = np.zeros([nx])
    rhs2 = np.zeros([nx])

    fpi2 = 4.0
    al = (45.0 * fpi2 * pi * pi - 272.0) / (2.0 * (45.0 * fpi2 * pi * pi - 208.0))
    a = (6.-9.*al) / 4. / (dx*dx)
    b = (-3.+24.*al) / 5. / (4.*dx*dx)
    c = (2.-11.*al) / 20. / (9.*dx*dx)

    if nc == 0:
        fbx[0] = 0.
        fcx[0] = 2.
        ffx[0] = al
        fbx[1] = al
        fcx[1] = 1.
        ffx[1] = al
        fbx[2] = al
        fcx[2] = 1.
        ffx[2] = al
        rhs1[0] = a*(f[1]-2.*f[0]+f[nx-1]) + b*(f[2]-2.*f[0]+f[nx-2]) + c*(f[3]-2.*f[0]+f[nx-3])
        rhs1[1] = a*(f[2]-2.*f[1]+f[0]) + b*(f[3]-2.*f[1]+f[nx-1]) + c*(f[4]-2.*f[1]+f[nx-2])
        rhs1[2] = a*(f[3]-2.*f[2]+f[1]) + b*(f[4]-2.*f[2]+f[0]) + c*(f[5]-2.*f[2]+f[nx-1])
        rhs2[0] = -1.
        rhs2[1] = 0.
        rhs2[2] = 0.
        for i in range(3,nx-3):
            fbx[i] = al
            fcx[i] = 1.
            ffx[i] = al
            rhs1[i] = a*(f[i+1]-2.*f[i]+f[i-1]) + b*(f[i+2]-2.*f[i]+f[i-2]) + c*(f[i+3]-2.*f[i]+f[i-3])
            rhs2[i] = 0.
        fbx[nx-3] = al
        fcx[nx-3] = 1.
        ffx[nx-3] = al
        fbx[nx-2] = al
        fcx[nx-2] = 1.
        ffx[nx-2] = al
        fbx[nx-1] = al
        fcx[nx-1] = 1. + al*al
        ffx[nx-1] = 0

        rhs1[nx-3] = a*(f[nx-2]-2.*f[nx-3]+f[nx-4]) + b*(f[nx-1]-2.*f[nx-3]+f[nx-5]) + c*(f[0]-2.*f[nx-3]+f[nx-6])
        rhs1[nx-2] = a*(f[nx-1]-2.*f[nx-2]+f[nx-3]) + b*(f[0]-2.*f[nx-2]+f[nx-4]) + c*(f[1]-2.*f[nx-2]+f[nx-5])
        rhs1[nx-1] = a*(f[0]-2.*f[nx-1]+f[nx-2]) + b*(f[1]-2.*f[nx-1]+f[nx-3]) + c*(f[2]-2.*f[nx-1]+f[nx-4])
        rhs2[nx-3] = 0.
        rhs2[nx-2] = 0.
        rhs2[nx-1] = al

        ans1 = TDMAsolver(fbx,fcx,ffx,rhs1)
        ans2 = TDMAsolver(fbx,fcx,ffx,rhs2)
        ans = ans1 - ans2*(ans1[0] -al*ans1[nx-1])/(1.+ans2[0]-al*ans2[nx-1])
        ans = ans.reshape([nx])
        return ans
    else:
        fbx[0] = 0.
        fcx[0] = 1.
        ffx[0] = 11.
        fbx[1] = 1./10.
        fcx[1] = 1.
        ffx[1] = 1./10.
        fbx[2] = 2./11.
        fcx[2] = 1.
        ffx[2] = 2./11.
        rhs1[0] = 1./(dx*dx) * (13.*f[0] - 27.*f[1] + 15.*f[2] - f[3])
        rhs1[1] = 6./5./(dx*dx) * (f[2]-2.*f[1]+f[0])
        rhs1[2] = 12./11./(dx*dx) * (f[3]-2.*f[2]+f[1]) + 3./11./(4.*dx*dx)*(f[4]-2.*f[2]+f[0])
        for i in range(3,nx-3):
            fbx[i] = al
            fcx[i] = 1.
            ffx[i] = al
            rhs[i] = a*(f[i+1]-2.*f[i]+f[i-1]) + b*(f[i+2]-2.*f[i]+f[i-2]) + c*(f[i+3]-2.*f[i]+f[i-3])
        fbx[nx-3] = 2./11.
        fcx[nx-3] = 1.
        ffx[nx-3] = 2./11.
        fbx[nx-2] = 1./10.
        fcx[nx-2] = 1.
        ffx[nx-2] = 1./10.
        fbx[nx-1] = 11.
        fcx[nx-1] = 1.
        ffx[nx-1] = 0.

        rhs1[nx-1] = 1./(dx*dx) * (13.*f[nx-1] - 27.*f[nx-2] + 15.*f[nx-3] - f[nx-4])
        rhs1[nx-2] = 6./5./(dx*dx) * (f[nx-3]-2.*f[nx-2]+f[nx-1])
        rhs1[nx-3] = 12./11./(dx*dx) * (f[nx-2]-2.*f[nx-3]+f[nx-4]) + 3./11./(4.*dx*dx)*(f[nx-1]-2.*f[nx-3]+f[nx-5])

        ans = TDMAsolver(fbx,fcx,ffx,rhs)
        ans = ans.reshape([nx])
        return ans

def TDMAsolver(a, b, c, d):
    n = len(a)
    c[0] /= b[0]
    d[0] /= b[0]
    for i in range(1,n):
        c[i] = c[i] / (b[i]-a[i]*c[i-1])
        d[i] = (d[i]-a[i]*d[i-1]) / (b[i]-a[i]*c[i-1])
    x = np.zeros([n,1])
    x[n-1] = d[n-1]
    for i in range(n-2,-1,-1):
        x[i] = d[i] - c[i]*x[i+1]

    return x


    return xc