import numpy as np
pi = np.pi

def refine_incompact3d_1d(yly, ncly, ny, istret, beta):
    yp = np.zeros([ny])
    yeta = np.zeros([ny])
    ypi = np.zeros([ny])
    yetai = np.zeros([ny])
    ppy = np.zeros([ny])
    pp2y = np.zeros([ny])
    pp4y = np.zeros([ny])
    ppyi = np.zeros([ny])
    pp2yi = np.zeros([ny])
    pp4yi = np.zeros([ny])
    


    yinf = -yly / 2. 
    den = 2. * beta * yinf 
    xnum = -yinf - np.sqrt(pi * pi * beta * beta + yinf * yinf) 
    alpha = abs(xnum / den) 
    xcx = 1. / beta / alpha 
    
    if alpha != 0.:
        if istret == 1:
            yp[0] = 0. 
     
        if istret == 2:
            yp[0] = 0. 
     
        if istret == 1:
            yeta[0] = 0. 
     
        if istret == 2:
            yeta[0] = -1. / 2.
     
        if istret == 3:
            yp[0] = 0.
     
        if istret == 3:
            yeta[0] = -1. / 2.
     
        for j in range(1,ny):
            if istret == 1:
                if ncly == 0:
                    yeta[j] = (j - 1.) * (1. / ny)
         
            if (ncly == 1) | (ncly == 2):
                yeta[j] = (j - 1.) * (1. / (ny - 1.))

            if istret == 2:
                if ncly == 0:
                    yeta[j] = (j - 1.) * (1. / ny) - 0.5
         
                if (ncly == 1) | (ncly == 2):
                    yeta[j] = (j - 1.) * (1. / (ny - 1.)) - 0.5

            if istret == 3:
                if ncly == 0:
                    yeta[j] = ((j - 1.) * (1. / 2. / ny) - 0.5)
         
                if (ncly == 1) | (ncly == 2):
                    yeta[j] = ((j - 1.) * (1. / 2. / (ny - 1.)) - 0.5)
         
     
            den1 = np.sqrt(alpha * beta + 1.)
            xnum = den1 / np.sqrt(alpha / pi) / np.sqrt(beta) / np.sqrt(pi)
            den = 2. * np.sqrt(alpha / pi) * np.sqrt(beta) * pi * np.sqrt(pi)
            den3 = ((np.sin(pi * yeta[j])) * (np.sin(pi * yeta[j])) / beta / pi) + alpha / pi
            den4 = 2. * alpha * beta - np.cos(2. * pi * yeta[j]) + 1.
            if (ncly != 0) & (j == ny) & (istret == 1):
                xnum1 = 0.

            else:
                xnum1 = (np.arctan(xnum * np.tan(pi * yeta[j]))) * den4 / den1 / den3 / den
 
            cst = np.sqrt(beta) * pi / (2. * np.sqrt(alpha) * np.sqrt(alpha * beta + 1.))
            if istret == 1:
                if yeta[j] < 0.5:
                    yp[j] = xnum1 - cst - yinf
     
                if yeta[j] == 0.5:
                    yp[j] = 0. - yinf
     
                if yeta[j] > 0.5:
                    yp[j] = xnum1 + cst - yinf

            if istret == 2:
                if yeta[j] < 0.5:
                    yp[j] = xnum1 - cst + yly

                if yeta[j] == 0.5:
                    yp[j] = 0. + yly

                if yeta[j] > 0.5:
                    yp[j] = xnum1 + cst + yly

            if istret == 3:
                if yeta[j] < 0.5:
                    yp[j] = (xnum1 - cst + yly) * 2.

                if yeta[j] == 0.5:
                    yp[j] = (0. + yly) * 2.

                if yeta[j] > 0.5:
                    yp[j] = (xnum1 + cst + yly) * 2.

    if alpha == 0.:
        yp[0] = -1.e10
        for j in range(1,ny):
            yeta[j] = (j - 1.) * (1. / ny)
            yp[j] = -beta * np.cos(pi * yeta[j]) / np.sin(yeta[j] * pi)
 
 
    if alpha!=0.:
        for j in range(0,ny):
            if istret == 1:
                if ncly == 0:
                    yetai[j] = (j - 0.5) * (1. / ny)

                if (ncly == 1) | (ncly == 2):
                    yetai[j] = (j - 0.5) * (1. / (ny - 1.))

            if istret == 2:
                if ncly == 0:
                    yetai[j] = (j - 0.5) * (1. / ny) - 0.5

                if (ncly == 1) | (ncly == 2):
                    yetai[j] = (j - 0.5) * (1. / (ny - 1.)) - 0.5

            if istret == 3:
                if ncly == 0:
                    yetai[j] = (j - 0.5) * (1. / 2. / ny) - 0.5

                if (ncly == 1) | (ncly == 2):
                    yetai[j] = (j - 0.5) * (1. / 2. / (ny - 1.)) - 0.5

            den1 = np.sqrt(alpha * beta + 1.)
            xnum = den1 / np.sqrt(alpha / pi) / np.sqrt(beta) / np.sqrt(pi)
            den = 2. * np.sqrt(alpha / pi) * np.sqrt(beta) * pi * np.sqrt(pi)
            den3 = ((np.sin(pi * yetai[j])) * (np.sin(pi * yetai[j])) / beta / pi) + alpha / pi
            den4 = 2. * alpha * beta - np.cos(2. * pi * yetai[j]) + 1.
            xnum1 = (np.arctan(xnum * np.tan(pi * yetai[j]))) * den4 / den1 / den3 / den
            cst = np.sqrt(beta) * pi / (2. * np.sqrt(alpha) * np.sqrt(alpha * beta + 1.))
            if istret == 1:
                if yetai[j] < 0.5:
                    ypi[j] = xnum1 - cst - yinf

                if yetai[j] == 0.5:
                    ypi[j] = 0. - yinf

                if yetai[j] > 0.5:
                    ypi[j] = xnum1 + cst - yinf
 
 
            if istret == 2:
                if yetai[j] < 0.5:
                    ypi[j] = xnum1 - cst + yly

                if yetai[j] == 0.5:
                    ypi[j] = 0. + yly

                if yetai[j] > 0.5:
                    ypi[j] = xnum1 + cst + yly
     
 
            if istret == 3:
                if yetai[j] < 0.5:
                    ypi[j] = (xnum1 - cst + yly) * 2.

                if yetai[j] == 0.5:
                    ypi[j] = (0. + yly) * 2.

                if yetai[j] > 0.5:
                    ypi[j] = (xnum1 + cst + yly) * 2.

    if alpha == 0.:
        ypi[0] = -1.e10
        for j in range(1,ny):
            yetai[j] = (j - 1.) * (1. / ny)
            ypi[j] = -beta * np.cos(pi * yetai[j]) / np.sin(yetai[j] * pi)
 
 

    for j in range(0,ny):
        ppy[j] = yly * (alpha / pi + (1. / pi / beta) * np.sin(pi * yeta[j]) * np.sin(pi * yeta[j]))
        pp2y[j] = ppy[j] * ppy[j]
        pp4y[j] = (-2. / beta * np.cos(pi * yeta[j]) * np.sin(pi * yeta[j]))
 
    for j in range(0,ny):
        ppyi[j] = yly * (alpha / pi + (1. / pi / beta) * np.sin(pi * yetai[j]) * np.sin(pi * yetai[j]))
        pp2yi[j] = ppyi[j] * ppyi[j]
        pp4yi[j] = (-2. / beta * np.cos(pi * yetai[j]) * np.sin(pi * yetai[j]))
 
    return yp,ypi,ppy,pp2y,pp4y,ppyi,pp2yi,pp4yi