from util import IO_util, util_geom, util_grid, util_f
import matplotlib.pyplot as plt
import numpy as np
from .. import top

def gridgen(dims, body, shock, filename_out, check=False):
    ## define body
    x_right = body[:,0]
    y_right = body[:,1]

    ## define shock
    x_left = shock[:,0]
    y_left = shock[:,1]

    ## define outlet
    p_l1,v_l1 = top.get_pv_for_top(x_left, y_left)
    p_r1,v_r1 = top.get_pv_for_top(x_right, y_right)
    x_top, y_top = top.decide_top(p_l1, v_l1, p_r1, v_r1, 2**16)

    ## define central line
    x_bot = np.linspace(x_left[0], x_right[0], 2**16)
    y_bot = np.zeros_like(x_bot)

    ## for a look
    if check:
        plt.plot(x_bot,y_bot)
        plt.plot(x_top,y_top)
        plt.plot(x_left,y_left)
        plt.plot(x_right,y_right)
        plt.gca().set_aspect(1.)
        plt.show()

    ## define params for grid gen
    u = np.zeros(dims)
    v = np.zeros(dims)

    choice_left = np.stack([x_left,y_left], axis=1)
    choice_right = np.stack([x_right,y_right], axis=1)
    choice_bottom = np.stack([x_bot,y_bot], axis=1)
    choice_top = np.stack([x_top,y_top], axis=1)

    ## refinement
    choice_left = util_f.cm_gridgen.refine_boundary(choice_left, 2**16/choice_left.shape[0])
    choice_right = util_f.cm_gridgen.refine_boundary(choice_right, 2**16/choice_right.shape[0])


    ## computation
    u,v = util_grid.gridgen_orth(dims, choice_left, choice_right, choice_bottom, choice_top, \
                                    filename_out, varname=['x','z'], \
                                    tol_in=1.e-8, tol_bdry=0, nsave=64, \
                                    fixed_boundary=[0,0,0,0], bend_boundary=[0,0,1,0])
    
    return u,v
