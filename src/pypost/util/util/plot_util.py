import matplotlib.pyplot as plt
import matplotlib.transforms as trans

def get_perfect_extent(fig, ax, yxratio, dpi=80, fig_full_size=8):
    ax.axis('off')
    fig.set_dpi(dpi)
    fig.set_size_inches( (fig_full_size, fig_full_size) )

    if yxratio <= 1.:
        ax.set_position( trans.Bbox.from_bounds(0., 0., 1., yxratio ))
    else:
        ax.set_position( trans.Bbox.from_bounds(0., 0., 1./yxratio, 1. ))
    extent = ax.get_window_extent().transformed( fig.dpi_scale_trans.inverted() )

    return extent