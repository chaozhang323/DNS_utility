import h5py
import matplotlib.pyplot as plt

def labels(axe, xlabel, ylabel, fontsize=16):
    axe.set_xlabel(xlabel, fontdict={'size':fontsize})
    axe.set_ylabel(ylabel, fontdict={'size':fontsize})
    axe.tick_params(axis='both', labelsize=fontsize)
    return

def title(axe, title, fontsize=16):
    axe.set_title(title, fontdict={'size':fontsize})
    return