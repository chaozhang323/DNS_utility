import numpy as np

def tinker(a, nmlst):
    if nmlst is not None:
        nmlst_out = nmlst
    else:
        keys = a.keys()
        nmlst_out = {key:None for key in keys}
    return nmlst_out



def slice_(a, slc, namelist=None):
    nmlst_work = tinker(a, namelist)
    
    out = {}
    for key in nmlst_work.keys():
        a1 = a[key]
        if type(a1) is dict:
            out[key] = slice_(a1, slc, nmlst_work[key])
        else:
            if a1.ndim==len(slc):
                out[key] = a1[slc]
            else:
                out[key] = a1
    return out

def add(a, b, namelist=None):
    nmlst_work = tinker(a, namelist)

    out = {}
    for key in nmlst_work.keys():
        a1 = a[key]
        b1 = b[key]
        if type(a1) is dict:
            out[key] = add(a1, b1, nmlst_work[key])
        else:
            out[key] = a1 + b1
    return out    

def divide(a, b, namelist=None):
    nmlst_work = tinker(a, namelist)

    out = {}
    for key in nmlst_work.keys():
        a1 = a[key]
        b1 = b[key]
        if type(a1) is dict:
            out[key] = divide(a1, b1, nmlst_work[key])
        else:
            out[key] = a1 / b1
    return out    

def divide_scalar(a, b, namelist=None):
    nmlst_work = tinker(a, namelist)

    out = {}
    for key in nmlst_work.keys():
        a1 = a[key]
        if type(a1) is dict:
            out[key] = divide_scalar(a1, b, nmlst_work[key])
        else:
            out[key] = a1 / b
    return out    
