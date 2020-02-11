import numpy as np

def correlation_grid(Fs, Fi, idx0):
    domain = Fi.shape
    box_size = Fs.shape
    idx = idx0
    corr_p = np.abs( np.corrcoef( Fs.ravel(), Fi[idx[0]:idx[0]+box_size[0], idx[1]:idx[1]+box_size[1]].ravel() )[0,1] )

    corr_nb = np.zeros([4,1])
    def get_corr_nb(Fs, Fi, dir, idx, box_size, domain):
        if dir==0 and idx[0]-1 >= 0:
            corr = np.abs( np.corrcoef( Fs.ravel(), Fi[idx[0]-1:idx[0]-1+box_size[0], idx[1]:idx[1]+box_size[1]].ravel() )[0,1] )
        elif dir==1 and idx[0]+1+box_size[0] <= domain[0]:
            corr = np.abs( np.corrcoef( Fs.ravel(), Fi[idx[0]+1:idx[0]+1+box_size[0], idx[1]:idx[1]+box_size[1]].ravel() )[0,1] )
        elif dir==2 and idx[1]-1 >= 0:
            corr = np.abs( np.corrcoef( Fs.ravel(), Fi[idx[0]:idx[0]+box_size[0], idx[1]-1:idx[1]-1+box_size[1]].ravel() )[0,1] )
        elif dir==3 and idx[1]+1+box_size[1] <= domain[1]:
            corr = np.abs( np.corrcoef( Fs.ravel(), Fi[idx[0]:idx[0]+box_size[0], idx[1]+1:idx[1]+1+box_size[1]].ravel() )[0,1] )
        else:
            corr = 0.
        return corr

    for i in range(4):
        corr_nb[i] = get_corr_nb(Fs, Fi, i, idx, box_size, domain)

    idx_corrmax = np.argmax(corr_nb)
    corrmax = corr_nb[idx_corrmax]

    def step(idx, idx_corrmax, corrmax, corr_p):
        idx = idx.copy()
        if idx_corrmax==0:
            idx[0] -= 1
            idx_back = 1
        if idx_corrmax==1:
            idx[0] += 1
            idx_back = 0
        if idx_corrmax==2:
            idx[1] -= 1
            idx_back = 3
        if idx_corrmax==3:
            idx[1] += 1
            idx_back = 2
        corr_nb[idx_back] = corr_p
        corr_p = corrmax
        return idx, idx_back, corr_p


    if corrmax > corr_p:
        idx, idx_back, corr_p = step(idx, idx_corrmax, corrmax, corr_p)
    else:
        return idx-idx0, corr_p

    ## iterating
    while True:
        for i in range(4):
            if i!=idx_back:
                corr_nb[i] = get_corr_nb(Fs, Fi, i, idx, box_size, domain)
        idx_corrmax = np.argmax(corr_nb)
        corrmax = corr_nb[idx_corrmax]
        if corrmax > corr_p:
            idx, idx_back, corr_p = step(idx, idx_corrmax, corrmax, corr_p)
        else:
            return idx-idx0, corr_p


def correlation_patch_max(Fs, Fi, idx0, patch_size=5):
    idx = idx0.copy()
    corr_patch = np.zeros([patch_size*2+1, patch_size*2+1])

    idx1 = idx.copy()
    while True:
        for i in range(patch_size*2+1):
            for j in range(patch_size*2+1):
                slc = (slice(idx[0]+i-patch_size,idx[0]+Fs.shape[0]+i-patch_size), \
                       slice(idx[1]+j-patch_size,idx[1]+Fs.shape[1]+j-patch_size))
                if if_slice_in_domain(slc, Fi.shape):
                    corr_patch[i,j] = np.abs( np.corrcoef( Fs.ravel(), Fi[slc].ravel() )[0,1] )
                    if np.isnan(corr_patch[i,j]):
                        corr_patch[i,j] = 0.
                else:
                    corr_patch[i,j] = 0.
        corr_max = np.max(corr_patch)
        corr_max_idx = np.unravel_index( np.argmax(corr_patch), corr_patch.shape )

        corr_p = corr_patch[patch_size,patch_size]
        if corr_max > corr_p:
            idx = idx1 + np.array(corr_max_idx) - patch_size
            idx1[:] = idx
        else:
            return idx-idx0




def if_slice_in_domain(slc, domain):
    if slc[0].start>=0 and slc[0].stop<=domain[0] and slc[1].start>=0 and slc[1].stop<=domain[1]:
        result = True
    else:
        result = False
    return result





