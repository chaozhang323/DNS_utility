import numpy as np
import scipy.interpolate as interpolate
import matplotlib.pyplot as plt

def main():
    filename_data_big = "/usr/local/home/yl8bc/duannas/duanl/Acoustics/TestFluentConversion/WallProfile_fromFluentStructured.dat"
    filename_data = "/usr/local/home/yl8bc/duannas/duanl/Acoustics/TestFluentConversion/WallProfile_fromFluentStructured.dat"
    #filename_data = "/usr/local/home/yl8bc/duannas/duanl/Acoustics/TestFluentConversion/NozzleProfile_fromFluentStructured_pointwise4.dat"
    #filename_data = "/usr/local/home/yl8bc/duannas/duanl/Acoustics/TestFluentConversion/NozzleContour_fromExcel_raw.dat"
    
    num_skiprows = 12
    
    with open(filename_data, 'r') as fid:
        header = [next(fid) for x in range(num_skiprows)]
    print(header)
    
    data_big = np.loadtxt(filename_data_big, skiprows=num_skiprows)
    

    data = np.loadtxt(filename_data, skiprows=num_skiprows)
    x = data[:,0]
    y = data[:,1]
    
    ##
    m = len(x)
    tck = interpolate.splrep(x, y, s=100.e-9, k=5)
    yi = interpolate.splev(x, tck)

    ## 
    fig, ax = plt.subplots(3)
    #ax.plot(x, y, '.')
    #ax.plot(x, yi, '-')
    xp, yp = cd(x, y)
    xpp, ypp = cd(xp, yp)
    xip, yip = cd(x, yi)
    xipp, yipp = cd(xip, yip)
    
    ## 
    print(np.linalg.norm(y-yi, np.inf))
    print(np.linalg.norm(yp-yip))
    print(np.linalg.norm(ypp-yipp))
    
    ##
    ax[0].plot(x, yi, '.-')
    ax[0].plot(x, y, '-')
    ax[1].plot(xip, yip, '.-')
    ax[1].plot(xp, yp, '-')
    ax[2].plot(xipp, yipp, '.-')
    ax[2].plot(xpp, ypp, '-')
    
    #ax[0].plot(data_big[:,0], data_big[:,1], 'r') 
    plt.show()
  
    np.savetxt('original_profile.dat', data, header=''.join(header), comments='')
    data[:,1] = yi
    np.savetxt('interpolated_profile.dat', data, header=''.join(header), comments='')

def cd(x, y):
    yp = (y[2:] - y[:-2]) / (x[2:] - x[:-2])
    xp = x[1:-1]
    return xp, yp

def fd(x, y):
    yp = (y[1:] - y[:-1]) / (x[1:] - x[:-1])
    xp = x[0:-1]
    return xp, yp


if __name__=="__main__":
    main()

