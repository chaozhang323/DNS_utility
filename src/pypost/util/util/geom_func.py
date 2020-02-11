import numpy as np
from scipy.spatial.distance import *

def get_2points_with_largest_distance(geom):
    nx, ny = geom.shape
    points = list()
    for i in range(nx):
        for j in range(ny):
            if geom[i,j] != 0:
                points.append([i,j])

    distance = cdist(points, points)
    p_max = np.argmax(distance)
    p_max = np.unravel_index(p_max, distance.shape)
    p1_max = points[p_max[0]]
    p2_max = points[p_max[1]]
    return p1_max, p2_max

def get_angle_from_2points(idx1, idx2):
    dx = idx2[0] - idx1[0]
    dy = idx2[1] - idx1[1]
    dx = float(dx)
    dy = float(dy)
    if dx < 0.:
        dy = -dy
        dx = -dx

    if dx == 0.:
        angle = np.pi/2.
    else:
        angle = np.arctan2(dy, dx)

    return angle