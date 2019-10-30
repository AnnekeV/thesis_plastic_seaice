import numpy as np 


def area(latitude, degree_grid):

    rad = (latitude * np.pi) / (180)
    dist_lon = 111.320*np.cos(rad) * 1000 * degree_grid # Distance of longitude cell in m
    dist_lat = 110.574 * 1000 * degree_grid # Distance of latitude cell in m

    return dist_lon*dist_lat
