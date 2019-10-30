import numpy as np


def grid_per_x_km(xkm = 1e4, ylims = [50,90]):
    "Makes a grid per x number of m, evenly spaced per latitude. Default is 10 km between 50 and 90 lat. Returns lons, lats "
    R     = 6371.0e3              #m Radius earth
    xdeg  = xkm/(R*np.pi/180.)    # per degree
    lat   = np.arange(ylims[0],ylims[1], xdeg)
    n     = np.array((np.cos(lat/180.*np.pi)*R/xkm), dtype='int')
    
    lats, lons = [], []
    for i in range(len(n)):
        lats.extend(np.repeat(lat[i], n[i]))
        lons.extend(np.linspace(-180,180,n[i], endpoint=False))
        
    return lons, lats