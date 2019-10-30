import sys
import xarray as xr
import glob
import cartopy.crs as ccrs
import cmocean
import matplotlib.pyplot as plt
import re
import datetime
from functions.grid_to_km import area
import cartopy
import cartopy.crs as ccrs
import matplotlib as mpl
import numpy as np
from scipy.sparse import csr_matrix
import matplotlib.animation as animation
from matplotlib import rc


import advect_ice
import lay_out



fig_dir = "/home/students/6252699/thesis/parcels2/figures/"
out_dir = "/scratch/AnnekeV/output/"

plt.style.use(['default'])


def set_axes():
    plt.figure(figsize = [8,8])
    ax = plt.axes(projection=ccrs.NorthPolarStereo())
    ax.gridlines(xlocs = np.arange(-180,185,20), ylocs = np.arange(-90, 95,5), color='black', alpha=0.5, linestyle=':')
    ax.add_feature(cartopy.feature.OCEAN)
    ax.add_feature(cartopy.feature.LAND)
    ax.set_extent([-180,180,50,90],  ccrs.PlateCarree())
    ax.stock_img()
    return ax

output_name = advect_ice.run_experiment(start_date='2014-02-01' ,  custom_kernel = "AdvectionRK4_ocean", simdays=61) 

ds = xr.open_dataset(output_name + '.nc', decode_times = False)

print("Ran advect ice")
print(ds)



# NEW WAY OF CLASSIFYING

def remove_sinks(M):
    '''Remove sinks once'''
    Mdia = np.diag(M)[:-1]
    sinks = np.where(Mdia[:-1]==1)[0]
    Mnosinks = M.copy()
    Mnosinks[sinks, sinks] = 0
    Mnosinks_sum = np.sum(Mnosinks, axis=0)
    Mnosinks_sum[Mnosinks_sum==0] =1
    Mnosinks/=Mnosinks_sum
    return Mnosinks

res = .5
N   = 360//res

lons_init = ds.lon.isel(obs= 0).values
lats_init = ds.lat.isel(obs= 0).values
in_ice_i  = np.floor(ds.in_ice.isel(obs=0).values+0.85)     # fill with the correct values!!!


nparts = len(lons_init)
particle_grid_index = np.zeros((2, nparts), dtype=int)
ncels = int ((360//res) * (40//res))

print("ncels = {}")


for i in range(nparts):
    if np.isnan(lons_init[i]):
        particle_grid_index[0,i] = -1
    else:
        particle_grid_index[0,i] = int((lats_init[i]-50)//res)*N  + int((lons_init[i]+ 180)//res) + int(in_ice_i[i]) * ncels
        


lons_f   = ds.lon.isel(obs= -1).values
lats_f   = ds.lat.isel(obs= -1).values
in_ice_f = ds.in_ice.isel(obs=-1).values


for i in range(nparts):
    if np.isnan(lons_f[i]):
        particle_grid_index[1,i] = -1
    else:
        particle_grid_index[1,i] = int((lats_f[i]-50)//res)*N  + int((lons_f[i]+ 180)//res) + int(in_ice_f[i]) * ncels
        


T = np.zeros((ncels*2+1,ncels*2+1))


for i in range(nparts):
    T[particle_grid_index[1,i], particle_grid_index[0,i]] += 1

sumT =  np.sum(T, axis=0)
sumT[sumT ==0] = 1
T /=sumT


max_sink = (np.diag(T)[:-1]).max()
count    = 0
while (max_sink == 1 and count <10):
    T        = remove_sinks(T)
    count   += 1
    max_sink = (np.diag(T)[:-1]).max()

print ("Count is {}" .format(count))
T[:,-1]  = 0
T[-1,-1] = 1

latbins, lonbins= np.arange(50,90+res, res),  np.arange(-180,180+res, res)

times = 5
nlat     = len(latbins)-1
nlon     = len(lonbins)-1

R        = np.zeros((times, ncels*2+1))
R[0,:]   = 1.
M_sparse = csr_matrix(T)

for t in range(1,times):
     R[t,:] = M_sparse.dot(R[t-1,:])

print("Ready with transition matrix")
        
line = sorted(R[R!=0])
vmax = line[int(len(line)*0.95)]  # 95 percentil
vmin = line[int(len(line)*0.01)]  # 1 percentile 
print (vmin, vmax)
nlat, nlon = len(latbins)-1, len(lonbins)-1

times = np.size(R, 0)

for t in range(times):
    ax   = set_axes()
    p    = ax.pcolormesh(lonbins, latbins, R[t,0:ncels].reshape(nlat,nlon), transform = ccrs.PlateCarree(),  
                         norm = mpl.colors.LogNorm(vmin = vmin, vmax = vmax), cmap = 'viridis')
    plt.colorbar(p)
    ax.set_title("Ocean \nt = {}".format(t) )
#     plt.show()
    plt.savefig(fig_dir + "/transitionplots/ocean_separate_ocean_res_{}_t_{}.png".format(res, t))
#     else: plt.show

print("Ready with half the pictures")

for t in range(times):
    ax   = set_axes()
    p    = ax.pcolormesh(lonbins, latbins, R[t,ncels:-1].reshape(nlat,nlon), transform = ccrs.PlateCarree(),  
                         norm = mpl.colors.LogNorm(vmin = vmin, vmax = vmax), cmap = 'viridis')
    plt.colorbar(p)
    ax.set_title("Ice\nt = {}".format(t) )
    plt.savefig(fig_dir + "/transitionplots/ocean_separate_ice_res_{}_t_{}.png".format(res,t))

print("Done!")

#     plt.show()