import matplotlib as mpl
mpl.use('Agg')
import plot_functions as pf

import numpy as np
from parcels import FieldSet, ParticleSet, ScipyParticle, JITParticle, ErrorCode, AdvectionRK4, Variable
from argparse import ArgumentParser
from datetime import timedelta
from datetime import datetime
from glob import glob
import netCDF4 as netcdf
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
import xarray as xr
import cmocean
from timeit import default_timer as timer
import re
import pandas as pd
import sys

# -------- FUNCTIONS -----------------------------



def ax_features(ax):
    ax.gridlines(xlocs = np.arange(-180,185,20), ylocs = np.arange(-90, 95,5), color='black', alpha=0.5, linestyle=':')
#     ax.add_feature(pf.ocean_50m)
    ax.add_feature(pf.land_50m, facecolor = "lightgrey", edgecolor=None)
    pf.set_circular_boundary(ax)
    ax.set_extent([-180,180, 50,90],  ccrs.PlateCarree())
    return ax


# -------- IMPORT FILES -----------------------------

fnames    = sorted(glob("/data/oceanparcels/CMEMS-GLORYS12V1-Arctic/*.nc"))
ds_glorys = xr.open_mfdataset(fnames)

latitude  = ds_glorys.latitude
longitude = ds_glorys.longitude
ds_glorys = ds_glorys.isel(depth=0)
ds_season = ds_glorys.groupby('time.season').mean('time')


print("Imported files...")



# -------- PLOT ICE THICKNESS PER SEASON -----------------------------



count  = 0

variable   = ds_glorys.sithick
var_season = ds_season.sithick
levels     = np.arange(0,4.5,.5)

fig        = plt.figure(figsize=[12,12])


for season in var_season.season:
    count+=1
   
    ax     = fig.add_subplot(2,2, count, projection=ccrs.NorthPolarStereo())
    ax     = ax_features(ax)

    p      = ax.pcolormesh(longitude, latitude, var_season.sel(season=season),  
                           vmin = 0, 
                           vmax= np.max(levels),
                           transform = ccrs.PlateCarree())
    tit = '{:s}'.format(season.values)
    print tit

    ax.set_title(tit)
    
    
fig.suptitle(variable.long_name)
fig.subplots_adjust(bottom = 0.1)
cbar_ax = fig.add_axes([0.15, 0.06, 0.7, .03])   #left, bottom, width , height

cbar = fig.colorbar(p, cax=cbar_ax, orientation = 'horizontal', extend = 'max')
cbar.set_label(variable.units)
cbar.set_ticks(levels)

plt.savefig("/home/students/6252699/thesis/parcels2/figures/variables/{}_per_season.png".format(variable.name), bbox_inches="tight", pad_inches=0.1)

print("Plotted per season...")


# -------- PLOT ICE THICKNESS PER YEAR -----------------------------

count           = 0
var_year        = variable.groupby('time.year').mean('time')

fig    = plt.figure(figsize=[12,7])

for year in [2001,2008,2016]:
    count+=1
    ax     = fig.add_subplot(1,3, count, projection=ccrs.NorthPolarStereo())
    ax     = ax_features(ax)
    
    p      = ax.pcolormesh(longitude, latitude, var_year.sel(year=year),  
                           vmin = 0, 
                           vmax=np.max(levels),
#                            levels = levels,
                           transform = ccrs.PlateCarree())

    ax.set_title('{}'.format(year))
    
fig.suptitle(variable.long_name)
fig.subplots_adjust(bottom = 0.1)
cbar_ax = fig.add_axes([0.15, 0.06, 0.7, .03])   #left, bottom, width , height

cbar = fig.colorbar(p, cax=cbar_ax, orientation = 'horizontal', extend = 'max')
cbar.set_label(variable.units)
cbar.set_ticks(levels)

plt.savefig("/home/students/6252699/thesis/parcels2/figures/variables/{}_per_year.png".format(variable.name), bbox_inches="tight", pad_inches=0.1)

print("Done...")




