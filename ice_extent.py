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
    ax.add_feature(cartopy.feature.OCEAN, facecolor = "lightseagreen" )
    ax.add_feature(pf.land_50m)
    pf.set_circular_boundary(ax)
    ax.add_feature(cartopy.feature.LAND, zorder=1, facecolor = "lightgrey")
    ax.set_extent([-180,180, 50,90],  ccrs.PlateCarree())
    return ax


for year in [2001, 2016]:

    # -------- IMPORT FILES -----------------------------

    fnames    = sorted( glob("/data/oceanparcels/CMEMS-GLORYS12V1-Arctic/*_{}*.nc".format(year)))
    ds_glorys = xr.open_mfdataset(fnames)

    print("Imported files...")

    latitude  = ds_glorys.latitude
    longitude = ds_glorys.longitude

    sic_season = ds_glorys.siconc.fillna(0.).groupby('time.season').mean('time')
    sic_month  = ds_glorys.siconc.fillna(0.).groupby('time.month').mean('time')
    sic_year   = ds_glorys.siconc.fillna(0.).groupby('time.year').mean('time')

    print("Calculated means...")






    # -------- SAVE FILES -----------------------------

    data_dir = "/scratch/AnnekeV/data/" 

    # print("\nOcean velocities")

    # sic_season.to_netcdf(path=data_dir + "sic_no_nan_seasonal_mean_2001-2016")

    # print("Exported seasonal means...")

    # sic_month.to_netcdf(path=data_dir + "sic_no_nan_monthly_mean_2001-2016")
    # print("Exported monthly means...")

    # sic_year.to_netcdf(path=data_dir + "sic_no_nan_yearly_mean_2001-2016")
    # print("Exported yearly means...")

    sic_month.to_netcdf(path=data_dir + "sic_no_nan_monthly_mean_{}".format(year))

'''
# -------- PLOT EASTWARD VELOCITIES -----------------------------


maxv =.5
count = 0

for month in [6,7,8,9,10]:
    count+=1
   
    fig    = plt.figure(figsize=[15,15])
    ax     = plt.axes(projection=ccrs.NorthPolarStereo())
    ax     = ax_features(ax)

    p      = ax.pcolormesh(longitude, latitude, uo_month.sel(month=month), transform = ccrs.PlateCarree(), 
                           vmin = -maxv, vmax=maxv,
                           cmap = "seismic")

    ax.set_title('uo  {:02d}'.format(month))
    fig.subplots_adjust(bottom = 0.1)
    cbar_ax = fig.add_axes([0.15, 0.06, 0.7, .03])   #left, bottom, width , height
    cbar = fig.colorbar(p, cax=cbar_ax, orientation = 'horizontal', extend = 'both')
    cbar.set_label("m/s")
    cbar.set_ticks([0.5, 0.4, 0.2, 0, -0.2, -0.4, -0.5])

    plt.savefig("/home/students/6252699/thesis/parcels2/figures/variables/ocean_speed/uo_month_{}.png".format(month), bbox_inches="tight", pad_inches=0.1)


    print("Made figure {}... ".format(count))


'''

'''

# -------- PLOT EASTWARD VELOCITIES -----------------------------


maxv =.5

for season in uo.season:
   
    fig    = plt.figure(figsize=[12,12])
    ax     = plt.axes(projection=ccrs.NorthPolarStereo())
    ax     = ax_features(ax)

    p      = ax.pcolormesh(longitude, latitude, uo.sel(season=season), transform = ccrs.PlateCarree(), 
                           vmin = -maxv, vmax=maxv,
                           cmap = "seismic")

    ax.set_title('uo')
    fig.subplots_adjust(bottom = 0.1)
    cbar_ax = fig.add_axes([0.15, 0.06, 0.7, .03])   #left, bottom, width , height
    cbar = fig.colorbar(p, cax=cbar_ax, orientation = 'horizontal', extend = 'both')
    cbar.set_label("m/s")
    cbar.set_ticks([0.5, 0.4, 0.2, 0, -0.2, -0.4, -0.5])

    plt.savefig("/home/students/6252699/thesis/parcels2/figures/variables/ocean_speed/uo_{:s}.png".format(season.values), bbox_inches="tight", pad_inches=0.1)


    print("Made figure 1... ")


    # -------- PLOT NORTHWARD VELOCITIES -----------------------------


    fig = plt.figure(figsize=[12,12])
    ax  = plt.axes(projection=ccrs.NorthPolarStereo())
    ax  = ax_features(ax)

    p   = ax.pcolormesh(longitude, latitude, vo.sel(season=season), transform = ccrs.PlateCarree(), 
                        vmin = -maxv, vmax=maxv, 
                        cmap = "seismic")

    ax.set_title('vo')
    
    fig.subplots_adjust(bottom = 0.1)
    cbar_ax = fig.add_axes([0.15, 0.06, 0.7, .03])   #left, bottom, width , height
    cbar = fig.colorbar(p, cax=cbar_ax, orientation = 'horizontal', extend = 'both')
    cbar.set_label("m/s")
    cbar.set_ticks([0.5, 0.4, 0.2, 0, -0.2, -0.4, -0.5])

    plt.savefig("/home/students/6252699/thesis/parcels2/figures/variables/ocean_speed/vo_{:s}.png".format(season.values), bbox_inches="tight", pad_inches=0.1)



    print("Made figure 2...")



'''

# -------- PLOT QUIVER VELOCITIES -----------------------------

'''

spacing = 1
regrid_shape = 50

x=longitude[::spacing]
y=latitude[::spacing]
u=uo[::spacing, ::spacing]
v=vo[::spacing, ::spacing]

# V=u**2+v**2

X,Y = np.meshgrid(x,y)

print np.shape(u)

start = timer()

fig = plt.figure(figsize=[12,12])
ax = plt.axes(projection=ccrs.NorthPolarStereo())
ax.gridlines(xlocs = np.arange(-180,185,20), ylocs = np.arange(-90, 95,5), color='black', alpha=0.5, linestyle=':')
ax.add_feature(cartopy.feature.OCEAN, facecolor = "lightseagreen", zorder=1 )
ax.add_feature(pf.land_50m, zorder=3)
pf.set_circular_boundary(ax)
ax.add_feature(cartopy.feature.LAND, zorder=1)
ax.set_extent([-180,180, 50,90],  ccrs.PlateCarree())


p   = ax.pcolormesh(longitude, latitude, std, transform = ccrs.PlateCarree(),zorder=2)
                    
q = ax.quiver(X,Y, u.values, v.values,  
#               V.values, 
              transform = ccrs.PlateCarree(),  regrid_shape=regrid_shape, zorder =4)
              

ax.quiverkey(q, X=0.5, Y=1.1, U=0.5,
             label='0.5 m/s', labelpos='E')
             
plt.colorbar(p)      
             
plt.savefig("/home/students/6252699/thesis/parcels2/figures/quiver_mean_whole_north_onecolor.png", bbox_inches="tight", pad_inches=0.1)

print("{:.1f} s".format(timer()-start))
'''
'''
# -------- PLOT QUIVER VELOCITIES LOCAL -----------------------------

spacing = 1
regrid_shape = 50

x=longitude[::spacing]
y=latitude[::spacing]
u=uo[::spacing, ::spacing]
v=vo[::spacing, ::spacing]

# V=u**2+v**2

X,Y = np.meshgrid(x,y)

print np.shape(u)

start = timer()

fig = plt.figure(figsize=[12,12])
ax = plt.axes(projection=ccrs.PlateCarree())
ax.gridlines(xlocs = np.arange(-180,185,20), ylocs = np.arange(-90, 95,5), color='black', alpha=0.5, linestyle=':', draw_labels=True)
ax.add_feature(pf.ocean_50m, facecolor = "lightseagreen" )
ax.add_feature(pf.land_50m)
# pf.set_circular_boundary(ax)
ax.add_feature(cartopy.feature.LAND, zorder=1)
ax.set_extent([0,20, 60,75],  ccrs.PlateCarree())


q = ax.quiver(X,Y, u.values, v.values,  
              v.values, cmap= 'seismic', 
              transform = ccrs.PlateCarree(), 
              regrid_shape=regrid_shape, zorder =2)

ax.quiverkey(q, X=0.5, Y=1.1, U=0.5,
             label='0.5 m/s', labelpos='E')

lon = np.linspace(10,10.5, 10)
lat = np.repeat(67, 10)
                

ax.scatter(lon, lat, color = 'black', transform = ccrs.PlateCarree())



plt.savefig("/home/students/6252699/thesis/parcels2/figures/quiver_mean_norway_with_scatter.png", bbox_inches="tight", pad_inches=0.1)

print("{:.1f} s".format(timer()-start))


'''