import matplotlib as mpl
mpl.use('Agg')
import plot_functions
import glob
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np



ten_per_day  = sorted(glob.glob("/scratch/AnnekeV/output/single_particles/ten_per_day/*time.nc"))
sorted_names = ["bering_ocean", "bering_ice", "norway_ocean", "norway_ice"]

for i in range(4):
    ds = xr.open_dataset(ten_per_day[i], decode_times=False)

    fig = plt.figure(figsize = [12,12])
    ax  = fig.add_subplot(1,1,1, projection=ccrs.NorthPolarStereo())
    ax.add_feature(plot_functions.land_50m, facecolor ='lightgrey')
    ax.set_extent([-180,180,50,90],  ccrs.PlateCarree())
    ax.gridlines(xlocs = np.arange(-180,185,30), ylocs = np.arange(0,95,10), color='black', alpha=0.5, linestyle=':')
    plot_functions.set_circular_boundary(ax)


    npart = np.size(ds.lon, 0)
    p   = np.arange(0,npart)
    t   = np.arange(1700)

    sc       = ax.scatter(ds.lon.isel(traj=p, obs=t), ds.lat.isel(traj=p, obs=t), c=ds.time.isel(traj=p, obs=t)-ds.time.isel(traj=p, obs=0), s=.005, cmap='jet', transform=ccrs.Geodetic(), zorder=100, alpha=.5)
    sc_start = ax.scatter(ds.lon.isel(traj=0, obs=0), ds.lat.isel(traj=0, obs=0), s=300, marker='*', transform=ccrs.Geodetic(), zorder=101, color = 'white')


    ax_cbar = fig.add_axes([0.1, 0.05, 0.8, 0.01])
    cb      = fig.colorbar(sc, ax_cbar,  orientation='horizontal')
    ticks   = np.arange(0,1700,250)

    cb.set_ticks(ticks*86400)
    cb.set_ticklabels(ticks)
    cb.set_label("Time (days)")
  

    plt.savefig("/home/students/6252699/thesis/parcels2/figures/maps/traj_map/ten_per_day_4_years_{}.png".format(sorted_names[i]), bbox_inches="tight", pad_inches=0.1)
    
    



