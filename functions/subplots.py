import matplotlib.pyplot as plt
import numpy as np
import cartopy
import cartopy.crs as ccrs
import lay_out


def set_subaxes(ax):
    ax.gridlines(xlocs = np.arange(-180,185,20), ylocs = np.arange(-90, 95,5), color='black', alpha=0.5, linestyle=':')
    ax.add_feature(cartopy.feature.OCEAN)
    ax.add_feature(cartopy.feature.LAND, facecolor = (0,0,0,.8), edgecolor = (0,0,0,0))
    ax.set_extent([-180,180,50,90],  ccrs.PlateCarree())
    return ax


lons, lats = np.linspace(0,360,1000), np.linspace(80,90,100)
X,Y        = np.meshgrid(lons,lats)
fake_data  = X + Y**2
vmin, vmax = 1, 1000

fig = plt.figure(figsize = [10,6])
fig.suptitle("Example plots", fontsize=20)

for i in range(6):
    ax = fig.add_subplot(2, 3, i+1, projection=ccrs.NorthPolarStereo())
    ax = set_subaxes(ax)
    
    ax.set_title("i = {}".format(i) )

    p    = ax.pcolormesh(lons, lats, fake_data*i,   
                         cmap = 'viridis',
                         vmin = 0, vmax = 40e3,
                         transform = ccrs.PlateCarree())

fig.subplots_adjust(bottom = 0.10)
cbar_ax = fig.add_axes([0.15, 0.06, 0.7, .03])
cbar = fig.colorbar(p, cax=cbar_ax, orientation = 'horizontal', extend = 'max')

plt.savefig("6subplots.png", bbox_inches="tight", pad_inches=0.1)
    