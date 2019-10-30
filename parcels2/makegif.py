import matplotlib as mpl
mpl.use('Agg')
import plot_functions
import glob



from timeit import default_timer as timer


start = timer()

x=14
print("x = {}".format(x))


fn_bering_ice =  '/scratch/AnnekeV/output/single_particles/ten_per_day/bering-everyday-05-30_npart_10950_start_2009-02-01_simdays_2800_kernel_AdvectionRK4_prob_dtdays_1_10_at_a_time_shifted.nc'
plot_functions.plot_animation_new(fn_bering_ice, save_gif=True,specific_savename="7_be_ice_days_yes_ice_in_plot_shifted_range_x_{}".format(x), 
                                  ice = True,  every_x_day=x, scatter_size=1, 
                                  title = "date", fps=5)

print("Done with be ice")

fn_norway_ice  = "/scratch/AnnekeV/output/single_particles/ten_per_day/norway-everyday-05-30_npart_10950_start_2009-02-01_simdays_2800_kernel_AdvectionRK4_prob_dtdays_1_10_at_a_time_shifted.nc"
plot_functions.plot_animation_new(fn_norway_ice, save_gif=True,specific_savename="7_nor_ice_days_yes_ice_in_plot_shifted_range_x_{}".format(x), 
                                  ice = True,  every_x_day=x, scatter_size=1, 
                                  color='deeppink', title = "date", fps=5)


print("{:.0f} minutes".format((timer()-start)/60.))

# for simdays in [2800]:
#     fn_shifted =  glob.glob('/scratch/AnnekeV/output/single_particles/bering*simdays_{:d}*prob*shifted.nc'.format(simdays))[-1]
#     plot_functions.plot_animation_new(fn_shifted, save_gif=True,specific_savename="{:d}_days_ice".format(simdays), ice = True,  every_x_day=2)


# fn_shifted = "/scratch/AnnekeV/output/single_particles/bering-everyday-05-06_npart_1095_start_2009-02-01_simdays_2800_kernel_AdvectionRK4_prob_dtdays_1shifted.nc"
# fn_shifted =  glob.glob('/scratch/AnnekeV/output/single_particles/bering*2800*ocean*shifted.nc')[0]
# plot_functions.plot_animation_new(fn_shifted, save_gif=True,specific_savename="ocean", ice=False)

'''



import cartopy 
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr 



ymin, ymax = 50,80
xmin, xmax = -180,180
dx, dy = 10, 5



fn = fn_shifted
ds = xr.open_dataset(fn, decode_times=False)
npart = len(ds.traj)

ymin, ymax = 50,900
xmin, xmax = -180,180
dx, dy = 20, 5


plt.figure(figsize = [8,8])
ax = plt.axes(projection=ccrs.NearsidePerspective(central_latitude=90, satellite_height=3e6))
ax.gridlines(xlocs = np.arange(-180,185,dx), ylocs = np.arange(0,95,dy), color='black', alpha=0.5, linestyle=':')
ax.add_feature(cartopy.feature.OCEAN)
ax.add_feature(cartopy.feature.LAND)
ax.set_extent([xmin,xmax,ymin,ymax],  ccrs.PlateCarree())
# ax.coastlines()


for i in range(npart):
    ax.plot(ds.lon.isel(traj=i), ds.lat.isel(traj=i), transform=ccrs.Geodetic(), linewidth=.5)

plt.savefig("/home/students/6252699/thesis/parcels2/figures/trajectories/bering_daily_2800_days_ocean.png",bbox_inches="tight", pad_inches=0.1)
# plt.show()

'''