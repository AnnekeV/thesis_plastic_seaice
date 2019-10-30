import numpy as np
import netCDF4 as netcdf
import matplotlib as mpl
# mpl.use('Agg')
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
import cmocean
import gc
import os
import imageio
import xarray as xr
from matplotlib import rc
import re
import datetime
from functions.grid_to_km import area
import matplotlib.animation as animation
import cartopy.feature as cfeature
import os.path

fig_dir       = "/home/students/6252699/thesis/parcels2/figures/"
# data_dir      = "/science/projects/oceanparcels/CMEMS-GLORYS12V1-Arctic/"
# data_dir      = "/scratch/AnnekeV/reanalysis_data/"
data_dir      = "/data/oceanparcels/CMEMS-GLORYS12V1-Arctic/"

import matplotlib as mpl

################################################################
plt.style.use(['default'])

# plotting parameters
mpl.rcParams['font.size'] = 14
mpl.rcParams['font.weight'] = 'normal'

mpl.rcParams["figure.figsize"] = [8 , 4]
mpl.rcParams['legend.fontsize'] = 12      #legend size
#title params
mpl.rcParams['axes.titlesize'] = 20
mpl.rcParams['axes.titleweight'] = 'bold'
mpl.rcParams['font.style'] = 'italic'    #all fonts italic
#axes params
mpl.rcParams['axes.labelsize'] = 16
mpl.rcParams["xtick.labelsize"] = 13
mpl.rcParams["ytick.labelsize"] = 13
# line width and grid
mpl.rcParams['lines.linewidth'] = 2
mpl.rcParams["axes.grid"] = 1
# mpl.rcParams["figure.subplot.left"] = 0.05
# mpl.rcParams["figure.subplot.right"] = 0.1
# mpl.rcParams["figure.subplot.bottom"] = 0.11
mpl.rcParams["savefig.bbox"] = 'tight'
mpl.rcParams["savefig.pad_inches"] = 0.1
mpl.rcParams['legend.loc'] =  'best'


land_50m = cfeature.NaturalEarthFeature('physical', 'land', '50m',
                                        edgecolor='face',
                                        facecolor=cfeature.COLORS['land'])

ocean_50m = cfeature.NaturalEarthFeature('physical', 'ocean', '50m',
                                        edgecolor='face',
                                        facecolor="lightblue")
land_10m = cfeature.NaturalEarthFeature('physical', 'land', '10m',
                                        edgecolor='face',
                                        facecolor='lightgrey')


def set_fig():
    fig = plt.figure(figsize = [8,8])
    ax = plt.axes(projection=ccrs.NorthPolarStereo())
    ax.gridlines(xlocs = np.arange(-180,185,20), ylocs = np.arange(-90, 95,5), color='black', alpha=0.5, linestyle=':')
    ax.add_feature(cartopy.feature.OCEAN, facecolor = "black" )
    ax.add_feature(land_50m, facecolor='lightgrey')
#     ax.stock_img()
    ax.set_extent([-180,180,50,90],  ccrs.PlateCarree())
#     ax.coastlines()
  
    return fig,ax 

def set_axes():
    '''Returns ax with no coastlines'''
    plt.figure(figsize = [8,8])
    ax = plt.axes(projection=ccrs.NorthPolarStereo())
    ax.gridlines(xlocs = np.arange(-180,185,20), ylocs = np.arange(-90, 95,5), color='black', alpha=0.5, linestyle=':')
    ax.add_feature(cartopy.feature.OCEAN)
    ax.add_feature(land_50m)
    ax.set_extent([-180,180,50,90],  ccrs.PlateCarree())
#     ax.coastlines()
    ax.stock_img()
    return ax
    
    

def set_circular_boundary(ax):
    '''
    Compute a circle in axes coordinates, which we can use as a boundary
    for the map. We can pan/zoom as much as we like - the boundary will be
    permanently circular.
    '''
    import matplotlib.path as mpath
  
    theta = np.linspace(0, 2*np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(verts * radius + center)

    ax.set_boundary(circle, transform=ax.transAxes)     
    
    
    
def plot_extent(date, ax, no_xarray , cmap = cmocean.cm.ice_r): 
    '''
    Plots the filled  outline of the 15% extent ice on top of an existing ax.
    DATE in datetimeformat
    AX  axes object    
    '''
    ds     = xr.open_dataset('/scratch/AnnekeV/reanalysis_data/GLOBAL_REANALYSIS_PHY_001_030-TDS_{}{:02d}{:02d}_uv_uvice_con_thick.nc'.format(date.year, date.month, date.day))
    siconc = ds.siconc
    
    if no_xarray: 
        lats = ds.latitude.values
        lons = ds.longitude.values
        sicsic = siconc.where(siconc > 0.15).values
        ax.imshow(lons,lats, sicsic, transform = ccrs.PlateCarree(), cmap = cmap, 
                                                             levels = [0.15,1.1])
    else:
        siconc.where(siconc > 0.15).isel(time=0).plot.imshow(ax = ax, transform = ccrs.PlateCarree(), cmap = cmap, 
                                                             levels = [0.15,1.1], add_colorbar=False)

def plot_extent_contour(date, ax, cmap = "cool"): 
    '''
    Plots the contour outline of the 15% extent ice on top of an existing ax.
    DATE in datetimeformat
    AX  axes object    
    '''
    ds     = xr.open_dataset(data_dir+'GLOBAL_REANALYSIS_PHY_001_030-TDS_{}{:02d}{:02d}_uv_uvice_con_thick.nc'.format(date.year, date.month, date.day))
    siconc = ds.siconc
    sic    = siconc.isel(time=0).plot.contour(ax = ax, transform = ccrs.PlateCarree(), cmap = cmap, levels = [0.15,1.1], add_colorbar=False, linewidths=1.5)
    

def plot_trajectory_map(file_dir, dt_days = 5, extra_title = "", y_extent = 60, scatter=False, savefig = False, choose_less = False):
   
    '''
    plot the trajectories of particles on a map,  file_dir = the file with the trajectories   y_extent = the extent of the map, choose_less, if there are very many particles, reduce nr. trajectories plotted
    '''    
    test = netcdf.Dataset(file_dir)
    lat = test['lat'][:,:]
    lon  = test['lon'][:,:]
    if scatter: sic = test['sip'][:,:]
    match = re.search(r'\d{4}-\d{2}-\d{2}', test['time'].units)                     #find something date format
    start_date  = datetime.datetime.strptime(match.group(), '%Y-%m-%d').date()       #find start date from netcdf file
    test.close()   
        
    t      = np.size(lat,1)
    n_days = t*dt_days-dt_days
    n_part = np.size(lat,0)
 
    if choose_less:
        lat    = lat[::(n_part//1000)]
        lon    = lon[::(n_part//1000)]
        n_part = np.size(lat,0)
        linewidth = 1
        print("Number of particles plotted: {}".format(n_part))
    else:
        linewidth = .1

    
    fig, ax = set_fig()
    ax.add_feature(cartopy.feature.OCEAN)
    ax.add_feature(cartopy.feature.LAND)
    ax.set_extent([-180,180,60,90],  ccrs.PlateCarree())
    colors = mpl.cm.hsv(np.linspace(0, 1 , n_part))    # colors npart

    
    plot_extent(start_date, ax)



#     plt.scatter(lon[:,0], lat[:,0],transform = ccrs.Geodetic(), color = "blue")
    for i in range(n_part):
        ax.plot(lon[i,:], lat[i,:],
                 color=colors[i], 
                 linewidth = linewidth, 
                 transform=ccrs.Geodetic()
                 )
        if scatter:
            for t in range(np.size(lat,1)):
                ice_color = "black"
                if sic[i,t] > 0.15: ice_color = "snow"
                plt.scatter(lon[i,t], lat[i,t], color = ice_color, transform=ccrs.Geodetic() , zorder =10, s=1)
#         if (i%(n_part/20)==0): print("{:.0f}%".format(float(i)/float(n_part)*100))



    plt.title(extra_title + "\nSimulation for {} days \nand {} particles".format(n_days, n_part))
    if savefig: plt.savefig( fig_dir + "trajectory_simdays_{}_npart_{}_latmin_{:.0f}_{}_{}.png".format(n_days,n_part, np.min(lat[:,0]), start_date, extra_title))
    plt.show()

def plot_output_scatter(file_dir, dt_days = 5, extra_title = " ", satellite_height = 2.0e6, only_final=False):
    '''
    plot the trajectories of particles on a map
    file_dir = the file with the trajectories
    y_extent = the extent of the map
    '''

    test = netcdf.Dataset(file_dir + ".nc")
    lat = test['lat'][:,:]
    lon  = test['lon'][:,:]
    test.close()
    
    n_days = int((np.size(lat,1))*dt_days-dt_days)
    n_part = np.size(lat,0)
    n_timesteps = np.size(lat,1)
    colors = mpl.cm.gist_rainbow(np.linspace(0, 1 , n_timesteps))    # colors npart
    
    plt.figure(figsize = [8,8])

#     ax = plt.axes(projection=ccrs.NearsidePerspective(
#         central_latitude=90, satellite_height = satellite_height))
    ax = plt.axes(projection=ccrs.NorthPolarStereo())
#     ax.set_extent([-180,180,70,90],ccrs.PlateCarree())
    ax.gridlines(xlocs = np.arange(-180,185,20), ylocs = np.arange(0, 95,5))
    ax.set_global()

    if only_final: n_timesteps=[-1]
    for i in range(n_timesteps):
        plt.scatter(lon[:,i], lat[:,i],
#                  marker='o',
#                  linewidth = 1,
                 transform=ccrs.Geodetic(),
                 s = 1,
                 color = colors[i]
                 )
    
#         plt.plot(lon[i,0], lat[i,0],
#                  color=colors[i], 
#                  transform=ccrs.PlateCarree(),
#                  marker= 'X',
#                  label = "{:.2f},{:.2f}".format(lat[i,0], lon[i,0]),
#                  markersize = 17,
#                  linestyle = 'None'
#                  )
    plt.legend()
    ax.add_feature(cartopy.feature.OCEAN)
    ax.add_feature(cartopy.feature.LAND)

    ax.set_extent([-180,180,70,90],ccrs.PlateCarree())
    plt.title("Inital position \nSimulation for {} days ".format(n_days)+ extra_title)
    plt.savefig( fig_dir + "initial_position_simdays_{}_npart_{}_latmin_{:.0f}.png".format(n_days,n_part, np.min(lat[:,0])))
    plt.show()
    
def plot_single_variable(f_part, variable_name, dt_days = 5 , ylimit = False):
    '''Plot a variable against time
    f_part = file with the trajectories of particles
    variable_name = the desired variable'''
    
    DAY = 24*3600

    f_adv  = netcdf.Dataset(f_part + ".nc")    # file containing advected particles
    lat    = f_adv['lat'][:,:]
    lon    = f_adv['lon'][:,:]
    var    = f_adv[variable_name][:,:]
    time   = f_adv['time'][:,:]
    f_adv.close()
#     lat    = lat[::10,:]
#     lon    = lon[::10,:]
#     var    = var[::10,:]
#     time   = time[::10,:]
        
    n_part = int(np.size(var,0))        # number of particles
    n_sampling_points = np.size(var,1)  # number of output points
    n_days = n_sampling_points*dt_days -dt_days  # number of days
    colors = mpl.cm.gist_rainbow(np.linspace(0, 1 , n_part))    # colors npart

    print(n_part)
    print("size var", np.size(var))
    print("size time", np.size(time))
    
    plt.figure()
    plt.xlabel("Time (days)")
    plt.ylabel(variable_name) 
    plt.title(variable_name)
    plt.xticks(rotation=45)
    for i in range(n_part):  
        if lat[i,0] != lat[i,-1]:
            plt.plot(time[i,:]/DAY,var[i,:], color = colors[i] )
    
    if (ylimit == True): plt.ylim(-1e-16,2e-16)
#     plt.savefig( fig_dir + variable_name + "_simdays_{}_npart_{}_latmin_{:.0f}.png".format(n_days,n_part, np.min(lat[:,0])))

    plt.show()    
    


def make_background(datadir, variable, t, ax, fig):
    "make a background for the gif of sit, with changing time so changing files, datadir is where a background file is, variable is the background scalar feeld and t the background date. Ax and fig the figure that they should be added too " 
    
    counter = 0 
    while counter <20:
        background_filename = datadir + "ORCA0083-N06_{:d}{:02d}{:02d}d05I.nc".format(t.year,t.month, t.day)
        if os.path.isfile(background_filename): 
            counter +=20
        else:
            t -= datetime.timedelta(days = 1)
            counter +=1
            if counter == 19: print("is there an ice file missing? the spacing seems to be at least 19 days")
            
    ds = xr.open_dataset(background_filename)        # file containing ice fraction information    
    dsvariable = "ds." + variable
    
    
    min_variable = (eval(dsvariable).min())
    assert(min_variable > -999)
    
    cax, kw = mpl.colorbar.make_axes(ax,location='bottom',pad=0.03,shrink=0.8)

    '''Plot background ice fraction'''
    im = ax.pcolormesh(ds.nav_lon,
                       ds.nav_lat.fillna(-999),
                       eval(dsvariable)[0,:,:].fillna(-999),
                       cmap=cmocean.cm.matter,
                       vmin= min_variable, 
                       vmax= 5., 
                       transform=ccrs.PlateCarree(),
                       label = "{:d}/{:02d}/{:02d}".format(t.year,t.month, t.day), 
                       zorder = 2)
     
    '''Colorbar'''
    c_label = "Sea ice thickness (m)"   #colorbar lable
    im.cmap.set_under('white')
    cbar = fig.colorbar(im, cax=cax, extend='both', **kw)
    cbar.ax.tick_params(labelsize=14)
    cbar.set_label(c_label, size=16)
    
#     def box(text , x , y , ax , color):
    props = dict(boxstyle='round', facecolor='black', alpha=0.5)
    plt.text( 0.01, 0.95  , "{:d}/{:02d}/{:02d}".format(t.year,t.month, t.day), 
             transform = ax.transAxes, bbox = props ,size = 14 , style = "italic" , weight = "bold")    

    
    
    
# ========================================================================================================================================================================
# ANIMATIONS
# ========================================================================================================================================================================



    
def plot_gif(file_dir, dt_days = 5, extra_title = " ", with_background = False):
    '''
    plot the trajectories of particles on a map
    file_dir = the file with the trajectories
    y_extent = the extent of the map
    '''

    test = netcdf.Dataset(file_dir + ".nc")
    lat = test['lat'][:,:]
    lon  = test['lon'][:,:]
    time = test['time'][:,:]
    match = re.search(r'\d{4}-\d{2}-\d{2}', test['time'].units)                     #find something date format
    start_date  = datetime.datetime.strptime(match.group(), '%Y-%m-%d').date()       #find start date from netcdf file
    test.close()
    

    
    timesteps   = np.size(lat,1)
    n_days      = int(timesteps*dt_days-dt_days)
    n_part      = np.size(lat,0)
    colors      = mpl.cm.gist_rainbow(np.linspace(0, 1 , n_part))    # colors npart
    
    
    '''add all non stationary particles to a list'''
    non_stationary = []
    for i in range(n_part):
        if lat[i,0] != lat[i,-1]: non_stationary.append(i)
    print("number not staionary: {}".format(len(non_stationary)))
    
    files    = []
    
    for t in range(timesteps):
        fig = plt.figure(figsize = [8,8])
        ax = plt.axes(projection=ccrs.NorthPolarStereo())
        ax.gridlines(xlocs = np.arange(-180,185,20), ylocs = np.arange(0, 95,5), zorder = 4)
        ax.set_global()
        ax.add_feature(cartopy.feature.OCEAN, zorder = 1)
        ax.add_feature(cartopy.feature.LAND, zorder =3)
        ax.set_extent([-180,180,60,90],ccrs.PlateCarree())
        plt.title(extra_title + "\nt = {} days ".format(t*dt_days))
        
        for i in non_stationary:
            plt.scatter(lon[i,t], lat[i,t],
                     marker='o',
                     transform=ccrs.Geodetic(),
                     color = "black",
#                      color = colors[i], 
                     zorder= 5)      
            
        if with_background:
            background_time = start_date + datetime.timedelta(seconds=time[0,t])
            make_background('/data2/imau/oceanparcels/hydrodynamic_data/NEMO-MEDUSA/ORCA0083-N006/means/',"sit", background_time,  ax = ax, fig = fig)
        
        filename = (fig_dir + "gif/simdays_{}_npart_{}_latmin_{:.0f}_t_{:04d}.png".format(n_days,n_part, 
                                                                               np.min(lat[:,0]), t*dt_days) )
        files.append(filename)
        plt.savefig(filename)
        plt.close()
        print("{:.0f} %".format(float(t)/float(timesteps)*100.))
    
    print("Done with plotting, making gif")
    images = []
    for filename in files:
        images.append(imageio.imread(filename))
    imageio.mimsave(fig_dir + 'movement_particles_{}.gif'.format(extra_title), images , duration = 0.002)
    gc.collect() 
    print("Gif is ready in "+ fig_dir + 'movement_particles_{}.gif'.format(extra_title))
    
    for delete_filename in files:     # delete all images gifs made
        if os.path.isfile(delete_filename): os.unlink(delete_filename)

    


      
def plot_animation(fname, color_var="sic", label_name = "Sea ice concentration", save_gif = False, specific_savename = ""):
    '''
    Plot an animation for a certain output file FNAME from parcels. COLOR_VAR define variable on which color should be based. LABEL_NAME the label of the colorbar. In order to save the output, use save_gif = TRUE, SPECIFIC_SAVENAME is a  string to specify an extra variable          for the name when saved. 
    NEW VERSION OF PLOT GIF
    '''
    
    ds = xr.open_dataset(fname, decode_times=False)
    N  = np.size(ds.lat,1)
    start_date = find_start_date(ds.time.units)
    time       = start_date
    
    # Initialize figure 
    fig, ax = set_fig()
#     scat1 = ax.scatter(ds.lon.values[:,-1], ds.lat.values[:,-1], s =3, c=ds[color_var].values[:,-1] , cmap = cmocean.cm.ice)
    scat = ax.scatter(ds.lon.values[:,0], ds.lat.values[:,0], s =3, c=ds[color_var].values[:,0], transform = ccrs.PlateCarree() , cmap = cmocean.cm.ice)
    ttl  = ax.set_title("{:d}/{:02d}/{:02d}".format(time.year,time.month, time.day) )
    
    cbar = plt.colorbar(scat)
    cbar.set_label(label_name)
    
    print("Background is made")
    
    def update_loc(t):
        scat.set_offsets(np.vstack((ds.lon.values[:,t], ds.lat.values[:,t])).transpose())
        scat.set_array(ds[color_var].values[:,t])
        time_values = ds.time.values[:,t]   # in case there infinite time values
        time = start_date + datetime.timedelta(seconds=time_values[np.isfinite(time_values)][0])        
        ttl.set_text("{:d}/{:02d}/{:02d}".format(time.year,time.month, time.day) )
        return scat
    rc('animation', html='html5')
    anim = animation.FuncAnimation(fig, update_loc, frames=range(1,N), interval=300, blit=False)
    print("Animation is done")
    
    
    if save_gif: anim.save(fig_dir + 'animations/animation_startdate_{:d}{:02d}{:02d}_timesteps_{}_{}_{}.gif'.format(
        start_date.year,start_date.month, start_date.day, N, color_var, specific_savename), 
                           writer='imagemagick', fps=20)
    print("Done")
    plt.close()
    
    return anim


def find_glorys_ds(time, incl_coords=False):
    '''
    Takes time datetime object and boolean inc_coords (default is False)
    Returns siconc larger than 15%, and when incl_coords also lons and lats
    If specific file doesn't exist, uses first date before that does exist
    ''' 
       
    count = 0
    while count<100:
        fname = data_dir + 'GLOBAL_REANALYSIS_PHY_001_030-TDS_{}{:02d}{:02d}_uv_uvice_con_thick.nc'.format(time.year, time.month, time.day)
        if os.path.isfile(fname):
            dss = xr.open_dataset(fname)
            print(time)
            return dss
            count+=1000
        count+=1
        time-=datetime.timedelta(days=1.)
#     if incl_coords: return dss.siconc.where(dss.siconc>0.15).values, dss.longitude.values, dss.latitude.values
    return 0





def plot_animation_new(fname, save_gif = False, specific_savename = "", ice=True, every_x_day=10, color='gold', scatter_size=3, title="date", fps=40):
    '''
    Plot an animation for a certain output file FNAME from parcels only one color and with extent. COLOR_VAR define variable on which color should be based. LABEL_NAME the label of the colorbar. In order to save the output, use save_gif = TRUE, SPECIFIC_SAVENAME is a  string to specify an extra variable          for the name when saved. 
    NEW VERSION OF PLOT GIF
    '''
    
    ds         = xr.open_dataset(fname, decode_times=False).isel(obs=np.arange(0,2800,every_x_day))
    N          = np.size(ds.lat,1)
#     N=2800
    start_date = find_start_date(ds.time.units)
    time       = start_date
    cmap       = cmocean.cm.ice
    
    '''Initialize figure '''
    fig, ax = set_fig()
    scat   = ax.scatter(ds.lon.values[:,0], ds.lat.values[:,0], s =scatter_size, c=color, transform = ccrs.PlateCarree())
    if title == "date":      ttl    = ax.set_title("{:d}/{:02d}/{:02d}".format(time.year,time.month, time.day))
    if title == "timedelta": ttl = ax.set_title("0 days")
    set_circular_boundary(ax)
    '''Sea ice extent'''
    if ice:
        dssi   = find_glorys_ds(time)
        sicsic = dssi.where(dssi.siconc > 0.0, other=0.0).siconc.isel(time=0).plot.imshow(ax = ax, transform = ccrs.PlateCarree(), cmap = cmap, add_colorbar=False, alpha=1)    
    print("Background is made")
    
    def update_loc(t):
        '''Title'''
        time_values = ds.time.values[:,t]   # in case there infinite time values
        time_passed = datetime.timedelta(seconds=time_values[np.isfinite(time_values)][0]) 
        time        = start_date + time_passed
        if title == "date":
            ttl.set_text("{:d}/{:02d}/{:02d}".format(time.year,time.month, time.day) )
        if title == "timedelta":
            ttl.set_text("{} days".format(time_passed.days))
            
        '''Update location particles'''
        scat.set_offsets(np.vstack((ds.lon.values[:,t], ds.lat.values[:,t])).transpose())
        '''Update sea ice extent'''
        if ice:
            dssi = find_glorys_ds(time)
            sicsic = dssi.where(dssi.siconc > 0.0, 0.0).siconc.isel(time=0).plot.imshow(ax = ax, transform = ccrs.PlateCarree(), cmap = cmap, add_colorbar=False, alpha=1)    
        return scat
    rc('animation', html='html5')
    anim = animation.FuncAnimation(fig, update_loc, frames=range(1,N), interval=100, blit=False)
    print("Animation is done")
    
    
    if save_gif: anim.save(fig_dir + 'animations/animation_startdate_{:d}{:02d}{:02d}_timesteps_{}_{}_fps_{}.gif'.format(
        start_date.year,start_date.month, start_date.day, N, specific_savename, fps), 
                           writer='imagemagick', fps=fps)
    print("Done")
    plt.close()
    
    return anim



def plot_animation_10_per_day(fname, save_gif = False, specific_savename = "", ice=True, every_x_day=10, color='gold', scatter_size=3, title="date"):
    '''
    Plot an animation for a certain output file FNAME from parcels only one color and with extent. COLOR_VAR define variable on which color should be based. LABEL_NAME the label of the colorbar. In order to save the output, use save_gif = TRUE, SPECIFIC_SAVENAME is a  string to specify an extra variable          for the name when saved. 
    NEW VERSION OF PLOT GIF
    '''
    def shift_nans(ds):
        '''Shift the thing backward'''

        def shift_for_10(key, n=10, ndays=2800):
            '''
            Specifically for n releases at a time
            Default is n=10 particle release at a time
            For a simulation time of ndays=2800 days
            '''
            super_mat    = np.zeros([npart, npart//n+nobs])
            super_mat[:] = np.nan

            for i in range(npart):
                data                    = ds[key].values[i,:].copy()
                super_mat[i,i//n:i//n+nobs]   = data

            return super_mat.copy()[:,:ndays]

#         ds = xr.open_dataset(fn, decode_times=False)

        npart  = np.size(ds.time, 0)
        nobs   = np.size(ds.time, 1)

        time1 = shift_for_10(key='time')
        lon  = shift_for_10(key='lon')
        lat  = shift_for_10(key='lat')

        return time1, lon, lat    
    
    
    ds         = xr.open_dataset(fname, decode_times=False)
    start_date = find_start_date(ds.time.units)
    ds['time'].values[:], ds['lon'].values[:], ds['lat'].values[:] = shift_nans(ds)
    
    print start_date
    ds         = ds.isel(obs=np.arange(0,2800,every_x_day))
    N          = np.size(ds.lat,1)
    time       = start_date
    cmap       = cmocean.cm.ice_r
    
    
    '''Initialize figure '''
    fig, ax = set_fig()
    scat   = ax.scatter(ds.lon.values[:,0], ds.lat.values[:,0], s =scatter_size, c=color, transform = ccrs.PlateCarree())
    if title == "date":      ttl    = ax.set_title("{:d}/{:02d}/{:02d}".format(time.year,time.month, time.day))
    if title == "timedelta": ttl = ax.set_title("0 days")
    set_circular_boundary(ax)
    '''Sea ice extent'''
    if ice:
        dssi   = find_glorys_ds(time)
        sicsic = dssi.siconc.where(dssi.siconc > 0.15).isel(time=0)
        sicplot = ax.imshow(sicsic, transform = ccrs.PlateCarree(), cmap = cmap, levels = [0.15,1.1], add_colorbar=False, alpha=0.5)    
    print("Background is made")
    
    def update_loc(t):
        '''Title'''
        time_values = ds.time.values[:,t]   # in case there infinite time values
        time_passed = datetime.timedelta(seconds=time_values[np.isfinite(time_values)][0]) 
        time        = start_date + time_passed
        if title == "date":
            ttl.set_text("{:d}/{:02d}/{:02d}".format(time.year,time.month, time.day) )
        if title == "timedelta":
            ttl.set_text("{} days".format(time_passed.days))
            
        '''Update location particles'''
        scat.set_offsets(np.vstack((ds.lon.values[:,t], ds.lat.values[:,t])).transpose())
        '''Update sea ice extent'''
        if ice:
            dssi = find_glorys_ds(time)
            sicsic = dssi.siconc.where(dssi.siconc > 0.15).isel(time=0)
            sic_plot.set_array(sicsic)
#             ax = ax, transform = ccrs.PlateCarree(), cmap = cmap, levels = [0.15,1.1], add_colorbar=False, alpha=0.5)    
        return scat
    rc('animation', html='html5')
    anim = animation.FuncAnimation(fig, update_loc, frames=range(1,N), interval=100, blit=False)
    print("Animation is done")
    
    
    if save_gif: anim.save(fig_dir + 'animations/animation_startdate_10_per_day_{:d}{:02d}{:02d}_timesteps_{}_{}_{}.gif'.format(start_date.year,start_date.month, start_date.day, 
                                                                                                                                N, specific_savename, title), 
                           writer='imagemagick', fps=40)
    print("Done")
    plt.close()
    
    return anim

# ========================================================================================================================================================================
# HISTOGRAMS
# ========================================================================================================================================================================


def plot_histogram(file_dir, y_extent = 60, res =1, save = True, specific_savename = "", dt_days =5, final = False, tracer = True):
    '''plot a histogram, FILE_DIR should include .nc,  RES is resolution of the pcolormesh, Y_EXTENT is lower limit of figure in latitude, FINAL is a boolean: whether to make only a histogram of final timestep or whole period
    TRACER tells you whether the start was tracers or nr of particles, SPECIFIC_SAVENAME = name for saving'''
    '''Import data'''
    test        = netcdf.Dataset(file_dir)
    lat         = test['lat'][:,:]  
    lon         = test['lon'][:,:]
    match       = re.search(r'\d{4}-\d{2}-\d{2}', test['time'].units)                     #find something date format
    start_date  = datetime.datetime.strptime(match.group(), '%Y-%m-%d').date()       #find start date from netcdf file

    test.close()
    
    n_days      = int(np.size(lat,1)*dt_days-dt_days)
    n_part      = np.size(lat,0)

    if final:  H, xedges, yedges = np.histogram2d(lon[:,-1], lat[:,-1], bins=(np.arange(-180,360,res),np.arange(30,90,res)))
    else: H, xedges, yedges = np.histogram2d(lon.flatten(), lat.flatten(), bins=(np.arange(-180,360,res),np.arange(30,90,res)))

    
    X, Y     = np.meshgrid(xedges, yedges)
    lat_area = area((Y[:-1,:-1]+Y[1:,1:])/2., res)/1e6
    H        = H.T
    H        = H/lat_area # Let each row list bins with common y range.

    plt.figure(figsize = [8,8])
    ax = plt.axes(projection=ccrs.NorthPolarStereo())
    ax.gridlines(xlocs = np.arange(-180,185,20), ylocs = np.arange(-90, 95,5))
    ax.add_feature(cartopy.feature.OCEAN)
    ax.add_feature(cartopy.feature.LAND)
    ax.set_extent([-180,180,y_extent,90],  ccrs.PlateCarree())
    vmin, vmax = 1./np.max(lat_area), np.max(H)

    im = ax.pcolormesh(X, Y, H, transform = ccrs.PlateCarree(), 
                       norm = mpl.colors.LogNorm(vmin=0.5/np.max(lat_area),vmax=np.max(H)), 
                       cmap = 'spring_r'
                      )
    cbar = ax.figure.colorbar(im)
    if not tracer:
        print("Assume that start was with particles")
        ticks = [1e-4,1e-3,1e-2,0.1, 1, 5, 10]
        cbar.set_ticks(ticks)
        cbar.set_ticklabels(ticks)
        cbar.set_label(r"No. particles / $km^2$")

    plt.savefig( fig_dir + "histogram_glorys_simdays_{}_npart_{}_latmin_{:.0f}_{}_{}.png".format(n_days,n_part, np.min(lat[:,0]), start_date, specific_savename))

def plot_extent(date, ax): 
    "Use set_axes() for the ax and date use find_start_date(str(2014-09-01)"
    ds     = xr.open_dataset('/scratch/AnnekeV/reanalysis_data/GLOBAL_REANALYSIS_PHY_001_030-TDS_{}{:02d}{:02d}_uv_uvice_con_thick.nc'.format(date.year, date.month, date.day))
    siconc = ds.siconc
    siconc.where(siconc > 0.15).isel(time=0).plot.imshow(ax = ax, transform = ccrs.PlateCarree(), cmap = cmocean.cm.ice_r, levels = [0.15,1.1])
    
def find_start_date(time_units):
    match = re.search(r'\d{4}-\d{2}-\d{2}',time_units)                     #find something date format
    start_date  = datetime.datetime.strptime(match.group(), '%Y-%m-%d').date() 
    return start_date


    

          
                
def create_discrete_colorbar(lowerbound, upperbound, number_bins, sequential=True, divergent=False):
    '''
    Creates a colormap with discreet colors from lowerbound to upperbound for number_bins. Choose between sequential or divergent
    Returns cmap, bounds, norm
    call for plotfunction(cmap=cmap,norm=norm)
    plt.colorbar(cmap=cmap,norm=norm,ticks=bounds,boundaries=bounds)
    '''
    
    if sequential:
        cmap = plt.cm.viridis # define the colormap
    if divergent:
        cmap = plt.cm.PiYG
                
    # extract all colors from the cmap
    cmaplist = [cmap(i) for i in range(cmap.N)]

    # create the new map
    cmap = mpl.colors.LinearSegmentedColormap.from_list(
        'Custom cmap', cmaplist, cmap.N)

    # define the bins and normalize
    bounds = np.linspace(lowerbound, upperbound, number_bins, endpoint=True)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    return cmap, bounds, norm

def plot_amount_within_bin(T, latbins,lonbins):
    '''Plots the fraction of T that stays in the same bin'''
    plt.figure(figsize = [8,8])
    ax = plt.axes(projection=ccrs.NorthPolarStereo())
    ax.gridlines(xlocs = np.arange(-180,185,30), ylocs = np.arange(0,95,10), color='black', alpha=0.5, linestyle=':')
    ax.add_feature(cartopy.feature.OCEAN)
    ax.add_feature(cartopy.feature.LAND)
    ax.set_extent([-180,180,latbins.min(),latbins.max()],  ccrs.PlateCarree())
    ax.coastlines()
    plt.title("Amount not leaving bin\n res = {} deg\n".format(lonbins[2]-lonbins[1]))

    cmap, bounds, norm = create_discrete_colorbar(0.2, 1, 17)

    p = ax.pcolormesh(lonbins,latbins, np.diag(T)[:-1].reshape(len(latbins)-1,len(lonbins)-1), transform = ccrs.PlateCarree(),
                      cmap = cmap, norm = norm
                     )
    cb =plt.colorbar(p, cmap=cmap, norm=norm,
        spacing='proportional', ticks=bounds, boundaries=bounds, format='%.2f', extend = 'max')
    cb.set_label("Fraction")

    add_grid_labels(ax,-180,180,50,90,30,10, latlabels=False)
    set_circular_boundary(ax)
    
'''    
if __name__ == "__main__":
    output_name =  "/scratch/AnnekeV/output/" + "run_02_07_bering_npart_1500_start_2001-08-01_lat_66.0_66.5_simdays_2000_AdvectionRK4_ice"
    plot_gif(output_name, dt_days = 5, extra_title = "bering2000", with_background = True)

def add_grid_labels(ax, x0, x1, y0, y1, dx, dy,
                     latlabels=True, fontsize=10, alpha=.8):
    """Add grid line labels manually for projections that aren't supported

    Args:
        ax (geoaxes.GeoAxesSubplot) 
        x0 (scalar)
        x1 (scalar)
        y0 (scalar)
        y1 (scalar)
        dx (scalar)
        dy (scalar)
        latlabels (Boolean)
    """
    dtype = int
    
    if latlabels:
        for lat in np.arange(y0, y1, dy, dtype=dtype):
            text = ax.text(180, lat, r" {}$^o$N".format(lat), ha='left', va = 'bottom', transform = ccrs.PlateCarree(),fontsize=fontsize,alpha=alpha)

        
    for lon in np.arange(x0, x1, dx, dtype=dtype):
            if abs(lon) < 90:
                va = 'top'
            elif abs(lon) > 90:
                va = 'bottom'
            elif abs(lon) == 90:
                va = 'center'
                
            if (lon) < 0:
                ha = 'right'
            elif (lon) > 0:
                ha = 'left'
            elif (lon == 0 or abs(lon) == 180):
                ha = 'center'
            
            text = ax.text(lon,y0, r"{}$^o$E".format(lon), ha = ha, va = va, transform = ccrs.PlateCarree(), fontsize=fontsize, alpha=alpha)
                '''