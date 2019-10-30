import numpy as np
import matplotlib as mpl
# mpl.use('Agg')
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
import xarray as xr
import lay_out
import sys
import xarray as xr
import glob
import dask
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
import scipy.sparse.linalg as sla
import matplotlib.animation as animation
from matplotlib import rc
import functions.lay_out
import plot_functions as pf
import os


fig_dir = "/home/students/6252699/thesis/parcels2/figures/"
out_dir = "/scratch/AnnekeV/output/"


plt.style.use(['default'])

def set_axes(fig):
    ax = plt.axes(projection=ccrs.NorthPolarStereo())
    ax.gridlines(xlocs = np.arange(-180,185,20), ylocs = np.arange(-90, 95,5), color='black', alpha=0.5, linestyle=':')
    ax.add_feature(cartopy.feature.OCEAN)
    ax.add_feature(cartopy.feature.LAND)
    ax.add_feature(cartopy.feature.RIVERS)

    ax.set_extent([-180,180,50,90],  ccrs.PlateCarree())
    ax.coastlines(linewidth = 0.)
    ax.stock_img()
    return ax

def set_subaxes(ax):
    ax.gridlines(xlocs = np.arange(-180,185,20), ylocs = np.arange(-90, 95,5), color='black', alpha=0.5, linestyle=':')
    ax.add_feature(cartopy.feature.OCEAN)
    ax.add_feature(cartopy.feature.LAND)
    ax.add_feature(cartopy.feature.RIVERS, linewidth = 0.2, zorder = 0)
    ax.set_extent([-180,180,50,90],  ccrs.PlateCarree())
#     ax.coastlines(alpha =.2 , linewidth = .5)
    return ax



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

'''Make sure the grid looks nice''' 
def grid_subplots(nr_subplots):
    '''
    Return nr of rows (r) and columns (c) in the grid
    In a more optimal way
    '''
    
    r = int(np.sqrt(nr_subplots))
    c = int(np.sqrt(nr_subplots))
    
    while r*c < nr_subplots: 
        c+=1
        
    return r, c
  


def trans_mat_digit(ds_name, res, final_index =-1, in_ice = False, void = True, sink = True):
    '''
    Makes a transition matrix.
    ds_name = the name of the file (incl. direcotry)
    res     = resolution
    final_index in case you don't want it to be the last value
    ===============================================
    Returns Transition matrix, latbins, lonbins, sparse matrix
    '''
    
    ds = xr.open_dataset(ds_name, decode_times = False)

    N   = 360//res

    lons_i   = ds.lon.isel(obs= 0).values
    lats_i   = ds.lat.isel(obs= 0).values
    lons_f   = ds.lon.isel(obs= final_index).values
    lats_f   = ds.lat.isel(obs= final_index).values
    
    print("Opened dataset")
    
    nparts              = len(lons_i)
    particle_grid_index = np.zeros((2, nparts), dtype=int)
    ncels               = int ((360//res) * (40//res))

    if in_ice:
        in_ice_i = ds.in_ice.isel(obs=0).values  
        in_ice_f = ds.in_ice.isel(obs=final_index).values
        T        = np.zeros((ncels*2+1,ncels*2+1))
        
    else:
        in_ice_i = np.zeros(np.shape(lons_i))
        in_ice_f = np.zeros(np.shape(lons_i))
        T        = np.zeros((ncels+1,ncels+1))

    
    '''Find initial indices'''
    for i in range(nparts):
        if np.isnan(lons_i[i]) or np.isnan(lats_i[i]) :
            particle_grid_index[0,i] = -1
        else:
            particle_grid_index[0,i] = int((lats_i[i]-50)//res)*N  + int((lons_i[i]+ 180)//res) + int(in_ice_i[i]) * ncels

    '''Find final indices'''
    for i in range(nparts):
        if np.isnan(lons_f[i]) or np.isnan(lats_f[i]):
            particle_grid_index[1,i] = -1
        else:
            particle_grid_index[1,i] = int((lats_f[i]-50)//res)*N  + int((lons_f[i]+ 180)//res) + int(in_ice_f[i]) * ncels

    for i in range(nparts):
        T[particle_grid_index[1,i], particle_grid_index[0,i]] += 1

    
    print("Made transition matrix")
    '''Normalize matrix'''
    sumT =  np.sum(T, axis=0)
    sumT[sumT ==0] = 1
    T /=sumT

    '''Remove sinks'''
    if sink:
        max_sink = (np.diag(T)[:-1]).max()
        count    = 0
        while (max_sink == 1 and count <10):
            T        = remove_sinks(T)
            count   += 1
            max_sink = (np.diag(T)[:-1]).max()
        print ("Count removal for sinks is {}" .format(count))

    T[:,-1]  = 0   # all nan's go to nan
    T[-1,-1] = 1

    
    
    if void: 
        zeros = np.where(T.sum(axis=0)==0)[0]
        Tkopy = T[:,:]
        Tkopy[zeros,:] = 0 
        sumTkopy =  np.sum(Tkopy, axis=0)
        sumTkopy[sumTkopy ==0] = 1.
        Tkopy = Tkopy/sumTkopy
        T = Tkopy[:,:]
        
    latbins, lonbins= np.arange(50,90+res, res),  np.arange(-180,180+res, res)          
    M_sparse = csr_matrix(T)
    print("Done with T")

    return T, latbins , lonbins, M_sparse




def trans_mat_digit_mf(ds_names, res, final_index =-1, in_ice = False, void = True, sink = True):
    '''
    Makes a transition matrix.
    ds_names = the name of the multiple files (incl. direcotry)
    res     = resolution
    final_index in case you don't want it to be the last value
    ===============================================
    Returns Transition matrix, latbins, lonbins, sparse matrix
    '''
    
    nfiles = len(ds_names)
    ncels  = int ((360//res) * (40//res))
    superT = np.zeros((ncels+1,ncels+1))

    
    for f in range(nfiles):
        ds = xr.open_mfdataset(ds_names[f], decode_times = False)

        N   = 360//res

        lons_i   = ds.lon.isel(obs= 0).values
        lats_i   = ds.lat.isel(obs= 0).values
        lons_f   = ds.lon.isel(obs= final_index).values
        lats_f   = ds.lat.isel(obs= final_index).values

        print("Opened dataset")

        nparts              = len(lons_i)
        particle_grid_index = np.zeros((2, nparts), dtype=int)

        if in_ice:
            in_ice_i = ds.in_ice.isel(obs=0).values  
            in_ice_f = ds.in_ice.isel(obs=final_index).values
            T        = np.zeros((ncels*2+1,ncels*2+1))

        else:
            in_ice_i = np.zeros(np.shape(lons_i))
            in_ice_f = np.zeros(np.shape(lons_i))
            T        = np.zeros((ncels+1,ncels+1))


        '''Find initial indices'''
        for i in range(nparts):
            if np.isnan(lons_i[i]) or np.isnan(lats_i[i]) :
                particle_grid_index[0,i] = -1
            else:
                particle_grid_index[0,i] = int((lats_i[i]-50)//res)*N  + int((lons_i[i]+ 180)//res) + int(in_ice_i[i]) * ncels

        '''Find final indices'''
        for i in range(nparts):
            if np.isnan(lons_f[i]) or np.isnan(lats_f[i]):
                particle_grid_index[1,i] = -1
            else:
                particle_grid_index[1,i] = int((lats_f[i]-50)//res)*N  + int((lons_f[i]+ 180)//res) + int(in_ice_f[i]) * ncels

        for i in range(nparts):
            T[particle_grid_index[1,i], particle_grid_index[0,i]] += 1

        superT += T
        
    print("Made transition matrix")
    '''Normalize matrix'''
    sumT =  np.sum(superT, axis=0)
    sumT[sumT ==0] = 1
    superT /=sumT

    '''Remove sinks'''
    if sink:
        max_sink = (np.diag(superT)[:-1]).max()
        count    = 0
        while (max_sink == 1 and count <10):
            superT        = remove_sinks(superT)
            count   += 1
            max_sink = (np.diag(T)[:-1]).max()
        print ("Count removal for sinks is {}" .format(count))

    superT[:,-1]  = 0   # all nan's go to nan
    superT[-1,-1] = 1

    
    
    if void: 
        zeros = np.where(superT.sum(axis=0)==0)[0]
        Tkopy = superT.copy()[:,:]
        Tkopy[zeros,:] = 0 
        sumTkopy =  np.sum(Tkopy, axis=0)
        sumTkopy[sumTkopy ==0] = 1.
        Tkopy = Tkopy/sumTkopy
        superT = Tkopy.copy()[:,:]
        
    latbins, lonbins= np.arange(50,90+res, res),  np.arange(-180,180+res, res)          
    M_sparse = csr_matrix(superT)
    print("Done with T")

    return superT, latbins , lonbins, M_sparse


def movement(M,latbins, lonbins, grid_plot, times=5, in_ice = False, title = "",  savefig = False, save_name=""):
    '''
    Plot the results form the multiplication of the transitionmatir with the tracer
    T is the sparse eigenvector
    '''
    nlat     = len(latbins)-1
    nlon     = len(lonbins)-1
    ncels    = nlat*nlon
    res      = latbins[1]-latbins[0]

    R        = np.zeros((times, ncels+1))
    if in_ice: R = np.zeros((times, ncels*2+1))
    R[0,:]   = 1.

    for t in range(1,times):
         R[t,:] = M.dot(R[t-1,:])
            
        
    line  = sorted(R[R!=0])
    vmax  = line[int(len(line)*0.95)]  # 95 percentil
    vmax  = 2.
    vmin  = line[int(len(line)*0.01)]  # 1 percentile 
    vmin  = 1e-3
    ticks = np.linspace(vmin,vmax,7)



    

#     for t in show_plot_on_t:

    if grid_plot:
        '''Initialize figure raster of 6'''
        plot_count = 6
        count = 0
        fig   = plt.figure(figsize = [plot_count/2*5, 10])
        fig.suptitle(title, fontsize=20)

        for t in [1,times//5, times//5*2, times//5*3, times//5*4, times-1]:
            count += 1
            ax = fig.add_subplot(2, 3, count, projection=ccrs.NorthPolarStereo())
            ax = set_subaxes(ax) 

            p    = ax.pcolormesh(lonbins, latbins, R[t,0:ncels].reshape(nlat,nlon), transform = ccrs.PlateCarree(),  
                                 norm = mpl.colors.LogNorm(vmin = vmin, vmax = vmax), 
                                 cmap = 'viridis')
            ax.set_title("t = {}".format(t) )

            fig.subplots_adjust(bottom = 0.10)
            cbar_ax = fig.add_axes([0.15, 0.06, 0.7, .03])
            cbar = fig.colorbar(p, cax=cbar_ax, orientation = 'horizontal', extend = 'max')
            if savefig: plt.savefig(fig_dir + "/transitionplots/{}res_{}_t_{}.png".format(save_name,res, t))
                
    
    elif in_ice:
        
        '''Initialize figure  plot seperate ocean and ice'''
        fig   = plt.figure(figsize = [8, 4])
        fig.suptitle(title, fontsize=20)
        t     = times-1
        
        '''Ocean'''
        x_coor, y_coor = (lonbins[:-1]+lonbins[1:])/2., (latbins[:-1]+latbins[1:])/2.
        ax = fig.add_subplot(1,2,1, projection=ccrs.NorthPolarStereo())
        ax = set_subaxes(ax) 
        p    = ax.pcolormesh(lonbins, latbins, R[t,0:ncels].reshape(nlat,nlon), transform = ccrs.PlateCarree(),  
                             norm = mpl.colors.LogNorm(vmin = vmin, vmax = vmax), cmap = 'viridis')
        ax.contour(x_coor, y_coor, R[t,0:ncels].reshape(nlat,nlon), transform = ccrs.PlateCarree(), levels = [2], cmap = 'Greys_r', linewidths = 1)
        ax.set_title(" Ocean")
        
        '''Ice'''
        ax2 = fig.add_subplot(1,2,2, projection=ccrs.NorthPolarStereo())
        ax2 = set_subaxes(ax2) 
        p    = ax2.pcolormesh(lonbins, latbins, R[t,ncels:-1].reshape(nlat,nlon), transform = ccrs.PlateCarree(),  
                             norm = mpl.colors.LogNorm(vmin = vmin, vmax = vmax), cmap = 'viridis')
        ax2.contour(x_coor, y_coor, R[t,ncels:-1].reshape(nlat,nlon), transform = ccrs.PlateCarree(), levels = [2], cmap = 'Greys_r', linewidths = 1)
        ax2.set_title(" Ice")
        
        '''Colorbar'''
        fig.subplots_adjust(bottom = 0.10)
        cbar_ax = fig.add_axes([0.15, 0.06, 0.7, .03])
        cbar = fig.colorbar(p, cax=cbar_ax, orientation = 'horizontal', extend = 'max')
        
        
        
        
        

    else: 
        '''Initialize figure single plot'''
        fig   = plt.figure(figsize = [8, 9])
        t     = times-1
        x_coor, y_coor = (lonbins[:-1]+lonbins[1:])/2., (latbins[:-1]+latbins[1:])/2.
        ax = fig.add_subplot(1,1,1, projection=ccrs.NorthPolarStereo())
        ax = set_subaxes(ax) 
        p    = ax.pcolormesh(lonbins, latbins, R[t,0:ncels].reshape(nlat,nlon), transform = ccrs.PlateCarree(),  
                             norm = mpl.colors.LogNorm(vmin = vmin, vmax = vmax), cmap = 'viridis')
        ax.contour(x_coor, y_coor, R[t,0:ncels].reshape(nlat,nlon), transform = ccrs.PlateCarree(), levels = [2], cmap = 'Greys_r', linewidths = 1)
        fig.suptitle(title, fontsize = 20)
#             fig.subplots_adjust(bottom = 0.10)
#             cbar_ax = fig.add_axes([0.15, 0.06, 0.7, .03])
#             cbar = fig.colorbar(p, cax=cbar_ax, orientation = 'horizontal', extend = 'max')

    return ax

        


            

def box(text , x , y , ax , color= 'Grey'):
    '''Makes a box, anchored to the top left'''
    props = dict(boxstyle='round', facecolor=color, alpha=0.5)
    plt.text( x, y, text, 
             transform = ax.transAxes, bbox = props,
             size = 8 , style = "italic" , weight = "bold",
             horizontalalignment='left', 
             verticalalignment='top')   
        
       
    
def eigen_vector(M,  latbins, lonbins, title,  k=10,  savefig = False, save_name="", in_ice = False):
    '''
    Calculates from sparse transition matrix M the eigenvectors
    
    M = sparse matrix
    latbins, lonbins. The edges of the longitudes and latitudes
    title = title of the matrix
    k = number of eigenvalues calculated
    '''
    
    nlat, nlon = len(latbins)-1, len(lonbins)-1
    ncels      = nlat*nlon
    eig_val_l, eig_vec_l = sla.eigs(M.transpose(), k = k, which = 'LM')
    eig_val_r, eig_vec_r = sla.eigs(M, k = k, which = 'LM')

    max_color = np.round(np.max([abs(eig_vec_l), abs(eig_vec_r)])*10, decimals=5)
    max_color = 0.005
    labels    = np.round(np.linspace(-max_color, max_color, 5), decimals=6)

    for i in range(k):
        fig = plt.figure(figsize = [10,6])
        fig.suptitle("{}\n".format(title) + r"$\lambda = $ {} , {}".format(eig_val_l[i], eig_val_r[i]), fontsize=16)
        ax1 = fig.add_subplot(1, 2, 1, projection=ccrs.NorthPolarStereo())
        ax1 = set_subaxes(ax1)    
        p    = ax1.pcolormesh(lonbins, latbins, eig_vec_l.real[:ncels, i].reshape(nlat,nlon),  
                             vmin=-max_color, vmax = max_color,
                             cmap = 'PiYG',
                             transform = ccrs.PlateCarree())
#         cbar = plt.colorbar(p, extend = 'both')
#         cbar.set_ticks(labels)
#         cbar.set_ticklabels(labels)
        ax1.set_title("left")
        ax1.coastlines(linewidth=.5, alpha = 0.8)

        
        ax2 = fig.add_subplot(1, 2, 2, projection=ccrs.NorthPolarStereo())
        ax2 = set_subaxes(ax2) 
        ax2.pcolormesh(lonbins, latbins, eig_vec_r.real[:ncels, i].reshape(nlat,nlon),  
                             vmin=-max_color, vmax = max_color,
                             cmap = 'PiYG',
                             transform = ccrs.PlateCarree())

        ax2.set_title("right" )
        ax2.coastlines(linewidth=.5, alpha = 0.8)
        fig.subplots_adjust(right=0.9)
        cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
        cbar = fig.colorbar(p, cax=cbar_ax, extend = 'both')
        cbar.set_ticks(labels)
        cbar.set_ticklabels(labels)
#         cbar = plt.colorbar(p, extend = 'both')
        box("{:.2f}".format(eig_vec_r[-1,i] ), .1 , .9 , ax2 )

        
        if savefig: plt.savefig(fig_dir + "/transitionplots/eigen/{}_{}_res_{}.png".format(save_name, i,res))


        if in_ice: 
            ax   = set_axes()
            p    = ax.pcolormesh(lonbins, latbins, eig_vec.real[ncels:-1,i].reshape(nlat,nlon), 
        #                          norm = mpl.colors.LogNorm(vmin = vmin, vmax = vmax), cmap = 'viridis',
                                 vmin=-max_color, vmax = max_color,  
                                 cmap = 'PiYG',
                                 transform = ccrs.PlateCarree() )
            plt.colorbar(p, extend = 'both')
            ax.set_title(r"Ice: $\lambda = $ {}".format(eig_val[i]) )
            cbar.set_ticks(labels)
            cbar.set_ticklabels(labels)
            if save_fig: plt.savefig(fig_dir + "/transitionplots/eigen/{}_left_separate_ice_res_{}.png".format(i, res))
        plt.show()

        
        



def trans_mat(res, lon_init, lon_final, lat_init, lat_final, nosinks = True):
    '''Fill First vector and transition matrix. Res is the resolution, first x than y, or scalar. lon init the inital longitudes, lat init the initial latitudes etcetera
    Returns M, latbins , lonbins'''
    
    if np.isscalar(res): res = [res, res]
    
    '''Set parameters'''     
    nparticles = len(lat_init)

    latbins = np.arange(50,90+res[1], res[1])
    lonbins = np.arange(-180,180+res[0],res[0])  
    nlat = len(latbins)-1
    nlon = len(lonbins)-1
    ncel = nlat*nlon


    M    = np.zeros((ncel+1,ncel+1))

    for n in range(nparticles):
        if np.isfinite(lon_init[n]) and np.isfinite(lat_init[n]):
            j   = np.argmax(np.histogram2d([lon_init[n]],[lat_init[n]], bins = [lonbins, latbins])[0].T.flatten())

            if np.isfinite(lon_final[n]) and np.isfinite(lat_final[n]):   
                '''If first and final value are not nan'''
                i   = np.argmax(np.histogram2d([lon_final[n]],[lat_final[n]], bins = [lonbins, latbins])[0].T.flatten())
                M[i,j] += 1.

            else:
                M[ncel,j]+=1.  #              '''If only final value is a nan'''

    M[ncel, ncel] = 1    # all nans go to nan


    Msum = np.sum(M, axis=0)
    Msum[Msum==0] =1
    M/=Msum
    
    if nosinks:
        max_sink = (np.diag(M)[:-1]).max()
        count    = 0
        while (max_sink == 1 and count <10):
            M        = remove_sinks(M)
            count   += 1
            max_sink = (np.diag(M)[:-1]).max()
        
        print ("Count is ", count)
            
    return M, latbins , lonbins





def get_vectors_in_time(M,  latbins, lonbins, times = 5):
    '''
    Multiply tracer 1 with the transition matrix to find the distribution after x times.
    M = the transition matrix, LATBINS, LONBINS, arrays with the edges of the bins, TIMES how often you want to apply the transition matrix
    Returns R, vector with tracer distirbution
    '''
    
    nlat     = len(latbins)-1
    nlon     = len(lonbins)-1

    R        = np.zeros((times, nlat*nlon+1))
    R[0,:]   = 1.
    M_sparse = csr_matrix(M)

    for t in range(1,times):
         R[t,:] = M_sparse.dot(R[t-1,:])
            
    return R


def plot_R_trans(R, latbins, lonbins, savefig = False, savetitle = ""):
    '''
    Plots all the figures at every timestep
    '''
    
    line = sorted(R[R!=0])
    vmax = line[int(len(line)*0.95)]  # 95 percentil
    vmin = line[int(len(line)*0.01)]  # 1 percentile 
    print (vmin, vmax)
    nlat, nlon = len(latbins)-1, len(lonbins)-1
    
    times = np.size(R, 0)
    
    for t in range(times):
        ax   = set_axes()
        p    = ax.pcolormesh(lonbins, latbins, R[t,:-1].reshape(nlat,nlon), transform = ccrs.PlateCarree(),  
                             norm = mpl.colors.LogNorm(vmin = vmin, vmax = vmax), cmap = 'viridis')
        plt.colorbar(p)
        ax.set_title("t = {}".format(t) )
        if savefig:  plt.savefig(fig_dir + "/transitionplots/{}_t_.png".format(savetitle, t))
        else: plt.show()


def figure_where_leaving_domain(M, latbins, lonbins, savefig = False, savetitle = ""):           
    '''
    This functions where particles/tracers are leaving the domain (based on where to go from an actual lan/lot value to nan) and a plot with a line per longitude
    '''
    
    "MAP"
    ax   = set_axes()
    nlat, nlon = len(latbins)-1, len(lonbins)-1
    p    = ax.pcolormesh(lonbins, latbins, M[-1,:-1].reshape(nlat,nlon), transform = ccrs.PlateCarree())#,  norm = mpl.colors.LogNorm(vmin = vmin, vmax = 1))
    plt.title("Final row: \nFrom value to nan")
    plt.colorbar(p)
    if savefig: plt.savefig(fig_dir + "/transitionplots/{}_map.png".format(savetitle))
    else: plt.show()
    
    
    "GRAPH PER LONGITUDE"
    atlantic       = [np.argmin(abs(lonbins+100)), np.argmin(abs(lonbins-80))] #atlantic goes from from -100 to 80 E
    val_to_nan     =  M[-1,:-1].reshape(nlat,nlon)
    sum_per_lon    =  np.sum(val_to_nan/nlat/nlon*100, axis = 0)
    atlantic_total = np.sum(sum_per_lon[atlantic[0]:atlantic[1]])
    pacific_total  = np.sum(sum_per_lon)-atlantic_total
    
    plt.plot((lonbins[1:]+lonbins[:-1])/2.,sum_per_lon)
    plt.xlabel("Longitude (degrees east)")
    plt.ylabel("Percentage out of domain")
    plt.title("Atlantic total = {}\nPacific total  = {}".format( atlantic_total, pacific_total))
    
    if savefig: plt.savefig(fig_dir + "/transitionplots/{}_per_lon.png".format(savetitle))
    else:  plt.show()



if __name__ == "__main__":
    
    print("Running transition_matrices.py\n")
    print("------------Example------------\n")
    print("Kernel sic test 1000 days, but 60 days transition matrix\n")
    ds = xr.open_dataset("/scratch/AnnekeV/output/kernel_test_run_02_28_npart_11386_start_2014-02-01_simdays_1000_AdvectionRK4_ice_sic.nc", decode_times=False)

    lat_init  = ds.lat.isel(obs = 0).values
    lon_init  = ds.lon.isel(obs = 0).values
    lat_final = ds.lat.isel(obs = 13).values   # after two months
    lon_final = ds.lon.isel(obs = 13).values   # after two months

    res  = 2
    M, latbins, lonbins = trans_mat(res, lon_init, lon_final, lat_init, lat_final)
    figure_where_leaving_domain(M, latbins, lonbins, savefig = False, savetitle = "")
    R                   = get_vectors_in_time(M,  latbins, lonbins, times = 5)
    plot_R_trans(R, latbins, lonbins, savefig = False, savetitle = "kernel_sic_test_60_days")



def grid_subplots(nr_plots):
    '''Calculate best distribution for grid subplots'''
    sqrt    = np.sqrt(nr_plots)
    rows    = int(sqrt)
    columns = int(np.ceil(sqrt))
    while rows*columns<nr_plots:
        columns+=1
    return rows,columns


    
def movement_strait(M,latbins, lonbins, days, year, strait ="Bering", times=7, nr_subplots=6, ice_extent=True):
    '''
    Plot the results form the multiplication of the transition matrix with the tracer, for 6 times. for a particular entrance of the arctic ocean.
    Either choose the "Bering"  strait or the "Iceland" or "Norway" strait
    M is the sparse matrix as input
    latbins, lonbins also from trans_mat_digit func. 
    Times is number of times the multipliction needs to be done
    Bering strait 7 times is default
    grid_subplots = gives how many subplots should be given. 
    '''
    
    '''variables'''
    nlat     = len(latbins)-1
    nlon     = len(lonbins)-1
    ncels    = nlat*nlon
    res      = latbins[1]-latbins[0]
    N        = 360//res
    

    if strait == "Bering":
        lons = np.linspace(-172,-164,20)
        lats = np.linspace(65,66,10)
        
    elif strait == "Iceland":
        lons = np.linspace(-18,0,100)
        lats = np.linspace(65,66,10)
        
    elif strait == "Norway":
        lons = np.linspace(0,15,100)
        lats = np.linspace(60,61,10)
    else: 
        print("Choose a strait, either bering or Iceland to Norway")
        return 0
    
    lons, lats                 = np.meshgrid(lons,lats)
    particle_grid_index        =  ( (lats.flatten()-50)//res*N + (lons.flatten()+180)//res ).astype(int)   #digitize the lons and lats into position in vector
    R                          = np.zeros((times, ncels+1))
    R[0,particle_grid_index]   = 1.


    '''Integrate in time'''
    for t in range(1,times):
         R[t,:] = M.dot(R[t-1,:])
            
        
    vmax  = 1
    vmin  = 1e-6
     
        
    '''Make sure the grid looks nice''' 
    r,c = grid_subplots(nr_subplots)   #nr of rows (r) and columns (c) in the grid
    ax_width = 4
    
    '''Initialize figure raster of 6'''
    count = 0
    fig   = plt.figure(figsize = [c*ax_width, r*ax_width])
    fig.suptitle(strait + ", resolution = {} deg".format(res), fontsize=20)

    for t in np.linspace(0, times-1, nr_subplots, endpoint=True, dtype = int):
        count += 1
        ax = fig.add_subplot(r, c, count, projection=ccrs.NorthPolarStereo())
        ax = set_subaxes(ax) 
        ax.add_feature(cartopy.feature.LAND, facecolor = (0,0,0,.8), edgecolor = (0,0,0,0) )


        p    = ax.pcolormesh(lonbins, latbins, R[t,0:ncels].reshape(nlat,nlon),   
                             norm = mpl.colors.LogNorm(vmin = vmin, vmax = vmax), 
                             cmap = 'viridis',
                            transform = ccrs.PlateCarree())
        fig.subplots_adjust(bottom = 0.10)
        cbar_ax = fig.add_axes([0.15, 0.06, 0.7, .03])
        cbar = fig.colorbar(p, cax=cbar_ax, orientation = 'horizontal', extend = 'max')
        if ice_extent:        
            extentdate = datetime.datetime.strptime("{}-03-01".format(year), "%Y-%m-%d").date()
            pf.plot_extent_contour(extentdate, ax)
        ax.set_title("t = {} days".format(t*days) )



        
    
def final_R(M,latbins, lonbins, days, year, strait ="Bering", times=7, nr_subplots=6, ice_extent=True):
    '''
    Plot the FINAL results of the multiplication of the transition matrix with the tracer,for a particular entrance of the arctic ocean.
    Either choose the "Bering"  strait or the "Iceland" or "Norway" strait
    M is the sparse matrix as input
    latbins, lonbins also from trans_mat_digit func. 
    Times is number of times the multipliction needs to be done
    Bering strait 7 times is default
    grid_subplots = gives how many subplots should be given. 
    '''
    
    '''variables'''
    nlat     = len(latbins)-1
    nlon     = len(lonbins)-1
    ncels    = nlat*nlon
    res      = latbins[1]-latbins[0]
    N        = 360//res
    

    if strait == "Bering":
        lons = np.linspace(-172,-164,20)
        lats = np.linspace(65,66,10)
        
    elif strait == "Iceland":
        lons = np.linspace(-18,0,100)
        lats = np.linspace(65,66,10)
        
    elif strait == "Norway":
        lons = np.linspace(0,12,100)
        lats = np.linspace(60,61,10)    
        
    elif strait == "Everywhere":
        print("")
    else: 
        print("Choose a strait, either bering or Iceland to Norway")
        return 0
    
    
    "Fill initial R"
    R  = np.zeros((times, ncels+1))
    if strait == "Everywhere":
        R[0,:] = 1.
    else:
        lons, lats                 = np.meshgrid(lons,lats)
        particle_grid_index        =  ( (lats.flatten()-50)//res*N + (lons.flatten()+180)//res ).astype(int)   #digitize the lons and lats into position in vector
        R[0,particle_grid_index]   = 1.

        

        
    '''Integrate in time'''
    for t in range(1,times):
         R[t,:] = M.dot(R[t-1,:])
    
    return R


    
def open_sic(date):
    date   = datetime.datetime.strptime(date, "%Y-%m-%d").date()
    ds     = xr.open_dataset('/scratch/AnnekeV/reanalysis_data/GLOBAL_REANALYSIS_PHY_001_030-TDS_{}{:02d}{:02d}_uv_uvice_con_thick.nc'.format(date.year, date.month, date.day))
    siconc = ds.siconc.isel(time=0)
    return siconc

def right_eigen_vector_grid(M,  latbins, lonbins, title, date,  k=10,  savefig = False, save_name="", in_ice = False):
    '''
    Calculates from sparse transition matrix M the k biggest right eigenvectors
    in a grid
    
    M = sparse matrix
    latbins, lonbins. The edges of the longitudes and latitudes
    title = title of the matrix
    k = number of eigenvalues calculated
    date is a string '2014-02-01' to calculate the extent
    '''
    
    
    nlat, nlon   = len(latbins)-1, len(lonbins)-1
    ncels        = nlat*nlon
    eig_val_r, eig_vec_r = sla.eigs(M, k = k, which = 'LM')
    idx          = eig_val_r.argsort()[::-1]   
    eig_val_r    = eig_val_r[idx]
    eig_vec_r    = eig_vec_r[:,idx]
    
    print("Calculated right eigenvectors, currently plotting")

    max_color    = 0.0001
    labels       = np.round(np.linspace(-max_color, max_color, 5), decimals=6)
    axwidth      = 4
    r, c         = grid_subplots(k)
    
    count = 0
    fig   = plt.figure(figsize = [c*axwidth, r*axwidth])
    fig.suptitle(title, fontsize=20)
    sic   = open_sic(date=date)

    for i in range(k):
        count += 1
        ax = fig.add_subplot(r, c, count, projection=ccrs.NorthPolarStereo())
        ax = set_subaxes(ax) 
        ax.coastlines(linewidth=.5, alpha = 0.8)
      
        p    = ax.pcolormesh(lonbins, latbins, np.real(eig_vec_r[:ncels, i].reshape(nlat,nlon)),  
                         vmin=-max_color, vmax = max_color,
                         cmap = 'PiYG',
                         transform = ccrs.PlateCarree())
        sic.plot.contour(ax = ax, transform = ccrs.PlateCarree(), cmap ="binary_r", levels = [0.15,1.1], add_colorbar=False, linewidths=1)
        ax.set_title(r"$\lambda = $ {:.3f}".format( eig_val_r[i]), fontsize =10)
        box("{:.2f}".format(abs(eig_vec_r[-1,i]) ), .05 , .95 , ax)
        pf.set_circular_boundary(ax)
        ax.add_feature(cartopy.feature.LAND, facecolor = (0,0,0,.5), edgecolor = (0,0,0,0), zorder = 10)


    fig.subplots_adjust(bottom=0.1, wspace=.1, hspace=.15)
    cbar_ax = fig.add_axes([0.15, 0.06, 0.7, 0.03/r])   #left, bottom, width , height
    cbar = fig.colorbar(p, cax=cbar_ax, orientation = 'horizontal', extend = 'both')
    cbar.set_ticks(labels)
    cbar.set_ticklabels(labels)
    

def right_eigen_vector_single(M,  latbins, lonbins, date, res, year, days,  k=10,  savefig = False, colorbar=True):
    '''
    Calculates from sparse transition matrix M the k biggest right eigenvectors and plots them individually
    
    M = sparse matrix
    latbins, lonbins. The edges of the longitudes and latitudes
    title = title of the matrix
    k = number of eigenvalues calculated
    date is a string '2014-02-01' to calculate the extent
    '''
    
    
    nlat, nlon   = len(latbins)-1, len(lonbins)-1
    ncels        = nlat*nlon
    eig_val_r, eig_vec_r = sla.eigs(M, k = k, which = 'LM')
    idx          = eig_val_r.argsort()[::-1]   
    eig_val_r    = eig_val_r[idx]
    eig_vec_r    = eig_vec_r[:,idx]
    
    print("Calculated right eigenvectors, currently plotting")

    max_color    = 0.0001
    labels       = np.round(np.linspace(-max_color, max_color, 5), decimals=6)
    axwidth      = 4
    
    count = 0
    sic   = open_sic(date=date)

    for i in range(k):
        
        fig   = plt.figure(figsize = [axwidth, axwidth])
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.NorthPolarStereo())
        ax = set_subaxes(ax) 
#         ax.coastlines(linewidth=.5, alpha = 0.8)
      
        p    = ax.pcolormesh(lonbins, latbins, np.real(eig_vec_r[:ncels, i].reshape(nlat,nlon)),  
                         vmin=-max_color, vmax = max_color,
                         cmap = 'PiYG',
                         transform = ccrs.PlateCarree())
        sic.plot.contour(ax = ax, transform = ccrs.PlateCarree(), cmap ="cool", levels = [0.15,1.1], add_colorbar=False, linewidths=1)
        ax.set_title(r"$\lambda = $ {:.3f}".format( eig_val_r[i]), fontsize =10)
        box("{:.2f}".format(abs(eig_vec_r[-1,i]) ), .05 , .95 , ax)
        pf.set_circular_boundary(ax)
        ax.add_feature(cartopy.feature.LAND, facecolor = (0,0,0,.3), edgecolor = (0,0,0,0.3), zorder = 10)


        if colorbar: 
            fig.subplots_adjust(bottom=0.1, wspace=.1, hspace=.15)
            cbar_ax = fig.add_axes([0.15, 0.06, 0.7, 0.03])   #left, bottom, width , height
            cbar = fig.colorbar(p, cax=cbar_ax, orientation = 'horizontal', extend = 'both')
            cbar.set_ticks(labels)
            cbar.set_ticklabels(labels)

        if savefig: 
            single_fig_dir  = fig_dir + "transitionplots/04_21/winter_right_eigenvector_{}_days_res_{}/".format(days,res)
            if not os.path.exists(single_fig_dir): os.makedirs(single_fig_dir)       
            plt.savefig(single_fig_dir + "win_eig_ri_{}_days_year_{}_res_{}_k_{}.png".format(days,year,res,i+1), bbox_inches="tight", pad_inches=0.1)
#         plt.show()

    
    


      
        
'''
==============================================================
HOW TO CALL FOR THESE FUNCTIONS:
==============================================================

-------------------------------------------------------------
FOR SEPERATE ICE - WATER EIGENVECTORS
-------------------------------------------------------------
for i in range(12):
    run_prop = run(names_13_march[i])
    title = "{}, {} simdays, {}".format(run_prop.start_date,run_prop.simdays, run_prop.kernel)
    T, latbins, lonbins, M = transition_matrices.trans_mat_digit(names_13_march[i], res=1., final_index =-1, in_ice = False, void = True)
    transition_matrices.eigen_vector(M,latbins, lonbins, title = title, k=10)
    plt.show()
    
    
-------------------------------------------------------------
TRANSITIONS MULTIPLICATION PER FILE
-------------------------------------------------------------
    
res    = 1.   # degrees
length = 2400 # days

for i in [-1]:
    T, latbins, lonbins, M = transition_matrices.trans_mat_digit(names_13_march[i], res=res, final_index =-1, in_ice = False, void = True)
    date, days, kernel = prop_from_name(names_13_march[i])
    times = length//eval(days)
    title = "{}, {} simdays, {}".format(date,days, kernel)
    ax = transition_matrices.movement(M,latbins, lonbins, times=times,  grid_plot = False, title = title, savefig = False, save_name="")
    extentdate = datetime.datetime.strptime("2015-02-11", "%Y-%m-%d").date()
    plot_extent(extentdate, ax)
    plt.show()

    
'''