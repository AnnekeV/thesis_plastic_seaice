'''
File with function to backtrack different horizons
May 2
Anneke Vries
'''

import numpy as np
from parcels import FieldSet, ParticleSet, ScipyParticle, JITParticle, ErrorCode, AdvectionRK4, Variable
from argparse import ArgumentParser
from datetime import timedelta
from datetime import datetime
from glob import glob
import netCDF4 as netcdf
import matplotlib.pyplot as plt
import matplotlib as mpl
import cartopy
import cartopy.crs as ccrs
import xarray as xr
import cmocean
from timeit import default_timer as timer
import re
import pandas as pd
import sys
import plot_functions
from functions.calculations import round_down




def central_difference(data):
    '''
    Calculates the central difference in case of sea ice thickness
    Input dataset xarray with .sit and .obs
    Returns sea ice growth in array form
    '''
    sit_array  = data.sit
    obs_array  = data.obs
    sit_growth = np.zeros(len(sit_array.values))


    for i in range(1,len(sit_array.values)-1):
        dt            = obs_array[i+1]-obs_array[i-1]
        dh            = sit_array[i+1]-sit_array[i-1]
        sit_growth[i] = dh/dt

    sit_growth[0]  =  (sit_array[1]-sit_array[0])   /( obs_array[1]-obs_array[0])
    sit_growth[-1] =  (sit_array[-1]-sit_array[-2]) /( obs_array[-1]-obs_array[-2])
    
    return sit_growth

def roll_average(data, windowsize):
    backupdata = data.copy()[:]
    data       = data.rolling(obs=windowsize, center=True, min_periods=1).mean()
    return data.where(~np.isnan(data), backupdata)



def get_sign_derivative(data):
    data['sit'] = roll_average(data.sit,50)
    sig         = central_difference(data)    # sea ice growth
    return sig/abs(sig)



def find_origin_layers_core(particle_trajectory,
                            horizon_thickness = 0.1,
                            make_figure= True,
                            date = datetime(2014,3,1),
                            saveit = False,
                            savename = "sorted_per_layer_matrix",
                            out_dir = "/home/students/6252699/thesis/parcels2/output/melt_freeze/" ,
                            print_stuff = True):

    '''Finds the origin of layers in a core'''

    sitsit          = particle_trajectory.sit.values
    lonlon          = particle_trajectory.lon.values
    latlat          = particle_trajectory.lat.values
    loncore,latcore = lonlon[0], latlat[0]   
    core_thickness  = sitsit[0]
    core_layers     = np.arange(0,core_thickness, 0.01)   # for every centimeter find origin
    nr_pos_pairs    = len(core_layers)        # nr possible pairs
    ndays           = len(particle_trajectory.obs)

    print("Core thickness is {:.2f} m".format(core_thickness))

    matrix_for_core = np.ones([nr_pos_pairs,8])*-999
    counter = 0



    melt_or_freeze  = get_sign_derivative(particle_trajectory) # if larger than 0, then melt (since backtracking)

    '''FIND PAIRS'''
    print("\nFinding pairs...\n")

    if print_stuff: print("i\tj\tsit_i\tsit_j\tlon_i\tlat_i\tlon_j\tlat_j")
    if print_stuff: print("===============================================================")
    for i in range(nr_pos_pairs):

        if melt_or_freeze[i] ==1.:
            j = 0
            while j < ndays:
                if sitsit[j]<core_layers[i] and melt_or_freeze[j] == -1.:  # when is smaller and freezing 
                    if print_stuff: print("{}\t{}\t{:.2f}\t{:.2f}\t{:.1f}\t{:.1f}\t{:.1f}\t{:.1f}".format(i,j,core_layers[i],sitsit[j], loncore, latcore,lonlon[j],latlat[j]))
                    matrix_for_core[counter, :] = [i,j,core_layers[i],sitsit[j],loncore,latcore,lonlon[j],latlat[j]]
                    counter+=1
                    j+=1e9
                j+=1

    matrix_for_core[matrix_for_core==-999.] = np.nan  # replace -999 values

    # '''DATAFRAME'''
    # ds_core = pd.DataFrame(data=matrix_for_core, columns=['i', 'j', 'sit_i', 'sit_j', 'lon_i', 'lat_i', 'lon_j', 'lat_j'])
    # ds_core = ds_core.dropna()


    '''GROUP IN LAYERS OF 10 CM'''
    print("\nGrouping layers...\n")

    sit_i_discrete_core = matrix_for_core[:,2]//horizon_thickness   # 0.1 m groups
    layer_id             = sit_i_discrete_core[0]
    nr_pairs            = np.argwhere(np.isnan(sit_i_discrete_core))[0][0]
    i                   = 0
    all_layers_core     = []

    if print_stuff: print("Nr pairs is {}".format(nr_pairs))

    while i < range(nr_pairs):
        layer   = []
        while layer_id == sit_i_discrete_core[i]:
            layer.append(matrix_for_core[i,:])
            i+=1
            if (print_stuff and i%10==0):print ("{:.0f} %".format(i/np.float(nr_pairs)*100.))
        all_layers_core.append(layer)
        if (i == nr_pairs): break
        layer_id = sit_i_discrete_core[i]

    nr_layers = len(all_layers_core)


    if saveit:
        np.save(out_dir + savename, all_layers_core)

    if make_figure:
        print("\nMaking figure...\n")

        ymin, ymax = 60,90
        xmin, xmax = -180,180
        dx, dy = 30, 10

        cmap = plt.cm.jet  # define the colormap
        cmaplist = [cmap(i) for i in range(cmap.N)]     # extract all colors from the .jet map


        plt.figure(figsize = [5,5])
        ax = plt.axes(projection=ccrs.NorthPolarStereo())
        ax.gridlines(xlocs = np.arange(-180,185,dx), ylocs = np.arange(0,95,dy), color='black', alpha=0.5, linestyle=':')
        ax.add_feature(plot_functions.ocean_50m)
        ax.add_feature(plot_functions.land_50m)
        ax.set_extent([xmin,xmax,ymin,ymax],  ccrs.PlateCarree())
        plot_functions.set_circular_boundary(ax)

        stepsize = int(len(cmaplist)//10)
        for i in range(5): cmaplist += cmaplist
        spec_layer  = np.array(all_layers_core[0])
        scat = ax.scatter(spec_layer[:,4], spec_layer[:,5], transform=ccrs.PlateCarree(), marker = '^', color='black', label = "Core location", zorder=10)

        for i in range(nr_layers):
            spec_layer  = np.array(all_layers_core[i])
            layer_id    = int(spec_layer[0,2]*10)
            layer_label = "{}-{}".format(round_down(spec_layer[0,2]), round_down(spec_layer[0,2])+horizon_thickness)
            scat = ax.scatter(spec_layer[:,6], spec_layer[:,7], transform=ccrs.PlateCarree(), color=cmaplist[stepsize*layer_id], label = layer_label, zorder =11)


        sic = plot_functions.plot_extent_contour(date, ax)
        ax.set_title("")
        plt.legend(bbox_to_anchor=(1.1, 1))
    
    return all_layers_core
    

if __name__ == "__main__":

    fn = "/scratch/AnnekeV/output/single_particles/backtrack-04-26_npart_10_start_2016-12-01_simdays_1400_kernel_AdvectionRK4_prob_dtdays_1.nc"
    ds = xr.open_dataset(fn, decode_times=False)
    trajectory = ds.isel(traj=0)
    _ = find_origin_layers_core(particle_trajectory= trajectory)
