import numpy as np
import glob
import os
from parcels import FieldSet, ParticleSet, ScipyParticle, JITParticle, ErrorCode, AdvectionRK4, Variable,ParticleFile, Field, VectorField
from kernels import  DeleteParticle, SampleIce, AdvectionRK4_ice, AdvectionRK4_ocean, AdvectionRK4_ice_sic, periodicBCC, IceOrOcean, AdvectionRK4_prob
from datetime import timedelta
from datetime import datetime
from timeit import default_timer as timer
import math
import random

def IceOrOcean15(particle, fieldset, time):
    if random.random() < math.atan((particle.sic*100)-15)/math.pi +.5  and particle.sic >.01:
        particle.in_ice = 1
    else:
        particle.in_ice = 0
        
def IceOrOcean5(particle, fieldset, time):
    if random.random() < math.atan((particle.sic*100)-5)/math.pi +.5  and particle.sic >.01:
        particle.in_ice = 1
    else:
        particle.in_ice = 0
                
def IceOrOcean30(particle, fieldset, time):
    if random.random() < math.atan((particle.sic*100)-30)/math.pi +.5  and particle.sic >.01:
        particle.in_ice = 1
    else:
        particle.in_ice = 0
        
def IceOrOcean60(particle, fieldset, time):
    if random.random() < math.atan((particle.sic*100)-60)/math.pi +.5  and particle.sic >.01:
        particle.in_ice = 1
    else:
        particle.in_ice = 0
        
def IceOrOcean95(particle, fieldset, time):
    if random.random() < math.atan((particle.sic*100)-95)/math.pi +.5  and particle.sic >.01:
        particle.in_ice = 1
    else:
        particle.in_ice = 0
        
        
def grid_per_x_km(xkm = 1e4, ylims = [50,90]):
    '''Makes a grid with a particle per x number of meter, 
    evenly spaced per latitude. xkm is in m, not km. Default is 10 km between 50 and 90 N. Returns lons, lats op particles'''
    R     = 6371.0e3              #m Radius earth
    xdeg  = xkm/(R*np.pi/180.)    # per degree
    lat   = np.arange(ylims[0],ylims[1], xdeg)
    n     = np.array((np.cos(lat/180.*np.pi)*R/xkm), dtype='int')
    
    lats, lons = [], []
    for i in range(len(n)):
        lats.extend(np.repeat(lat[i], n[i]))
        lons.extend(np.linspace(-180,180,n[i], endpoint=False))
        
    return lons, lats

    
def make_glorys_fieldset(year):
    '''Make the fieldset with the glorys data for 2014-2016, sithickness and sea ice concentration with ice and ocean velocities'''
    
    data_dir = "/scratch/AnnekeV/reanalysis_data/"
    
    ifiles = sorted(glob.glob(data_dir + "GLOBAL_REANALYSIS_PHY_001_030-TDS_{}????_uv_uvice_con_thick.nc".format(year) )
#                   + sorted(glob.glob("/scratch/AnnekeV/reanalysis_data/" + "GLOBAL_REANALYSIS_PHY_001_030-TDS_{:d}*_uv_uvice_con_thick.nc".format(2015) ) )
#                   + sorted(glob.glob("/scratch/AnnekeV/reanalysis_data/" + "GLOBAL_REANALYSIS_PHY_001_030-TDS_{:d}*_uv_uvice_con_thick.nc".format(2016) ) )
                 )


    mesh_mask = data_dir + "GLOBAL_REANALYSIS_PHY_001_030-TDS_20150101_uv_uvice_con_thick.nc" #only works on an a grid! for c grid v and u are on a different grid
    filenames = {'U':   {'lon': mesh_mask, 'lat': mesh_mask, 'data': ifiles},
                 'V':   {'lon': mesh_mask, 'lat': mesh_mask, 'data': ifiles},
                 'sit': {'lon': mesh_mask, 'lat': mesh_mask, 'data': ifiles}, # sea ice thickness
                 'sic': {'lon': mesh_mask, 'lat': mesh_mask, 'data': ifiles}} # sea ice concentration
    variables = {'U'  : 'usi',
                 'V'  : 'vsi',
                 'sit': 'sithick',
                 'sic': 'siconc'
                }
    dimensions = {'lat'  : 'latitude',
                  'lon'  : 'longitude',
                  'time' : 'time'
                  }

    fieldset = FieldSet.from_netcdf(filenames, variables, dimensions, allow_time_extrapolation=False)

    '''Add the ocean velocities'''
    dimensionsU = {'data': 'uo', 'lon': 'longitude', 'lat': 'latitude', 'time': 'time'}
    dimensionsV = {'data': 'vo', 'lon': 'longitude', 'lat': 'latitude', 'time': 'time'}
    Uocean = Field.from_netcdf(ifiles, 'Uocean', dimensionsU, fieldtype='U', allow_time_extrapolation=False)
    Vocean = Field.from_netcdf(ifiles, 'Vocean', dimensionsV, fieldtype='V', allow_time_extrapolation=False) 
    fieldset.add_field(Uocean)
    fieldset.add_field(Vocean)
    uv_ocean = VectorField('UVocean', fieldset.Uocean, fieldset.Vocean)
    fieldset.add_vector_field(uv_ocean)     
    fieldset.add_periodic_halo(zonal=True)
    fieldset.add_constant('halo_west', fieldset.U.grid.lon[0])
    fieldset.add_constant('halo_east', fieldset.U.grid.lon[-1])
    
    return fieldset





def get_particle_set(fieldset, start_date, extent):
    '''Makes a particle set based on a fieldset and a fixed starting date'''
    
    class Sea_Ice_Particle(JITParticle):         # Define a new particle class to sample sea ice
        sit    = Variable('sit', initial=fieldset.sit)    # Sea ice thickness
        sic    = Variable('sic', initial=fieldset.sic)    # Sea ice concentration
        in_ice = Variable('in_ice', initial = 0.)

#     lon,lat = grid_per_x_km(xkm=10e4, ylims = [50,90])   # xkm = km between data points
    lon, lat = np.meshgrid(np.arange(-180,180, 2./15.), np.arange(50,90, 2./15.))
    lon, lat = lon.flatten(), lat.flatten()
    
    
    npart   = len(lat)
    time    = np.repeat(np.datetime64(start_date), npart)

    print("Number of particles before: {}".format(npart))

    '''Make temporary pset in order to determine where land is ''' 
    pset_try  = ParticleSet.from_list(fieldset, pclass = Sea_Ice_Particle, time=time, lon=lon, lat=lat)
    fieldset.computeTimeChunk(pset_try.particles[0].time,1)  # Load timechunk for 1 day
    latp,lonp = [],[]

    for n in range(npart):
        p = pset_try.particles[n]
        velo_p = np.max(np.abs(fieldset.UVocean[p.time,  p.depth, p.lat, p.lon]))   # Determine velocity at initial location
        if velo_p >0.0:
            latp.append(p.lat)
            lonp.append(p.lon)
    

                
    npart = len(latp)  #Determine new nr.particles (not on land)
    time = np.repeat(np.datetime64(start_date),npart)
    pset = ParticleSet.from_list(fieldset, pclass = Sea_Ice_Particle, time=time, lon=lonp, lat=latp) 

    for n in range(npart):
        p = pset.particles[n]            
        '''Include initial condition for particle in ice or not'''
        psic = fieldset.sic[p.time,  p.depth, p.lat, p.lon]
        if random.random() < math.atan((psic*100)-extent)/math.pi +.5  and psic >.01:
            p.in_ice = 1.
        else:
            p.in_ice = 0.
    
    return  pset, npart, latp





def run_experiment(start_date, simdays = 35, custom_kernel = "AdvectionRK4_prob", outputdir = '/home/students/6252699/thesis/parcels2/output/', year = 2014, extent=15):
    '''
    Run the experiment for simdays,  advection of ice or only ocean AdvectionRK4_ocean.
    Start_date should be start_date as '2014-02-01'
    Outputdir depends on whether run in terminal (T) or from jupyter notebook (JN).
    Terminal/qsub:   /scratch
    JupyterNotebook: /home
    '''    

    start     = timer()
    fieldset  = make_glorys_fieldset(year)
    print("\nFieldset is ready")  
    end = timer()
    print("Time elapsed = {:.1f} s ".format(end-start))
    
    
    pset, npart, latp  = get_particle_set(fieldset, start_date, extent)
    print("\nParticle set is ready")
    end = timer()
    print("Time elapsed = {:.1f} s ".format(end-start))
    
    ice_or_ocean_kernel = "IceOrOcean{}".format(extent)
    
    kernels     =  pset.Kernel(eval(custom_kernel)) + pset.Kernel(SampleIce) + pset.Kernel(periodicBCC)  + pset.Kernel(eval(ice_or_ocean_kernel))
    output_name = outputdir   + "nr2-04-05_npart_{}_start_{}_simdays_{}__extent_{}_kernel_".format(npart, start_date,  simdays, extent) + custom_kernel
    
    output_file = pset.ParticleFile(name=output_name, outputdt=timedelta(days=5))

    print("\nNo. particles is {}".format(npart))
    print("Start time is {} "    .format(start_date))
    print("Runtime is {}"        .format(simdays))
    print("Output name is "      +output_name)
    print("Latitude: {:.0f}-{:.0f}\n".format(latp[0], latp[-1]))

    pset.execute(kernels, 
                 runtime     = timedelta(days=simdays), 
                 dt          = timedelta(minutes=10), 
                 output_file = output_file,
                 recovery    = {ErrorCode.ErrorOutOfBounds: DeleteParticle})

    end = timer()
    print("Time elapsed = {:.1f} s ".format(end-start))
    return output_name




if __name__ == "__main__":
    
    
    
    kernel      = "AdvectionRK4_ocean"
    start_date  = "2015-02-01"
    simdays     = 60
    year        = 2015
    output_name = run_experiment(start_date = start_date, simdays = simdays, custom_kernel = kernel, outputdir = '/scratch/AnnekeV/output/sensitivity/', year = year)

    
        
    kernel      = "AdvectionRK4_prob"
    start_date  = "2015-08-01"
    simdays     = 60
    year        = 2015
    output_name = run_experiment(start_date = start_date, simdays = simdays, custom_kernel = kernel, outputdir = '/scratch/AnnekeV/output/sensitivity/', year = year)

    
            
    kernel      = "AdvectionRK4_prob"
    start_date  = "2015-02-01"
    simdays     = 60
    year        = 2015
    output_name = run_experiment(start_date = start_date, simdays = simdays, custom_kernel = kernel, outputdir = '/scratch/AnnekeV/output/sensitivity/', year = year, 
                                 extent = 5)

    
    
    
    kernel      = "AdvectionRK4_prob"
    start_date  = "2015-02-01"
    simdays     = 60
    year        = 2015
    output_name = run_experiment(start_date = start_date, simdays = simdays, custom_kernel = kernel, outputdir = '/scratch/AnnekeV/output/sensitivity/', year = year, 
                                 extent = 30)

    
    kernel      = "AdvectionRK4_prob"
    start_date  = "2015-02-01"
    simdays     = 60
    year        = 2015
    output_name = run_experiment(start_date = start_date, simdays = simdays, custom_kernel = kernel, outputdir = '/scratch/AnnekeV/output/sensitivity/', year = year, 
                                 extent = 60)    
    
    
    
    kernel      = "AdvectionRK4_prob"
    start_date  = "2015-02-01"
    simdays     = 60
    year        = 2015
    output_name = run_experiment(start_date = start_date, simdays = simdays, custom_kernel = kernel, outputdir = '/scratch/AnnekeV/output/sensitivity/', year = year, 
                                 extent = 95) 
#     kernels = ["AdvectionRK4_prob", "AdvectionRK4_ocean"]

#     for kernel in kernels:
#         for month in [2, 8]:
#             for year in [2014, 2015, 2016]:
#                 '''60 days'''
#                 simdays     = 60
#                 try:                 
#                     start_date  = '{}-{:02d}-01'.format( year, month)
#                     print("\n Running \nyear: {}\nkernel: {}\ndays: {}\n".format(year, kernel, simdays))
#                     output_name = run_experiment(start_date = start_date, simdays = simdays, custom_kernel = kernel, outputdir = '/scratch/AnnekeV/output/bigrun/', year = year)

#                 except:
#                     print("\n\n You're run didn't work for \nyear: {}\nkernel: {}\ndays: {}\n\n".format(year, kernel, simdays))

#                 '''40 days'''
#                 simdays     = 40
#                 try:                 
#                     start_date  = '{}-{:02d}-11'.format( year, month)   
                    

#                     print("\n Running \nyear: {}\nkernel: {}\ndays: {}\n".format(year, kernel, simdays))
#                     output_name = run_experiment(start_date = start_date, simdays = simdays, custom_kernel = kernel, outputdir = '/scratch/AnnekeV/output/bigrun/', year = year)

#                 except:
#                     print("\n\n You're run didn't work for \nyear: {}\nkernel: {}\ndays: {}\n\n".format(year, kernel, simdays))

#                 simdays     = 20
#                 '''20 days'''
#                 try:      
#                     start_date  = '{}-{:02d}-21'.format(year, month)
                    

#                     print("\n Running \nyear: {}\nkernel: {}\ndays: {}\n".format(year, kernel, simdays))
#                     output_name = run_experiment(start_date = start_date, simdays = simdays, custom_kernel = kernel, outputdir = '/scratch/AnnekeV/output/bigrun/', year = year)

#                 except:
#                     print("\n\n You're run didn't work for \nyear: {}\nkernel: {}\ndays: {}\n\n".format(year, kernel, simdays))

