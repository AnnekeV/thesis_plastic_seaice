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


datadir = "/scratch/AnnekeV/reanalysis_data/"
# datadir = "/science/projects/oceanparcels/CMEMS-GLORYS12V1-Arctic/"

def IceOrOcean(particle, fieldset, time):
    if random.random() < math.atan((particle.sic*100)-15)/math.pi +.5  and particle.sic >.01:
        particle.in_ice = 1
    else:
        particle.in_ice = 0



    
def make_glorys_fieldset(year):
    '''Make the fieldset with the glorys data for 2014-2016, sithickness and sea ice concentration with ice and ocean velocities'''
    
    ifiles = sorted( glob.glob(datadir + "GLOBAL_REANALYSIS_PHY_001_030-TDS_????????_uv_uvice_con_thick.nc") )[-3000:]


    mesh_mask = "/scratch/AnnekeV/reanalysis_data/" + "GLOBAL_REANALYSIS_PHY_001_030-TDS_20150101_uv_uvice_con_thick.nc" #only works on an a grid! for c grid v and u are on a different grid
    
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
    Uocean = Field.from_netcdf(ifiles, 'Uocean', dimensionsU, fieldtype='U', allow_time_extrapolation=True)
    Vocean = Field.from_netcdf(ifiles, 'Vocean', dimensionsV, fieldtype='V', allow_time_extrapolation=True) 
    fieldset.add_field(Uocean)
    fieldset.add_field(Vocean)
    uv_ocean = VectorField('UVocean', fieldset.Uocean, fieldset.Vocean)
    fieldset.add_vector_field(uv_ocean)     
    fieldset.add_periodic_halo(zonal=True)
    fieldset.add_constant('halo_west', fieldset.U.grid.lon[0])
    fieldset.add_constant('halo_east', fieldset.U.grid.lon[-1])
    
    return fieldset





def get_particle_set(fieldset, start_date, kernel):
    '''Makes a particle set based on a fieldset and a fixed starting date'''
    
    class Sea_Ice_Particle(JITParticle):         # Define a new particle class to sample sea ice
        sit    = Variable('sit', initial=fieldset.sit)    # Sea ice thickness
        sic    = Variable('sic', initial=fieldset.sic)    # Sea ice concentration
        in_ice = Variable('in_ice', initial = 0.)
    
    
    ndays = 365*3
    lon   = np.tile(np.linspace(-169,-168.5,10), (ndays,1)).flatten()
    lat   = np.repeat(66, ndays*10)
    time  = np.repeat( np.arange(np.datetime64(start_date), np.datetime64(start_date) + np.timedelta64(ndays, 'D')), 10)
    npart = ndays*10
                
    pset = ParticleSet.from_list(fieldset, pclass = Sea_Ice_Particle, time=time, lon=lon, lat=lat) 

    if( kernel == "AdvectionRK4_prob"):
        for n in range(npart):
            p = pset.particles[n]  
            fieldset.computeTimeChunk(p.time,100 )

            '''Include initial condition for particle in ice or not'''
            psic = fieldset.sic[p.time,  p.depth, p.lat, p.lon]
            if random.random() < math.atan((psic*100)-15)/math.pi +.5  and psic >.01:
                p.in_ice = 1.
            else:
                p.in_ice = 0.
    
    return  pset, npart, lat





def run_experiment(start_date, simdays = 35, custom_kernel = "AdvectionRK4_prob", outputdir = '/home/students/6252699/thesis/parcels2/output/', year = 2014, dt_days=5):
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
    
    
    pset, npart, latp  = get_particle_set(fieldset, start_date, kernel)
    print("\nParticle set is ready")
    end = timer()
    print("Time elapsed = {:.1f} s ".format(end-start))
    
    
    kernels     =  pset.Kernel(eval(custom_kernel)) + pset.Kernel(SampleIce) + pset.Kernel(periodicBCC)  + IceOrOcean
    output_name = outputdir   + "bering-everyday-05-30_npart_{}_start_{}_simdays_{}_kernel_{}_dtdays_{}_10_at_a_time".format(npart, start_date,  simdays, custom_kernel, dt_days)
    
    output_file = pset.ParticleFile(name=output_name, outputdt=timedelta(days=dt_days))

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
    
    
    kernels     = ["AdvectionRK4_prob","AdvectionRK4_ocean"]
    start_date  = "2009-02-01"
    simdays     = 2800
    year        = 2011
    dt_days     = 1
    for kernel in kernels: 
        output_name = run_experiment(start_date = start_date, simdays = simdays, custom_kernel = kernel, outputdir = '/scratch/AnnekeV/output/single_particles/', year = year, dt_days = dt_days)
