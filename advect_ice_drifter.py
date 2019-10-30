import numpy as np
import glob
import os
from parcels import FieldSet, ParticleSet, ScipyParticle, JITParticle, ErrorCode, AdvectionRK4, Variable,ParticleFile, Field, VectorField
from kernels import  DeleteParticle, SampleIce, AdvectionRK4_ice, AdvectionRK4_ocean, AdvectionRK4_ice_sic, periodicBCC, IceOrOcean, AdvectionRK4_prob, IceOrOcean2
from datetime import timedelta
from datetime import datetime
from timeit import default_timer as timer
import math
import random
import pandas as pd

data_dir = "/data/oceanparcels/CMEMS-GLORYS12V1-Arctic/"

def IceOrOcean(particle, fieldset, time):
    if random.random() < math.atan((particle.sic*100)-15)/math.pi +.5  and particle.sic >.01:
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
    
    ifiles = sorted( sorted(glob.glob(data_dir + "GLOBAL_REANALYSIS_PHY_001_030-TDS_*_uv_uvice_con_thick.nc" ) )
#                   + sorted(glob.glob(data_dir + "GLOBAL_REANALYSIS_PHY_001_030-TDS_{:d}*_uv_uvice_con_thick.nc".format(2010) ) )
#                   + sorted(glob.glob(data_dir + "GLOBAL_REANALYSIS_PHY_001_030-TDS_{:d}*_uv_uvice_con_thick.nc".format(2008) ) )[-10:]
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





def get_particle_set(fieldset, start_date, init_lon, init_lat, b_id):
    '''Makes a particle set based on a fieldset and a fixed starting date'''
    

    class Sea_Ice_Particle(JITParticle):         # Define a new particle class to sample sea ice
        sit        = Variable('sit', initial=fieldset.sit)    # Sea ice thickness
        sic        = Variable('sic', initial=fieldset.sic)    # Sea ice concentration
        in_ice     = Variable('in_ice', initial = 1.)
        prev_sit   = Variable('prev_sit', initial=fieldset.sit)
        prev_state = Variable('prev_state', initial=1.)
        buoy_id    = Variable('buoy_id', initial=b_id)
        
    def create_time_array(start_date, dt_hours=12):
        start = np.datetime64(start_date)-2*np.timedelta64(dt_hours, 'h') 
        end   = np.datetime64(start_date)+3*np.timedelta64(dt_hours, 'h')
        time = np.arange(start, end, np.timedelta64(dt_hours, 'h'))
        return time

    
    npart   = 5
    time    = create_time_array(start_date)
    lon     = np.repeat(init_lon, npart)
    lat     = np.repeat(init_lat, npart)


    print("Number of particles before: {}".format(npart))

    pset = ParticleSet.from_list(fieldset, pclass = Sea_Ice_Particle, time=time, lon=lon, lat=lat) 
    
    for n in range(npart):
        p = pset.particles[n]
        fieldset.computeTimeChunk(p.time,100 )
        '''Include initial condition for particle in ice or not'''
        psic = fieldset.sic[p.time,  p.depth, p.lat, p.lon]
        if random.random() < math.atan((psic*100)-15)/math.pi +.5  and psic >.01:
            p.in_ice = 1.
        else:
            p.in_ice = 0.

    
    return  pset, npart





def run_experiment(start_date, simdays = 35, custom_kernel = "AdvectionRK4_prob", outputdir = '/home/students/6252699/thesis/parcels2/output/', year = 2014):
    '''
    Run the experiment for simdays,  advection of ice or only ocean AdvectionRK4_ocean.
    Start_date should be start_date as '2014-02-01'
    Outputdir depends on whether run in terminal (T) or from jupyter notebook (JN).
    Terminal/qsub:   /scratch
    JupyterNotebook: /home
    '''  
    
    df = pd.read_pickle("/home/students/6252699/thesis/parcels2/output/df_init")

    start     = timer()
    fieldset  = make_glorys_fieldset(year)
    print("\nFieldset is ready")  
    end = timer()
    print("Time elapsed = {:.1f} s ".format(end-start))
    
    
    #     ============== buoy id =================
    
    for index, row in df.iterrows():
        try: 
            b_id        = row.buoy_id
            start_date  = row.date
            init_lon, init_lat = row.lon, row.lat
            simdays     = row.simdays
            end_date    = row.end_date

            if start_date.year <2001: continue
            if end_date.year   >2016: continue

            print("\n")
            print("b_id    : {} ".format(b_id))
            print("nstart  : {}".format(start_date))
            print("simdays : {}".format(simdays))
            print("{}, {} ".format(init_lon, init_lat))



        #     ============== NEW =====================
            pset, npart  = get_particle_set(fieldset, start_date,  init_lon, init_lat, b_id)
            print("\nParticle set is ready")
            end = timer()
            print("Time elapsed = {:.1f} s ".format(end-start))


            kernels_new =  pset.Kernel(eval(custom_kernel)) + pset.Kernel(SampleIce) + pset.Kernel(periodicBCC)  + IceOrOcean2
            output_name = outputdir   + "06-26_drifter{}_new_npart_{}_start_{:4d}-{:02d}-{:02d}_dt_12_hours_simdays_{}_kernel_".format(b_id, npart, start_date.year,
                                                                                                                            start_date.month, start_date.day,  simdays) + custom_kernel

            output_file = pset.ParticleFile(name=output_name, outputdt=timedelta(hours=12))

            print("\nNo. particles is {}".format(npart))
            print("Start time is {} "    .format(start_date))
            print("Runtime is {}"        .format(simdays))
            print("Output name is "      + output_name)


            pset.execute(kernels_new, 
                         runtime     = timedelta(days=simdays), 
                         dt          = timedelta(minutes=10), 
                         output_file = output_file,
                         recovery    = {ErrorCode.ErrorOutOfBounds: DeleteParticle})
            print("Done with new kernel")



        #     ============== OLD =====================
            pset, npart  = get_particle_set(fieldset, start_date,  init_lon, init_lat, b_id)
            print("\nParticle set is ready")
            end = timer()
            print("Time elapsed = {:.1f} s ".format(end-start))


            kernels_old =  pset.Kernel(eval(custom_kernel)) + pset.Kernel(SampleIce) + pset.Kernel(periodicBCC)  + IceOrOcean
            output_name = outputdir   + "06-14_drifter{}_old_npart_{}_start_{:4d}-{:02d}-{:02d}_dt_12_hours_simdays_{}_kernel_".format(b_id, npart, start_date.year, start_date.month, start_date.day,  simdays) + custom_kernel

            output_file = pset.ParticleFile(name=output_name, outputdt=timedelta(hours=12))

            print("\nBuoy id is {}"      .format(b_id))
            print("No. particles is {}"  .format(npart))
            print("Start time is {} "    .format(start_date))
            print("Runtime is {}"        .format(simdays))
            print("Output name is "      +output_name)


            pset.execute(kernels_old, 
                         runtime     = timedelta(days=simdays), 
                         dt          = timedelta(minutes=10), 
                         output_file = output_file,
                         recovery    = {ErrorCode.ErrorOutOfBounds: DeleteParticle})
            print("Done with old kernel")

            end = timer()
            print("Time elapsed = {:.1f} s ".format(end-start))
    
        except:
            print("Error for {}".format(b_id))
    return output_name




if __name__ == "__main__":
    
    
    
    kernel      = "AdvectionRK4_prob"
    start_date  = '2011-01-01T00:00'
    simdays     = 600
    year        = 2011
    output_name = run_experiment(start_date = start_date, simdays = simdays, custom_kernel = kernel, outputdir = '/scratch/AnnekeV/output/drifters/all/', year = year)

    
    
