import numpy as np
from parcels import FieldSet, ParticleSet, ScipyParticle, JITParticle, ErrorCode, AdvectionRK4, Variable,ParticleFile, Field, VectorField
from datetime import timedelta
from datetime import datetime
from timeit import default_timer as timer
from advect_collect_double_field import  DeleteParticle, Sample_sit, AdvectionRK4_ice, AdvectionRK4_ocean
import glob
import os
     

    
def periodicBCC(particle, fieldset, time):
    '''Make periodic boundary kernel with halo'''
    if particle.lon < fieldset.halo_west:
        particle.lon += fieldset.halo_east - fieldset.halo_west
    elif particle.lon > fieldset.halo_east:
        particle.lon -= fieldset.halo_east - fieldset.halo_west
     
    
def grid_per_x_km(xkm = 1e4, ylims = [50,90]):
    "Makes a grid per x number of m, evenly spaced per latitude. Default is 10 km between 50 and 90 lat. Returns lons, lats "
    R     = 6371.0e3              #m Radius earth
    xdeg  = xkm/(R*np.pi/180.)    # per degree
    lat   = np.arange(ylims[0],ylims[1], xdeg)
    n     = np.array((np.cos(lat/180.*np.pi)*R/xkm), dtype='int')
    
    lats, lons = [], []
    for i in range(len(n)):
        lats.extend(np.repeat(lat[i], n[i]))
        lons.extend(np.linspace(-180,180,n[i], endpoint=False))
        
    return lons, lats

def AdvectionRK4_ice_sic(particle, fieldset, time):
    """Advection of particles using fourth-order Runge-Kutta integration.
    Function needs to be converted to Kernel object before execution in ice or in water
    In ice is defined as 15% of ice """
    
    if particle.sip >= 0.15:
        (u1, v1) = fieldset.UV[time, particle.depth, particle.lat, particle.lon]
        lon1, lat1 = (particle.lon + u1*.5*particle.dt, particle.lat + v1*.5*particle.dt)
        (u2, v2) = fieldset.UV[time + .5 * particle.dt, particle.depth, lat1, lon1]
        lon2, lat2 = (particle.lon + u2*.5*particle.dt, particle.lat + v2*.5*particle.dt)
        (u3, v3) = fieldset.UV[time + .5 * particle.dt, particle.depth, lat2, lon2]
        lon3, lat3 = (particle.lon + u3*particle.dt, particle.lat + v3*particle.dt)
        (u4, v4) = fieldset.UV[time + particle.dt, particle.depth, lat3, lon3]
        particle.lon += (u1 + 2*u2 + 2*u3 + u4) / 6. * particle.dt
        particle.lat += (v1 + 2*v2 + 2*v3 + v4) / 6. * particle.dt
        
        
    else:
        (u1, v1) = fieldset.UVocean[time, particle.depth, particle.lat, particle.lon]
        lon1, lat1 = (particle.lon + u1*.5*particle.dt, particle.lat + v1*.5*particle.dt)
        (u2, v2) = fieldset.UVocean[time + .5 * particle.dt, particle.depth, lat1, lon1]
        lon2, lat2 = (particle.lon + u2*.5*particle.dt, particle.lat + v2*.5*particle.dt)
        (u3, v3) = fieldset.UVocean[time + .5 * particle.dt, particle.depth, lat2, lon2]
        lon3, lat3 = (particle.lon + u3*particle.dt, particle.lat + v3*particle.dt)
        (u4, v4) = fieldset.UVocean[time + particle.dt, particle.depth, lat3, lon3]
        particle.lon += (u1 + 2*u2 + 2*u3 + u4) / 6. * particle.dt
        particle.lat += (v1 + 2*v2 + 2*v3 + v4) / 6. * particle.dt     


    
def make_glorys_fieldset():
    '''Make the fieldset with the glorys data for 2014-2016, sithickness and sea ice concentration with ice and ocean velocities'''
    
    all_glorys = glob.glob("/scratch/AnnekeV/reanalysis_data/" + "GLOBAL_REANALYSIS_PHY_001_030-TDS_2014????_uv_uvice_con_thick.nc" )# + glob.glob("/scratch/AnnekeV/reanalysis_data/" + "GLOBAL_REANALYSIS_PHY_001_030-TDS_{:d}*_uv_uvice_con_thick.nc".format(2015) ) + glob.glob("/scratch/AnnekeV/reanalysis_data/" + "GLOBAL_REANALYSIS_PHY_001_030-TDS_{:d}*_uv_uvice_con_thick.nc".format(2016) ) 
    ifiles = sorted(all_glorys)


    mesh_mask = "/scratch/AnnekeV/reanalysis_data/" + "GLOBAL_REANALYSIS_PHY_001_030-TDS_20150101_uv_uvice_con_thick.nc" 
    filenames = {'U':   {'lon': mesh_mask, 'lat': mesh_mask, 'data': ifiles},
                 'V':   {'lon': mesh_mask, 'lat': mesh_mask, 'data': ifiles},
                 'sit': {'lon': mesh_mask, 'lat': mesh_mask, 'data': ifiles}, # sea ice thickness
                 'sip': {'lon': mesh_mask, 'lat': mesh_mask, 'data': ifiles}} # sea ice concentration
    variables = {'U': 'usi',
                 'V': 'vsi',
                 'sit': 'sithick',
                 'sip': 'siconc'
                }
    dimensions = {'lat'  : 'latitude',
                  'lon'  : 'longitude',
                  'time' : 'time'
                  }

    print(filenames)
    fieldset = FieldSet.from_netcdf(filenames, variables, dimensions, allow_time_extrapolation=False)

    '''Add the ocean velocities'''
    dimensionsU = {'data': 'uo', 'lon': 'longitude', 'lat': 'latitude', 'time': 'time'}
    dimensionsV = {'data': 'vo', 'lon': 'longitude', 'lat': 'latitude', 'time': 'time'}
    Uocean = Field.from_netcdf(ifiles, 'Uocean', dimensionsU, fieldtype='U', allow_time_extrapolation=False)
    Vocean = Field.from_netcdf(ifiles, 'Vocean', dimensionsV, fieldtype='V', allow_time_extrapolation=False) #,                                  grid=Uocean.grid, dataFiles=Uocean.dataFiles)
    fieldset.add_field(Uocean)
    fieldset.add_field(Vocean)
    uv_ocean = VectorField('UVocean', fieldset.Uocean, fieldset.Vocean)
    fieldset.add_vector_field(uv_ocean)     

    fieldset.add_periodic_halo(zonal=True)
    fieldset.add_constant('halo_west', fieldset.U.grid.lon[0])
    fieldset.add_constant('halo_east', fieldset.U.grid.lon[-1])
    
    return fieldset





def get_particle_set(fieldset, start_date):
    
    class Sea_Ice_Particle(JITParticle):         # Define a new particle class to sample sea ice
        sit = Variable('sit', initial=np.nan)  
        sip = Variable('sip', initial=np.nan) 

    lon,lat = grid_per_x_km(xkm=1e4, ylims = [50,90])
    npart   = len(lat)
    time    = np.repeat(np.datetime64(start_date), npart)

    '''Make temporary pset in order to determine where land is ''' 
    pset_try  = ParticleSet.from_list(fieldset, pclass = Sea_Ice_Particle, time=time, lon=lon, lat=lat)
    latp,lonp = [],[]
    fieldset.computeTimeChunk(pset_try.particles[0].time,1)  # load timechunk for 1 day

    for n in range(npart):
        p = pset_try.particles[n]
        velo_p = np.max(np.abs(fieldset.UVocean[p.time,  p.depth, p.lat, p.lon]))   #determine velocity at initial location
        if velo_p >0.0:
            latp.append(p.lat)
            lonp.append(p.lon)

    npart = len(latp)  #determine new nr.particles (not on land)
    time = np.repeat(np.datetime64(start_date),npart)
    pset = ParticleSet.from_list(fieldset, pclass = Sea_Ice_Particle, time=time, lon=lonp, lat=latp) # make a new particle set only with particles that are not on land

    return  pset, npart, latp





def run_experiment(start_date, simdays = 35, custom_kernel = "AdvectionRK4_ice", outputdir = '/home/students/6252699/thesis/parcels2/output/'):
    '''
    Run the experiment for simdays,  advection of ice or only ocean AdvectionRK4_ocean.
    Start_date should be start_date as '2014-02-01'
    Outputdir depends on whether run in terminal (T) or from jupyter notebook (JN).
    T:  /scratch
    JN: /home
    '''    

    start     = timer()
    fieldset  = make_glorys_fieldset()
    print("\nFieldset is ready")  
    end = timer()
    print("Time elapsed = {:.1f} s ".format(end-start))
    
    pset, npart, latp  = get_particle_set(fieldset, start_date)
    print("\nParticle set is ready")
    end = timer()
    print("Time elapsed = {:.1f} s ".format(end-start))
    
    kernels     =  pset.Kernel(eval(custom_kernel)) + pset.Kernel(Sample_sit) + pset.Kernel(periodicBCC)  

    output_name = outputdir   + "kernel_test_run_02_21_npart_{}_start_{}_simdays_{}_".format(npart, start_date,  simdays) + custom_kernel
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
#     for m in range(1,13):
    for kernel in ["AdvectionRK4_ice_sic"]:
        start_date = '2014-{:02d}-01'.format(2)
        output_name = run_experiment(start_date = start_date, simdays =100, custom_kernel = kernel, outputdir = '/scratch/AnnekeV/output/')

