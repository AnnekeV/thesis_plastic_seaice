import numpy as np
from parcels import FieldSet, ParticleSet, ScipyParticle, JITParticle, ErrorCode, AdvectionRK4, Variable,ParticleFile, Field, VectorField
from datetime import timedelta
from datetime import datetime
from timeit import default_timer as timer
# from plot_functions import plot_gif
import glob
import os




def periodicBC(particle, fieldset, time):
    """
    Kernel for periodic boundaries in longitude, not sure if this is correct
    """
    if particle.lon < 0.:
        particle.lon += 360.
    elif particle.lon > 360.:
        particle.lon -= 360.
        
def DeleteParticle(particle, fieldset, time):
    """Kernel for deleting particles if they are out of bounds."""
    particle.delete()

def Sample_sit(particle, fieldset, time):  # Custom function that samples fieldset.P at particle location
    """
    Kernel for attributing sea ice thickness to a particle
    """
    particle.sit = fieldset.sit[time, particle.depth, particle.lat, particle.lon]
    particle.sip = fieldset.sip[time, particle.depth, particle.lat, particle.lon]
    
    
def AdvectionRK4_ice(particle, fieldset, time):
    """Advection of particles using fourth-order Runge-Kutta integration.
    Function needs to be converted to Kernel object before execution in ice or in water"""
    
    if particle.sit > 1e-16:
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
        
        
def AdvectionRK4_ocean(particle, fieldset, time):
    """Advection of particles using fourth-order Runge-Kutta integration.
    Function needs to be converted to Kernel object before execution"""
    (u1, v1) = fieldset.UVocean[time, particle.depth, particle.lat, particle.lon]
    lon1, lat1 = (particle.lon + u1*.5*particle.dt, particle.lat + v1*.5*particle.dt)
    (u2, v2) = fieldset.UVocean[time + .5 * particle.dt, particle.depth, lat1, lon1]
    lon2, lat2 = (particle.lon + u2*.5*particle.dt, particle.lat + v2*.5*particle.dt)
    (u3, v3) = fieldset.UVocean[time + .5 * particle.dt, particle.depth, lat2, lon2]
    lon3, lat3 = (particle.lon + u3*particle.dt, particle.lat + v3*particle.dt)
    (u4, v4) = fieldset.UVocean[time + particle.dt, particle.depth, lat3, lon3]
    particle.lon += (u1 + 2*u2 + 2*u3 + u4) / 6. * particle.dt
    particle.lat += (v1 + 2*v2 + 2*v3 + v4) / 6. * particle.dt       
   
    
def run_experiment(simdays = 2000, lat_npart = 30, lon_npart=50, ocean_currents = True, custom_kernel = "AdvectionRK4_ice" ):
    start = timer()
    Peeken = False
    
    res       = "0083"  
    data_dir  = '/data2/imau/oceanparcels/hydrodynamic_data/NEMO-MEDUSA/ORCA%s-N006/means/' %res #Directory for nemo data
    outputdir = '/scratch/AnnekeV/output/' #Directory for output files
#     outputdir = '/home/students/6252699/thesis/parcels2/output/'
    
    if (simdays < 305):
        ifiles    = sorted(glob.glob(data_dir+'ORCA%s-N06_2001????d05I.nc' %res))
        ufiles    = sorted(glob.glob(data_dir+'ORCA%s-N06_2001????d05U.nc' %res))
        vfiles    = sorted(glob.glob(data_dir+'ORCA%s-N06_2001????d05V.nc' %res))
        
    else:
        ifiles    = sorted(glob.glob(data_dir+'ORCA%s-N06_200?????d05I.nc' %res))
        ufiles    = sorted(glob.glob(data_dir+'ORCA%s-N06_200?????d05U.nc' %res))
        vfiles    = sorted(glob.glob(data_dir+'ORCA%s-N06_200?????d05V.nc' %res))
        

    mesh_mask = data_dir + "ORCA%s-N06_20090813d05U.nc" %res
    filenames = {'U':   {'lon': mesh_mask, 'lat': mesh_mask, 'data': ifiles},
                 'V':   {'lon': mesh_mask, 'lat': mesh_mask, 'data': ifiles},
                 'sit': {'lon': mesh_mask, 'lat': mesh_mask, 'data': ifiles}, # sea ice thickness
                 'sip': {'lon': mesh_mask, 'lat': mesh_mask, 'data': ifiles}} # sea ice presence

    variables = {'U': 'uice_ipa',
                 'V': 'vice_ipa',
                 'sit': 'sit'  ,
                 'sip':'sip'
                }

    dimensions = {'lat'  : 'nav_lat',
                  'lon'  : 'nav_lon',
                  'time' : 'time_counter'
                  }
    
    fieldset = FieldSet.from_nemo(filenames, variables, dimensions, allow_time_extrapolation=False)

    if ocean_currents:
        dimensionsU = {'data': 'uo', 'lon': 'nav_lon', 'lat': 'nav_lat', 'time': 'time_counter'}
        dimensionsV = {'data': 'vo', 'lon': 'nav_lon', 'lat': 'nav_lat', 'time': 'time_counter'}
        Uocean = Field.from_netcdf(ufiles, 'Uocean', dimensionsU, fieldtype='U', allow_time_extrapolation=False)
        Vocean = Field.from_netcdf(vfiles, 'Vocean', dimensionsV, fieldtype='V', allow_time_extrapolation=False) 
        fieldset.add_field(Uocean)
        fieldset.add_field(Vocean)
        uv_ocean = VectorField('UVocean', fieldset.Uocean, fieldset.Vocean)
        fieldset.add_vector_field(uv_ocean)     
    
    lat  = np.linspace(66,66.5,lat_npart)
    lon  = np.linspace(360-169.5,360-168,lon_npart)
    x,y  = np.meshgrid(lat,lon)
    
    npart = lat_npart *lon_npart
    time = np.repeat(np.datetime64('2001-02-01'),npart)
    

    if Peeken:
        '''Latitude and longitude from peeken''' 
        path = '/home/students/6252699/thesis/data/'
        coordinates = np.loadtxt(path + "peeken_table1.txt", skiprows = 1)
        lat_peeken  = coordinates[:,0]
        lon_peeken  = coordinates[:,1]
        '''Pick last three since they were measured in 2005, hence in a time slot that we can calculate stuff with'''
        lat_peeken    = lat_peeken[-3:]
        lon_peeken    = lon_peeken[-3:]
        time = [np.datetime64('2016-08-29'), np.datetime64('2016-08-18'), np.datetime64('2016-09-19')]
    

    
    
    class Sea_Ice_Particle(JITParticle):         # Define a new particle class
        sit = Variable('sit', initial=np.nan)  # Variable 'p' initialised by sampling the pressure
        sip = Variable('sip', initial=np.nan) 

    pset = ParticleSet.from_list(fieldset, 
                                 pclass = Sea_Ice_Particle, 
                                 time=time,
                                 lon=y, 
                                 lat=x)

    kernels    = pset.Kernel(periodicBC) + pset.Kernel(Sample_sit) +  pset.Kernel(eval(custom_kernel))

    output_name = outputdir   + "run_02_09_bering_npart_{}_start_{}_lat_{}_{}_simdays_{}_".format(npart, time[0],lat[0],lat[-1],simdays) + custom_kernel
    output_file =  pset.ParticleFile(
                    name=output_name, 
                    outputdt=timedelta(days=5))
    
    print("No. particles is {}".format(npart))
    print("Start time is ", time[0])
    print("Output name is " + output_name)

    pset.execute(kernels, 
                 runtime     = timedelta(days=simdays), 
                 dt          = timedelta(minutes=10), 
                 output_file = output_file,
                 recovery    = {ErrorCode.ErrorOutOfBounds: DeleteParticle})

    end = timer()
    
    print("Time elapsed = {:.1f} s ".format(end-start))
    return output_name


    
if __name__ == "__main__":
    output_name = run_experiment()
    plot_gif(output_name, dt_days = 5, extra_title = "bering2000winter", with_background = True)

