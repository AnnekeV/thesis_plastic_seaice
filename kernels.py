# These are the different kernels used in the project by Anneke Vries for simulating plastic in sea ice
# Last kernel used to determine state of particle ('in ice' or 'in water') was IceOrOcean2

import glob
import os
from parcels import FieldSet, ParticleSet, ScipyParticle, JITParticle, ErrorCode, AdvectionRK4, Variable,ParticleFile, Field, VectorField
import random

random.seed()

     

    
def periodicBCC(particle, fieldset, time):
    '''Make periodic boundary kernel with halo'''
    if particle.lon < fieldset.halo_west:
        particle.lon += fieldset.halo_east - fieldset.halo_west
    elif particle.lon > fieldset.halo_east:
        particle.lon -= fieldset.halo_east - fieldset.halo_west

        
def DeleteParticle(particle, fieldset, time):
    """Kernel for deleting particles if they are out of bounds."""
    particle.delete()

def SampleIce(particle, fieldset, time):  # Custom function that samples fieldset.P at particle location
    """
    Kernel for attributing sea ice thickness and concentration to a particle
    """
    particle.sit = fieldset.sit[time, particle.depth, particle.lat, particle.lon]
    particle.sic = fieldset.sic[time, particle.depth, particle.lat, particle.lon]
    
    
def TotalDistance(particle, fieldset, time):
    # Calculate the distance in latitudinal direction (using 1.11e2 kilometer per degree latitude)
    lat_dist = (particle.lat - particle.prev_lat) * 1.11e2
    # Calculate the distance in longitudinal direction, using cosine(latitude) - spherical earth
    lon_dist = (particle.lon - particle.prev_lon) * 1.11e2 * math.cos(particle.lat * math.pi / 180)
    # Calculate the total Euclidean distance travelled by the particle
    particle.distance += math.sqrt(math.pow(lon_dist, 2) + math.pow(lat_dist, 2))

    particle.prev_lon = particle.lon  # Set the stored values for next iteration.
    particle.prev_lat = particle.lat
    

def IceOrOcean(particle, fieldset, time):
    """
    Determines state (in ice or in water) particle depending on sea ice concentration
    """
    if random.random() < math.atan((particle.sic*100)-15)/math.pi +.5  and particle.sic >.01:
        particle.in_ice = 1
    else:
        particle.in_ice = 0
    
    
    
    
def AdvectionRK4_ice_sic(particle, fieldset, time):
    """Advection of particles using fourth-order Runge-Kutta integration.
    Function needs to be converted to Kernel object before execution in ice or in water
    In ice is defined as 15% of ice """
    
    if particle.sic >= 0.15:
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
        
        
def AdvectionRK4_ice(particle, fieldset, time):
    """Advection of particles using fourth-order Runge-Kutta integration.
    Function needs to be converted to Kernel object before execution in ice or in water
    Ice defined as more than 0 thickness"""
    
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

        
def AdvectionRK4_prob(particle, fieldset, time):
    """Advection of particles using fourth-order Runge-Kutta integration.
    Function needs to be converted to Kernel object before execution in ice or in water
    Ice defined on a probability function arctan with sic >0.15"""
    
    if particle.in_ice ==1:
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
        

def IceOrOcean2(particle, fieldset, time):
    """
    This kernel determines if a praticle is in the ocean or in ice (in order to later be advected by one or the other). It first determines the state of the particle at timestep t and then for timestep t+1 it checks 
    if the ice thickness has in or decreased. If the particle was previously in ice, and the ice has grown, it stays in ice, if the ice has become thinner, it has a chance (depending on sea ice concentration) to 
    enter water. If the particle was previously in water, and the ice has become thinner, it will remain in water, but if the ice has increase, it has a chance (depending on sea ice concentration) to enter the ice. 
    """
        if math.isnan(particle.prev_sit):
            particle.in_ice == 0
        
        # Previous state was in ice
        elif particle.prev_state == 1:
            # If was in ice but now growing
            if particle.sit-particle.prev_sit >= 0:
                particle.in_ice = 1

            # If was in ice but now melting
            else:
                # If melting but still above arctan
                if random.random() < math.atan((particle.sic*100)-15)/math.pi +.5  and particle.sic >.01:
                    particle.in_ice = 1
                else:
                    particle.in_ice = 0        

        # Previous state was in water
        else: 
            # If was in water but now growing
            if particle.sit-particle.prev_sit >= 0:
                if random.random() < math.atan((particle.sic*100)-15)/math.pi +.5  and particle.sic >.01:
                    particle.in_ice = 1
                else:
                    particle.in_ice = 0
            # If was in water and now melting 
            else: 
                particle.in_ice = 0

        particle.prev_sit   = particle.sit
        particle.prev_state = particle.in_ice

