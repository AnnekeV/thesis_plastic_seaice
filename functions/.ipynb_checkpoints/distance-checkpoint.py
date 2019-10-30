from numpy import arccos, sin, cos, pi

def great_circle(lon1, lat1, lon2, lat2):
    '''Calculates the distance between two points on Earth in m'''
    R     = 6371.0e3 # m radius earth
    lon1  = lon1/180.*pi
    lat1  = lat1/180.*pi
    lon2  = lon2/180.*pi
    lat2  = lat2/180.*pi
    dlon  = abs(lon2-lon1)
    angle = arccos( sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2)*cos(dlon) )
    return angle * R