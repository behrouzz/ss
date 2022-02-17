import numpy as np
from numpy import array, pi, sqrt, sin, cos, arctan2
import pandas as pd
from datetime import datetime
from earth_oe import func
from utils import rd,dg, nrm,rev, obl_ecl, getE, datetime_to_jd, jd_to_datetime

def elements_to_ecliptic(N,i,w,a,e,M):
    """Get heliocentric ecplictic coordinates"""
    E = getE(e, M*dg, 15)
    E = E*rd
    v = 2*arctan2( sqrt(1+e)*sin(E/2), sqrt(1-e)*cos(E/2) )
    r = a * (1 - e*cos(E))
    ox = r * cos(v)
    oy = r * sin(v)
    # oz = 0
    x = ox*(cos(w)*cos(N)-sin(w)*cos(i)*sin(N)) - oy*(sin(w)*cos(N)+cos(w)*cos(i)*sin(N))
    y = ox*(cos(w)*sin(N)+sin(w)*cos(i)*cos(N)) + oy*(cos(w)*cos(i)*cos(N)-sin(w)*sin(N))
    z = ox*sin(w)*sin(i) + oy*cos(w)*sin(i)
    return x, y, z

au = 149597870700 # m

t = datetime(2025,8, 1)
d = datetime_to_jd(t)

N,i,w,a,e,M = func['N'](d), func['i'](d), func['w'](d), func['a'](d), func['e'](d), func['M'](d)
M = rev(M)
N,i,w,M = N*rd, i*rd, w*rd, M*rd
a = a*au
x, y, z = elements_to_ecliptic(N,i,w,a,e,M)
arr1 = np.array([x, y, z])
print('Khodam : ', arr1/au)


from astropy.coordinates import get_body, HeliocentricTrueEcliptic, HeliocentricMeanEcliptic
from astropy.coordinates import solar_system_ephemeris
from astropy.time import Time
tt = Time(t)
hte = HeliocentricTrueEcliptic(obstime=tt, equinox=tt)
hme = HeliocentricMeanEcliptic(obstime=tt, equinox=tt)
with solar_system_ephemeris.set('jpl'):
    mm = get_body('earth', tt)
mm_h_ec1 = mm.transform_to(hte)
mm_h_ec2 = mm.transform_to(hme)
print('Astropy: ', mm_h_ec1.cartesian.xyz.value/(0.001*au))
print('Astropy: ', mm_h_ec2.cartesian.xyz.value/(0.001*au))
