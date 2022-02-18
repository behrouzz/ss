# Ref:
# https://downloads.rene-schwarz.com/download/M001-Keplerian_Orbit_Elements_to_Cartesian_State_Vectors.pdf

from numpy import array, pi, sqrt, sin, cos, arctan2
from new_utils import rd,dg, nrm, obl_ecl, getE, datetime_to_jd, jd_to_datetime
from datetime import datetime
from new_oe import elements

au = 149597870700 # m

def elements_to_ecliptic(N,i,w,a,e,M):
    """Get heliocentric ecplictic coordinates"""
    E = getE(e, M*dg)
    E = E*rd
    v = 2 * arctan2( sqrt(1+e)*sin(E/2), sqrt(1-e)*cos(E/2) )
    r = a * (1 - e*cos(E))
    ox = r * cos(v)
    oy = r * sin(v)
    # oz = 0
    x = ox*(cos(w)*cos(N)-sin(w)*cos(i)*sin(N)) - oy*(sin(w)*cos(N)+cos(w)*cos(i)*sin(N))
    y = ox*(cos(w)*sin(N)+sin(w)*cos(i)*cos(N)) + oy*(cos(w)*cos(i)*cos(N)-sin(w)*sin(N))
    z = ox*sin(w)*sin(i) + oy*cos(w)*sin(i)
    return x, y, z


def obl_ecl(t):
    incl = 24.312603556024428 - 3.5623324732981026e-07 * t
    return incl*rd

def ecliptic_to_equatorial(x_ecl, y_ecl, z_ecl, t):
    ecl = obl_ecl(t)
    x_equ = x_ecl
    y_equ = y_ecl*cos(ecl) - z_ecl*sin(ecl)
    z_equ = y_ecl*sin(ecl) + z_ecl*cos(ecl)
    return x_equ, y_equ, z_equ


# 2022-02-13 10:27:50.400
t = datetime(2020, 1, 1)
jd = datetime_to_jd(t)

N,i,w,a,e,M = elements('venus', jd)

N,i,w,M = N*rd, i*rd, w*rd, M*rd
a = a*au

x, y, z = elements_to_ecliptic(N,i,w,a,e,M)
arr1 = array([x, y, z])
print('HelioEcliptic:', arr1/au)
x_equ, y_equ, z_equ = ecliptic_to_equatorial(x, y, z, jd)



from astropy.coordinates import get_body, HeliocentricTrueEcliptic, HeliocentricMeanEcliptic
from astropy.coordinates import solar_system_ephemeris
from astropy.time import Time

tt = Time(t)
hte = HeliocentricTrueEcliptic(obstime=tt)#, equinox=t)
hme = HeliocentricMeanEcliptic(obstime=tt)#, equinox=t)

with solar_system_ephemeris.set('jpl'):
    mm = get_body('venus', tt)

mm_h_ec1 = mm.transform_to(hte)
mm_h_ec2 = mm.transform_to(hme)

print('Astropy:')
print('hel_ecl_xyz: ', mm_h_ec1.cartesian.xyz.value/(0.001*au))
print('hel_ecl_xyz: ', mm_h_ec2.cartesian.xyz.value/(0.001*au))

