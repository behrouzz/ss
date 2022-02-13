from astropy.coordinates import get_sun
from astropy.coordinates import get_body, EarthLocation, AltAz, GeocentricTrueEcliptic, HeliocentricTrueEcliptic
from astropy.coordinates import get_body_barycentric, solar_system_ephemeris
from astropy.time import Time
import numpy as np

t = Time('2022-02-11 10:00:00')

gg = GeocentricTrueEcliptic(obstime=t)#, equinox=t)
hh = HeliocentricTrueEcliptic(obstime=t)#, equinox=t)

with solar_system_ephemeris.set('jpl'):
    mm = get_body('venus', t)

mm_g_ec = mm.transform_to(gg)
mm_h_ec = mm.transform_to(hh)

au = 149597870.700

print('Astropy:')
"""
print('ra : ', mm.ra.value)
print('dec: ', mm.dec.value)
print('r  : ', mm.distance.value/au)
print('equ_xyz: ', mm.cartesian.xyz.value/au)
print('geo_ecl_xyz: ', mm_g_ec.cartesian.xyz.value/au)
print('geo_ecl_lon: ', mm_g_ec.lon.value)
print('geo_ecl_lat: ', mm_g_ec.lat.value)
print('geo_ecl_r  : ', mm_g_ec.distance.value/au)
"""
print('hel_ecl_xyz: ', mm_h_ec.cartesian.xyz.value/au)
#print('hel_ecl_lon: ', mm_h_ec.lon.value)
#print('hel_ecl_lat: ', mm_h_ec.lat.value)
#print('hel_ecl_r  : ', mm_h_ec.distance.value/au)
print()

from conversions import elements_to_ecliptic
from new_oe import elements
from orbital_elements import *
from utils import day

name = 'venus'
d = day(2022, 2, 11, 10)

N,i,w,a,e,M = elements(name, d)
N_,i_,w_,a_,e_,M_ = venus_oe(d)

x_ecl, y_ecl, z_ecl = elements_to_ecliptic(name, N,i,w,a,e,M)
x_ecl_, y_ecl_, z_ecl_ = elements_to_ecliptic(name, N_,i_,w_,a_,e_,M_)

print('OLD xyz_h_ecl:')
print(x_ecl, y_ecl, z_ecl)

print()
print('NEW xyz_h_ecl:')
print(x_ecl_, y_ecl_, z_ecl_)
