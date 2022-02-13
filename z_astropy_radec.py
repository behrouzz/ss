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

print('ra : ', mm.ra.value)
print('dec: ', mm.dec.value)
print('r  : ', mm.distance.value/au)
print()

