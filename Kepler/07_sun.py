# Ref:
# https://downloads.rene-schwarz.com/download/M001-Keplerian_Orbit_Elements_to_Cartesian_State_Vectors.pdf

from numpy import array, pi, sqrt, sin, cos, arctan2
from new_utils import rd,dg, nrm,rev, obl_ecl, getE, datetime_to_jd, jd_to_datetime
from datetime import datetime
from new_oe import elements

au = 149597870700 # m

def sun_yaroo(y, m, D, UT):
    
    def day(y, m, D, UT):
        d = 367*y - 7 * ( y + (m+9)//12 ) // 4 + 275*m//9 + D - 730530
        d = d + UT/24.0
        return d
    
    def sun_oe(d):
        N = 0.0
        i = 0.0
        w = 282.9404 + 4.70935E-5 * d
        a = 1.000000  # (AU)
        e = 0.016709 - 1.151E-9 * d
        M = rev(356.0470 + 0.9856002585 * d)
        return N,i,w,a,e,M

    d = day(y, m, D, UT)
    N,i,w,a,e,M = sun_oe(d)
    print(N,i,w,a,e,M)
    # all in deg
    E = getE(e, M, 15)
    E = M + e*(180/pi) * sin(M*rd) * ( 1.0 + e * cos(M*rd) )
    xv = cos(E*rd) - e
    yv = sqrt(1-e**2) * sin(E*rd)
    r = sqrt(xv**2 + yv**2)
    v = rev(arctan2(yv, xv)*dg)
    lon = rev(w+v)
    x = r * cos(lon*rd)
    y = r * sin(lon*rd)
    z = 0.0
    return x, y, z, lon
    


def elements_to_ecliptic(N,i,w,a,e,M):
    """Get heliocentric ecplictic coordinates"""
    E = getE(e, M*dg, 15)
    E = E*rd
    v = 2*arctan2( sqrt(1+e)*sin(E/2), sqrt(1-e)*cos(E/2) )
    print('True Anomaly:', v*dg)
    r = a * (1 - e*cos(E))
    print('r =',r)
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
t = datetime(2025, 1, 1)
jd = datetime_to_jd(t)

N,i,w,a,e,M = elements('earth', jd)
print('N,i,w,a,e,M  (my fitted)')
print(N,i,w,a,e,M)
print()

# Let's check elements
import pandas as pd
ch = pd.read_csv('fitting/earth_10yr_365.csv')
ch = ch[ch['jd']==jd].iloc[0]
print('N,i,w,a,e,M  (NASA JPL)')
print(ch.N, ch.i, ch.w, ch.a, ch.e, ch.M)
#N,i,w,a,e,M = ch.N, ch.i, ch.w, ch.a, ch.e, ch.M
print()


N,i,w,M = N*rd, i*rd, w*rd, M*rd
a = a*au
x, y, z = elements_to_ecliptic(N,i,w,a,e,M)

#x, y, z, lon = sun_yaroo(2025, 1, 1, 0)
arr1 = array([x, y, z])
print('HelioEcliptic:', arr1/au)
x_equ, y_equ, z_equ = ecliptic_to_equatorial(x, y, z, jd)



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

print('Astropy:')
print('hel_ecl_xyz: ', mm_h_ec1.cartesian.xyz.value/(0.001*au))
print('hel_ecl_xyz: ', mm_h_ec2.cartesian.xyz.value/(0.001*au))

