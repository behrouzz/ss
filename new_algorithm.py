# Ref:
# https://downloads.rene-schwarz.com/download/M001-Keplerian_Orbit_Elements_to_Cartesian_State_Vectors.pdf

from numpy import array, pi, sqrt, sin, cos, arctan2
from utils import datetime_to_jd, jd_to_datetime
from datetime import datetime

au = 149597870700 # m
mu_sun = 1.32712440041e+20 # m3/s2
M0_mercury = 1.747958829169579E+02 *rd
M0_venus = 5.011477187404461E+01 *rd
M0_earth = 3.586172562416435E+02 *rd
M0_mars = 1.935648274483629E+01 *rd
M0_jupiter = 1.881562710723796E+01 *rd
M0_saturn = 3.204255572596800E+02 *rd
M0_uranus = 1.428889650363792E+02 *rd
M0_neptune = 2.661050265674106E+02 *rd
M0_moon = 1.466732747774804E+02 *rd #Geocentric ecliptic


rd = pi/180
dg = 180/pi



def nrm(x):
    if x >= pi:
        while x >= pi:
            x = x - pi
    elif x < 0:
        while x < 0:
            x = pi + x
    return x


def floor(x):
    if x<0:
        return int(x)-1
    else:
        return int(x)

def getE(ec, m, dp=5):
    # http://www.jgiesen.de/kepler/kepler.html
    K = pi/180
    maxIter=30
    i=0
    delta = 10**-dp
    m = m/360.0
    m = 2 * pi * (m-floor(m))
    if ec<0.8:
        E=m
    else:
        E=pi
    F = E - ec*sin(m) - m
    while ((abs(F)>delta) and (i<maxIter)):
        E = E - F/(1-ec*cos(E))
        F = E - ec*sin(E) - m
        i = i + 1
    E = E/K
    return round(E*(10**dp)) / (10**dp)

def elements_to_ecliptic(t, N,i,w,a,e, M0, mu=mu_sun, t0=2451545.0):
    """Get heliocentric ecplictic coordinates"""
    dt = (t - t0) * 86400
    M = nrm(M0 + dt*sqrt(mu/a**3))
    E = getE(e, M)
    v = 2 * arctan2( sqrt(1+e)*sin(E/2), sqrt(1-e)*cos(E/2) )
    r = a * (1 - e*cos(E))
    ox = r * cos(v)
    oy = r * sin(v)
    oz = 0
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

"""
N,i,w,a,e,M0 = ...
t = datetime(1990, 4, 19)
jd = datetime_to_jd(t)

x, y, z = elements_to_ecliptic(t=jd, N=N*rd, i=i*rd, w=w*rd, a=a*au, e=e, M0=M0_venus)
x_equ, y_equ, z_equ = ecliptic_to_equatorial(x, y, z, jd)
arr = array([x_equ, y_equ, z_equ])
print(arr/au)
"""
