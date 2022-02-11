import numpy as np
from numpy import sin, cos, sqrt, arctan2
from utils import rev, obl_ecl, getE, rd, dg

def cartesian_to_spherical(x, y ,z):
    lon = rev(np.arctan2(y, x)*dg)
    lat = np.arctan2(z, np.sqrt(x**2 + y**2))*dg
    r = np.sqrt(x**2 + y**2 + z**2)
    return lon, lat, r

def ecliptic_to_equatorial(x_ecl, y_ecl, z_ecl, d):
    ecl = obl_ecl(d)
    x_equ = x_ecl
    y_equ = y_ecl * np.cos(ecl*rd) - z_ecl * np.sin(ecl*rd)
    z_equ = y_ecl * np.sin(ecl*rd) + z_ecl * np.cos(ecl*rd)
    return x_equ, y_equ, z_equ

def spherical_to_cartesian(lon, lat, r):
    x = r * np.cos(lat*rd) * np.cos(lon*rd)
    y = r * np.cos(lat*rd) * np.sin(lon*rd)
    z = r * np.sin(lat*rd)
    return x, y, z

def elements_to_ecliptic(name, N,i,w,a,e,M):
    E = getE(e, M, dp=5)
    xv = a * (cos(E*rd) - e)
    yv = a * (sqrt(1 - e**2) * sin(E*rd))
    r = sqrt(xv**2 + yv**2) # distance
    v = rev(arctan2(yv, xv)*dg) # true anomaly
    if name=='sun':
        lon = rev(w+v)
        x_ecl = r * cos(lon*rd)
        y_ecl = r * sin(lon*rd)
        z_ecl = 0.0
        return x_ecl, y_ecl, z_ecl, lon
    else:
        x_ecl = r * ( cos(N*rd) * cos((v+w)*rd) - sin(N*rd) * sin((v+w)*rd) * cos(i*rd) )
        y_ecl = r * ( sin(N*rd) * cos((v+w)*rd) + cos(N*rd) * sin((v+w)*rd) * cos(i*rd) )
        z_ecl = r * ( sin((v+w)*rd) * sin(i*rd) )
        return x_ecl, y_ecl, z_ecl

