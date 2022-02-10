import numpy as np
from utils import rev, obl_ecl, rd, dg

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
