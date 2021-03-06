import numpy as np
from numpy import sin, cos, sqrt, arctan2
from utils import rev, obl_ecl, getE, rd, dg

def mag(x):
    return np.linalg.norm(np.array(x))

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

def state_to_element(r, v, mu=1.32712440018e+20):
    """
    N (rad) : longitude of the ascending node
    i (rad) : inclination to the ecliptic
    w (rad) : argument of perihelion
    a (m)   : semi-major axis, or mean distance from Sun
    e       : eccentricity (0=circle, 0-1=ellipse, 1=parabola)
    M (rad) : mean anomaly
    """
    h = np.cross(r, v)

    # eccentricity vector
    ev = (np.cross(v,h)/mu) - (r/mag(r))
    e = mag(ev) # orbit eccentricity

    # vector pointing towards the ascending node
    n = np.cross(np.array([0,0,1]), h)

    # true anomaly (in radian)
    tmp = np.inner(ev,r)/(e*mag(r))

    if np.inner(r,v) >= 0:
        vv = np.arccos(tmp)
    else:
        vv = 2*np.pi - np.arccos(tmp)

    # inclination (in radian)
    i = np.arccos(h[-1]/mag(h))

    # eccentric anomaly
    E = 2 * np.arctan(np.tan(vv/2)/np.sqrt((1+e)/(1-e)))

    # longitude of the ascending node (in radian)
    if n[1] >= 0:
        N = np.arccos(n[0]/mag(n))
    else:
        N = 2*pi - np.arccos(n[0]/mag(n))

    # argument of perihelion
    tmp = np.inner(n, ev) / (mag(n)*e)
    if ev[-1] >= 0:
        w = np.arccos(tmp)
    else:
        w = 2*np.pi - np.arccos(tmp)

    # mean anomaly
    M = E - e*np.sin(E)

    # semi-major axis
    a = 1 / ((2/mag(r)) - (mag(v)**2 / mu))

    return N,i,w,a,e,M
