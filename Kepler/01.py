from numpy import pi, sqrt, sin, cos, arctan2
from utils import rev, getE, datetime_to_jd, jd_to_datetime, day, rd,dg
from datetime import datetime

mu_sun = 1.32712440041e+20 # m3/s2
au = 149597870700 # m
"""
mu_sun = 0.00029591220819207774 AU3 / d2
dg = 180 / np.pi
"""


def nrm(x):
    if x >= pi:
        while x >= pi:
            x = x - pi
    elif x < 0:
        while x < 0:
            x = pi + x
    return x

def state_vector(t, N,i,w,a,e, M0, mu=mu_sun, t0=2451545.0):
    dt = (t - t0) * 86400
    M = nrm(M0 + dt*sqrt(mu/a**3))
    E = getE(e, M)
    v = 2 * arctan2( sqrt(1+e)*sin(E/2), sqrt(1-e)*cos(E/2) )
    r = a * (1 - e*cos(E))
    ox = r * cos(v)
    oy = r * sin(v)
    oz = 0
    xh = ox*(cos(w)*cos(N)-sin(w)*cos(i)*sin(N)) - oy*(sin(w)*cos(N)+cos(w)*cos(i)*sin(N))
    yh = ox*(cos(w)*sin(N)+sin(w)*cos(i)*cos(N)) + oy*(cos(w)*cos(i)*cos(N)-sin(w)*sin(N))
    zh = ox*sin(w)*sin(i) + oy*cos(w)*sin(i)
    return xh, yh, zh


def venus_oe(d):
    N =  76.6799 + 2.46590E-5 * d
    i = 3.3946 + 2.75E-8 * d
    w =  54.8910 + 1.38374E-5 * d
    a = 0.723330  # (AU)
    e = 0.006773 - 1.302E-9 * d
    M =  rev(48.0052 + 1.6021302244 * d)
    return N,i,w,a,e,M

def yaroo(N,i,w,a,e,M):    
    E = getE(e, M, dp=5)
    xv = a * (cos(E*rd) - e)
    yv = a * (sqrt(1 - e**2) * sin(E*rd))
    r = sqrt(xv**2 + yv**2) # distance
    v = rev(arctan2(yv, xv)*dg) # true anomaly
    x_ecl = r * ( cos(N*rd) * cos((v+w)*rd) - sin(N*rd) * sin((v+w)*rd) * cos(i*rd) )
    y_ecl = r * ( sin(N*rd) * cos((v+w)*rd) + cos(N*rd) * sin((v+w)*rd) * cos(i*rd) )
    z_ecl = r * ( sin((v+w)*rd) * sin(i*rd) )
    return x_ecl, y_ecl, z_ecl

t = datetime(1990, 4, 19)

d = day(1990, 4, 19, 0)
N,i,w,a,e,M = venus_oe(d)

# yaroo
x_ecl, y_ecl, z_ecl = yaroo(N,i,w,a,e,M)

# man
a = a*au
N,i,w = N*rd, i*rd, w*rd
MA= 5.011477187404461E+01 *rd # venus
jd = datetime_to_jd(t)
xh, yh, zh = state_vector(jd, N,i,w,a,e, MA)

print(xh/au, yh/au, zh/au)
print(x_ecl, y_ecl, z_ecl)
