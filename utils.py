from numpy import pi, sin, cos, tan, sqrt, arctan2, arcsin, arctan, arccos, log10

rd = pi/180
dg = 180/pi

def day(y, m, D, UT):
    d = 367*y - 7 * ( y + (m+9)//12 ) // 4 + 275*m//9 + D - 730530
    d = d + UT/24.0
    return d

def getUT(d):
    if d < 0:
        UT = (d - int(d)+1) * 24
    else:
        UT = (d - int(d)) * 24
    return UT

def floor(x):
    if x<0:
        return int(x)-1
    else:
        return int(x)
    
def rev(x):
    if x >= 360:
        while x >= 360:
            x = x - 360
    elif x < 0:
        while x < 0:
            x = 360 + x
    return x

def obl_ecl(d):
    """obliquity of the ecliptic (tilt of the Earth's axis of rotation)"""
    ecl = 23.4393 - 3.563E-7 * d
    return ecl


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
