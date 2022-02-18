import pandas as pd
from numpy import pi, sqrt, sin, cos, arctan2
from new_utils import rd,dg, nrm, obl_ecl, getE, datetime_to_jd, jd_to_datetime

au = 149597870700 # m
rd = pi/180
dg = 180/pi

dc = {'TBD': '2020-01-01 00:00:00.000',
      'JDTDB': 2458849.5,
      'EC': 0.0067446119596257,
      'QR': 0.7184447354793101,
      'IN': 3.394579443347534,
      'OM': 76.62474799040017,
      'W': 54.90669305494785,
      'Tp': 2458928.636045873,
      'N': 1.602160239452148,
      'MA': 233.2113737942087,
      'TA': 232.5955379827174,
      'A': 0.7233232702585717,
      'AD': 0.7282018050378334,
      'PR': 224.6966259274418}

# Nasa JPL
pos = (0.7232003000354648, 0.05254837843508311, -0.04101276226074742)
print(pos)

mu = 1.32712440041e+20 # sun (m3/s2)
M0 = 5.011477187404461E+01 *rd # venus
t0 = 2451545.0



t = dc['JDTDB']
N = dc['OM'] *rd
i = dc['IN'] *rd
w = dc['W'] *rd
a = dc['A'] *au
e = dc['EC']

#dt = (t - t0) * 86400
#M = nrm(M0 + dt*sqrt(mu/a**3))
M = dc['MA'] *rd

E = getE(e, M*dg, 20)
E = E*rd

v = 2 * arctan2( sqrt(1+e)*sin(E/2), sqrt(1-e)*cos(E/2) )
r = a * (1 - e*cos(E))
ox = r * cos(v)
oy = r * sin(v)
# oz = 0
x = ox*(cos(w)*cos(N)-sin(w)*cos(i)*sin(N)) - oy*(sin(w)*cos(N)+cos(w)*cos(i)*sin(N))
y = ox*(cos(w)*sin(N)+sin(w)*cos(i)*cos(N)) + oy*(cos(w)*cos(i)*cos(N)-sin(w)*sin(N))
z = ox*sin(w)*sin(i) + oy*cos(w)*sin(i)
print((x/au,y/au,z/au))
"""
# 2.9591293262935277E-04 au^3/d^2
"""
