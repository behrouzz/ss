import pandas as pd
from numpy import pi, sqrt, sin, cos, arctan2
from new_utils import rd,dg, nrm, obl_ecl, getE, datetime_to_jd, jd_to_datetime

au = 149597870700 # m
rd = pi/180
dg = 180/pi

dc = {'TBD': '2020-01-01 00:00:00.000',
      'JDTDB': 2458849.5,
      'EC': 0.01712127709447231,
      'QR': 0.9832461589959328,
      'IN': 0.002777009999140959,
      'OM': 159.6966862721255,
      'W': 304.3574361359509,
      'Tp': 2458853.7319454737,
      'N': 0.9850567216484872,
      'MA': 355.8312936655853,
      'TA': 355.685564731535,
      'A': 1.0003738366513,
      'AD': 1.017501514306667,
      'PR': 365.4611882628869}



# Nasa JPL
pos = (-0.1663457602264735, 0.9691203923795578, -4.12553307140163e-05)
print(pos)


N = dc['OM'] *rd
i = dc['IN'] *rd
w = dc['W'] *rd
a = dc['A'] *au
e = dc['EC']
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

# 2.9591293262935277E-04 au^3/d^2

