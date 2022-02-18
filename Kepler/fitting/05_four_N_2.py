import numpy as np
from numpy import array, pi, sin
import pandas as pd
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

df = pd.read_csv('earth_10yr_365.csv')

x = df['jd'].values
y = df['N'].values

# The main frequencies are:
fs = array([0.00547495, 0.0339447,  0.03941966])

# So the omegas should be:
w1,w2,w3 = 2*pi*fs



def fit_func(tt, yy):

    def func(t, A1,A2,A3, p1,p2,p3, c):
        return  A1*sin(w1*t+p1) + A2*sin(w2*t+p2) + A3*sin(w3*t+p3) + c
    
    p0 = [18,18,18, 0,0,0, 175]
    
    popt, pcov = curve_fit(func, tt, yy, maxfev=50000, p0=p0)
    return popt

popt = fit_func(x, y)
print(popt)
A1,A2,A3, p1,p2,p3, c = popt

def the_func(t):
    return A1*sin(w1*t+p1) + A2*sin(w2*t+p2) + A3*sin(w3*t+p3) + c

residuals = y - the_func(x)
ss_res = np.sum(residuals**2)
ss_tot = np.sum((y - np.mean(y))**2)
r_squared = 1 - (ss_res / ss_tot)
print('R2 =', r_squared)

fig, ax = plt.subplots()
ax.plot(x, y)
ax.plot(x, the_func(x), c='r', alpha=0.5)
plt.show()

"""
w1,w2,w3 = (0.03440012539754288, 0.2132808402966189, 0.24768102852601487)

A1,A2,A3, p1,p2,p3, c = 
[6.49677968e+00 1.87012233e+01 1.57589634e+01 6.98302551e-02
 6.12457893e-01 4.97318647e-01 1.74975476e+02]
 
R2 = 0.7644438239716355
"""
