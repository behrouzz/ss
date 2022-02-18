import numpy as np
from numpy import array, pi, sin
import pandas as pd
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

df = pd.read_csv('earth_10yr_365.csv')

x = df['jd'].values
y = df['i'].values

# The main frequencies are:
fs = array([0.00547495, 0.0339447 , 0.03941966, 0.07363811])

# So the omegas should be:
w1,w2,w3,w4 = 2*pi*fs



def fit_func(tt, yy):

    def func(t, A1,A2,A3,A4, p1,p2,p3,p4, m, c):
        return  A1*sin(w1*t+p1) + A2*sin(w2*t+p2) + A3*sin(w3*t+p3) + A4*sin(w4*t+p4) + m*t + c
    
    p0 = (0.00017268569217032112, -0.0010020496242858664, 0.0008679414004990366, 0.00011464311271332716, 1.764321438180201, 2.2390995811071677, 2.1224384925912814, 2.2353018109724756, 0.0034600873394944618, 3.617555106895709e-07)
    
    popt, pcov = curve_fit(func, tt, yy, maxfev=50000, p0=p0)
    return popt

popt = fit_func(x, y)
print(popt)
A1,A2,A3,A4, p1,p2,p3,p4, m, c = popt

def the_func(t):
    return A1*sin(w1*t+p1) + A2*sin(w2*t+p2) + A3*sin(w3*t+p3) + A4*sin(w4*t+p4) + m*t + c

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
R2 = 0.8351742114016133

w1,w2,w3,w4 = (0.03440012539754288, 0.2132808402966189, 0.24768102852601487, 0.4626818908004742)

A1,A2,A3,A4, p1,p2,p3,p4, m, c =
[0.00016035100526651472, -0.0010023693702712888, 0.0008662679851036481, 0.00011350346627919421, 1.6635895632970603, 2.242420255689696, 2.1197422160640023, 2.226495441466261, 3.5789304324961146e-07, -0.8771987347518402]

A1*sin(w1*t+p1) + A2*sin(w2*t+p2) + A3*sin(w3*t+p3) + A4*sin(w4*t+p4) + m*t + c

"""
