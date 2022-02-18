import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


# Carrier
A, w, p, c = [-20.569851652817754, 0.2136839635257106,
              -1515.9328827259324, 174.97910307507195]

def fit_func(tt, yy):
    
    def func(t, Am, Wm, Pm):
        return  (A + Am*np.sin(Wm*x+Pm)) * np.sin(w*x+p) + c
    
    p0 = [200, 0.03583347766207066, 0]
    
    popt, pcov = curve_fit(func, tt, yy, maxfev=50000, p0=p0)
    return popt


df = pd.read_csv('../new_current/earth.csv')
df['t'] = pd.to_datetime(df['t'])

x = df['jd'].values
y = df['N'].values



popt = fit_func(x, y)
Am, Wm, Pm = popt
coefs = Am, Wm, Pm
print(coefs)

def the_func(t):
    return (A + Am*np.sin(Wm*x+Pm)) * np.sin(w*x+p) + c

residuals = y - the_func(x)
ss_res = np.sum(residuals**2)
ss_tot = np.sum((y - np.mean(y))**2)
r_squared = 1 - (ss_res / ss_tot)
print('R2 =', r_squared)

fig, ax = plt.subplots()
ax.plot(df['t'], y)
ax.plot(df['t'], the_func(x), c='r', alpha=0.5)
plt.grid()
plt.show()
