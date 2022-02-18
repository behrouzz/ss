import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def fit_func(tt, yy):
    
    def func(t, A, w, p, c):
        return A * np.sin(w*t + p) + c

    ff = np.fft.fftfreq(len(tt), (tt[1]-tt[0]))
    Fyy = abs(np.fft.fft(yy))

    guess_freq = abs(ff[np.argmax(Fyy[1:])+1])
    guess_amp = np.std(yy) * 2**0.5
    guess_offset = np.mean(yy)

    guess = np.array([guess_amp, 2*np.pi*guess_freq, 0, guess_offset])
    
    popt, pcov = curve_fit(func, tt, yy, p0=guess, maxfev=5000)
    A, w, p, c = popt
    return A, w, p, c

df = pd.read_csv('../new_current/earth.csv')
df['t'] = pd.to_datetime(df['t'])

x = df['jd'].values
y = df['N'].values

# sin trend
A, w, p, c = fit_func(x, y)
coefs = A, w, p, c
print(coefs)

def sin_trend_func(x):
    return A * np.sin(w*x + p) + c

residuals = y - sin_trend_func(x)
ss_res = np.sum(residuals**2)
ss_tot = np.sum((y - np.mean(y))**2)
r_squared = 1 - (ss_res / ss_tot)
print('R2 =', r_squared)

fig, ax = plt.subplots()
ax.plot(df['t'], y)
ax.plot(df['t'], sin_trend_func(x), c='r', alpha=0.5)
plt.show()
