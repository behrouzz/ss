import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#========================================================
from scipy import optimize

def fit_sin(tt, yy):
    
    def func(t, A, w, p, c):
        return A * np.sin(w*t + p) + c

    ff = np.fft.fftfreq(len(tt), (tt[1]-tt[0]))
    Fyy = abs(np.fft.fft(yy))
    guess_freq = abs(ff[np.argmax(Fyy[1:])+1])
    guess_amp = np.std(yy) * 2**0.5
    guess_offset = np.mean(yy)
    guess = np.array([guess_amp, 2*np.pi*guess_freq, 0, guess_offset])
    popt, pcov = optimize.curve_fit(func, tt, yy, p0=guess)
    A, w, p, c = popt
    return A, w, p, c
#========================================================

planets = ['mercury', 'venus', 'earth', 'mars',
           'jupiter', 'saturn', 'uranus', 'neptune']

ls = []

for i in planets:
    df = pd.read_csv('current/'+i+'.csv')
    df['time'] = pd.to_datetime(df['time'])
    df.set_index('time', inplace=True)

    x = df['day'].values
    y = df['a'].values

    A, w, p, c = fit_sin(x, y)

    dc = {}
    dc['planet'] = i
    dc['A'] = A
    dc['w'] = w
    dc['p'] = p
    dc['c'] = c
    ls.append(dc)
"""
sin_func = lambda t: A * np.sin(w*t + p) + c
dc['planet'] = i
d


fig, ax = plt.subplots()
ax.plot(x, y)
ax.plot(x, sin_func(x), c='r', alpha=0.5)
plt.show()
"""

df_a = pd.DataFrame(ls)

"""
Fitting < a >
=============
Sin is the best fit for <a> of ALL planets.

    planet             A         w          p          c
0  mercury -7.302882e-07  0.015619  -0.166265   0.387098
1    venus -4.398727e-06  0.021536   7.387954   0.723330
2    earth  8.991877e-04  0.212770   6.440250   1.000002
3     mars  6.999730e-05  0.015287   1.132729   1.523677
4  jupiter  7.842231e-04  0.001335  11.201102   5.203064
5   saturn  3.270544e-02  0.000989   6.415190   9.549749
6   uranus  6.841860e-02  0.001430   1.266557  19.234973
7  neptune  1.489743e-01  0.001395   2.680258  30.151265

sin_func = lambda t: A * np.sin(w*t + p) + c
"""
