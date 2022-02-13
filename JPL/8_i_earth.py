import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def fit_func(tt, yy):
    
    def func(t, A, w, p, c, m):
        return A * np.sin(w*t + p) + c + m*t

    ff = np.fft.fftfreq(len(tt), (tt[1]-tt[0]))
    Fyy = abs(np.fft.fft(yy))

    guess_freq = abs(ff[np.argmax(Fyy[1:])+1])
    guess_amp = np.std(yy) * 2**0.5
    guess_offset = np.mean(yy)

    guess = np.array([guess_amp, 2*np.pi*guess_freq, 0, guess_offset, 0])
    
    popt, pcov = curve_fit(func, tt, yy, p0=guess)
    A, w, p, c, m = popt
    return A, w, p, c, m



planets = ['mercury', 'venus', 'earth', 'mars',
           'jupiter', 'saturn', 'uranus', 'neptune']
i = 2

df = pd.read_csv('current/'+planets[i]+'.csv')
df['time'] = pd.to_datetime(df['time'])
df.set_index('time', inplace=True)

x = df['day'].values
y = df['i'].values

# polynomial
my_coefs = np.polyfit(x, y, 9)
my_func = np.poly1d(my_coefs)
print(planets[i])
print(str(list(my_coefs)))

# sin trend
A, w, p, c, m = fit_func(x, y)
def sin_trend_func(x):
    return A * np.sin(w*x + p) + c + m*x


fig, ax = plt.subplots()
ax.plot(x, y)
ax.plot(x, my_func(x), c='y', alpha=0.5)
ax.plot(x, sin_trend_func(x), c='r', alpha=0.5)
plt.show()
