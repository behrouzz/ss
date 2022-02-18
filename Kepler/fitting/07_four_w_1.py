import numpy as np
import pandas as pd
from scipy.signal import find_peaks
import matplotlib.pyplot as plt

df = pd.read_csv('earth_10yr_365.csv')

t = df['jd'].values
s = df['w'].values

fft = np.fft.fft(s)
N = len(s)

#f = np.linspace(0.004, 0.04, N)

a = np.abs(fft)
a = a[:N//2]
f = np.linspace(0, 1, N)[:N//2]


# ignore zero frequency
f = f[1:]
a = a[1:]

peaks, _ = find_peaks(a, height=10000)

print(f[peaks])
# [0.00547495, 0.0339447,  0.03941966]

fig, ax = plt.subplots()
ax.plot(f, a)
ax.scatter(f[peaks], a[peaks], c='r')
ax.set_xlabel("Frequency [Hz]")
ax.set_ylabel("Amplitude")
plt.show()

"""
The main frequencies are:
freqs = [0.00547495, 0.0339447,  0.03941966]

So the omegas should be:
ws = 2*np.pi*np.array(freqs)
"""
