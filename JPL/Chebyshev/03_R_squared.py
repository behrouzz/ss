import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('data/earth.csv')

x = df['JDTDB'].values
y = df['OM'].values


degrees = []
r2s = []
for i in range(200, 1000, 50):
    f = np.polynomial.Chebyshev.fit(x, y, deg=i)

    # R2
    residuals = y - f(x)
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((y - np.mean(y))**2)
    r2 = 1 - (ss_res / ss_tot)

    degrees.append(i)
    r2s.append(r2)

fig, ax = plt.subplots()
ax.plot(degrees, r2s)
ax.scatter(degrees, r2s)
plt.grid()
plt.show()



