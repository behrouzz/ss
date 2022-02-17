import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('data/earth.csv')

x = df['JDTDB'].values
y = df['OM'].values

f = np.polynomial.Chebyshev.fit(x, y, deg=500)

print('Coef shape :', f.coef.shape)
print('Domain     :', f.domain)
print('Window     :', f.window)

# R2
residuals = y - f(x)
ss_res = np.sum(residuals**2)
ss_tot = np.sum((y - np.mean(y))**2)
r2 = 1 - (ss_res / ss_tot)
print('R-SQUARED  :', r2)


fig, ax = plt.subplots()
ax.plot(x, y)
ax.plot(x, f(x), c='r')
plt.show()


