import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('data/earth.csv')
cf = pd.read_csv('data/coefs.csv')

x = df['JDTDB'].values
y = df['OM'].values


OMcf = cf['OM'].values
OMdm = np.array([2458849.5, 2462502.5])

f = np.polynomial.chebyshev.Chebyshev(coef=OMcf, domain=OMdm)

fig, ax = plt.subplots()
ax.plot(x, y)
ax.plot(x, f(x), c='r', alpha=0.5)
plt.show()
