import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('../new_current/earth.csv')
df['t'] = pd.to_datetime(df['t'])

x = df['jd'].values
y = df['w'].values

coefs = np.polyfit(x, y, 1)
func = np.poly1d(coefs)

print('Mean  :', y.mean())
print('Coefs :', str(tuple(coefs)).replace('(','[').replace(')',']'))

residuals = y - func(x)
ss_res = np.sum(residuals**2)
ss_tot = np.sum((y - np.mean(y))**2)
r_squared = 1 - (ss_res / ss_tot)
print('R2 =', r_squared)

fig, ax = plt.subplots()
ax.plot(df['t'], y)
ax.plot(df['t'], func(x), c='r', alpha=0.5)
plt.show()
