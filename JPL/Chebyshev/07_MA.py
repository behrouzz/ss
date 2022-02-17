import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('data/earth.csv')

a = list(df['MA'])
b = []
n = 0
for i in range(1, len(a)):
    if abs(a[i]-a[i-1])>300:
        n = n + 360
    b.append(n+a[i])
b = [a[0]] + b
df['M_'] = b

x = df['JDTDB'].values
y = df['M_'].values

# polynomial
coefs = np.polyfit(x, y, 1)
f = np.poly1d(coefs)
print(coefs)

residuals = y - f(x)
ss_res = np.sum(residuals**2)
ss_tot = np.sum((y - np.mean(y))**2)
r2 = 1 - (ss_res / ss_tot)
print('R2:', r2)


fig, ax = plt.subplots()
ax.plot(x, y)
ax.plot(x, f(x), c='r', alpha=0.5)
plt.show()

"""
[ 9.85646217e-01 -2.42319894e+06]
R2: 0.9999975309895057
"""
