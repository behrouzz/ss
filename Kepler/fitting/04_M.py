import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('../new_current/earth.csv')
df.set_index('t', inplace=True)


a = list(df['M'])
b = []
n = 0
for i in range(1, len(a)):
    if abs(a[i]-a[i-1])>300:
        n = n + 360
    b.append(n+a[i])
b = [a[0]] + b
df['M_'] = b

x = df['jd'].values
y = df['M_'].values

coefs = np.polyfit(x, y, 1)
func = np.poly1d(coefs)

print('Mean  :', y.mean())
print('Coefs :', str(tuple(coefs)).replace('(','[').replace(')',']'))

fig, ax = plt.subplots()
ax.plot(x, y)
ax.plot(x, func(x), c='r', alpha=0.5)
plt.show()
