import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

planets = ['mercury', 'venus', 'earth', 'mars',
           'jupiter', 'saturn', 'uranus', 'neptune']

df = pd.read_csv('current/'+planets[1]+'.csv')
df['time'] = pd.to_datetime(df['time'])
df.set_index('time', inplace=True)

x = df['day'].values
y = df['a'].values

# Linrear Regression
m, b = np.polyfit(x, y, 1)
print('m:', m)
print('b:', b)

# Polynomial
coefs = np.polyfit(x, y, 15)
func = np.poly1d(coefs)
print(coefs)

fig, ax = plt.subplots()
ax.plot(x, y)
ax.plot(x, m*x+b, c='y', alpha=0.5)
ax.plot(x, func(x), c='r', alpha=0.5)
plt.show()
