import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

planets = ['mercury', 'venus', 'earth', 'mars',
           'jupiter', 'saturn', 'uranus', 'neptune']

df = pd.read_csv('current/'+planets[4]+'.csv')
df['time'] = pd.to_datetime(df['time'])
df.set_index('time', inplace=True)

a = list(df['M'])
b = []
n = 0
for i in range(1, len(a)):
    if abs(a[i]-a[i-1])>300:
        n = n + 360
    b.append(n+a[i])
b = [a[0]] + b
df['M_'] = b

x = df['day'].values
y = df['M'].values

# polynomial
my_coefs = np.polyfit(x, y, 1)
my_func = np.poly1d(my_coefs)
print(my_coefs)


fig, ax = plt.subplots()
ax.plot(x, y)
ax.plot(x, my_func(x), c='r', alpha=0.5)
plt.show()
