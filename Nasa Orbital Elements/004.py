import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

cols = ['N', 'i', 'w', 'a', 'e', 'M']

df = pd.read_csv('data/uranus.csv')

x = df['d'].values
y = df['a'].values

m,b = np.polyfit(x, y, 1)
print('m:', m)
print('b:', b)

fig, ax = plt.subplots()

ax.plot(x, y)
ax.plot(x, m*x+b, c='r', alpha=0.5)


plt.show()


