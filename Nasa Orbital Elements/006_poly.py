import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

cols = ['N', 'i', 'w', 'a', 'e', 'M']

df = pd.read_csv('data/saturn.csv')

x = df['d'].values
y = df['M'].values

coefs = np.polyfit(x, y, 15)
func = np.poly1d(coefs)

print(coefs)

fig, ax = plt.subplots()
ax.plot(x, y)
ax.plot(x, func(x), c='r', alpha=0.5)
plt.show()

# for M I sould use another fit method (cyclic)


