import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

planets = ['mercury', 'venus', 'earth', 'mars',
           'jupiter', 'saturn', 'uranus', 'neptune']

df = pd.read_csv('current/jupiter.csv')
df['time'] = pd.to_datetime(df['time'])
df.set_index('time', inplace=True)

a = list(df['M'])

b = []
n = 0
for i in range(1, len(a)):
    if abs(a[i]-a[i-1])>300:
        print(i, '=>', a[i-1], a[i])
        n = n + 360
    b.append(n+a[i])
b = [a[0]] + b
df['M_'] = b

"""
5570 => 359.9633796622275 0.0095682163356004
5574 => 0.010415702658915 359.9633580124951
5579 => 359.978438992675 0.039901110184885
"""

print('-'*50)
for i in range(5560,5590):
	print(i, ':', a[i])

"""
5560 : 359.7603612470963
5561 : 359.8388477569227
5562 : 359.9084092576795
5563 : 359.9497317206942
5564 : 359.9568604106043
5565 : 359.9387457454839
5566 : 359.9147090304264
5567 : 359.9054831221474
5568 : 359.9229162423763
5569 : 359.9633796622275
5570 : 0.0095682163356004
5571 : 0.0400512936539174
5572 : 0.0405038867675263
5573 : 0.010415702658915
5574 : 359.9633580124951
5575 : 359.9214438724177
5576 : 359.9059469658282
5577 : 359.9271081374841
5578 : 359.978438992675
5579 : 0.039901110184885
5580 : 0.0883610435959414
5581 : 0.1087222182320326
5582 : 0.1003015357893088
5583 : 0.0762890169630462
5584 : 0.0577245042380699
5585 : 0.0636023546346355
5586 : 0.1010958696881474
5587 : 0.1612251692210549
5588 : 0.2235385271182447
5589 : 0.2671495853147123
"""
c = a[:5570] + [a[5570]+360] + [a[5570]+360] + [a[5572]+360] +\
    [a[5573]+360] + a[5574:5579] + [i+360 for i in a[5579:]]

df['M2'] = c

	

x = df['day'].values
y = df['M2'].values

# polynomial
my_coefs = np.polyfit(x, y, 1)
my_func = np.poly1d(my_coefs)
print(my_coefs)


fig, ax = plt.subplots()
ax.plot(x, y)
ax.plot(x, my_func(x), c='r', alpha=0.5)
plt.show()
