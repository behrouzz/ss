import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('data/earth.csv')

x = df['JDTDB'].values


for deg in range(100, 1001, 100):
    print('DEG =', deg)
    print('-------------')
    for col in ['OM', 'IN', 'W', 'A', 'EC', 'MA']:
        y = df[col].values
        f = 0
        f = np.polynomial.Chebyshev.fit(x, y, deg=deg)

        # R2
        residuals = y - f(x)
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((y - np.mean(y))**2)
        r2 = 1 - (ss_res / ss_tot)
        print(col, ':', r2)
    print('='*50)
