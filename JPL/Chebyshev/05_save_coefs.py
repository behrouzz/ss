import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('data/earth.csv')

x = df['JDTDB'].values

dc = {}
for col in ['OM', 'IN', 'W', 'A', 'EC']:
    y = df[col].values
    f = np.polynomial.Chebyshev.fit(x, y, deg=500)
    print(col)
    print('Domain:', f.domain)
    print('Window:', f.window)
    dc[col] = f.coef

df = pd.DataFrame(dc)
df.set_index('OM').to_csv('data/coefs.csv')

# ATTENTION:
# ---------
# Don't forget to save domain (and perhaps window)
