import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def rev(x):
    if x >= 360:
        while x >= 360:
            x = x - 360
    elif x < 0:
        while x < 0:
            x = 360 + x
    return x

df = pd.read_csv('data/earth.csv')
cf = pd.read_csv('data/coefs.csv')

x = df['JDTDB'].values

cf_N = cf['OM'].values
cf_i = cf['IN'].values
cf_w = cf['W'].values
cf_a = cf['A'].values
cf_e = cf['EC'].values

dm = np.array([2458849.5, 2462502.5])

f_N = np.polynomial.chebyshev.Chebyshev(coef=cf_N, domain=dm)
f_i = np.polynomial.chebyshev.Chebyshev(coef=cf_i, domain=dm)
f_w = np.polynomial.chebyshev.Chebyshev(coef=cf_w, domain=dm)
f_a = np.polynomial.chebyshev.Chebyshev(coef=cf_a, domain=dm)
f_e = np.polynomial.chebyshev.Chebyshev(coef=cf_e, domain=dm)
f_M = np.poly1d([0.9856462168694512, -2423198.943616535]) # rev

funcs = [f_N, f_i, f_w, f_a, f_e, f_M]
cols  = ['OM','IN','W', 'A','EC','MA']

for i in range(6):
    f = funcs[i]
    y = df[cols[i]].values
    if cols[i]=='MA':
        y_pred = np.array([rev(i) for i in f(x)])
    else:
        y_pred = f(x)
    residuals = y - y_pred
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((y - np.mean(y))**2)
    r2 = 1 - (ss_res / ss_tot)
    print(cols[i], ':', r2)
    
