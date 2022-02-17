import numpy as np
import pandas as pd
from utils import rev

cf = pd.read_csv('../JPL/Chebyshev/data/coefs.csv')

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
f_M = np.poly1d([0.9856462168694512, -2423198.943616535])
#f_M = lambda t: rev(f_M(t))

func = {'N':f_N, 'i':f_i, 'w':f_w, 'a':f_a, 'e':f_e, 'M':f_M}
