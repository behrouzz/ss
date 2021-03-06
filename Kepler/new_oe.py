import numpy as np
from numpy import sin

def rev(x):
    if x >= 360:
        while x >= 360:
            x = x - 360
    elif x < 0:
        while x < 0:
            x = 360 + x
    return x

def _get_sin(i, dc):
    A, w, p, c = dc[i]
    func = lambda t: A * np.sin(w*t + p) + c
    return func


def _get_sin_trend(i, dc):
    A, w, p, c, m = dc[i]
    func = lambda t: A * np.sin(w*t + p) + c + m*t
    return func


def _get_fsunN(i, dc):
    w1,w2,w3, A1,A2,A3, p1,p2,p3, c = dc[i][1:]
    func = lambda t: A1*np.sin(w1*t+p1) + A2*np.sin(w2*t+p2) + A3*np.sin(w3*t+p3) + c
    return func

def _get_fsuni(i, dc):
    w1,w2,w3,w4, A1,A2,A3,A4, p1,p2,p3,p4, m, c = dc[i][1:]
    func = lambda t: A1*sin(w1*t+p1) + A2*sin(w2*t+p2) + A3*sin(w3*t+p3) + A4*sin(w4*t+p4) + m*t + c
    return func

def _get_fsunw(i, dc):
    w1,w2,w3, A1,A2,A3, p1,p2,p3, m, c = dc[i][1:]
    func = lambda t: A1*sin(w1*t+p1) + A2*sin(w2*t+p2) + A3*sin(w3*t+p3) + m*t + c
    return func

def get_func(i, dc):
    if isinstance(dc[i], tuple):
        if len(dc[i])==4:
            func = _get_sin(i, dc)
        elif len(dc[i])==5:
            func = _get_sin_trend(i, dc)
        else:
            raise Exception('There is a problem!')
    elif isinstance(dc[i], list):
        if isinstance(dc[i][0], str):
            if dc[i][0]=='fsunN':
                func = _get_fsunN(i, dc)
            elif dc[i][0]=='fsuni':
                func = _get_fsuni(i, dc)
            elif dc[i][0]=='fsunw':
                func = _get_fsunw(i, dc)
        else:
            func = np.poly1d(dc[i])
    elif isinstance(dc[i], float):
        func = lambda t: dc[i]
    else:
        raise Exception('Only list or tuple are accepted!')
    return func


def elements(name, d):
    N = get_func('N', dc_all[name])(d)
    i = get_func('i', dc_all[name])(d)
    w = get_func('w', dc_all[name])(d)
    a = get_func('a', dc_all[name])(d)
    e = get_func('e', dc_all[name])(d)
    M = rev(get_func('M', dc_all[name])(d))
    return N,i,w,a,e,M

# w1,w2,w3, A1,A2,A3, p1,p2,p3, c 
ear_dc = {
    'N' : ['fsunN', 0.03440012539754288, 0.2132808402966189, 0.24768102852601487, 6.49677967578373, 18.70122332050655, 15.758963369927942, 0.06983025510013369, 0.6124578925927411, 0.4973186473867709, 174.97547589015696],
    'i' : ['fsuni', 0.03440012539754288, 0.2132808402966189, 0.24768102852601487, 0.4626818908004742, 0.00016035100526651472, -0.0010023693702712888, 0.0008662679851036481, 0.00011350346627919421, 1.6635895632970603, 2.242420255689696, 2.1197422160640023, 2.226495441466261, 3.5789304324961146e-07, -0.8771987347518402],
    'w' : ['fsunw', 0.03440012539754288, 0.2132808402966189, 0.24768102852601487, 6.491471798065003, 18.700216981300503, 15.758233422831537, 0.0703434311306654, 0.6124602149854734, 0.4973483382104187, -0.00010780596359698544, 440.2510231700193],
    #'N' : (-20.56854255511579, 0.21368417754323393, -1516.4594773979945, 927.1665019656947, -0.00030568323546088826),
    #'i' : (-0.0011012553615547334, 0.2136982233893679, -1552.5909728787956, -0.8867004309185633, 3.617555106895709e-07),
    #'w' : (20.648043760837105, 0.21368479213803637, -1517.971405278454, -368.5384387195423, 0.0002668264147489876),
    'a' : (0.000899192533198006, 0.21276989638906812, 719.2605279990571, 1.0013470585676556, -5.465201294386334e-10),
    'e' : (-0.0006247262437435298, 0.2299832617334431, 657.9560596124007, 0.032764047157516546, -6.528090430012576e-09),
    'M' : [0.9856471837850469, -2423201.3228826406]
    }

ven_dc = {
    'N' : [2.1395300025119835e-26, -7.369171808595139e-20, -1.0363834526970677e-13, 2.5490473565048103e-07, 1.0978900467125792, -1929226.925138333],
    'i' : [-4.820198757058084e-28, 1.6623693399032919e-21, 2.3348689782529537e-15, -5.760217892244847e-09, -0.024776900952577906, 43600.06428482244],
    'w' : (-0.17398271400047768, 0.0043633578358246855, 1953.9118925950233, 7.534066003433867, 1.9283341129067293e-05),
    'a' : (4.411807607291706e-06, 0.021538100001284543, 1966.87926011159, 0.7235088858787811, -7.279699941769788e-11),
    'e' : (2.3584919367315106e-05, 0.0042364680421193, -1966.8931624260942, 0.01611569670848598, -3.802527763317829e-09),
    'M' : [1.6021214166378048, -3939142.3007621677]
    }

dc_all = {'earth':ear_dc, 'venus':ven_dc}
