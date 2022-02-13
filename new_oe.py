from utils import day, rev
from orbital_elements import venus_oe
import numpy as np

def get_func(i, dc):
    if isinstance(dc[i], tuple):
        if len(dc[i])==4:
            func = _get_sin(i, dc)
        elif len(dc[i])==5:
            func = _get_sin_trend(i, dc)
        else:
            raise Exception('There is a problem!')
    elif isinstance(dc[i], list):
        func = np.poly1d(dc[i])
    else:
        raise Exception('Only list or tuple are accepted!')
    return func

def _get_sin(i, dc):
    A, w, p, c = dc[i]
    func = lambda t: A * np.sin(w*t + p) + c
    return func

def _get_sin_trend(i, dc):
    A, w, p, c, m = dc[i]
    func = lambda t: A * np.sin(w*t + p) + c + m*t
    return func


d = day(2022, 2, 11, 10)

N,i,w,a,e,M = venus_oe(d)
M = rev(M)

print('N:', N)
print('i:', i)
print('w:', w)
print('a:', a)
print('e:', e)
print('M:', M)
print()

sun_dc = {
    'N' : (20.552775873722602, 0.2136832880120676, -3.580978347278808, 174.97389301955806),
    'i' : (0.001101954782290859, 0.21369711394557486, -5.275796592097575, 0.00015489408592731231, 3.621325113811721e-07),
    'w' : (-20.630936449162093, 0.21368386463476155, -3.5858347681354386, 288.0401174521317),
    'a' : (0.0008991876561645148, 0.21276978132643024, 6.440250190274306, 1.000002237215457),
    'e' : (0.000624492722824836, 0.2299825556661901, 6.272769979941623, 0.01670036421395417),
    'M' : [ 9.85645883e-01, -6.84436431e+03]
    }

ven_dc = {
    'N' : [8.883608712132133e-33, -5.6784062266386335e-28, 1.4587320685450877e-23, -1.759389614126735e-19, 5.272692894172828e-16, 1.2541053620125329e-11, -1.8200822892364454e-07, 0.0011218340122082697, -3.4777646043631467, 4507.8193294953935],
    'i' : [1.4064221204108598e-34, -8.521753182907468e-30, 2.0605705954325883e-25, -2.2962854986367836e-21, 5.114135925769306e-18, 1.6976569276153283e-13, -2.191256592521741e-09, 1.2437865744731987e-05, -0.03561566031280168, 45.19853435473836],
    'w' : [3.885724519806261e-73, -1.6805644141761093e-68, 1.929439129881774e-64, 7.812718352683747e-61, -1.5921442042476676e-56, -1.4821551680656297e-52, 6.610922236196345e-49, 2.0023949729735412e-44, 8.917924029651088e-41, -1.3918385119037133e-36, -2.099452830930911e-32, -1.7767295899649145e-29, 2.147006898584575e-24, 1.657208323216468e-20, -1.3988074666227706e-16, -2.418704223238284e-12, 1.0269141962488598e-08, 0.00026047171826611663, -2.7316389611311918, 10215.801090948815, -13901043.735505505],
    'a' : (-4.398726985167627e-06, 0.021535750705983677, 7.38795406021702, 0.7233297518634283),
    'e' : [-1.0381777086897013e-76, 4.655154633928052e-72, -5.689854175005972e-68, -1.8959279377799044e-64, 4.9285125355206404e-60, 4.0370676752167875e-56, -2.5458926933730127e-52, -6.076432248109977e-48, -2.064776467667157e-44, 4.905821239612119e-40, 6.32215778587972e-36, -3.328436567169255e-33, -7.314701325241423e-28, -4.705581172143132e-24, 5.474258549178694e-20, 7.777640567573351e-16, -4.398865291026659e-12, -8.695678283345249e-08, 0.0010099689841294244, -4.01899261968326, 5774.214226892281],
    'M' : [ 1.60212135e+00, -1.14719551e+04]
    }

print('N:', get_func('N', ven_dc)(d))
print('i:', get_func('i', ven_dc)(d))
print('w:', get_func('w', ven_dc)(d))
print('a:', get_func('a', ven_dc)(d))
print('e:', get_func('e', ven_dc)(d))
print('M:', rev(get_func('M', ven_dc)(d)))


