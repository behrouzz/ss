"""
!$$SOF
MAKE_EPHEM=YES
COMMAND=299
EPHEM_TYPE=ELEMENTS
CENTER='500@10'
START_TIME='1990-04-19 00:00:00'
STOP_TIME='1990-04-19 01:00:00'
STEP_SIZE='1'
REF_SYSTEM='ICRF'
REF_PLANE='ECLIPTIC'
OUT_UNITS='AU-D'
ELM_LABELS='YES'
TP_TYPE='ABSOLUTE'
CSV_FORMAT='YES'
OBJ_DATA='NO'



    JDTDB    Julian Day Number, Barycentric Dynamical Time
      EC     Eccentricity, e
      QR     Periapsis distance, q (au)
      IN     Inclination w.r.t X-Y plane, i (degrees)
      OM     Longitude of Ascending Node, OMEGA, (degrees)
      W      Argument of Perifocus, w (degrees)
      Tp     Time of periapsis (Julian Day Number)
      N      Mean motion, n (degrees/day)
      MA     Mean anomaly, M (degrees)
      TA     True anomaly, nu (degrees)
      A      Semi-major axis, a (au)
      AD     Apoapsis distance (au)
      PR     Sidereal orbit period (day)

OM : N
IN : i
W  : w
A  : a
EC : e
MA : M
"""
from urllib.request import urlopen
from datetime import datetime
import pandas as pd

BASE_URL = 'https://ssd.jpl.nasa.gov/horizons_batch.cgi?batch=1&'

def vector_url(target, t1, t2, step, center='500@10', ref_plane='ECLIPTIC'):
    params=f"""MAKE_EPHEM=YES
    COMMAND='{target}'
    EPHEM_TYPE=ELEMENTS
    CENTER='{center}'
    START_TIME='{t1}'
    STOP_TIME='{t2}'
    STEP_SIZE='{step}'
    REF_SYSTEM='ICRF'
    REF_PLANE='{ref_plane}'
    OUT_UNITS='AU-D'
    ELM_LABELS='YES'
    TP_TYPE='ABSOLUTE'
    CSV_FORMAT='YES'
    OBJ_DATA='NO'"""
    params = params.replace('\n', '&').replace(' ', '')
    url = BASE_URL + params
    return url



def get_request(target, t1, t2, step):
    error_msg = ''
    url = vector_url(target, t1, t2, step)
    
    with urlopen(url) as r:
        text = r.read().decode('utf-8')
    if ('$$SOE' not in text) or ('$$EOE' not in text):
        error_msg = text[:text.find('$$SOF')]
        return [error_msg, np.array([]), np.array([])]
    mark1 = text.find('$$SOE')
    text = text[mark1+6:]
    mark2 = text.find('$$EOE')
    text = text[:mark2]
    rows = text.split('\n')[:-1]
    n_columns = len(rows[0].split(',')[:-1])
    
    times = [i.split(',')[1].strip()[5:] for i in rows]
    times = [datetime.strptime(i, '%Y-%b-%d %H:%M:%S.%f') for i in times]

    OM = [i.split(',')[5].strip() for i in rows] # N
    IN = [i.split(',')[4].strip() for i in rows] # i
    W  = [i.split(',')[6].strip() for i in rows] # w
    A  = [i.split(',')[11].strip() for i in rows]# a
    EC = [i.split(',')[2].strip() for i in rows] # e
    MA = [i.split(',')[9].strip() for i in rows] # M

    df = pd.DataFrame({'N':OM, 'i':IN, 'w':W, 'a':A, 'e':EC, 'M':MA})
    df.index = times
    
    return df

target=299
t1 = '1990-04-19T00:00:00'
t2 = '1990-04-19T01:00:00'
step = 1

df = get_request(target, t1, t2, step)
df = df.astype(float)
print(df.iloc[0])
