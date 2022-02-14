from urllib.request import urlopen
from datetime import datetime
import pandas as pd

BASE_URL = 'https://ssd.jpl.nasa.gov/horizons_batch.cgi?batch=1&'

def download_jpl_elements(target, t1, t2, step, center='500@10', ref_plane='ECLIPTIC'):
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
    print(url)
    error_msg = ''
    with urlopen(url) as r:
        text = r.read().decode('utf-8')
    if ('$$SOE' not in text) or ('$$EOE' not in text):
        error_msg = text[:text.find('$$SOF')]
        return error_msg
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

    df = pd.DataFrame({'time': times, 'N':OM, 'i':IN, 'w':W, 'a':A, 'e':EC, 'M':MA})
    df['time'] = pd.to_datetime(df['time'])
    df[['N','i','w','a','e','M']] = df[['N','i','w','a','e','M']].astype(float)

    df['_y'] = df['time'].dt.year
    df['_m'] = df['time'].dt.month
    df['_d'] = df['time'].dt.day
    df['_H'] = df['time'].dt.hour
    df['_M'] = df['time'].dt.minute
    df['_S'] = df['time'].dt.second
    df['_MS'] = df['time'].dt.microsecond

    df['day'] = 367*df['_y'] - 7 * ( df['_y'] + (df['_m']+9)//12 ) // 4 + 275*df['_m']//9 + df['_d'] - 730530
    df['UT'] = df['_H'] + df['_M']/60 + df['_S']/3600 + df['_MS']/3600000000
    df['day'] = df['day'] + df['UT']/24

    df = df[['time', 'day', 'N', 'i', 'w', 'a', 'e', 'M']]
    df = df.set_index('time')
    
    return df, text, rows, n_columns, times

t1 = '2022-02-11T10:00:00'
t2 = '2022-02-12T10:00:00'

df, text, raw_rows, n_columns, times = download_jpl_elements(399, t1, t2, step=5)

cols = ['JDTDB', 'Calendar Date (TDB)', 'EC', 'QR', 'IN', 'OM', 'W', 'Tp', 'N', 'MA', 'TA', 'A', 'AD', 'PR']

raw_rows = [i.split(',')[:-1] for i in raw_rows]
rows = []
for r in raw_rows:
    row = [i.strip() for i in r]
    rows.append(row)

df = pd.DataFrame(data=rows, columns=cols)
