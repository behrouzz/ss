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
    raw_rows = text.split('\n')[:-1]
    n_columns = len(raw_rows[0].split(',')[:-1])

    cols = ['JDTDB', 'Calendar Date (TDB)', 'EC', 'QR', 'IN', 'OM', 'W', 'Tp', 'N', 'MA', 'TA', 'A', 'AD', 'PR']

    raw_rows = [i.split(',')[:-1] for i in raw_rows]
    rows = []
    for r in raw_rows:
        row = [i.strip() for i in r]
        rows.append(row)

    df = pd.DataFrame(data=rows, columns=cols)
    
    return df

def modify_jpl_df(df):
    df['Calendar Date (TDB)'] = \
                 df['Calendar Date (TDB)'].str[5:]\
                 .apply(lambda i: datetime.strptime(i, '%Y-%b-%d %H:%M:%S.%f'))
    df = df.set_index('Calendar Date (TDB)')
    df.index.name = 'TBD'
    df = df.astype(float)
    return df


t1 = '2020-01-01T00:00:00'
t2 = '2030-01-01T00:00:00'
df = download_jpl_elements(399, t1, t2, step='1 DAYS')
df = modify_jpl_df(df)
df.to_csv('data/earth.csv')
