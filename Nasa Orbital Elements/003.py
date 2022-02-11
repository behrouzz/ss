import pandas as pd

def day(y, m, D, UT):
    d = 367*y - 7 * ( y + (m+9)//12 ) // 4 + 275*m//9 + D - 730530
    d = d + UT/24.0
    return d

t1 = '2000-01-01T00:00:00'
t2 = '2050-01-01T00:00:00'

name = 'neptune'

df = pd.read_csv('data/'+name+'.csv')
df['time'] = pd.to_datetime(df['time'])

df['_y'] = df['time'].dt.year
df['_m'] = df['time'].dt.month
df['_d'] = df['time'].dt.day
df['_H'] = df['time'].dt.hour
df['_M'] = df['time'].dt.minute
df['_S'] = df['time'].dt.second
df['_MS'] = df['time'].dt.microsecond

df['d'] = 367*df['_y'] - 7 * ( df['_y'] + (df['_m']+9)//12 ) // 4 + 275*df['_m']//9 + df['_d'] - 730530

df['UT'] = df['_H'] + df['_M']/60 + df['_S']/3600 + df['_MS']/3600000000
df['d'] = df['d'] + df['UT']/24


df = df[['time', 'd', 'N', 'i', 'w', 'a', 'e', 'M']]

df = df.set_index('time')
df.to_csv('data/'+name+'.csv')
