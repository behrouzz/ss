from nasajpl import get_request

target=899
t1 = '2000-01-01T00:00:00'
t2 = '2050-01-01T00:00:00'
step = 18250

df = get_request(target=target, t1=t1, t2=t2, step=step)
df = df.astype(float)
df.index.name = 'time'

df.to_csv('data/neptune.csv')
