import pandas as pd

for i in ['earth', 'mercury', 'mars', 'jupiter', 'saturn', 'uranus', 'neptune']:
    df = pd.read_csv(i+'_all.csv')

    df = df[['TBD', 'JDTDB', 'EC', 'IN', 'OM', 'W', 'A','MA']]
    df.columns = ['t', 'jd', 'e', 'i', 'N', 'w', 'a','M']

    df = df[['t','jd','N','i','w','a','e','M']]
    df.set_index('t').to_csv(i+'.csv')

