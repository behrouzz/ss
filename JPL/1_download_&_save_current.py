from nasajpl import download_jpl_elements

dc = {'mercury':199, 'venus':299, 'earth':399, 'mars':499,
      'jupiter':599, 'saturn':699, 'uranus':799, 'neptune':899}

t1 = '2020-01-01T00:00:00'
t2 = '2030-01-01T00:00:00'
step = 10*365*5

for k in dc.keys():
    df = download_jpl_elements(target=dc[k], t1=t1, t2=t2, step=step)
    _ = input(f'Press ENTER to save <{k}>...')
    df.to_csv('current/'+k+'.csv')

