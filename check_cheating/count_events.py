import numpy as np
import pandas as pd

if __name__ == '__main__':

    cheated = pd.read_csv('events_cheated_writer.txt', sep=' ', header=None)
    uncheated = pd.read_csv('events_non_cheated_writer.txt', sep=' ', header=None)

    bad_runs = pd.concat([~uncheated[2].isin(cheated[2]), uncheated[2]], axis=1)

    count = 0
    for is_bad, run in bad_runs.itertuples(index=False, name=None):

        if is_bad:
            print(f'{run}', end=', ')
            count += 1
            if count % 50 == 0:
                print()
    print(f'\n\n{count = }')
