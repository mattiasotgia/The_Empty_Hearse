import numpy as np
import pandas as pd

if __name__ == '__main__':

    cheated = pd.read_csv('events_cheated_writer.txt', sep=' ', header=None)
    uncheated = pd.read_csv('events_non_cheated_writer.txt', sep=' ', header=None)

    bad_runs_1 = pd.concat([~uncheated[2].isin(cheated[2]), uncheated[2]], axis=1) # these runs are in the uncheated, but not in the cheated
    bad_runs_2 = pd.concat([~cheated[2].isin(uncheated[2]), cheated[2]], axis=1) # these runs are in the cheated, but not in the uncheated

    count = 0

    print('These are bad runs, present in the un-cheated events but not in the cheated events')
    for is_bad, run in bad_runs_1.itertuples(index=False, name=None):
        if is_bad:
            print(f'{run}', end=', ')
            count += 1
            if count % 50 == 0:
                print()
    print(f'\n\n{count = }')

    count = 0

    print('These are bad runs, present in the cheated events but not in the un-cheated events')
    for is_bad, run in bad_runs_2.itertuples(index=False, name=None):
        if is_bad:
            print(f'{run}', end=', ')
            count += 1
            if count % 50 == 0:
                print()
    print(f'\n\n{count = }')
