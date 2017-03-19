import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np

genomelist_all = [line[:-1] for line in open('/export/data/kchukreev/1_input/all_ids.txt').readlines()]
shifts = []

for name in genomelist_all:
    lines = []
    try:
        lines = [line[:-1] for line in open('/export/data/kchukreev/5_mutation_statistics/frameshifts/'+name+'.txt')]
    except:
        lines = []

    if len(lines) != 0:
        shifts.append([name, len(lines)/2])

number_shifts = []

for el in shifts:
    number_shifts.append(el[1])

number_shifts = np.asarray(number_shifts)

fig = plt.figure()
plt.hist(number_shifts, bins=300)
plt.xlim(0, 50)
fig.savefig('shifts.png', format='png')

