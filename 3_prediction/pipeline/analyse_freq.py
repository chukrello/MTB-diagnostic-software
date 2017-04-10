#Script analyses mutations and offer new dictionary

import sys

mutations = [line[:-1].split('\t') for line in open('frequency.tsv').readlines()]

bad_mutations = []

for mut in mutations:
    R = int(mut[1])
    S = int(mut[2])

    if R == 0:
        bad_mutations.append(mut)

mutations = [mut for mut in mutations if mut not in bad_mutations]

def keysort(mut):
    return float(mut[1])/(int(mut[1]) + int(mut[2]))

mutations.sort(key=keysort)

file = open('new_dictionary.txt', 'w')

for mut in mutations:
    file.write(mut[0] + '\t' + mut[1] + '\t' + mut[2] + '\n')

file.close()