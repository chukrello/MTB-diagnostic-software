# Script for generating table of prediction
# 1 - resistant, -1 - senstive, 0 - idk
# generate_predict_file.py *FILENAME_OF_FULL_FILE* *output*

import sys

FILENAME = sys.argv[1]
DICTIONARY = sys.argv[2]
SUBSET = sys.argv[3]
OUTPUT = sys.argv[4]

available_mutations = [line[:-1] for line in open(DICTIONARY).readlines()]

#tabs in dictionary into spacec (full file)
for i in range(len(available_mutations)):
    available_mutations[i] = available_mutations[i].replace('\t', ' ')

subset = [line[:-1] for line in open(SUBSET).readlines()]
data = [line[:-1] for line in open(FILENAME).readlines()]

drug_list = ['INH', 'RIF', 'EMB', 'PZA', 'STR', 'CIP', 'MOX', 'OFX', 'AMI', 'CAP', 'KAN', 'PRO', 'ETH']

prediction = []

for sample_info in data:
    
    name = sample_info.split(',')[0]
    # if there is no mutations, only name
    if name[:-1] == '\n':
        name = name[:-1]

    if name in subset:

        #using new dictionary
        all_mutations = sample_info.split(',')[1:]
        mutations_filtered = [mutation for mutation in all_mutations if mutation in available_mutations]
        

        #creating prediction
        prediction.append([name])

        #searching for drug in mutations
        for i in range(len(drug_list)):
            isResistant = False
            for mutation in mutations_filtered:
                if drug_list[i] in mutation:
                    prediction[-1] += [1]
                    isResistant = True
                    break

            if isResistant == False:
                prediction[-1] += [-1]


file = open(OUTPUT, 'w')
for line in prediction:
    for el in line[:-1]:
        file.write(str(el) + '\t')
    file.write(str(el) + '\n')

file.close()