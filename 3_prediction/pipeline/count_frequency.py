# counting frequency of mutations from full file
# python count_frequency.py *full_file* *otput*

import sys, os, argparse, subprocess 

parser = argparse.ArgumentParser(description='Counts frequency of mutations for FULL file. FULL file is the file where contatins every drug-resistant mutation for each sample.')

parser.add_argument('-f','--full_file', help='FULL file', default='MEGAFULL_40.txt')
parser.add_argument('-o','--output', help='Output file', default='MEGAFULL_40_freq.txt')

args = parser.parse_args()

FILE = args.full_file
OUTPUT = args.output

lines = [line[:-1].split(',') for line in open(FILE).readlines()]


#make dictionary of mutations
dict_counter = {}
for line in lines:
    for i in range(1, len(line)):
        if line[i] != '' and line[i] not in dict_counter:
            dict_counter[line[i]] = [line[0]]
        elif line[i] != '':
            dict_counter[line[i]].append(line[0])

#write dictionary in array
mutations_freq = []
for mut in dict_counter:
    mutations_freq.append([mut, len(dict_counter[mut])])
    mutations_freq[-1] += dict_counter[mut]

#get phenotype list for statistics
phenotype_list = [line[:-1].split('\t') for line in open('../../1_input/phenotype.tsv').readlines()]

def GetResistanceByNameAndDrug(name, drug, phenotype_list):
    for el in phenotype_list:
        if el[0] == name:
            if drug == 'INH':
                drug_index = 1
                return el[drug_index]
            elif drug == 'RIF':
                drug_index = 2
                return el[drug_index]
            elif drug == 'EMB':
                drug_index = 3
                return el[drug_index]
            elif drug == 'PZA':
                drug_index = 4
                return el[drug_index]
            elif drug == 'STR':
                drug_index = 5
                return el[drug_index]
            elif drug == 'CIP':
                drug_index = 6
                return el[drug_index]
            elif drug == 'MOX':
                drug_index = 7
                return el[drug_index]
            elif drug == 'OFX':
                drug_index = 8
                return el[drug_index]
            elif drug == 'AMI':
                drug_index = 9
                return el[drug_index]
            elif drug == 'CAP':
                drug_index = 10
                return el[drug_index]
            elif drug == 'KAN':
                drug_index = 11
                return el[drug_index]
            elif drug == 'WTF': #Prothionamide
                drug_index = 12
                return el[drug_index]
            elif drug == 'ETH':
                drug_index = 13
                return el[drug_index]

mutations_statistics = []

#calculate statistics
for i in range(len(mutations_freq)):
    mutations_statistics.append([mutations_freq[i][0]])
    drug = mutations_freq[i][0].split(' ')[-1]

    R = 0
    S = 0
    for sample in mutations_freq[i][2:]:
        phen_result = GetResistanceByNameAndDrug(sample, drug, phenotype_list)
        if phen_result == 'S':
            S += 1
        elif phen_result == 'R':
            R += 1
    mutations_statistics[-1].append(str(R))
    mutations_statistics[-1].append(str(S))

#extracting zeros
zeros = [mut for mut in mutations_statistics if mut[1]==mut[2]=='0']
mutations_statistics = [mut for mut in mutations_statistics if mut not in zeros]

mutations_statistics.sort(key = lambda mut: float(int(mut[1]))/(int(mut[1]) + int(mut[2])))
mutations_statistics += zeros

file = open(OUTPUT, 'w')

for mut in mutations_statistics:
    file.write('\t'.join(mut) + '\n')
file.close()
