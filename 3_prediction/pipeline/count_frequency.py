#open file with counting dict mutations
FILE = 'full_all.txt'
lines = [line[:-1].split('\t') for line in open(FILE).readlines()]


#make dictionary of mutations
dict_counter = {}
for line in lines:
    for i in range(1, len(line)):
        if 'SAM' in line[i]:
            print(line)
        if line[i] != '' and line[i] not in dict_counter:
            dict_counter[line[i]] = [line[0]]
        elif line[i] != '':
            dict_counter[line[i]].append(line[0])

print(dict_counter)

#write dictionary in array
mutations_freq = []
for mut in dict_counter:
    mutations_freq.append([mut, len(dict_counter[mut])])
    mutations_freq[-1] += dict_counter[mut]

#get phenotype list for statistics
phenotype_list = [line[:-1].split('\t') for line in open('phenotype_db.tsv').readlines()]

def GetResistanceByNameAndDrug(name, drug, phenotype_list):
    for el in phenotype_list:
        if el[0] == name:
            if drug == 'INH':
                drug_index = 1
                return el[drug_index]
            elif drug == 'RMP':
                drug_index = 2
                return el[drug_index]
            elif drug == 'EMB':
                drug_index = 3
                return el[drug_index]
            elif drug == 'PZA':
                drug_index = 4
                return el[drug_index]
            elif drug == 'SM':
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
            elif drug == 'WTF':
                drug_index = 12
                return el[drug_index]
            elif drug == 'ETH':
                drug_index = 13
                return el[drug_index]

#calculate statistics
for i in range(len(mutations_freq)):
    drug = mutations_freq[i][0].split(' ')[-1]
    R = 0
    S = 0
    for sample in mutations_freq[i][2:]:
        phen_result = GetResistanceByNameAndDrug(sample, drug, phenotype_list)
        if phen_result == 'S':
            S += 1
        elif phen_result == 'R':
            R += 1
    mutations_freq[i].append(R)
    mutations_freq[i].append(S)

file = open('mutation_statistic.tsv', 'w')

for mut in mutations_freq:
    for el in mut[:-1]:
        file.write(str(el) + '\t')
    file.write(str(mut[-1])+'\n')

file.close()


#'rpoB H445Y RIF'
#SAMEA1018864