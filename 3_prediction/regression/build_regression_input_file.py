ids = [line[:-1] for line in open('list_names.txt').readlines()]

target = 'Gene\tkatG\t315.3333333333333\tS\tT'

mutations = []
counter = {}

for id in ids:
    data = [line[:-1] for line in open('/home/chukreev/Desktop/regression/'+id).readlines()]
    for el in data:
        if el not in mutations:
            mutations.append(el)
            counter[el] = 1
        else:
            counter[el] += 1

singletons = []

for el in counter:
    if counter[el] == 1:
        singletons.append(el)

synonymos = []

for mut in mutations:
    info = mut.split('\t')
    if info[-1] == info[-2]:
        synonymos.append(mut)

mutations = [element for element in mutations if element not in singletons]
mutations = [element for element in mutations if element not in synonymos]

print(len(mutations))



regression_table = []

for id in ids:
    new_regr = [id[:-4]]
    data = [line[:-1] for line in open('/home/chukreev/Desktop/regression/'+id).readlines()]
    for i in range(len(mutations)):
        if mutations[i] in data:
            new_regr.append(1)
        else:
            new_regr.append(0)
    regression_table.append(new_regr)

phenotype = [line[:-1].split('\t') for line in open('phenotype_db.tsv').readlines()][1:]

for i in range(len(regression_table)):
    for phen in phenotype:
        if regression_table[i][0] == phen[0]:
            regression_table[i].append(phen[4])

noneffective = []

for i in range(1, len(regression_table[0])-1):
    snp_index = i

    isEffective = False

    for j in range(len(regression_table)):
        if regression_table[j][snp_index] == 1 and regression_table[j][-1] == 'R':
            isEffective = True
            break

    if isEffective == False:
        noneffective.append(snp_index)

print(len(noneffective))

new_mutations = []

for i in range(len(mutations)):
    if i not in noneffective:
        new_mutations.append(mutations[i])

mutations = new_mutations

print(len(mutations))

regression_table = []

for id in ids:
    new_regr = [id[:-4]]
    data = [line[:-1] for line in open('/home/chukreev/Desktop/regression/'+id).readlines()]
    for i in range(len(mutations)):
        if mutations[i] in data:
            new_regr.append(1)
        else:
            new_regr.append(0)
    regression_table.append(new_regr)

phenotype = [line[:-1].split('\t') for line in open('phenotype_db.tsv').readlines()][1:]

for i in range(len(regression_table)):
    for phen in phenotype:
        if regression_table[i][0] == phen[0]:
            regression_table[i].append(phen[4])

file = open('regr_input_PZA.txt', 'w')

for line in regression_table:
    if line[-1] in ['R', 'S']:
        for el in line[:-1]:
            file.write(str(el)+'\t')
        file.write(str(line[-1]) + '\n')

file.close()