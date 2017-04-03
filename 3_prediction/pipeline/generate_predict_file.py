data = [line[:-1] for line in open('new_full_4.txt').readlines()]
file = open('dict_iteration_5.txt', 'w')


#FILE_AVAIL_MUT = ''

drug_list = ['INH', 'RMP', 'EMB', 'PZA', 'STR', 'CIP', 'MOX', 'OFX', 'AMI', 'CAP', 'KAN', 'PRO', 'ETH']

#available_mutations = [line[:-1] for line in open(FILE_AVAIL_MUT).readlines()]
prediction = []

for sample_info in data:
    name = sample_info.split('\t')[0]
    all_mutations = sample_info.split('\t')[1:]
    #filtered_mutations = [mutation for mutation in all_mutations if mutation in available_mutations]

    filtered_mutations = all_mutations

    prediction.append([name])

    for i in range(len(drug_list)):
        isResistant = False
        for mutation in filtered_mutations:
            if drug_list[i] in mutation:
                prediction[-1] += [1]
                isResistant = True
                break

        if isResistant == False:
            prediction[-1] += [-1]



for line in prediction:
    for el in line[:-1]:
        file.write(str(el) + '\t')
    file.write(str(el) + '\n')

file.close()



