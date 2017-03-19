PATH = '/export/data/kchukreev/5_mutation_statistics/'

lines = [line[:-1] for line in open(PATH + 'results_nonwalker.txt').readlines()]

results = []

for line in lines:
    if line[0] == 'S':
        results.append([line])
    else:
        results[-1].append(line)

lines = [line[:-1] for line in open(PATH + 'results_vphaser_nonwalker.txt').readlines()]

results_vphaser = []

for line in lines:
    if line[0] == 'S':
        results_vphaser.append([line])
    else:
        results_vphaser[-1].append(line)

common_results = results

for el in results_vphaser:
	for i in range(len(common_results)):
		if common_results[i][0] == el[0]:
			common_results[i] += el[1:]

file = open('results_common_nonwalker.txt', 'w')

for sample in common_results:
	for el in sample:
		file.write(el+'\n')

file.close()
