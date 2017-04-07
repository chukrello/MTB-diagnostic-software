SUBSET_NAMES = ['Farhat.txt', 'Pankhurst.txt']

for name in SUBSET_NAMES:
	file = open(name[:-4] + '_subset.txt', 'w')
	lines = [line.split('\t')[0] for line in open(name).readlines()]
	for line in lines:
		file.write(line + '\n')
	file.close()