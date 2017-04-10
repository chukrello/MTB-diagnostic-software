names = ['Pankhurst']

for name in names:
	
	lines = [line[:-1].split(' ') for line in open('../' + name + '.txt').readlines()]
	new_dictionary = []

	for line in lines:
		coord = int(line[1][1:-1])

		if coord < 4000:
			new_dictionary.append('Gene\t' + line[0] + '\t' + str(coord) + '\t' + line[1][0] + '\t' + line[1][-1] + '\t' + line[2] + '\n')
		else:
			new_dictionary.append('NotGene\t' + line[0] + '\t' + str(coord) + '\t' + line[1][0] + '\t' + line[1][-1] + '\t' + line[2] + '\n')

	file = open(name + '.txt', 'w')

	for el in new_dictionary:
		file.write(el)

	file.close()