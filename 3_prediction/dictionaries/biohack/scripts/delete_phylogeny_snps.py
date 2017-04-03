phylogeny = [line[:-1] for line in open('../phylogeny_dict.txt').readlines()]

def delete_phyl(name):

	dict = [line[:-5] for line in open('../' + name + '.txt').readlines()]
	new_dict_indexes = []
	for i in range(len(dict)):
		if dict[i] not in phylogeny:
			new_dict_indexes.append(i)

	new_dict = []
	lines = open('../' + name + '.txt').readlines()
	for i in new_dict_indexes:
		new_dict.append(lines[i])

	file = open(name + '.txt', 'w')

	for el in new_dict:
		file.write(el)

	file.close()

delete_phyl('Coll')
