def delete_duplicate(name):

	lines = list(set(open('../'+name+'.txt').readlines()))
	lines.sort()
	file = open(name+'.txt', 'w')

	for el in lines:
		file.write(el)


delete_duplicate('Walker')
delete_duplicate('Coll')
delete_duplicate('Tbdb')
delete_duplicate('Farhat')