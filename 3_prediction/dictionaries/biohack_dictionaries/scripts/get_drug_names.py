# Getting names of the drugs
# It could be different in different dictionaries

def get_drug_names(filename):
	return list(set([line.split(' ')[2][:-1] for line in open(filename).readlines()]))

filenames = ['Coll_right.txt', 'Farhat_right.txt', 'Pankhurst_right.txt', 'Tbdb_right.txt', 'Walker_right.txt']

for filename in filenames:
	print(filename + '    ' + str(get_drug_names(filename)))