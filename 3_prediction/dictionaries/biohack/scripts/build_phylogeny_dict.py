# dict builds phylogeny file by coll table

lines = [line[:-1].split(' ') for line in open('../phylogeny.txt').readlines()]

phyl = []

for line in lines:
	gene = line[0]
	codon = line[4]
	from_l = line[5][0]
	if len(line[5]) == 3:
		to_l = line[5][2]
	else:
		to_l = '-'
	phyl.append(gene + ' ' + from_l + codon + to_l)

for el in phyl:
	print(el)
