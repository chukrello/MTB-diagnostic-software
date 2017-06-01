tree = open('RAxML_myc_lables.nwk').read()

glossary = [line[:-1].split(' ') for line in open('../1_input/other/farhat glossary.txt').readlines()]

for el in glossary:
	tree = tree.replace(el[0], el[1])

file = open('new_tree.nwk', 'w')
file.write(tree)
file.close()