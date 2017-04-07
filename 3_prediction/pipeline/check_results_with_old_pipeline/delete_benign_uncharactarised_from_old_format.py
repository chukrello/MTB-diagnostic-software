lines = open('results_walker_walker_dict.txt').readlines()

file = open('old_results_walker_walker.txt', 'w')

for line in lines:
	if 'Benign' not in line and 'Uncharacterised' not in line:
		file.write(line)

file.close()