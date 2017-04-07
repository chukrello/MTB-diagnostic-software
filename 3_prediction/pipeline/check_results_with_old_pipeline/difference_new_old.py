lines1 = [line[:-1] for line in open('old_results_walker_walker.txt').readlines()]
lines2 = [line[:-1] for line in open('results_old_format.txt').readlines()]

def make_results(lines):
	res = []

	for line in lines:
		if line[:3] == 'SAM':
			res.append([line])
		else:
			res[-1].append(line)

	return res


results1 = make_results(lines1)
results2 = make_results(lines2)


for i in range(len(results1)):
	if len(results1[i]) != len(results2[i]):
		print(results1[i])
		print(results2[i])


def calculate_coord(gene, pos):
	if genes[gene][2] == '+':
		print(genes[gene][0] + 3*pos)
	else:
		print(genes[gene][1] - 3*pos)