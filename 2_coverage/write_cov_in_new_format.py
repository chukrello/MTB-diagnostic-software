# writes coverage in simple reduced format
# python write_cov_in_new_format.py *subset* *good_depth*


import sys
from joblib import Parallel, delayed
import multiprocessing

SUBSET = sys.argv[1]
GOOD_DEPTH = int(sys.argv[2])
PATH = 'new_format/'

subset = [line[:-1] for line in open(SUBSET).readlines()]

def write_cov_in_new_format(id, PATH):

	print(id)
	
	cov = [[int(line[:-1].split(' ')[0]), int(line[:-1].split(' ')[1])] for line in open('../../data/coverage_locuses/' + id + '.txt').readlines()]
	if cov[0][1] >= GOOD_DEPTH:
		new_cov = [['GOOD', cov[0][0], cov[0][0]]]
	else:
		new_cov = [cov[0]]

	for line in cov[1:]:
		if line[1] >= GOOD_DEPTH:
			if new_cov[-1][0] == 'GOOD' and line[0] == new_cov[-1][2] + 1:
				new_cov[-1][2] = line[0]
			else:
				new_cov.append(['GOOD', line[0], line[0]])
		else:
			new_cov.append(line)

	for i in range(len(new_cov)):
		for j in range(len(new_cov[i])):
			new_cov[i][j] = str(new_cov[i][j])

	file = open(PATH + id + '.txt', 'w')

	for el in new_cov:
		file.write(' '.join(el) + '\n')

	file.close()
			


Parallel(n_jobs=-1)(delayed(write_cov_in_new_format)(id, PATH) for id in subset)

