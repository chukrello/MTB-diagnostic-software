lines = [line[:-1].split(',') for line in open('out.txt').readlines()]

file = open('results_old_format.txt' ,'w')

for line in lines:
	for name in line:
		file.write(name+'\n')

file.close()