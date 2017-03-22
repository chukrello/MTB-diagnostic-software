headers = ['Coll', 'Walker', 'Farhat', 'Tbdb']
drugs = ['INH', 'EMB', 'PZA', 'RIF', 'AMI']

for head in headers:
	for drug in drugs:
		lines = open(head+'.txt').readlines()
		file = open(head + '_' + drug + '.txt', 'w')
		for line in lines:
			if drug in line:
				file.write(line)
		file.close()