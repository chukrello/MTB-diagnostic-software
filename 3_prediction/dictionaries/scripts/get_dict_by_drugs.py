# script divides dictionary by drugs
# now it has group AMG (aminoglycosides) = AMI + CAP + KAN

headers = ['Coll', 'Walker', 'Farhat', 'Tbdb']
drugs = ['INH', 'EMB', 'PZA', 'RIF']
AMI_drugs = ['AMI', 'KAN', 'CAP']

for head in headers:
	for drug in drugs:
		lines = open(head+'.txt').readlines()
		file = open(head + '_' + drug + '.txt', 'w')
		for line in lines:
			if drug in line:
				file.write(line)
		file.close()

	#for AMG
	drug = 'AMG'
	lines = open(head+'.txt').readlines()
	file = open(head + '_' + drug + '.txt', 'w')
	AMG_lines = []
	for line in lines:
		if 'AMI' in line or 'KAN' in line or 'CAP' in line:
			line = line[:-4] + 'AMG\n'
			AMG_lines.append(line)
	AMG_lines = set(list(AMG_lines))
	for el in AMG_lines:
		file.write(el)
	file.close()