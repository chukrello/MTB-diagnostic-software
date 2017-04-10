PATH_to_input = '../1_input/'

ids = [line[:-1] for line in open(PATH_to_input + 'all_ids.txt').readlines()]

for id in ids:
	try:
		data = [line[:-1].split('\t') for line in open('snps_samples_common_format/' + id + '.txt').readlines()]
	except Exception:
		continue
	
	new_data = []
	for el in data:
		if el[0] == 'Gene':
			new_data.append(el)
			new_data[-1][2] = str(int(float(new_data[-1][2])))
		else:
			new_data.append(el)
	
	file = open('tmp/'+id+'.txt', 'w')
	for el in new_data:
		file.write('\t'.join(el)+'\n')
	file.close()