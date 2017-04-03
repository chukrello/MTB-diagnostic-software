# Fixes drug names in dictionaries
# Convert to ['AMI', 'PAS', 'EMB', 'FLQ', 'CAP', 'KAN', 'OFX', 'MOX', 'INH', 'CIP', 'STR', 'RIF', 'ETH', 'PZA', 'LZD']

filenames = ['Coll.txt', 'Farhat.txt', 'Pankhurst.txt', 'Tbdb.txt', 'Walker.txt']
ideal_list = ['AMI', 'PAS', 'EMB', 'FLQ', 'CAP', 'KAN', 'OFX', 'MOX', 'INH', 'CIP', 'STR', 'RIF', 'ETH', 'PZA', 'LZD']
bad_list = ['AMK', 'AK', 'SM', 'RMP']

def fix_names(filename, ideal_list, bad_list):
	lines = open(filename).readlines()

	for i in range(len(lines)):
		tmp = lines[i].split(' ')
		drug_name = lines[i].split(' ')[-1][:-1]
		if drug_name not in ideal_list:
			if drug_name == 'AMK':
				lines[i] = tmp[0] + ' ' + tmp[1] + ' ' + 'AMI\n'
			elif drug_name == 'AK':
				lines[i] = tmp[0] + ' ' + tmp[1] + ' ' + 'AMI\n'
			elif drug_name == 'SM':
				lines[i] = tmp[0] + ' ' + tmp[1] + ' ' + 'STR\n'
			elif drug_name ==  'RMP':
				lines[i] = tmp[0] + ' ' + tmp[1] + ' ' + 'RIF\n'

	file = open(filename[:-4]+'_right.txt', 'w')

	for line in lines:
		file.write(line)

	file.close()


for filename in filenames:
	fix_names(filename, ideal_list, bad_list)