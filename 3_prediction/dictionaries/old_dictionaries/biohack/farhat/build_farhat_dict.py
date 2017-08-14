lines = [line[:-1] for line in open('farhat.txt').readlines()]

ideal_list = ['AMI', 'PAS', 'EMB', 'FLQ', 'CAP', 'KAN', 'OFX', 'MOX', 'INH', 'CIP', 'STR', 'RIF', 'ETH', 'PZA', 'LZD', 'LEVO']

dict = []

cur_drug = ''

for line in lines:
	if line in ideal_list:
		cur_drug = line
	else:
		if 'SNP_CN' in line or 'SNP_N' in line:
			gene = line.split('_')[-2]
			change = line.split('_')[-1]

			if change[-1] == '.':
				change = change[:-1] + '*'

			dict.append(gene + ' ' + change + ' ' + cur_drug)

file = open('farhat_dict.txt', 'w')

for el in dict:
	file.write(el + '\n')

file.close()