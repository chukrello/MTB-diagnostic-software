# Script takes predict file, take phenotype (not in params) file and count PPV, NPV, Spec, Sens, R_R, S_R, S_S, R_S
# count_statistics.py *predict_file* *output*

import sys
import numpy as np

PREDICT_FILE = sys.argv[1]
OUTPUT = sys.argv[2]

phenotype = [line[:-1].split('\t') for line in open('../../1_input/phenotype.tsv').readlines()]
predict_results = [line[:-1].split('\t') for line in open(PREDICT_FILE).readlines()]
drug_list = ['INH',
	'RIF',
	'EMB',
	'PZA',
	'STR',
	'CIP',
	'MOX',
	'OFX',
	'AMI',
	'CAP',
	'KAN',
	'PRO',
	'ETH']
file = open(OUTPUT, 'w')


def predict(predict_results, phenotype):
	phenotype_results = []

	# shorten phenotype to the columns we are interested in
	for el in phenotype[1:]:
		phenotype_results.append(el[:14])
	
	# turning letters in numbers
	for i in range(len(phenotype_results)):
		for j in range(1,len(phenotype_results[i])):
			if phenotype_results[i][j] == 'R':
				phenotype_results[i][j] = 1
			elif phenotype_results[i][j] == 'S':
				phenotype_results[i][j] = -1
			else:
				phenotype_results[i][j] = 0

	# filling genotype results
	genotype_results = []
	for prediction_samn in predict_results:
		genotype_results.append([prediction_samn[0]])
		for el in prediction_samn[1:]:
			genotype_results[-1] += [int(el)]

	current_phenotype = []
	for el in genotype_results:
		name = el[0]
		for phen in phenotype_results:
			if name == phen[0]:
				current_phenotype.append(phen)
	phenotype_results = current_phenotype


	# analysing of prediction
	predict_results = [[0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0]]

	for i in range(1,14):
		cnt = 0
		for phen in phenotype_results:
			for gen in genotype_results:
				if phen[0] == gen[0]:
					if phen[i] == 1:
						if gen[i] == 1:
							predict_results[i-1][3] += 1
						elif gen[i] == -1:
							predict_results[i - 1][5] += 1
						else:
							predict_results[i - 1][4] += 1
					elif phen[i] == -1:
						if gen[i] == 1:
							predict_results[i - 1][0] += 1
							cnt += 1
						elif gen[i] == -1:
							predict_results[i - 1][2] += 1
						else:
							predict_results[i - 1][1] += 1
					else:
						if gen[i] == 1:
							predict_results[i - 1][6] += 1
						elif gen[i] == -1:
							predict_results[i - 1][7] += 1
						else:
							predict_results[i - 1][8] += 1
					break
	return predict_results


def write_statistics_values_by_drugs(predict_results, drug_indexes):
	
	#Output
	predict_data = np.array(predict_results).T

	results = []
	for el in predict_data:
	    results.append(list(el))

	def calculate_parameters(results, index):
		return [int(results[0][index]),int(results[2][index]),int(results[3][index]),int(results[5][index])]


	def print_data(index):
		array = calculate_parameters(results, index)
		S_R = array[0]
		S_S = array[1]
		R_R = array[2]
		R_S = array[3]

		try:
			PPV = float(R_R)/(R_R + S_R)
			PPV *= 100
			PPV = int(round(PPV))
		except Exception:
			PPV = 'inf'

		try:
			NPV = float(S_S)/(S_S + R_S)
			NPV *= 100
			NPV = int(round(NPV))
		except Exception:
			NPV = 'inf'

		try:
			sens = float(R_R)/(R_R + R_S)
			sens *= 100
			sens = int(round(sens))
		except Exception:
			sens = 'inf'
		
		try:
			spec = float(S_S)/(S_S + S_R)
			spec *= 100
			spec = int(round(spec))
		except Exception:
			spec = 'inf'

		file.write(drug_list[i] + '\t' + str(PPV) + ' ' + str(R_R) + '/' + str(R_R + S_R) + '\t' + str(NPV) + ' ' + str(S_S) + '/' + str(S_S + R_S) + '\t' + str(sens) + ' ' + str(R_R) + '/' + str(R_R + R_S)  + '\t' + str(spec) + ' ' + str(S_S) + '/' + str(S_S + S_R) + '\n')
	
	for i in drug_indexes:
		print_data(i)

def write_common_statistic_values(predict_results, drug_indexes):

	#Output
	predict_data = np.array(predict_results).T

	results = []
	for el in predict_data:
		results.append(list(el))

	S_R = S_S = R_R = R_S = 0

	for index in drug_indexes:
		S_R += results[0][index]
		S_S += results[2][index]
		R_R += results[3][index]
		R_S += results[5][index]


	try:
		PPV = float(R_R)/(R_R + S_R)
		PPV *= 100
		PPV = int(round(PPV))
	except Exception:
		PPV = 'inf'

	try:
		NPV = float(S_S)/(S_S + R_S)
		NPV *= 100
		NPV = int(round(NPV))
	except Exception:
		NPV = 'inf'

	try:
		sens = float(R_R)/(R_R + R_S)
		sens *= 100
		sens = int(round(sens))
	except Exception:
		sens = 'inf'
	
	try:
		spec = float(S_S)/(S_S + S_R)
		spec *= 100
		spec = int(round(spec))
	except Exception:
		spec = 'inf'

	file.write('Common statistics\t' + str(PPV) + ' ' + str(R_R) + '/' + str(R_R + S_R) + '\t' + str(NPV) + ' ' + str(S_S) + '/' + str(S_S + R_S) + '\t' + str(sens) + ' ' + str(R_R) + '/' + str(R_R + R_S)  + '\t' + str(spec) + ' ' + str(S_S) + '/' + str(S_S + S_R) + '\n\n')

def write_numbers(predict_results):

	#Output
	predict_data = np.array(predict_results).T

	results = []
	for el in predict_data:
		results.append(list(el))

	for i in range(len(results)):
		for j in range(len(results[i])):
			results[i][j] = str(results[i][j])

	file.write('S\tR\t'+' '.join(results[0])+ '\n')
	file.write('S\tN\t'+' '.join(results[1])+ '\n')
	file.write('S\tS\t'+' '.join(results[2])+ '\n')
	file.write('R\tR\t'+' '.join(results[3])+ '\n')
	file.write('R\tN\t'+' '.join(results[4])+ '\n')
	file.write('R\tS\t'+' '.join(results[5])+ '\n')
	file.write('N\tR\t'+' '.join(results[6])+ '\n')
	file.write('N\tS\t'+' '.join(results[7])+ '\n')
	file.write('N\tN\t'+' '.join(results[8])+ '\n')


predict_results = predict(predict_results, phenotype)
write_statistics_values_by_drugs(predict_results, range(11))
write_common_statistic_values(predict_results, range(11))
write_numbers(predict_results)