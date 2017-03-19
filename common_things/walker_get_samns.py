#
# Сonverts wakler ids (ers and err) into samn with the help of Kuleshov files
# Works good with cluster and phenotype results (extracting, comparing)
# explaining in 18.01

PATH = '/home/chukreev/PycharmProjects/tuberculosis/resources/phenotype/'
kuleshov_info = [[line[:-1].split('\t')[0],line[:-1].split('\t')[3],line[:-1].split('\t')[4]] for line in open(PATH+'kuleshov_info.txt').readlines()]

def get_samns(kuleshov_info, ids):
	new = []
	for id in ids:
		for el in kuleshov_info:
			if id == el[0] or id == el[1] or el[2] in id:
				new.append(el[0])
	return new

def get_genotype_results():
	lines = [line[:-1].split('\t') for line in open(PATH + 'phenotype_database.tsv').readlines()]
	data = []
	for el in lines:
		data.append([el[0], el[1], el[2]] + el[22:-1])

	return data

def get_results_of_drug(data, drug_index):
	new = []
	for el in data:
		new.append([el[0],el[1],el[2],el[drug_index]])
	return new

def get_resistant_result_ids(drug_data):
	new = []
	for el in drug_data:
		if el[-1] == 'R':
			if el[0] != '-':
				new.append(el[0])
			elif el[1] != '-':
				new.append(el[1])
			else:
				new.append(el[2])
	return new

def get_cluster_data():
	lines = [line[:-1] for line in open('/home/chukreev/PycharmProjects/tuberculosis/results.txt').readlines()]
	cluster_data = []

	for i in range(len(lines)):
	    if 'SAM' in lines[i]:
	        cluster_data.append([lines[i]])
	    elif 'Resistant' in lines[i]:
	    	cluster_data[-1].append(lines[i])

	return cluster_data

def find_resistant_in_cluster_data(cluster_data, drug):
	resistant_drug_ids = []
	for sample in cluster_data:
		id = sample[0]
		for el in sample[1:]:
			if 'Resistant' in el and drug in el:
				resistant_drug_ids.append(id)

	return resistant_drug_ids


data = get_genotype_results()
results = get_results_of_drug(data,4)
res_RIF_ids = get_resistant_result_ids(results)
walker_res_RIF_samns = list(set(get_samns(kuleshov_info, res_RIF_ids)))


cluster_data = get_cluster_data()
cluster_res_RIF_ids = find_resistant_in_cluster_data(cluster_data, 'RIF')

print(len(walker_res_RIF_samns))
print(len(cluster_res_RIF_ids))

for el in walker_res_RIF_samns:
	if el not in cluster_res_RIF_ids:
		print(el)