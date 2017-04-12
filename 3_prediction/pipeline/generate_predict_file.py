# Script for generating table of prediction
# 1 - resistant, -1 - senstive, 0 - idk
# generate_predict_file.py *FILENAME_OF_FULL_FILE* *output*

# import sys

# FILENAME = sys.argv[1]
# DICTIONARY = sys.argv[2]
# SUBSET = sys.argv[3]
# MIN_DEPTH = sys.argv[4]
# GENE_PERCENTAGE = sys.argv[5]
# OUTPUT = sys.argv[6]

FILENAME = 'MEGAFULL_40.txt'
DICTIONARY = '../dictionaries/Walker_dictionary.txt'
SUBSET = '../../1_input/subsets/Walker_subset.txt'
MIN_DEPTH = 5
GENE_PERCENTAGE = 70
OUTPUT = 'test.txt'

drug_list = ['INH', 'RIF', 'EMB', 'PZA', 'STR', 'CIP', 'MOX', 'OFX', 'AMI', 'CAP', 'KAN', 'PRO', 'ETH']
available_mutations = [line[:-1].replace('\t', ' ') for line in open(DICTIONARY).readlines()]
subset = [line[:-1] for line in open(SUBSET).readlines()]
data = [line[:-1] for line in open(FILENAME).readlines()]


def get_gene_intervals(dictionary, genes):
    list_genes = list(set([mut.split(' ')[1] for mut in dictionary]))
    gene_intervals = []

    for gene in list_genes:
        gene_intervals.append([genes[gene][0], genes[gene][1], gene])

    return gene_intervals
# generate coverage_dictionary
# coverage_dict[sample_name][coord]
def get_coverage_dictionary(dictionary, subset):
    coverage_dict = {}

    for id in subset:
        coverage_dict[id] = {}
        coverage_coords = [line[:-1].split(' ') for line  in open('../../../data/coverage_locuses/'+id+'.txt').readlines()]
        for coord in coverage_coords:
            coverage_dict[id][int(coord[0])] = int(coord[1])

    return coverage_dict

# check every mutation whether it is covered or not
# if not, we cannot say S
def check_dictionary_position_for_coverage(dictionary, coverage_dict, genes, id):

    uncovered_mutations = []

    for mutation in dictionary:
        if mutation[0] == 'Gene':
            gene = mutation[1]
            pos = int(mutation[2])

            if genes[gene][2] == '+':
                coords = range(genes[gene][0] + (pos-1)*3,genes[gene][0] + (pos-1)*3 + 3)
            else:
                coords = reange(genes[gene][1] - (pos-1)*3-2,genes[gene][1] - (pos - 1)*3+1)

            for coord in coords:
                if coverage_dict[id][coord] <= MIN_DEPTH:
                    uncovered_mutations.append(mutation)

        else:
            coord = int(mutation[2])

            if coverage_dict[id][coord] <= MIN_DEPTH:
                    uncovered_mutations.append(mutation)


    uncovered_drugs = list(set([mut[-1] for mut in uncovered_mutations]))

    return uncovered_drugs

# check every gene whether it is covered or not
# if not, we cannot say S
def check_genes_for_coverage(gene_intervals, coverage_dict, drug_association, id):
    
    uncovered_genes = []
    uncovered_drugs = []

    for gene_ineterval in gene_intervals:

        cnt_bad_cov = 0

        for i in range(gene_interval[0], gene_interval[1] + 1):
            if coverage_dict[id][i] <= MIN_DEPTH:
                cnt_bad_cov += 1

        percentage = float(cnt_bad_cov)/(gene_interval[1] - gene_interval[0] + 1)

        if percentage <= GENE_PERCENTAGE:
            uncovered_genes.append(gene_interval[2])

    for gene in uncovered_genes:
        uncovered_drugs += drug_association[gene]

    uncovered_drugs = list(set(uncovered_drugs))

    return uncovered_drugs

# takes list of mutations and find which genes is important for drugs
def get_drug_gene_association_from_dictionary(dictionary, drug_list):
    drug_association = {}

    for line in dictionary:
        mut = line.split(' ')
        gene = mut[1]
        drug = mut[-1]

        if gene not in drug_association:
            drug_association[gene] = [drug]
        else:
            if drug not in drug_association[gene]:
                drug_association[gene] += [drug]

    return drug_association

def get_genes_info():
    import re
    lines = open('../../1_input/AL123456_rev.gff').readlines()
    gene_reg = re.compile('\tgene \w*')
    genes = []
    for line in lines:
        m = gene_reg.search(line)
        if m:
            genes.append(line)
    genes = [[line[:-1].split('\t')[3],line[:-1].split('\t')[4], line[:-1].split('\t')[8].split(' ')[1], line[:-1].split('\t')[6]] for line in genes]
    list_genes = ['ahpC', 'eis', 'embA', 'embB', 'embC', 'embR', 'fabG1', 'gid', 'gyrA', 'gyrB', 'inhA', 'iniA', 'iniC',
             'katG', 'manB', 'ndh', 'pncA', 'rmlD', 'rpoB', 'rpsA', 'rpsL', 'rrs', 'tlyA']
    def get_sequence(genes, gene):
        for element in genes:
            if element[2] == gene:
                return [int(element[0]), int(element[1]), element[3]]
    genes_seqs = dict()
    for gene in list_genes:
        genes_seqs[gene] = get_sequence(genes,gene)

    genes_seqs['gidB'] = genes_seqs['gid']
    del(genes_seqs['gid'])

    return genes_seqs



#start end strand seq
genes = get_genes_info()
# coverage_dict[sample_name][coord]
coverage_dict = get_coverage_dictionary(available_mutations, subset)
# coord1 coord2 gene_name
gene_intervals = get_gene_intervals(available_mutations, genes)
drug_association = get_drug_gene_association_from_dictionary(available_mutations, drug_list)


print(genes)

prediction = []

for sample_info in data:
    
    name = sample_info.split(',')[0]
    # if there is no mutations, only name
    if name[:-1] == '\n':
        name = name[:-1]

    if name in subset:

        uncovered_drugs = check_dictionary_position_for_coverage(available_mutations, coverage_dict, genes, name)
        uncovered_drugs += check_genes_for_coverage(gene_intervals, coverage_dict, drug_association, name)

        #using new dictionary
        all_mutations = sample_info.split(',')[1:]
        mutations_filtered = [mutation for mutation in all_mutations if mutation in available_mutations]

        #creating prediction
        prediction.append([name])

        #searching for drug in mutations
        for i in range(len(drug_list)):
            isResistant = False
            for mutation in mutations_filtered:
                if drug_list[i] in mutation:
                    prediction[-1] += [1]
                    isResistant = True
                    break

            if isResistant == False:
                prediction[-1] += [-1]


file = open(OUTPUT, 'w')
for line in prediction:
    for el in line[:-1]:
        file.write(str(el) + '\t')
    file.write(str(el) + '\n')

file.close()