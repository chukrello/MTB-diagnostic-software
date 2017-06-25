# Script for generating table of prediction
# 1 - resistant, -1 - senstive, 0 - idk
# generate_predict_file.py *FILENAME_OF_FULL_FILE* *output*

import sys
from joblib import Parallel, delayed
import multiprocessing

FILENAME = sys.argv[1]
DICTIONARY = sys.argv[2]
SUBSET = sys.argv[3]
MIN_DEPTH = int(sys.argv[4])
GENE_PERCENTAGE = int(sys.argv[5])
OUTPUT = sys.argv[6]


# FOR DEBUGGING
# FILENAME = 'MEGAFULL_40.txt'
# DICTIONARY = '../dictionaries/Walker_dictionary.txt'
# SUBSET = '../../1_input/subsets/Walker_subset.txt'
# MIN_DEPTH = 5
# GENE_PERCENTAGE = 70
# OUTPUT = 'test2.txt'


TEST_NAMES = []



drug_list = ['INH', 'RIF', 'EMB', 'PZA', 'STR', 'CIP', 'MOX', 'OFX', 'AMI', 'CAP', 'KAN', 'PRO', 'ETH']
available_mutations = [line[:-1].replace('\t', ' ') for line in open(DICTIONARY).readlines()]

#SUBSET CORRECTION
subset = [line[:-1] for line in open(SUBSET).readlines()]
data = [line[:-1] for line in open(FILENAME).readlines()]


def get_gene_intervals(dictionary, genes):
    list_genes = list(set([mut.split(' ')[1] for mut in dictionary]))
    gene_intervals = []

    for gene in list_genes:
        gene_intervals.append([genes[gene][0], genes[gene][1], gene])

    return gene_intervals


# function for parallel
def calculate_dict(id):
    coverage_dict = {}
    coverage_intervals = []

    coverage_coords = [line[:-1].split(' ') for line  in open('../../../data/coverage_new_format/'+id+'.txt').readlines()]

    for line in coverage_coords:

        if line[0] == 'GOOD':
            coverage_intervals.append([int(line[1]), int(line[2])])
        else:
            coverage_dict[int(line[0])] = int(line[1])

    return [coverage_intervals, coverage_dict]


# generate coverage_dictionary
# coverage_dict[sample_name][coord]
def get_coverage_dictionary(subset):
    
    coverage_dict = {}
    coverage_intervals = {}

    results = Parallel(n_jobs=-1)(delayed(calculate_dict)(id) for id in subset)

    for i in range(len(subset)):
       coverage_dict[subset[i]] = results[i][1]
       coverage_intervals[subset[i]] = results[i][0]
        
    return [coverage_intervals, coverage_dict]

# check every mutation whether it is covered or not
# if not, we cannot say S
def check_dictionary_position_for_coverage(dictionary, coverage_dict, coverage_intervals, genes, id):

    uncovered_mutations = []

    for mutation in dictionary:

        mutation = mutation.split(' ')

        if mutation[0] == 'Gene':
            gene = mutation[1]
            pos = int(mutation[2])

            if genes[gene][2] == '+':
                coords = range(genes[gene][0] + (pos-1)*3,genes[gene][0] + (pos-1)*3 + 3)
            else:
                if gene == 'rrs':
                    coords = [genes[gene][1] - (pos-1)]
                else:
                    coords = range(genes[gene][1] - (pos-1)*3-2,genes[gene][1] - (pos - 1)*3+1)

            for coord in coords:

                isGood = False

                for interval in coverage_intervals[id]:
                    if interval[0] <= coord <= interval[1]:
                        isGood = True
                        break 

                if isGood == True:
                    continue

                if coord not in coverage_dict[id]:
                    uncovered_mutations.append(mutation)
                elif coverage_dict[id][coord] <= MIN_DEPTH:
                    uncovered_mutations.append(mutation)

        else:
            coord = int(mutation[2])
            isGood = False

            for interval in coverage_intervals[id]:
                    if interval[0] <= coord <= interval[1]:
                        isGood = True
                        break 

            if isGood == True:
                continue

            if coord not in coverage_dict[id]:
                uncovered_mutations.append(mutation)
            elif coverage_dict[id][coord] <= MIN_DEPTH:
                uncovered_mutations.append(mutation)

    # if id in TEST_NAMES:
    #     print(id)
    #     print(uncovered_mutations)

    uncovered_drugs = list(set([mut[-1] for mut in uncovered_mutations]))

    return uncovered_drugs

# check every gene whether it is covered or not
# if not, we cannot say S
def check_genes_for_coverage(gene_intervals, coverage_dict, coverage_intervals, drug_association, id):
    
    uncovered_genes = []
    uncovered_drugs = []

    # if id in TEST_NAMES:
    #     print(id)

    for gene_interval in gene_intervals:

        cnt_bad_cov = 0

        for i in range(gene_interval[0], gene_interval[1] + 1):

            isGood = False

            for interval in coverage_intervals[id]:
                if interval[0] <= i <= interval[1]:
                    isGood = True
                    break

            if isGood == True:
                continue

            if i not in coverage_dict[id]:
                cnt_bad_cov += 1
            elif coverage_dict[id][i] <= MIN_DEPTH:
                cnt_bad_cov += 1

        bad_percentage = float(cnt_bad_cov)/(gene_interval[1] - gene_interval[0] + 1)


        # if id in TEST_NAMES:
        #     print(gene_interval[2])
        #     print(bad_percentage)

        if bad_percentage >= (100 - GENE_PERCENTAGE)/100.0:
            uncovered_genes.append(gene_interval[2])

    for gene in uncovered_genes:
        uncovered_drugs += drug_association[gene]

    # if id in TEST_NAMES:
    #     print(uncovered_genes)

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
             'katG', 'manB', 'ndh', 'pncA', 'rmlD', 'rpoB', 'rpsA', 'rpsL', 'rrs', 'tlyA', 'ethR', 'fpbC', 'iniB', 'kasA', 'ethA', 'fabD', 'efpA', 'thyA', 'panD', 'accD6', 'fbpC', 'nat', 'folC', 'rrl', 'rpoC', 'ribD', 'rplC']
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


print('genes')
#start end strand seq
genes = get_genes_info()
print('coverage dict')
# coverage_dict[sample_name][coord]
coverage = get_coverage_dictionary(subset)
coverage_intervals = coverage[0]
coverage_dict = coverage[1]
print('gene intervals')
# coord1 coord2 gene_name
gene_intervals = get_gene_intervals(available_mutations, genes)
print('drug asscociaton')
drug_association = get_drug_gene_association_from_dictionary(available_mutations, drug_list)

prediction = []

for sample_info in data:
    
    name = sample_info.split(',')[0]
    # if there is no mutations, only name
    if name[:-1] == '\n':
        name = name[:-1]

    if name in subset:

        if coverage_dict[name] != {}:
            uncovered_drugs = check_dictionary_position_for_coverage(available_mutations, coverage_dict, coverage_intervals, genes, name)
            uncovered_drugs += check_genes_for_coverage(gene_intervals, coverage_dict, coverage_intervals, drug_association, name)
        else:
            uncovered_drugs = []

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
                if drug_list[i] in uncovered_drugs:
                    prediction[-1] += [0]
                else:
                    prediction[-1] += [-1]    
                


file = open(OUTPUT, 'w')
for line in prediction:
    for el in line[:-1]:
        file.write(str(el) + '\t')
    file.write(str(el) + '\n')

file.close()
