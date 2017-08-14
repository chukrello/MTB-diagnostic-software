# Script build file with mutations identified in samns
#
# Output format:
# SAM___, mut1, mut2, ....
#
# Input line:
# build_full_samn_mutation_list.py list_samns dictionary output_file


import os
from Bio import SeqIO
from Bio.Seq import Seq
import sys

LIST_SAMNS = sys.argv[1]
DICTIONARY = sys.argv[2]
HEADER = sys.argv[3]
OUTPUT = sys.argv[4]

PATH_to_mutations = '../../../data/' + HEADER + '/'
PATH_to_fasta = '../../1_input/'

list_samns = [line[:-1]+'.txt' for line in open(LIST_SAMNS).readlines()]
dictionary = [line[:-1].split('\t') for line in open(DICTIONARY).readlines()]
file = open(OUTPUT, 'w')


#Getting reference ang .gff

name,sequence = '',''
fasta_sequences = SeqIO.parse(open(PATH_to_fasta + 'AL123456_rev.fa'),'fasta')
for fasta in fasta_sequences:
    name, sequence = fasta.id, fasta.seq.tostring()

def get_genes_info(sequence):
    import re
    lines = open(PATH_to_fasta + 'AL123456_rev.gff').readlines()
    gene_reg = re.compile('\tgene \w*')
    genes = []
    for line in lines:
        m = gene_reg.search(line)
        if m:
            genes.append(line)
    genes = [[line[:-1].split('\t')[3],line[:-1].split('\t')[4], line[:-1].split('\t')[8].split(' ')[1], line[:-1].split('\t')[6]] for line in genes]
    list_genes = ['ahpC', 'eis', 'embA', 'embB', 'embC', 'embR', 'fabG1', 'gid', 'gyrA', 'gyrB', 'inhA', 'iniA', 'iniC',
             'katG', 'manB', 'ndh', 'pncA', 'rmlD', 'rpoB', 'rpsA', 'rpsL', 'rrs', 'tlyA']
    def get_sequence(genes, gene, sequence):
        for element in genes:
            if element[2] == gene:
                return [int(element[0]), int(element[1]), element[3], sequence[int(element[0])-1:int(element[1])]]
    genes_seqs = dict()
    for gene in list_genes:
        genes_seqs[gene] = get_sequence(genes,gene,sequence)
    return genes_seqs
genes = get_genes_info(sequence) #start end strand seq
genes['gidB'] = genes['gid']
del(genes['gid'])


def find_mutations_in_file_by_dict(name, data, dictionary, file):
    results = []

    for line in data:
        for mutation in dictionary:
            if line[0] == mutation[0] == 'Gene':
                if line[1] == 'rrs':
                    #BUG FIXING:
                    if line[2] == mutation[2] and line[4] == mutation[4]:
                        results.append(mutation)
                else:
                    if line[1] == mutation[1] and line[2] == mutation[2] and line[3] == mutation[3] and line[4] == mutation[4]:
                        results.append(mutation)

            elif line[0] == mutation[0] == 'NotGene':
                if line[1] == mutation[2] and line[2] == mutation[3] and line[3] == mutation[4]:
                    results.append(mutation)

    if len(results) == 0:
        file.write(name[:-4] + '\n')
    else:
        file.write(name[:-4] + ',')
        for mutation in results[:-1]:
            file.write(' '.join(mutation) + ',')
        file.write(' '.join(results[-1]) + '\n')


# Script run

for name in list_samns:
    try:
        data = [line[:-1].split('\t') for line in open(PATH_to_mutations+name).readlines()]
    except Exception:
        continue

    for i in range(len(data)):
        if data[i][0] == 'Gene':
            if data[i][1] == 'gid':
                data[i][1] = 'gidB'
    find_mutations_in_file_by_dict(name, data, dictionary, file)

file.close()

