# Script build file with mutations identified in samns
#
# Output format:
# SAM___: mut1 \t mut2 \t ....
#
# Input line:
# build_full_samn_mutation_list.py list_samns dictionary output_file


import os
from Bio import SeqIO
from Bio.Seq import Seq
import sys


PATH_to_mutations = '/media/chukreev/06368dc5-6760-4d76-9ff4-c6f01b0ffc7e/Untitled Folder/mutations/'
PATH_to_fasta = '../../1_input/'


list_samns = [line[:-1].split('\t')[0]+'.txt' for line in open(sys.argv[1]).readlines()]
dictionary = [line[:-1].split(' ') for line in open(sys.argv[2]).readlines()]
file = open(sys.argv[3], 'w')

name,sequence = '',''
fasta_sequences = SeqIO.parse(open('AL123456_rev.fa'),'fasta')
for fasta in fasta_sequences:
    name, sequence = fasta.id, fasta.seq.tostring()

def get_genes_info(sequence):

    import re
    lines = open('AL123456_rev.gff').readlines()
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


def find_mutations_in_file_by_dict(data, dictionary, file):
    results = []

    for el in data:
        if el[0] == 'Gene':

            gene_data = el[1]
            pos_data = el[2]
            alt_data = el[4]

            for mutation in dictionary:

                gene = mutation[0]
                pos = int(mutation[1][1:-1])
                alt = mutation[1][-1]

                # print(gene)
                # print(pos)
                # print(ref)
                # print(alt)

                if gene == gene_data and int(pos) == int(pos_data) and alt == alt_data:
                    results.append(mutation)
        else:

            pos = el[1]
            alt = el[3]

            for mutation in dictionary:

                gene = mutation[0]

                pos_data = int(mutation[1][1:-1])
                alt_data = mutation[1][-1]

                if int(pos_data) > 3000:

                    if genes[gene][-1] == '+':

                        if int(pos) == int(pos_data) and alt == alt_data:
                            results.append(mutation)

                    else:

                        alt_data = str(Seq(alt_data).reverse_complement())

                        if int(pos) == int(pos_data) and alt == alt_data:
                            results.append(mutation)
                            print(mutation)


    for res in results[:-1]:
        file.write(res[0] + ' ' + res[1] + ' ' + res[2] + '\t')
    if len(results) > 0:
        res = results[-1]
        file.write(res[0] + ' ' + res[1] + ' ' + res[2] + '\n')
    else:
        file.write('\n')

for name in list_samns:
    file.write(name[:-4]+'\t')
    print(name)
    try:
        data = [line[:-1].split('\t') for line in open(PATH+name).readlines()]
    except Exception:
        continue

    for i in range(len(data)):
        if data[i][0] == 'Gene':
            data[i][2] = int(float(data[i][2]))
            if data[i][1] == 'gid':
                data[i][1] = 'gidB'

    find_mutations_in_file_by_dict(data, dictionary, file)

file.close()

