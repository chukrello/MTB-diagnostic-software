import os
import re
from Bio import SeqIO
from Bio.Seq import Seq

#
# Script helps you to extract sequence from reference
#

name,sequence = '',''
PATH_input = '/home/chukreev/PycharmProjects/tuberculosis/resources/reference/'
# PATH_input = '/home/chukreev/Downloads/sequence.fasta'
fasta_sequences = SeqIO.parse(open(PATH_input+'AL123456_rev.fa'),'fasta')
# fasta_sequences = SeqIO.parse(open(PATH_input),'fasta')
for fasta in fasta_sequences:
    name, sequence = fasta.id, fasta.seq.tostring()

# genes = [[2938152,2939687]]

def get_genes_info(sequence):
    lines = open(PATH_input+'AL123456_rev.gff').readlines()
    gene_reg = re.compile('\tgene \w*')
    genes = []
    for line in lines:
        m = gene_reg.search(line)
        if m:
            genes.append(line)
    genes = [[line[:-1].split('\t')[3],line[:-1].split('\t')[4], line[:-1].split('\t')[8].split(' ')[1], line[:-1].split('\t')[6]] for line in genes]
    list_genes = list(set(['ahpC', 'eis', 'embA', 'embB', 'embC', 'embR', 'fabG1', 'gid', 'gyrA', 'gyrB', 'inhA', 'iniA', 'iniC',
             'katG', 'manB', 'ndh', 'pncA', 'rmlD', 'rpoB', 'rpsA', 'rpsL', 'rrs', 'tlyA', 'rrl', 'rpsA', 'rpsL', 'pncA', 'rrs', 'rplC', 'thyA', 'Rv0678', 'ahpC', 'folC', 'panD', 'tlyA', 'katG', 'inhA', 'gyrB', 'gyrA', 'rpoC', 'rpoB', 'ribD', 'eis', 'embA', 'embB', 'embC', 'ethA', 'fabG1', 'embR', 'kasA', 'ethR']))
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

print(sequence[4076845-1])

# for el in genes:
#     print(el)
#     print(sequence[el[0]-1:el[1]])
