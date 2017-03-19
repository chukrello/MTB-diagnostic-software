#
# Scripts helps to generate comand which cuts .bam file
# in 30.01
# Original command: samtools view -b SAMN03647317_h37rv.bam ch1:1560000-1560900 > test.bam
# 

from Bio import SeqIO
from Bio.Seq import Seq

PATH = '/export/data/kkuleshov/myc/sra/'
PATH_input = '/home/chukreev/PycharmProjects/tuberculosis/resources/reference/'

name,sequence = '',''
fasta_sequences = SeqIO.parse(open(PATH_input+'AL123456_rev.fa'),'fasta')
for fasta in fasta_sequences:
    name, sequence = fasta.id, fasta.seq.tostring()

def cut_bam(data, header):
    for i in range(len(data)):
        print('samtools view -b ' + PATH + data[i][0] + '/' + data[i][0] + '_h37rv.bam ch1:'+str(data[i][1])+'-'+str(data[i][2])+' > ' + header + data[i][0]+'.bam')
    print(' ')

def make_data(samns, coords):
    data = []
    for i in range(len(samns)):
        data.append([samns[i], coords[i][0], coords[i][1]])
    return data

def index_bam(data, header):
    for i in range(len(data)):
        print('samtools index ' + header + data[i][0]+'.bam')
    print(' ')

def get_coord(gene, coord, isAminoAcidCoord):
    def get_genes_info(sequence):

        import re
        lines = open(PATH_input+'AL123456_rev.gff').readlines()
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

    if isAminoAcidCoord == True:
        if genes[gene][2] == '+':
            abs_coord = genes[gene][0] + 1 + (coord - 1)*3
        else:
            abs_coord = genes[gene][1] - 1 - (coord - 1)*3
    else:
        if genes[gene][2] == '+':
            abs_coord = genes[gene][0] + (coord - 1)
        else:
            abs_coord = genes[gene][1] - (coord - 1)

    return [abs_coord -5, abs_coord+5]


samns = [line[:-1] for line in open('bams_for_cutting.txt').readlines()]
coords = [get_coord('katG', 315, True)]*len(samns)
header = 'nonwalker_EMB_test_'


cut_bam(make_data(samns, coords), header)
index_bam(make_data(samns, coords), header)