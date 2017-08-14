import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description='Script compares .fasta files by genes. Change list of genes in the code')

parser.add_argument('-gff', '--gff_file', help = '.gff file for extracting coordinates', default = '/export/data/kchukreev/MTB-diagnostic-software/1_input/AL123456_rev.gff')
parser.add_argument('-fasta1', '--fasta_number_1', help = 'Path to first reference', default = '/export/data/kchukreev/MTB-diagnostic-software/1_input/AL123456_rev.fa')
parser.add_argument('-fasta2', '--fasta_number_2', help = 'Path to second reference', default = '/export/data/kchukreev/MTB-diagnostic-software/1_input/h37rv.fasta')

args = parser.parse_args()

gff_file = args.gff_file
ref1 = args.fasta_number_1
ref2 = args.fasta_number_2

list_genes = ['embB', 'rrs', 'katG', 'ethA', 'ethR', 'thyA', 'fabD', 'kasA', 'gyrA', 'tlyA', 'inhA', 'folC', 'pncA', 'rpsA', 'eis', 'efpA', 'ribD', 'embC', 'panD', 'embR', 'ndh', 'iniA', 'rplC', 'accD6', 'iniB', 'iniC', 'gyrB', 'nat', 'ahpC', 'rpoB', 'rpoC', 'rrl', 'gid', 'fabG1', 'embA', 'fbpC', 'rpsL']

def get_genes_info():
    import re
    lines = open(gff_file).readlines()
    gene_reg = re.compile('\tgene \w*')
    genes = []
    for line in lines:
        m = gene_reg.search(line)
        if m:
            genes.append(line)
    genes = [[line[:-1].split('\t')[3],line[:-1].split('\t')[4], line[:-1].split('\t')[8].split(' ')[1], line[:-1].split('\t')[6]] for line in genes]
    def get_sequence(genes, gene):
        for element in genes:
            if element[2] == gene:
                return [int(element[0]), int(element[1]), element[3]]
    genes_seqs = dict()
    for gene in list_genes:
        genes_seqs[gene] = get_sequence(genes,gene)

    return genes_seqs

def get_fasta_file(path_fasta):
    fasta_sequences = SeqIO.parse(open(path_fasta),'fasta')
    for fasta in fasta_sequences:
        name, sequence = fasta.id, fasta.seq.tostring()
    
    return sequence

fasta1 = get_fasta_file(ref1)
fasta2 = get_fasta_file(ref2)
genes = get_genes_info()

print(genes)

for gene in genes:
    str1 = fasta1[genes[gene][0]-1:genes[gene][1]-1]
    str2 = fasta2[genes[gene][0]-1:genes[gene][1]-1]

    if str1 == str2:
        print('OK!')
    else:
        print('NOT OK!')

if fasta1 == fasta2:
    print('OK!')
else:
    print('NOT OK!')

print(len(fasta1))
print(len(fasta2))