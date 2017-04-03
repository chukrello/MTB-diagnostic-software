# script checks if reference aminoacid is right

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_protein

PATH_input = '/home/chukreev/PycharmProjects/tuberculosis/resources/reference/'
name, sequence = '', ''
fasta_sequences = SeqIO.parse(open(PATH_input+'AL123456_rev.fa'),'fasta')
for fasta in fasta_sequences:
    name, sequence = fasta.id, fasta.seq.tostring()

NAME = 'Tbdb'


def make_dict_database(name):
    PATH = '/home/chukreev/PycharmProjects/tuberculosis/'
    mutations = [line[:-1].split(' ') for line in open('../'+name+'.txt').readlines()]
    
    return mutations
database_mutations = make_dict_database(NAME)

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
    list_genes = list(set(['ahpC', 'eis', 'embA', 'embB', 'embC', 'embR', 'fabG1', 'gid', 'gyrA', 'gyrB', 'inhA', 'iniA', 'iniC',
             'katG', 'manB', 'ndh', 'pncA', 'rmlD', 'rpoB', 'rpsA', 'rpsL', 'rrs', 'tlyA', 'rrl', 'rpsA', 'rpsL', 'pncA', 'rrs', 'rplC', 'thyA', 'Rv0678', 'ahpC', 'folC', 'panD', 'tlyA', 'katG', 'inhA', 'gyrB', 'gyrA', 'rpoC',
              'rpoB', 'ribD', 'eis', 'embA', 'embB', 'embC', 'ethA', 'fabG1', 'embR', 'kasA', 'ethR', 'iniB', 'accD6', 'efpA', 'fabD', 
              'fbpC', 'furA', 'nat', 'srmR']))
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



def search_aminoacid(prot_pos, gene, sequence):
    start = 0
    if genes[gene][2] == '+':
        pos = genes[gene][0] - 1 + (prot_pos - 1)*3
        triplet = sequence[pos:pos+3]
        protein = str(Seq(triplet).translate())
    else:

        if gene in ['rrs', 'rrl']:
            pos = genes[gene][1] - 1 - (prot_pos - 1)
            return str(Seq(sequence[pos]).reverse_complement())            

        pos = genes[gene][1] - 1 - (prot_pos - 1) * 3
        triplet = sequence[pos-2:pos+1]
        protein = str(Seq(triplet).reverse_complement().translate())
    return protein

def search_upstream_snp(nucl,gene,sequence):
    return sequence[nucl-1]

    
bad_elements = []


for m in database_mutations:
    if int(m[1][1:-1]) < 4000:
        if search_aminoacid(int(m[1][1:-1]), m[0], sequence) != m[1][0]:
            bad_elements.append(m)
    else:
        if search_upstream_snp(int(m[1][1:-1]), m[0], sequence) != m[1][0]:
            bad_elements.append(m)

for el in bad_elements:
    print(el)
