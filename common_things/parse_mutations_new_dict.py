#
# Scripts compares REF letter from Walker dictionary and letter in our reference
#

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_protein

PATH_input = '/home/chukreev/PycharmProjects/tuberculosis/resources/reference/'
name, sequence = '', ''
fasta_sequences = SeqIO.parse(open(PATH_input+'AL123456_rev.fa'),'fasta')
for fasta in fasta_sequences:
    name, sequence = fasta.id, fasta.seq.tostring()

def make_walker_database():
    PATH = '/home/chukreev/PycharmProjects/tuberculosis/'
    lines = [line[:-1].split('\t') for line in open(PATH+'new_database.txt').readlines()]
    mutations = []
    for mut in lines:
        isIndel = False
        if mut[-1] == '-':
            changes = mut[4].split('/')
            for change in changes:
                if len(change) > 1 or change == '-':
                    isIndel = True
                    break

        if mut[1] == 'rrs' or mut[1] == 'rrl':
            if isIndel == False:
                #gene mutation, with codone change
                mutations.append(['snp', mut[1], mut[3], mut[4].split('/')[0], mut[4].split('/')[1], mut[0]])

            else:
                if len(changes[0]) > 1:
                    mutations.append(['del', mut[1], mut[3], mut[4].split('/')[0], mut[4].split('/')[1], mut[0]])
                else:
                    mutations.append(['ins', mut[1], mut[3], mut[4].split('/')[0], mut[4].split('/')[1], mut[0]])
        else:
            if isIndel == False:
                if mut[-1] != '-':
                    #gene mutation, with codone change
                    mutations.append(['snp', mut[1], mut[5], mut[-1].split('/')[0], mut[-1].split('/')[1], mut[0]])
                else:
                    if mut[1].split('_')[1] != 'promoter':
                        print mut
                        exit(0)
                    #upstream change, no codone
                    mutations.append(['snp', mut[1].split('_')[0], mut[3], mut[4].split('/')[0], mut[4].split('/')[1], mut[0]])
            else:
                if len(changes[0]) > 1:
                    mutations.append(['del', mut[1], mut[3], mut[4].split('/')[0], mut[4].split('/')[1], mut[0]])
                else:
                    mutations.append(['ins', mut[1], mut[3], mut[4].split('/')[0], mut[4].split('/')[1], mut[0]])
    return mutations
database_mutations = make_walker_database()

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
# searching for the aminoacid

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
    #nucl < 0
    if genes[gene][2] == '+':
        pos = genes[gene][0] - 1 + nucl
        return sequence[pos]
    else:
        pos = genes[gene][1] - 1 - nucl
        return Seq(sequence[pos], generic_dna).reverse_complement().tostring()
    

    
bad_genes = []


for m in database_mutations:
    if m[1] == 'gyrB':
        m[2] = str(int(m[2]) + 39)

    if m[1] in ['eis','pncA', 'katG'] and m[2][0] == '-':
        m[3] = str(Seq(m[3]).reverse_complement())

    if m[2][0] != '-' and m[0] == 'snp':
        if search_aminoacid(int(m[2]), m[1], sequence) != m[3]:
            bad_genes.append(m)

    if m[2][0] == '-' and m[0] == 'snp':
        if search_upstream_snp(int(m[2]), m[1], sequence) != m[3]:
            bad_genes.append(m)
            

for el in bad_genes:
    print(el)

# list_genes = ['ahpC', 'eis', 'embA', 'embB', 'embC', 'embR', 'fabG1', 'gid', 'gyrA', 'gyrB', 'inhA', 'iniA', 'iniC',
#              'katG', 'manB', 'ndh', 'pncA', 'rmlD', 'rpoB', 'rpsA', 'rpsL', 'rrs', 'tlyA']
# list_genes = [[el] for el in list_genes]

# for i in range(len(list_genes)):
#     for gene in genes:
#         if gene[-2] == list_genes[i][0]:
#             list_genes[i] = [gene[0], gene[1], list_genes[i][0]]

# def get_mutations(gene, sequence):
#     print(gene[-1])
#     print(sequence[int(gene[0])-1:int(gene[1])] + '\n')