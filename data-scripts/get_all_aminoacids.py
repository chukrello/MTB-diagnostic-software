# IMPORTANT SCRIPT
# Converting .vcf format in my format
# python get_all_aminoacids.py *list_samns* *qval* *n_jobs* *header*

import os
from Bio import SeqIO
from Bio.Seq import Seq
import sys

LIST_VCFS = sys.argv[1]
QVAL = int(sys.argv[2])
N_JOBS = int(sys.argv[3])
HEADER = sys.argv[4]

list_vcfs = [line[:-1] for line in open(LIST_VCFS).readlines()]

# getting kuleshov fasta reference
name,sequence = '',''
PATH_input = '../1_input/'
fasta_sequences = SeqIO.parse(open(PATH_input+'AL123456_rev.fa'),'fasta')
for fasta in fasta_sequences:
    name, sequence = fasta.id, fasta.seq.tostring()

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
    list_genes = [el[2] for el in genes]
    def get_sequence(genes, gene, sequence):
        for element in genes:
            if element[2] == gene:
                return [int(element[0]), int(element[1]), element[3], sequence[int(element[0])-1:int(element[1])]]
    genes_seqs = dict()
    for gene in list_genes:
        genes_seqs[gene] = get_sequence(genes,gene,sequence)
    return genes_seqs
genes = get_genes_info(sequence) #start end strand seq

def BuildAminoAcids(sequence, strand, start_pos, isFirstModified, isSecondModified, isThirdModified, alt1, alt2, alt3):
    old_triplet = ''
    new_triplet = ''

    if strand == '+':

        old_triplet = sequence[start_pos:start_pos+3]
        if isFirstModified:
            new_triplet += alt1
        else:
            new_triplet += sequence[start_pos]

        if isSecondModified:
            new_triplet += alt2
        else:
            new_triplet += sequence[start_pos+1]


        if isThirdModified:
            new_triplet += alt3
        else:
            new_triplet += sequence[start_pos+2]

    else:

        old_triplet = str(Seq(sequence[start_pos-2:start_pos+1]).reverse_complement())

        if isFirstModified:
            new_triplet += alt1
        else:
            new_triplet += sequence[start_pos-2]


        if isSecondModified:
            new_triplet += alt2
        else:
            new_triplet += sequence[start_pos-1]

        if isThirdModified:
            new_triplet += alt3
        else:
            new_triplet += sequence[start_pos]

        new_triplet = str(Seq(new_triplet).reverse_complement())

    old_aminoacid = str(Seq(old_triplet).translate())
    new_aminoacid = str(Seq(new_triplet).translate())

    return [old_aminoacid, new_aminoacid]

def look_vcf_file(filename, name):
    #for SNP
    # open .vcf file

    print(name)

    file = open('../../data/' + HEADER + '/' + name +'.txt', 'w')

    info = []
    raw_data_2 = [[int(line.split('\t')[1]), line.split('\t')[3], line.split('\t')[4]] for line in open(filename).readlines() if line[0] != '#' and float(line.split('\t')[5]) >= QVAL]

    data = []
    for el in raw_data_2:
        if el not in data:
            data.append(el)

    for i in range(len(data)):
        pos = data[i][0]
        ref = data[i][1]
        alt = data[i][2]

        cur_gene = ''

        if len(ref) == len(alt):
            #SNP ANALYSIS

            for key in genes:
                if int(genes[key][0]) <= int(pos) and int(pos) <= int(genes[key][1]):
                    cur_gene = key
                    break

            if cur_gene != '':
                if genes[cur_gene][2] == '+':
                    protein_pos = (pos - genes[cur_gene][0])/3 + 1
                    nucleotide_pos = (pos - genes[cur_gene][0])%3

                    #NEIGHBOURS ANALYSIS
                    if nucleotide_pos == 0:
                        if i!=len(data)-1 and data[i+1][0] == pos+1:
                            if i!=len(data)-2 and data[i+2][0] == pos+2:
                                aminoacids = BuildAminoAcids(sequence, "+", pos-1, True, True, True, alt, data[i+1][2], data[i+2][2])
                            else:
                                aminoacids = BuildAminoAcids(sequence, "+", pos-1, True, True, False, alt, data[i+1][2], '-')
                        elif i!=len(data)-1 and data[i+1][0] == pos+2:
                            aminoacids = BuildAminoAcids(sequence, "+", pos-1, True, False, True, alt, '-', data[i+1][2])
                        else:
                            aminoacids = BuildAminoAcids(sequence, "+", pos-1, True, False, False, alt, '-', '-')

                    elif nucleotide_pos == 1:
                        if i!=0 and data[i-1][0] == pos-1:
                            if i!=len(data)-1 and data[i+1][0] == pos+1:
                                aminoacids = BuildAminoAcids(sequence, "+", pos-2, True, True, True, data[i-1][2], alt, data[i+1][2])
                            else:
                                aminoacids = BuildAminoAcids(sequence, "+", pos-2, True, True, False, data[i-1][2], alt, '-')
                        elif i!=len(data)-1 and data[i+1] == pos+1:
                            aminoacids = BuildAminoAcids(sequence, "+", pos-2, False, True, True, '-', alt, data[i+1][2])
                        else:
                            aminoacids = BuildAminoAcids(sequence, "+", pos-2, False, True, False, '-', alt, '-')
                    else:
                        if i != 0 and data[i-1][0] == pos-1:
                            if i != 1 and data[i-2][0] == pos-2:
                                aminoacids = BuildAminoAcids(sequence, "+", pos-3, True, True, True, data[i-2][2], data[i-1][2], alt)
                            else:
                                aminoacids = BuildAminoAcids(sequence, "+", pos-3, False, True, True, '-', data[i-1][2], alt)
                        elif i != 0 and data[i-1][0] == pos-2:
                            aminoacids = BuildAminoAcids(sequence, "+", pos-3, True, False, True, data[i-2][2], '-', alt)
                        else:
                            aminoacids = BuildAminoAcids(sequence, "+", pos-3, False, False, True, '-', '-', alt)

                    old_aminoacid = aminoacids[0]
                    new_aminoacid = aminoacids[1]

                    file.write('Gene\t'+cur_gene+'\t'+str(protein_pos)+'\t'+old_aminoacid+'\t'+new_aminoacid+'\n')

                else: # if strand is '-'
                    if cur_gene == 'rrs':
                        coord = genes[cur_gene][1] - (pos-1)
                        reference_nucl = Seq(sequence[coord-1]).reverse_complement()
                        file.write('Gene\t'+cur_gene+'\t'+str(scoord)+'\t'+str(reference_nucl)+'\t'+str(Seq(alt).reverse_complement())+'\n')
                    else:

                        protein_pos = (genes[cur_gene][1] - pos)/3 + 1
                        nucleotide_pos = (genes[cur_gene][1] - pos)%3

                        #NEIGHBOURS ANALYSIS
                        if nucleotide_pos == 2:
                            if i!=len(data)-1 and data[i+1][0] == pos+1:
                                if i!=len(data)-2 and data[i+2][0] == pos+2:
                                    aminoacids = BuildAminoAcids(sequence, "-", pos+1, True, True, True, alt, data[i+1][2], data[i+2][2])
                                else:
                                    aminoacids = BuildAminoAcids(sequence, "-", pos+1, True, True, False, alt, data[i+1][2], '-')
                            elif i!=len(data)-1 and data[i+1][0] == pos+2:
                                aminoacids = BuildAminoAcids(sequence, "-", pos+1, True, False, True, alt, '-', data[i+1][2])
                            else:
                                aminoacids = BuildAminoAcids(sequence, "-", pos+1, True, False, False, alt, '-', '-')

                        elif nucleotide_pos == 1:
                            if i!=0 and data[i-1][0] == pos-1:
                                if i!=len(data)-1 and data[i+1][0] == pos+1:
                                    aminoacids = BuildAminoAcids(sequence, "-", pos, True, True, True, data[i-1][2], alt, data[i+1][2])
                                else:
                                    aminoacids = BuildAminoAcids(sequence, "-", pos, True, True, False, data[i-1][2], alt, '-')
                            elif i!=len(data)-1 and data[i+1] == pos+1:
                                aminoacids = BuildAminoAcids(sequence, "-", pos, False, True, True, '-', alt, data[i+1][2])
                            else:
                                aminoacids = BuildAminoAcids(sequence, "-", pos, False, True, False, '-', alt, '-')
                        else:
                            if i != 0 and data[i-1][0] == pos-1:
                                if i != 1 and data[i-2][0] == pos-2:
                                    aminoacids = BuildAminoAcids(sequence, "-", pos-1, True, True, True, data[i-2][2], data[i-1][2], alt)
                                else:
                                    aminoacids = BuildAminoAcids(sequence, "-", pos-1, False, True, True, '-', data[i-1][2], alt)
                            elif i != 0 and data[i-1][0] == pos-2:
                                aminoacids = BuildAminoAcids(sequence, "-", pos-1, True, False, True, data[i-2][2], '-', alt)
                            else:
                                aminoacids = BuildAminoAcids(sequence, "-", pos-1, False, False, True, '-', '-', alt)

                        old_aminoacid = aminoacids[0]
                        new_aminoacid = aminoacids[1]

                        file.write('Gene\t'+cur_gene+'\t'+str(protein_pos)+'\t'+old_aminoacid+'\t'+new_aminoacid+'\n')
            else:
                file.write('NotGene\t'+str(pos)+'\t'+ref+'\t'+alt+'\n')

    file.close()
    return info

def get_aminoacids(file):
    name = file.split('_')[0]
    filename = file
    look_vcf_file(filename, name)

os.system('mkdir ../../data/' + HEADER)
from joblib import Parallel, delayed
Parallel(n_jobs=N_JOBS)(delayed(get_aminoacids)(file) for file in list_vcfs)

print('DONE')