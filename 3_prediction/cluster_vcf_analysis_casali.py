import os
from Bio import SeqIO
from Bio.Seq import Seq

name,sequence = '',''
PATH_input = '/export/data/kchukreev/1_input/'
fasta_sequences = SeqIO.parse(open(PATH_input+'AL123456_rev.fa'),'fasta')
for fasta in fasta_sequences:
    name, sequence = fasta.id, fasta.seq.tostring()

def make_walker_database():
    lines = [line[:-1].split('\t') for line in open('walker_database.txt').readlines()]
    mutations = []
    for mut in lines:
        if 'del' in mut[0]:
            gene = mut[0].split('_')[0]
            pos = mut[0].split('_')[1]
            letters = mut[0].split('_')[2][3:]
            mutations.append(['del', gene, pos, letters, mut[1], mut[2]])
        elif 'ins' in mut[0]:
            gene = mut[0].split('_')[0]
            pos = mut[0].split('_')[1]
            letters = mut[0].split('_')[2][3:]
            mutations.append(['ins', gene, pos, letters, mut[1], mut[2]])
        else:
            try:
                gene = mut[0].split('_')[0]
                change = mut[0].split('_')[1]
                from_let = change[0]
                to_let = change[-1]
                pos = change[1:-1]
                mutations.append(['snp', gene, pos, from_let, to_let, mut[1], mut[2]])
            except Exception:
                print(mut)

    return mutations
database_genes = make_walker_database()

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

def getting_walker_samns(names, kul_list):
    names = [line[:-1] for line in open(names).readlines()]
    kul_list = [line[:-1].split('\t') for line in open(kul_list).readlines()]

    for i in range(len(names)):
        if names[i][:3] != 'SAM':
            name = names[i].split(',')[0]
            for el in kul_list:
                if name == el[4] or name == el[5]:
                    names[i] = el[1]

    return names
casali_names = [line[:-1] for line in open(PATH_input+'casali_sample.txt').readlines()]

def make_upstream_snps(filename, genes):
    upstream_snps = []
    lines = [line[:-1].split('\t') for line in open(filename).readlines()]

    for line in lines:
        if '-' in line[0]:
            if 'del' not in line[0]:
                coord_otn = int(line[0].split('_')[1][2:-1])
                gene = line[0].split('_')[0]
                nucleotide = line[0][-1]
                if genes[gene][2] == '+':
                    start_gene = genes[gene][0]
                    coord_abs = sequence[start_gene-1-coord_otn]
                else:
                    start_gene = genes[gene][1]
                    coord_abs = sequence[start_gene-1+coord_otn]

                upstream_snps.append([gene,coord_abs,nucleotide])

    return upstream_snps
upstream_snps = make_upstream_snps('walker_database.txt', genes)

def look_vcf_file(filename, name, upstream_snp, file):
    #for SNP
    print('Looking at ' + name)
    info = []
    file.write(name+'\n')
    # open vcf file
    data = [[int(line.split('\t')[1]), line.split('\t')[3], line.split('\t')[4]] for line in open(filename).readlines() if line[0] != '#']

    for mutation in data:

        pos = mutation[0]
        ref = mutation[1]
        alt = mutation[2]

        ups = []

        # checking for upstream snps
        for snp in upstream_snp:
            if snp[0] == pos and snp[2] == alt:
                ups.append(snp)

        for el in ups:
            el[0] = 'upstream_snp'
            file.write(' '.join(el) + '\n')
            info.append(ups)

        cur_gene = ''

        if len(ref) == len(alt):
            #SNP

            for key in genes:
                if int(genes[key][0]) <= int(pos) and int(pos) <= int(genes[key][1]):
                    cur_gene = key
                    # print('Gene is ' + key)
                    break
            if cur_gene != '':
                if genes[cur_gene][2] == '+':
                    protein_pos = (pos - genes[cur_gene][0])/3 + 1
                    nucleotide_pos = (pos - genes[cur_gene][0])%3

                    for walker_mut in database_genes:
                        if walker_mut[1] == cur_gene and walker_mut[2][0] != '-' and int(walker_mut[2]) == protein_pos:
                            # print('HERE!')
                            if nucleotide_pos == 0:
                                old_triplet = sequence[pos-1:pos+2]
                                new_triplet = alt + sequence[pos] + sequence[pos+1]
                            elif nucleotide_pos == 1:
                                old_triplet = sequence[pos-2:pos+1]
                                new_triplet = sequence[pos-2] + alt + sequence[pos]
                            elif nucleotide_pos == 2:
                                old_triplet = sequence[pos-3:pos]
                                new_triplet = sequence[pos-3] + sequence[pos-2] + alt

                            # if (Seq(old_triplet).translate() == walker_mut[3] or walker_mut[3] == '*') and (Seq(new_triplet).translate() == walker_mut[4] or walker_mut[4] == '*'):
                            if (Seq(old_triplet).translate() == walker_mut[3]) and (Seq(new_triplet).translate() == walker_mut[4]):
                                # print('WOW')
                                file.write(' '.join(walker_mut) + '\n')
                                info.append(walker_mut)

                else: # if strand is '-'

                    if cur_gene == 'rrs':
                        walker_coord_start = genes[cur_gene][1] - (pos-1)
                        reference_nucl = Seq(sequence[walker_coord_start-1]).reverse_complement()

                        for walker_mut in database_genes:
                            if walker_mut[0] == 'snp' and walker_mut[1] == 'rrs':
                                if int(walker_mut[2]) == int(walker_coord_start) and str(Seq(alt).reverse_complement()) == walker_mut[4]:
                                    print('yes')
                                    print(' '.join(walker_mut))
                                    file.write(' '.join(walker_mut) + '\n')
                                    info.append(walker_mut)
                    else:

                        protein_pos = (genes[cur_gene][1] - pos)/3 + 1
                        nucleotide_pos = (genes[cur_gene][1] - pos)%3

                        for walker_mut in database_genes:
                            if walker_mut[1] == cur_gene and walker_mut[2][0] != '-' and int(walker_mut[2]) == protein_pos:
                                # print('HERE!')
                                if nucleotide_pos == 0:
                                    old_triplet = str(Seq(sequence[pos-3:pos]).reverse_complement())
                                    new_triplet = str(Seq(sequence[pos-3] + sequence[pos-2] + alt).reverse_complement())
                                elif nucleotide_pos == 1:
                                    old_triplet = str(Seq(sequence[pos-2:pos+1]).reverse_complement())
                                    new_triplet = str(Seq(sequence[pos-2] + alt + sequence[pos]).reverse_complement())
                                elif nucleotide_pos == 2:
                                    old_triplet = str(Seq(sequence[pos-1:pos+2]).reverse_complement())
                                    new_triplet = str(Seq(alt + sequence[pos] + sequence[pos+1]).reverse_complement())

                                # if (Seq(old_triplet).translate() == walker_mut[3] or walker_mut[3] == '*') and (Seq(new_triplet).translate() == walker_mut[4] or walker_mut[4] == '*'):
                                if (Seq(old_triplet).translate() == walker_mut[3]) and (Seq(new_triplet).translate() == walker_mut[4]):
                                    # print('WOW')
                                    file.write(' '.join(walker_mut) + '\n')
                                    info.append(walker_mut)
        elif len(ref) > len(alt):
            # deletion
            delta = len(ref) - len(alt)
            special = ref[1:delta+1]
            pos = pos + 1

            cur_gene = ''

            for key in genes:
                if int(genes[key][0]) <= int(pos) and int(pos) <= int(genes[key][1]):
                    cur_gene = key
                    print('Gene is ' + key)
                    break

            if cur_gene != '':

                if genes[cur_gene][2] == '+':
                    coord_otn = pos - genes[cur_gene][0] + 1
                else:
                    coord_otn = genes[cur_gene][1] - pos + 1

                for walker_mut in database_genes:
                    if walker_mut[0] == 'del' and cur_gene == walker_mut[1] and coord_otn == int(walker_mut[2]):
                        print('WOW')
                        print(walker_mut)
                        print(str(coord_otn) + ' ' + str(pos) + ' ' + str(special))
                        file.write(' '.join(walker_mut) + '\n')
                        info.append(walker_mut)
                    if walker_mut[0] == 'del' and cur_gene == walker_mut[1] and abs(int(coord_otn) - int(walker_mut[2])) < 10:
                        print(walker_mut)
                        print(str(coord_otn) + ' ' + str(pos) + ' ' + str(special))

        elif len(alt) > len(ref):
            # insertion
            delta = len(alt) - len(ref)
            special = alt[1:delta+1]
            pos = pos + 1

            cur_gene = ''

            for key in genes:
                if int(genes[key][0]) <= int(pos) and int(pos) <= int(genes[key][1]):
                    cur_gene = key
                    print('Gene is ' + key)
                    break

            if cur_gene != '':

                if genes[cur_gene][2] == '+':
                    coord_otn = pos - genes[cur_gene][0] + 1
                else:
                    coord_otn = genes[cur_gene][1] - pos + 1

                for walker_mut in database_genes:
                    if walker_mut[0] == 'ins' and cur_gene == walker_mut[1] and coord_otn == int(walker_mut[2]):
                        print('WOW')
                        print(walker_mut)
                        print(str(coord_otn) + ' ' + str(pos) + ' ' + str(special))
                        file.write(' '.join(walker_mut) + '\n')
                        info.append(walker_mut)
                    if walker_mut[0] == 'ins' and cur_gene == walker_mut[1] and abs(int(coord_otn) - int(walker_mut[2])) < 10:
                        print(walker_mut)
                        print(str(coord_otn) + ' ' + str(pos) + ' ' + str(special))

    return info

file = open('results_casali.txt', 'w')
mutation_info = [] # [name, [list of mutations]]

# from joblib import Parallel, delayed
#
# Parallel(n_jobs=100)(delayed(look_vcf_file)('/export/data/kkuleshov/myc/sra/' + name + '/' + name + '_h37rv.vcf',name, upstream_snps, file) for name in walker_names)

for name in casali_names[:-1]:
    mutation_info.append([name, look_vcf_file('/export/data/kkuleshov/myc/sra/' + name + '/' + name + '_h37rv.vcf.unf',name, upstream_snps, file)])

file.close()

# statistics

statistics_genes_snp = []
statistics_upstream_snp = []

for snp in database_genes:
    statistics_genes_snp.append([0, snp])
    for sample in mutation_info:
        for el in sample[1]:
            if el == snp:
                statistics_genes_snp[-1][0] += 1

for snp in upstream_snps:
    statistics_upstream_snp.append([0, snp])
    for sample in mutation_info:
        for el in sample[1]:
            if el == snp:
                statistics_genes_snp[-1][0] += 1

file = open('statistics_snp_casali.txt', 'w')

statistics_genes_snp = sorted(statistics_genes_snp, key = lambda snp: snp[0])[::-1]
statistics_upstream_snp = sorted(statistics_upstream_snp, key = lambda snp: snp[0])[::-1]

for snp in statistics_genes_snp:
    file.write(str(snp[0]) + ' ' + str(snp[1]) + '\n')
for snp in statistics_upstream_snp:
    file.write(str(snp[0]) + ' ' + str(snp[1]) + '\n')

file.close()