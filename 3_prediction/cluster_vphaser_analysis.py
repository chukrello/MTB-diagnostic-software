import os
from Bio import SeqIO
from Bio.Seq import Seq

# IMPORTANT SCRIPT
# Script gets all the VCFs from samples and compares mutations from VCF and from Walker dictionaries
# Works good only with SNP

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
all_names = [line[:-1] for line in open(PATH_input+'all_ids.txt').readlines()]
walker_names = [line[:-1] for line in open(PATH_input+'walker_ids.txt').readlines()]
nonwalker_names = [name for name in all_names if name not in walker_names]


def make_upstream_snps(database_mutations, genes):
    upstream_snps = []
    for mutation in database_mutations:
        if '-' in mutation[2] and mutation[0] == 'snp':
            
            #nucleotide shift, nucl < 0
            nucl = int(mutation[2])
            pos = 0
            alt = ''

            if genes[mutation[1]][2] == '+':
                pos = genes[mutation[1]][0] + nucl
                alt = mutation[4] 
            else:
                pos = genes[mutation[1]][1] - nucl
                alt = Seq(mutation[4]).reverse_complement().tostring()

            upstream_snps.append([pos, alt, mutation])

    return upstream_snps
upstream_snps = make_upstream_snps(database_mutations, genes)

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

def look_vcf_file(filename, name, upstream_snp, file):
    #for SNP
    print('Looking at ' + name)
    info = []
    file.write(name+'\n')
    # open .vcf file
    raw_data_2 = [[int(line.split('\t')[1]), line.split('\t')[3], line.split('\t')[4]] for line in open(filename[:-4]).readlines() if line[0] != '#' and float(line.split('\t')[5]) >= 600]
    data = []
    for el in raw_data_2:
        if el not in data:
            data.append(el)

    for i in range(len(data)):

        pos = data[i][0]
        ref = data[i][1]
        alt = data[i][2]

        #UPSTREAM SNP ANALYSIS
        for snp in upstream_snp:
            if snp[0] == pos and snp[1] == alt:
                # print(snp[2])
                file.write(' '.join(snp[2]) + '\n')
                info.append(snp[2])    


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

                    for walker_mut in database_mutations:
                        if walker_mut[0] == 'snp':
                            if walker_mut[1] == cur_gene and walker_mut[2][0] != '-' and int(walker_mut[2]) == protein_pos:
                                
                                if old_aminoacid != walker_mut[3]:
                                    print(walker_mut)
                                    print(data[i])

                                if (old_aminoacid == walker_mut[3]) and (new_aminoacid == walker_mut[4]):
                                    file.write(' '.join(walker_mut) + '\n')
                                    info.append(walker_mut)

                else: # if strand is '-'
                    if cur_gene == 'rrs':
                        walker_coord_start = genes[cur_gene][1] - (pos-1)
                        reference_nucl = Seq(sequence[walker_coord_start-1]).reverse_complement()

                        for walker_mut in database_mutations:
                            if walker_mut[0] == 'snp':
                                if walker_mut[0] == 'snp' and walker_mut[1] == 'rrs':
                                    if int(walker_mut[2]) == int(walker_coord_start) and str(Seq(alt).reverse_complement()) == walker_mut[4]:
                                        file.write(' '.join(walker_mut) + '\n')
                                        info.append(walker_mut)
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


                        for walker_mut in database_mutations:
                            if walker_mut[0] == 'snp':
                                if walker_mut[1] == cur_gene and walker_mut[2][0] != '-' and int(walker_mut[2]) == protein_pos:
                                    if old_aminoacid != walker_mut[3]:
                                        print(old_aminoacid)
                                        print(walker_mut)
                                        print(data[i])

                                    if old_aminoacid == walker_mut[3] and new_aminoacid == walker_mut[4]:
                                        file.write(' '.join(walker_mut) + '\n')
                                        info.append(walker_mut)
        elif len(ref) > len(alt):
            #DELETION ANALYSIS
            delta = len(ref) - len(alt)
            special = ref[1:delta+1]
            pos = pos + 1

            cur_gene = ''

            for key in genes:
                if int(genes[key][0]) <= int(pos) and int(pos) <= int(genes[key][1]):
                    cur_gene = key
                    break

            if cur_gene != '':

                if genes[cur_gene][2] == '+':
                    coord_otn = pos - genes[cur_gene][0] + 1
                else:
                    coord_otn = genes[cur_gene][1] - pos + 1

                for walker_mut in database_mutations:
                    if walker_mut[0] == 'del' and cur_gene == walker_mut[1] and coord_otn == int(walker_mut[2]):
                        # print('WOW')
                        # print(walker_mut)
                        # print(str(coord_otn) + ' ' + str(pos) + ' ' + str(special))
                        file.write(' '.join(walker_mut) + '\n')
                        info.append(walker_mut)
                    # if walker_mut[0] == 'del' and cur_gene == walker_mut[1] and abs(int(coord_otn) - int(walker_mut[2])) < 10:
                        # print(walker_mut)
                        # print(str(coord_otn) + ' ' + str(pos) + ' ' + str(special))

        elif len(alt) > len(ref):
            #INSERTION ANALYSIS
            delta = len(alt) - len(ref)
            special = alt[1:delta+1]
            pos = pos + 1

            cur_gene = ''

            for key in genes:
                if int(genes[key][0]) <= int(pos) and int(pos) <= int(genes[key][1]):
                    cur_gene = key
                    break

            if cur_gene != '':

                if genes[cur_gene][2] == '+':
                    coord_otn = pos - genes[cur_gene][0] + 1
                else:
                    coord_otn = genes[cur_gene][1] - pos + 1

                for walker_mut in database_mutations:
                    if walker_mut[0] == 'ins' and cur_gene == walker_mut[1] and coord_otn == int(walker_mut[2]):
                        # print('WOW')
                        # print(walker_mut)
                        # print(str(coord_otn) + ' ' + str(pos) + ' ' + str(special))
                        file.write(' '.join(walker_mut) + '\n')
                        info.append(walker_mut)
                    # if walker_mut[0] == 'ins' and cur_gene == walker_mut[1] and abs(int(coord_otn) - int(walker_mut[2])) < 10:
                        # print(walker_mut)
                        # print(str(coord_otn) + ' ' + str(pos) + ' ' + str(special))

    return info

def look_vphaser_file(filename, name, upstream_snp, file):
    #for SNP
    print('Looking at ' + name)
    info = []
    file.write(name+'\n')
    # open vphaser file
    try:
        raw_data_2 = [[int(line.split('\t')[0]), line.split('\t')[1], line.split('\t')[2][:-1]] for line in open(filename).readlines()]
    except Exception:
        return 0
    data = []
    for el in raw_data_2:
        if el not in data:
            data.append(el)

    for i in range(len(data)):

        pos = data[i][0]

        if sequence[pos-1] == data[i][2]:
            ref = data[i][2]
            alt = data[i][1]
        elif sequence[pos-1] == data[i][1]:
            ref = data[i][1]
            alt = data[i][2]
        else:
            continue

        #UPSTREAM SNP ANALYSIS
        for snp in upstream_snp:
            if int(snp[0]) == int(pos) and snp[1] == alt:
                # print(snp[2])
                file.write(' '.join(snp[2]) + '\n')
                print(snp)
                info.append(snp[2])    


        cur_gene = ''

        if len(ref) == len(alt):
            #SNP ANALYSIS

            for key in genes:
                if int(genes[key][0]) <= pos and pos <= int(genes[key][1]):
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

                    for walker_mut in database_mutations:
                        if walker_mut[0] == 'snp':
                            if walker_mut[1] == cur_gene and walker_mut[2][0] != '-' and int(walker_mut[2]) == protein_pos:
                                
                                if old_aminoacid != walker_mut[3]:
                                    print(walker_mut)
                                    print(data[i])

                                if (old_aminoacid == walker_mut[3]) and (new_aminoacid == walker_mut[4]):
                                    file.write(' '.join(walker_mut) + '\n')
                                    print(walker_mut)
                                    info.append(walker_mut)

                else: # if strand is '-'
                    if cur_gene == 'rrs':
                        walker_coord_start = genes[cur_gene][1] - (pos-1)
                        reference_nucl = Seq(sequence[walker_coord_start-1]).reverse_complement()

                        for walker_mut in database_mutations:
                            if walker_mut[0] == 'snp':
                                if walker_mut[0] == 'snp' and walker_mut[1] == 'rrs':
                                    if int(walker_mut[2]) == int(walker_coord_start) and str(Seq(alt).reverse_complement()) == walker_mut[4]:
                                        file.write(' '.join(walker_mut) + '\n')
                                        print(walker_mut)
                                        info.append(walker_mut)
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


                        for walker_mut in database_mutations:
                            if walker_mut[0] == 'snp':
                                if walker_mut[1] == cur_gene and walker_mut[2][0] != '-' and int(walker_mut[2]) == protein_pos:
                                    if old_aminoacid != walker_mut[3]:
                                        print(old_aminoacid)
                                        print(walker_mut)
                                        print(data[i])

                                    if old_aminoacid == walker_mut[3] and new_aminoacid == walker_mut[4]:
                                        file.write(' '.join(walker_mut) + '\n')
                                        print(walker_mut)
                                        info.append(walker_mut)
    return info

file = open('results_vphaser_nonwalker.txt', 'w')
mutation_info = [] # [name, [list of mutations]]
for name in nonwalker_names[:-1]:
    mutation_info.append([name, look_vphaser_file('/export/data/kchukreev/9_minor_variants/vcfs/' + name + '_vphaser.vcf',name, upstream_snps, file)])

file.close()

print('DONE')


#STATISTICS ANALYSIS
statistics_genes_snp = []

for snp in database_mutations:
    statistics_genes_snp.append([0, snp])
    for sample in mutation_info:
        for el in sample[1:]:
            if el == snp:
                statistics_genes_snp[-1][0] += 1

file = open('statistics_snp_vphaser_nonwalker.txt', 'w')

statistics_genes_snp = sorted(statistics_genes_snp, key = lambda snp: snp[0])[::-1]

for snp in statistics_genes_snp:
    file.write(str(snp[0]) + ' ' + str(snp[1]) + '\n')

file.close()