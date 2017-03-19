from joblib import Parallel, delayed
import re
PATH_input = '/export/data/kchukreev/1_input/'
lines = open(PATH_input+'AL123456_rev.gff').readlines()

gene_reg = re.compile('\tgene \w*')

genes = []
for line in lines:
    m = gene_reg.search(line)
    if m:
        genes.append(line)
genes = [[int(line[:-1].split('\t')[3]),int(line[:-1].split('\t')[4]), line[:-1].split('\t')[8].split(' ')[1], line[:-1].split('\t')[6]] for line in genes]

genomelist_all = [line[:-1] for line in open('/export/data/kchukreev/1_input/all_ids.txt').readlines()]

def FindFrameShift(name, genes):
    # open vcf
    try:
        data = [[int(line.split('\t')[1]), line.split('\t')[3], line.split('\t')[4]] for line in open('/export/data/kkuleshov/myc/sra/'+name+'/'+name+'_h37rv.vcf.unf').readlines() if line[0] != '#']
    except Exception:
        return 0

    file = open('/export/data/kchukreev/5_mutation_statistics/frameshifts/'+name+'.txt', 'w')

    for mut in data:
        if len(mut[1]) != len(mut[2]) and len(mut[1]) - len(mut[2]) %3 != 0:
            for gene in genes:
                if gene[0] < mut[0] < gene[1]:
                    #frameshift!
                    mut[0] = str(mut[0])
                    file.write(' '.join(mut) + '\n')
                    gene[0] = str(gene[0])
                    gene[1] = str(gene[1])
                    file.write(' '.join(gene) + '\n')

    file.close()

Parallel(n_jobs=140)(delayed(FindFrameShift)(name, genes) for name in genomelist_all)