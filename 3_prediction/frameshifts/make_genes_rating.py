genomelist_all = [line[:-1] for line in open('/export/data/kchukreev/1_input/all_ids.txt').readlines()]
shifts = []
genes = []
genes_names = []

for name in genomelist_all:
    lines = []
    try:
        lines = [line[:-1] for line in open('/export/data/kchukreev/5_mutation_statistics/frameshifts/'+name+'.txt')]
    except:
        lines = []

    for i in range(len(lines)):
        if i % 2 == 1:
            gene = lines[i].split(' ')[-2]
            if gene not in genes_names:
                genes_names.append(gene)
                genes.append([gene, 1])
            else:
                for j in range(len(genes)):
                    if genes[j][0] == gene:
                        genes[j][1] += 1



genes = sorted(genes, key = lambda gene: gene[1])[::-1]

for gene in genes:
    print(gene[0] + ' ' + str(gene[1]))