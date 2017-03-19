genomelist_walker = [line[:-1] for line in open('/export/data/kchukreev/1_input/walker_ids.txt').readlines()]
genomelist_all = [line[:-1] for line in open('/export/data/kchukreev/1_input/all_ids.txt').readlines()]
genomelist_nonwalker = [el for el in genomelist_all if el not in genomelist_walker]

def ShowIndels(name, insertion, deletion):

    def WriteIndel(pos,ref,alt, file):
        file.write(name + '\t' + pos + '\t' + ref + '\t' + alt + '\n')

     # open vcf file
    try:
        data = [[line.split('\t')[1], line.split('\t')[3], line.split('\t')[4]] for line in open('/export/data/kkuleshov/myc/sra/' + name + '/' + name + '_h37rv.vcf.unf').readlines() if line[0] != '#']
    except Exception:
        return ['']

    for mutation in data:
        pos = mutation[0]
        ref = mutation[1]
        alt = mutation[2]
        if len(ref) != len(alt):
            if len(ref) > len(alt):
                WriteIndel(pos, ref, alt, deletion)
            else:
                WriteIndel(pos, ref, alt, insertion)



file = open('indels.txt', 'w')
insertion = open('ins.txt', 'w')
deletion = open('dels.txt', 'w')

for name in genomelist_nonwalker:
    ShowIndels(name, insertion, deletion)

file.close()