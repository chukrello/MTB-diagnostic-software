def BinSearchVirt(depth, x):
    i = 0
    j = len(depth)-1
    while i < j:
        m = int((i+j)/2)
        if x > depth[m][0]:
            i = m+1
        else:
            j = m
    if depth[j][0] == x:
        return depth[j][1]
    else:
        return None

# trace
# gff = [[int(line.split('\t')[3]),int(line.split('\t')[4])] for line in open('/export/data/kchukreev/1_input/AL123456_rev.gff').readlines()]
# open .detph file

#new trace
list_genes = ['ahpC', 'eis', 'embA', 'embB', 'embC', 'embR', 'fabG1', 'gid', 'gyrA', 'gyrB', 'inhA', 'iniA', 'iniC',
             'katG', 'manB', 'ndh', 'pncA', 'rmlD', 'rpoB', 'rpsA', 'rpsL', 'rrs', 'tlyA']

gff_all = []

for line in open('/export/data/kchukreev/1_input/AL123456_rev.gff').readlines():
    for element in list_genes:
        if 'gene ' + element in line:
            gff_all.append([int(line.split('\t')[3]),int(line.split('\t')[4]), element])


gff = []

for el in gff_all:
    if el not in gff:
        gff.append(el)





import matplotlib.pyplot as plt
import numpy as np
import threading


def calculate_percentage(genome):

    PATH = '/export/data/kkuleshov/myc/sra/'
    more5 = 0
    common_length = 0
    depth = [[int(line.split('\t')[1]), int(line.split('\t')[2])] for line in open(PATH+genome+'/'+genome+'_h37rv.depth').readlines()]
    # print('Genome ' + genome + ' is processing')
    for el in gff:

        more5 = 0
        common_length = 0

        # searching for index
        for i in range(el[0], el[1]+1):
            cov = BinSearchVirt(depth,i)
            if cov != None and cov > 5:
                more5 += 1
            common_length += 1

        percentage = (more5*1.0)/common_length
        if percentage < 0.88:
            print(genome + ' ' + str(el[2]))
        else:
            print(genome + ' OK ' + str(percentage))



def calculate_average_coverage(genome):
    if genome!='':
        PATH = '/export/data/kkuleshov/myc/sra/'
        depthfile = open(PATH+genome+'/'+genome+'_h37rv.depth')
        cov = 0

        line = depthfile.readline()
        lines = 0
        while line != '':
            cov += int(line.split('\t')[2][:-1])
            lines += 1
            line = depthfile.readline()
        # global cnt
        # cnt += 1
        # print(cnt)

        file = open('/export/data/kchukreev/4_coverage/coverages/' + genome + '.txt', 'w')
        file.write(str(cov*1.0/lines)+'\n')
        file.close()

        # lock = threading.Lock()
        # lock.acquire()
        # file.write(str(cov*1.0/lines)+'\n')
        # lock.release()

genome_list = [line[:-1] for line in open('/export/data/kchukreev/1_input/walker_ids.txt').readlines()]

# for genome in genome_list:
#     calculate_percentage(genome)

from joblib import Parallel, delayed

cnt = 0
file = open('coverage_av.txt', 'w')

# Parallel(n_jobs=140)(delayed(calculate_percentage)(genome) for genome in genome_list)
Parallel(n_jobs=140)(delayed(calculate_percentage)(genome) for genome in genome_list)

# errors = [line[:-1].split(' ') for line in open('/home/chukreev/PycharmProjects/tuberculosis/errors.txt').readlines()]





