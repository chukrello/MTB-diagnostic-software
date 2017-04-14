# open .detph file and cut coordinates of the genes we nees
#


import matplotlib.pyplot as plt
import numpy as np
import threading
import os

list_genes = ['ahpC', 'eis', 'embA', 'embB', 'embC', 'embR', 'fabG1', 'gid', 'gyrA', 'gyrB', 'inhA', 'iniA', 'iniC',
             'katG', 'nat' 'ndh', 'pncA', 'rpoB', 'rpsA', 'rpsL', 'rrs', 'tlyA']
all_ids = [line[:-1] for line in open('../1_input/subsets/Full_subset.txt').readlines()]
os.system('mkdir coverage_output')


def get_genes_info():
    import re
    lines = open('../1_input/AL123456_rev.gff').readlines()
    gene_reg = re.compile('\tgene \w*')
    genes = []
    for line in lines:
        m = gene_reg.search(line)
        if m:
            genes.append(line)
    genes = [[line[:-1].split('\t')[3],line[:-1].split('\t')[4], line[:-1].split('\t')[8].split(' ')[1], line[:-1].split('\t')[6]] for line in genes]
    list_genes = ['ahpC', 'eis', 'embA', 'embB', 'embC', 'embR', 'fabG1', 'gid', 'gyrA', 'gyrB', 'inhA', 'iniA', 'iniC',
             'katG', 'manB', 'ndh', 'pncA', 'rmlD', 'rpoB', 'rpsA', 'rpsL', 'rrs', 'tlyA']
    def get_sequence(genes, gene):
        for element in genes:
            if element[2] == gene:
                return [int(element[0]), int(element[1]), element[3]]
    genes_seqs = dict()
    for gene in list_genes:
        genes_seqs[gene] = get_sequence(genes,gene)

    return genes_seqs

# trace
def calculate_trace(list_genes, genes):

    lines = [[int(line.split('\t')[3]),int(line.split('\t')[4]), line] for line in open('../1_input/AL123456_rev.gff').readlines()]

    trace = []

    for gene in list_genes:
        for line in lines:
            if gene in line[2]:
                if genes[gene][2] == '+':
                    trace.append([line[0]-150, line[1]])
                    break
                else:
                    trace.append([line[0], line[1] + 150])
                    break

    return trace

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

def extract_coverage(id, trace):
    print(id)
    PATH = '/export/data/kkuleshov/myc/sra/'
    depth = [[int(line.split('\t')[1]), int(line.split('\t')[2])] for line in open(PATH+ id +'/'+ id +'_h37rv.depth').readlines()]

    extracted_coverage = []

    for track in trace:
        for i in range(track[0], track[1]+1):
            cov = BinSearchVirt(depth,i)
            if cov != None:
                extracted_coverage.append([i, cov])
            else:
                extracted_coverage.append([i, 0])

    file = open('coverage_output/' + id + '.txt', 'w')
    for el in extracted_coverage:
        file.write(str(el[0]) + ' ' + str(el[1]) + '\n')

    file.close()

genes = get_genes_info()
trace = calculate_trace(list_genes, genes)
from joblib import Parallel, delayed
Parallel(n_jobs=-1)(delayed(extract_coverage)(id, trace) for id in all_ids)
