# open .detph file
import matplotlib.pyplot as plt
import numpy as np
import threading
import os

list_genes = ['ahpC', 'eis', 'embA', 'embB', 'embC', 'embR', 'fabG1', 'gid', 'gyrA', 'gyrB', 'inhA', 'iniA', 'iniC',
             'katG', 'nat' 'ndh', 'pncA', 'rpoB', 'rpsA', 'rpsL', 'rrs', 'tlyA']

# trace
def calculate_trace(list_genes):

	lines = [[int(line.split('\t')[3]),int(line.split('\t')[4]), line] for line in open('../1_input/AL123456_rev.gff').readlines()]

	trace = []

	for gene in list_genes:
        for line in lines:
			if gene in line[2]:
				trace.append([line[0], line[1]])
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

def extract_coverage(id, trace, id):
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


all_ids = [line[:-1] for line in open('../1_input/susbsets/Full_subset.txt').readlines()]
trace = calculate_trace(list_genes)
os.system('mkdir coverage_output')

from joblib import Parallel, delayed
Parallel(n_jobs=140)(delayed(calculate_average_coverage)(genome) for genome in genome_list)
