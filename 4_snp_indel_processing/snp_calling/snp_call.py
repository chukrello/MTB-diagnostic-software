#!/export/home/kkuleshov/.pyenv/versions/2.7.11/bin/python2.7

import sys, os, argparse, subprocess 
from joblib import Parallel, delayed

parser = argparse.ArgumentParser(description='SNP calling with the help of Freebayes')

parser.add_argument('-l','--list', help='List of sample ids (only ids!)', required=True)
parser.add_argument('-p','--path', help='Path to the bams with format ID_h37rv.bam', required=True)
parser.add_argument('-j','--jobs', help='Number of jobs', default='1')

args = parser.parse_args()

list_ids = args.list
path_to_bams = args.path
jobs = int(args.jobs)

def freebayes_in_parrallel(id):

    print(id)

    bam = path_to_bams + id + '_h37rv.bam'

    c = 'python launch_freebayes.py -f ' + bam
    print(c)
    subprocess.call(c, shell=True)
    # c = 'python launch_depth.py -f ' + bam
    # print c
    # subprocess.call(c, shell=True)        
    
    return


ids = [line[:-1] for line in open(list_ids).readlines()]

Parallel(n_jobs=jobs)(delayed(freebayes_in_parrallel)(id) for id in ids)

