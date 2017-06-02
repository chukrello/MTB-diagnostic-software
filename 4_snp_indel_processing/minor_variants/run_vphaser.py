import sys, os, argparse, subprocess 
from joblib import Parallel, delayed

parser = argparse.ArgumentParser(description='Runs VPhaser for ids')

parser.add_argument('-l','--list', help='List of sample ids (only ids!)', required=True)
parser.add_argument('-p','--path', help='Path to the bams with format ID_h37rv.bam', required=True)
parser.add_argument('-j','--jobs', help='Number of jobs', default='1')

args = parser.parse_args()

list = args.list
path_to_bams = args.path
jobs = int(args.jobs)

def run_vphaser(id):
	comand_mkdir = 'mkdir SAMNs_results/'+id

	comand_copy = 'cp ' + path_to_bams + '/' + id + '_h37rv.bam SAMNs_results/' + id + '/'

	comand_vphaser = 'OMP_NUM_THREADS=15 /export/home/kchukreev/soft/VPhaser-2-02112013/bin/variant_caller -a 0.01 -i SAMNs_results/' + id + '/' + id + '_h37rv.bam -o SAMNs_results/' + id

	comand_delete_region = 'rm SAMNs_results/' + id + '/*.region'
	
	comand_delete_bam = 'rm SAMNs_results/' + id + '/*.bam'

	try:
		os.system(comand_mkdir)
		os.system(comand_copy)
		os.system(comand_vphaser)
		os.system(comand_delete_region)
		os.system(comand_delete_bam)
	except Exception:
		a = 1



ids = [line[:-1] for line in open(list).readlines()]
Parallel(n_jobs=jobs)(delayed(run_vphaser)(id) for id in ids)
