# Full pipeline for predicting
# pipeline.py *dict* *subset*

import os, sys, argparse, subprocess

parser = argparse.ArgumentParser(description='Pipeline for generating prediction without precalculating FULL file.')

parser.add_argument('-d','--dictionary', help='Mutations dictionary.', required=True)
parser.add_argument('-s','--subset', help='Subset of samples for generating prediction.', default='../../1_input/subsets/Full_subset.txt')
parser.add_argument('-fh','--full_header', help='Header for the full file.', required=True)
parser.add_argument('-h','--header', help='Header for the filename.', required=True)

args = parser.parse_args()


DICTIONARY = args.dictionary
SUBSET = args.subset
DATA = args.full_header
HEADER = args.header


print('Building full file')
os.system('python build_full_file_samn_mutation_list.py ' + SUBSET + ' ' + DICTIONARY + ' ' + DATA + ' ' + HEADER + '_full_file_output.txt')
print('Counting frequency of mutations')
os.system('python count_frequency.py ' + HEADER + '_full_file_output.txt ' + HEADER + '_mutations_statistics.tsv')
print('Generating prediction')
os.system('python generate_predict_file.py MEGAFULL.txt ' + DICTIONARY + ' ' + SUBSET + ' ' +  HEADER + '_predict_file.txt')
print('Calculating statistics')
os.system('python count_statistics.py ' + HEADER + '_predict_file.txt ' + HEADER + '_out_stats.txt')
