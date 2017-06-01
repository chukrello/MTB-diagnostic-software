# Pipeline with full precalculation (full file on full dictionary and full dataset)
# python pipeline_2.0.py *dict* *subset* *header*

import os, sys, argparse, subprocess

parser = argparse.ArgumentParser(description='Pipeline for generating prediction with precalculating FULL file.')

parser.add_argument('-d','--dictionary', help='Mutations dictionary.', required=True)
parser.add_argument('-s','--subset', help='Subset of samples for generating prediction.', default='../../1_input/subsets/Full_subset.txt')
parser.add_argument('-head','--header', help='Header for the filename.', required=True)
parser.add_argument('-md','--min_depth', help='Minimal read depth for predicting.', default='5')
parser.add_argument('-gp','--gene_percentage', help='Minimal read depth for predicting.', default='70')



args = parser.parse_args()

DICTIONARY = args.dictionary
SUBSET = args.subset
HEADER = args.header
MIN_DEPTH = args.min_depth
GENE_PERCENTAGE = args.gene_percentage



os.system('mkdir ../output/pipeline/'+HEADER)
print('Generating prediction')
os.system('python generate_predict_file.py MEGAFULL.txt ' + DICTIONARY + ' ' + SUBSET + ' ' + MIN_DEPTH + ' ' + GENE_PERCENTAGE + ' ../output/pipeline/' + HEADER + '/' + HEADER + '_predict_file.txt')
print('Calculating statistics')
os.system('python count_statistics.py ../output/pipeline/'+HEADER+'/'+ HEADER + '_predict_file.txt ../output/pipeline/' +HEADER+'/' + HEADER + '_out_stats.tsv')

