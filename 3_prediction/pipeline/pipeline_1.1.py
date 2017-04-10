# Full pipeline for predicting
# pipeline.py *dict* *subset*

import os
import sys

DICTIONARY = sys.argv[1]
SUBSET = sys.argv[2]
HEADER = sys.argv[3]

print('Building full file')
os.system('python build_full_file_samn_mutation_list.py ' + SUBSET + ' ' + DICTIONARY + ' ' + HEADER + '_full_file_output.txt')
print('Counting frequency of mutations')
os.system('python count_frequency.py ' + HEADER + '_full_file_output.txt ' + HEADER + '_mutations_statistics.tsv')
print('Generating prediction')
os.system('python generate_predict_file.py MEGAFULL.txt ' + DICTIONARY + ' ' + SUBSET + ' ' +  HEADER + '_predict_file.txt')
print('Calculating statistics')
os.system('python count_statistics.py ' + HEADER + '_predict_file.txt ' + HEADER + '_out_stats.txt')
