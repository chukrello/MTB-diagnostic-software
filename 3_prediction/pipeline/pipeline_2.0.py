# Pipeline with full precalculation (full file on full dictionary and full dataset)
# python pipeline_2.0.py *dict* *subset* *header*

import os
import sys

DICTIONARY = sys.argv[1]
SUBSET = sys.argv[2]
HEADER = sys.argv[3]
MIN_DEPTH = sys.argv[4]
GENE_PERCENTAGE = sys.argv[5]

os.system('mkdir ../output/pipeline/'+HEADER)
print('Generating prediction')
os.system('python generate_predict_file.py MEGAFULL.txt ' + DICTIONARY + ' ' + SUBSET + ' ' + MIN_DEPTH + ' ' + GENE_PERCENTAGE + ' ../output/pipeline/'+HEADER+'/'+ HEADER + '_predict_file.txt')
print('Calculating statistics')
os.system('python count_statistics.py ../output/pipeline/'+HEADER+'/'+ HEADER + '_predict_file.txt ../output/pipeline/' +HEADER+'/' + HEADER + '_out_stats.tsv')

