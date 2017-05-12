#!/export/home/kkuleshov/.pyenv/versions/2.7.11/bin/python2.7
#-*- coding: utf-8 -*-

import sys, os, argparse, subprocess
import csv
import logging
from collections import defaultdict

'''
Настройки подключения к базе данных
'''

parser = argparse.ArgumentParser(description='Frebayes launch and procces files')
parser.add_argument('-f','--bam_file', help='Path to file of bam format', required=True)
parser.add_argument('-s','--samtools_path', help='Path to samtools', default='/export/data/kchukreev/.samtools-1.3.1/samtools')

args = parser.parse_args()
bam_file_path = args.bam_file
samtools = args.samtools_path

depth_file_name = os.path.basename(bam_file_path).split('.')[0]
file_path_to_sample = os.path.dirname(bam_file_path)       
   
def main_():
    c = samtools + ' depth ' + bam_file_path + ' > ' + os.path.join(file_path_to_sample, depth_file_name +'.depth') 
    subprocess.call(c, shell=True)
    return

main_()

