
import sys, os, argparse, subprocess 
import csv
import logging
from collections import defaultdict

parser = argparse.ArgumentParser(description='Frebayes launch and procces files')
parser.add_argument('-f','--bam_file', help='Path to file of bam format', required=True)
parser.add_argument('-b','--freebayes', help='Path to freebayes tool', default='/export/home/kchukreev/soft/freebayes/bin/freebayes')
parser.add_argument('-r','--ref_fasta', help='Path to reference fasta file', default='/export/data/kchukreev/MTB-diagnostic-software/1_input/h37rv.fasta')
parser.add_argument('-o','--output_folder', help='Path to output folder to save files', default='out/')
parser.add_argument('-altf','--minAltFrac', help='minAltFrac', default='0.75')
parser.add_argument('-cov','--minCoverage', help='minCoverage', default='10')

args = parser.parse_args()
freebayes = args.freebayes
bam_file_path = args.bam_file
ref_fasta = args.ref_fasta
output_folder = args.output_folder
minAltFrac = args.minAltFrac
minCoverage = args.minCoverage

vcf_file_name = os.path.basename(bam_file_path).split('.')[0]
file_path_to_sample = os.path.dirname(bam_file_path)    
   
def main_():
    subprocess.call('mkdir processing', shell=True)
    os.chdir('./processing')
    
    c = (freebayes + ' --bam ' + bam_file_path +
                ' --vcf ' + vcf_file_name + '.vcf.unf' + 
                ' --fasta-reference ' + ref_fasta + 
                ' --pvar 0.0001' +
                ' --ploidy 1' + 
                #' --no-mnps' +
                ' --min-mapping-quality 0' + 
                ' --min-base-quality 20' +
                ' --min-alternate-fraction ' + minAltFrac +
                ' --min-coverage ' + minCoverage)

    print(c)

    subprocess.call(c, shell=True)
    c = '../filterVcf.pl ' + vcf_file_name + '.vcf.unf --noindels -d ' + minCoverage + ' -o ' + vcf_file_name + '.vcf -b ' + vcf_file_name + '.vcf.badsites.txt'
    subprocess.call(c, shell=True)
    c = 'mv ' + vcf_file_name + '.vcf.unf ' + file_path_to_sample
    subprocess.call(c, shell=True)
    c = 'mv ' + vcf_file_name + '.vcf ' + file_path_to_sample
    subprocess.call(c, shell=True)
    c = 'mv ' + vcf_file_name + '.vcf.badsites.txt ' + file_path_to_sample
    subprocess.call(c, shell=True)
    os.chdir('../')
    # sql = "UPDATE `sample_read` SET `Snp_call` = '%s' WHERE `sample_read`.`Path_to_sample` = '%s'" %('done', file_path_to_sample)
    # exec_sql(sql)
    return

main_()

