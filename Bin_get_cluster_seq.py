#!/usr/bin/env python
# jinh 20200316
'''
This is a downstream script for Bin_cluster.py to obtain the fasta of each culster.
Here we set 200kbp as a min length of a cluster to outfile. The short cluster will write to unbinned.fasta 
'''
import os, sys
from Bio import SeqIO
import shutil
from collections import defaultdict
import argparse

# if len(sys.argv) != 5:
#     print "Usage: Bin_get_cluster_seq.py seq.fasta cluster_file outfolder prefix"
#     sys.exit()

parser = argparse.ArgumentParser(description='This is a downstream script for Bin_cluster.py to obtain the fasta of each culster. \n Here we set 200kbp as a min length of a cluster to outfile. The short cluster will write to unbinned.fasta ')
parser.add_argument("-i", '--input', required=True, help='input contigs fasta')
parser.add_argument("-b", '--bin_file', required=True, help='scaffold to bin files')
parser.add_argument("-o", '--outdir', required=True, help='output files directory')
parser.add_argument("-p", '--prefix', required=True, help='output file prefix')
parser.add_argument("-f", '--force', action = 'store_true', help='overwriting the existing files')

args = parser.parse_args()

input = args.input
cluster_file = open(args.bin_file, "r")
cluter_fa_outfolder = args.outdir
prefix = args.prefix

if not os.path.exists(cluter_fa_outfolder):
    os.makedirs(cluter_fa_outfolder)
else:
    if args.force == True:
        shutil.rmtree(cluter_fa_outfolder)
        os.makedirs(cluter_fa_outfolder)
    else:
        print("Outdir already exist, please rename or use --force")
        sys.exit()

cls_dic = {}
for line in cluster_file:
    line = line.strip().split("\t")
    cls_dic[line[1]] = line[0]

cls_out_dic = defaultdict(list)
cls_length_dic = defaultdict(int)
for seq in SeqIO.parse(input, "fasta"):
    if seq.description in cls_dic:
        cls_name = cls_dic[seq.description]
        cls_out_dic[cls_name].append(">" + seq.description + "\n" + str(seq.seq) + "\n")
        cls_length_dic[cls_name] += len(seq.seq)


out_ubinned = cluter_fa_outfolder + '/' + prefix + '.unbinned.fa'
out_ubinned_wt = open(out_ubinned, 'w')
for k in cls_out_dic:
    out_format = ''.join(cls_out_dic[k])
    if cls_length_dic[k] > 200000:
        out_cls = cluter_fa_outfolder + '/' + prefix + '.' + k + '.fa'
        if not os.path.exists(out_cls) :
            out_cls_wt = open(out_cls, 'w')
        out_cls_wt.write(out_format)
    else :
        out_ubinned_wt.write(out_format)

out_ubinned_wt.close()


    



