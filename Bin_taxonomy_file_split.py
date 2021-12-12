#!/usr/bin/env python

import os
import sys
import argparse
import shutil

parser = argparse.ArgumentParser(description='Split faa')
parser.add_argument("-i", '--input', required=True, help='input file')
parser.add_argument("-f", '--f2s', required=True, help='fasta2sample file')
parser.add_argument("-o", '--out_path', required=True, help='path of output file')

args = parser.parse_args()

input = args.input
f2s = args.f2s
out_path = args.out_path

if not os.path.exists(out_path):
    os.makedirs(out_path)

binfile_read = open(f2s, 'r')
bin_dict  = {}
for line in binfile_read:
    line_sp = line.strip().split("\t")
    bin_dict[line_sp[0]] = line_sp[1]

infile = open(input)

infile = open(input)
for line in infile:
    line_sp = line.strip().split("\t")
    seq_name = line_sp[0]
    outfile = out_path + "/" + bin_dict[seq_name] + '_con.tax'
    outfile_wt = open(outfile, 'a')
    outfile_wt.write(line)