#!/usr/bin/env python 
#jinh 20200315
import os,sys
from itertools import islice
import argparse
from itertools import islice
from scipy.stats.stats import pearsonr
from Bio import SeqIO
import shutil
import time
from collections import defaultdict
import pandas as pd
from itertools import islice


parser = argparse.ArgumentParser(description='Metagenome Binning based on Abundance, Tetranucleotide frequency, GC content, Taxnoomy and Single copy genes')
parser.add_argument("-a", '--fasta', required=True, help='path to contigs fasta file')
parser.add_argument("-j", '--jgi', required=True, help='path to output of jgi_summarize_bam_contig_depths')
parser.add_argument("-s", '--sgc', required=True, help='path to output of SGC hummer out')
parser.add_argument("-t", '--tax', required=True, help='path to output of contig taxonomic files')
parser.add_argument("-k", '--kmer', required=True, help='path to output of contig kmer freuence files')
parser.add_argument("-o", '--outdir', required=True, help='output directory')
parser.add_argument("-l", "--min_length", required=False, default=2000, help="Minimum size of a contig for binning (should be >=1500)")
parser.add_argument("-f", '--force', action = 'store_true', help='overwriting the existing files')

args = parser.parse_args()

# scg93_db_file=open("/nvmessdnode3/opt/database/123_SCG_hmm_db/bacteria_123_CSCG","r")
scg93_db_file=open("/mnt1/database/SGC_db/bacteria_123_CSCG","r")

fasta = args.fasta
jgi = open(args.jgi, "r")
tax = open(args.tax, "r")
kmer = args.kmer
sgc = args.sgc
outdir = args.outdir
fl = args.min_length
tnffile = open(kmer, "r")
scgfile = open(sgc, "r")

#id means identity, c_tax is a array, including tax, split with;
def compare_tax(c1_list, c2_list) :
    tax_score=0
    if c1_list[5] != "None" and  "unclassified" not in c1_list[5] and c1_list[5] == c2_list[5]:
            tax_score+= 6
    elif c1_list[4] != "None" and  "unclassified" not in c1_list[4] and c1_list[4] == c2_list[4]:
            tax_score+= 5
    elif c1_list[3] != "None" and  "unclassified" not in c1_list[3] and c1_list[3] == c2_list[3]:
            tax_score+= 4
    elif c1_list[2] != "None" and  "unclassified" not in c1_list[2] and c1_list[2] == c2_list[2]:
            tax_score+= 3
    elif c1_list[1] != "None" and  "unclassified" not in c1_list[1] and c1_list[1] == c2_list[1]:
            tax_score+= 2
    elif c1_list[0] != "None" and  "unclassified" not in c1_list[0] and c1_list[0] == c2_list[0]:
            tax_score+= 1
    else : 
        tax_score=0
    return tax_score

#g1 means cluster GC, g2 means new scaffold GC.
def compare_gc(g1,g2):
    gc_result=""
    g1=float(g1)
    g2=float(g2)
    if g2 >= g1*0.9 or g2 <= g1*1.1 : 
        gc_result="T"
    else: 
        gc_result="F"
    # print g1,g2
    return gc_result

#d1 cluster depth d2 scaffold depth
def compare_depth(d1,d2):
    depth_score=0
    d1=float(d1)
    d2=float(d2)
    if d1 <= 100 and d2 <= 100 :
        if abs(d1-d2) <= 10:
            depth_score=4
        elif abs(d1-d2) <= 20:
            depth_score=3
        elif abs(d1-d2) <= 30:
            depth_score=2
        elif abs(d1-d2) <= 40:
            depth_score=1
        else:
            depth_score=0
    else :
        if d2 >= d1*0.8 and d2 <= d1*1.2:
            depth_score=4
        elif d2 >= d1*0.75 and d2 <= d1*1.25:
            depth_score=3
        elif d2 >= d1*0.7 and d2 <= d1*1.3:
            depth_score=2
        elif d2 >= d1*0.65 and d2 <= d1*1.35:
            depth_score=1
        else:
            depth_score=0
    return depth_score

def compare_tnf(t1,t2):
    tnf_score=0
    tnf_corr=pearsonr(list(map(float, t1)), list(map(float, t2)))[0]
    if tnf_corr >= 0.98:
        tnf_score=6
    elif tnf_corr >= 0.97:
        tnf_score=5
    elif tnf_corr >= 0.95:
        tnf_score=4
    elif tnf_corr >= 0.93:
        tnf_score=3
    elif tnf_corr >= 0.91:
        tnf_score=2
    elif tnf_corr >= 0.7:
        tnf_score=1
    else: 
        tnf_score=0
    return tnf_score

def judge_dia_tax(total_hits, c_tax):
    c_list = []
    if c_tax[0] != 'nan' and float(c_tax[1])/float(total_hits) > 0.5 :  # phylum
        p = c_tax[0]
    else: 
        p = "None"
    if c_tax[3] != 'nan' and float(c_tax[4])/float(total_hits) > 0.5 :  # class
        c = c_tax[3]
    else: 
        c = "None"
    if c_tax[6] != 'nan' and float(c_tax[7])/float(total_hits) > 0.5 : # order
        o = c_tax[6]
    else: 
        o = "None"
    if c_tax[9] != 'nan' and float(c_tax[10])/float(total_hits) > 0.5 and float(c_tax[11]) >= 60 : # family
        f = c_tax[9]
    else: 
        f = "None"
    if c_tax[12] != 'nan' and float(c_tax[13])/float(total_hits) > 0.5 and float(c_tax[14]) >= 70 : # genus
        g = c_tax[12]
    else: 
        g = "None"
    if c_tax[15] != 'nan' and float(c_tax[16])/float(total_hits) > 0.6 and float(c_tax[17]) >= 90 : # species
        s = c_tax[15]
    else: 
        s = "None"
    return [p, c, o, f, g, s]

# prepare dict
# 1. get GC coutent and Length of contigs
scaffold = defaultdict(dict)
for rec in SeqIO.parse(fasta, 'fasta'):
    GC = (rec.seq.count("C") + rec.seq.count("G")) / float(len(rec.seq))
    Seq_len = len(rec.seq)
    scaffold[rec.description]["length"] = Seq_len
    scaffold[rec.description]["GC"] = GC

# 2. get contigs depth
for line in islice(jgi, 1 , None):
    lines = line.strip().split("\t")
    seq_name = lines[0]
    depth = lines[2]
    scaffold[seq_name]["depth"] = depth
    
# 3. get contig tax
for line in tax:
    lines = line.strip().split("\t")
    seq_name = lines[0]
    # print (line)
    scaffold[seq_name]["tax"] = judge_dia_tax(lines[3], lines[7:])
    # print (scaffold[seq_name]["tax"]) 

# 4. make tnf dict
for tnf in islice(tnffile, 1, None):
    tnf_nm=tnf.strip().split("\t")[0]
    tnf_tnf=tnf.strip().split("\t")[1:]
    if tnf_nm in scaffold:
        scaffold[tnf_nm]["TNF"]=tnf_tnf
    else:
        print ("Miss scaffold name, exit. (tnf)")
        sys.exit()

# 5. scg file
# SCG file ; if no hmmer out, use 0 instead
db_dic={}
scg_list=[]
scg_dic = defaultdict(defaultdict)

for line in scg93_db_file : 
    line_splt=line.strip().split("\t")
    db_dic[line_splt[0]]=line_splt[1]
    scg_list.append(line_splt[0])

for scaf in scaffold.keys():
    for sg in scg_list:
        scg_dic[scaf][sg]=0

for line_scg in scgfile:
    if not line_scg.startswith("#") :
        scg_line_splt=line_scg.strip().split(" ")
        while "" in scg_line_splt:
            scg_line_splt.remove("")
        scf_name="_".join(scg_line_splt[2].split("_")[0:-1])
        scg_name=scg_line_splt[1].split(".")[0]
        if float(scg_line_splt[5]) >= float(db_dic[scg_name]) and scf_name in scg_dic:
            scg_dic[scf_name][scg_name]+=1

# 6. filter scaffold by length
scaffold_fitered={}
for kf in scaffold:
    if float(scaffold[kf]["length"]) >= int(fl) :
        scaffold_fitered[kf]=scaffold[kf]
print ("%s > %s " %(len(scaffold_fitered.keys()),fl))

#s1 and s2 are dict
def plus_scg(s1,s2):
    tmp={}
    for k in s1:
        tmp[k]=float(s1[k])+float(s2[k])
    return tmp

def summerize_scg(tmp):
    mgc_count=0
    for i in tmp:
        if float(tmp[i]) > 1 : 
            mgc_count+=1
    return mgc_count

print("----------------------------------------Clustering--------------------------------------")
def cluster_tnf_summarize(cluster_tnf, cluster_len, scaffold_tnf, sacfd_len):
    cluster_tnf_summarize_out=[]
    for i in range(0,len(cluster_tnf)):
        tnf_cls=float(cluster_tnf[i])
        len_cls=float(cluster_len)
        tnf_scf=float(scaffold_tnf[i])
        len_scf=float(sacfd_len)
        ave=(tnf_cls*len_cls+tnf_scf*len_scf)/(len_cls+len_scf)
        cluster_tnf_summarize_out.append(ave)
    return cluster_tnf_summarize_out

def score_result(tax_s,tnf_s,depth_s):
    if tax_s >= 6 and depth_s >=1 and tnf_s>=1:
        result="T"
    elif tax_s >= 5 and depth_s >=2 and tnf_s>=2:
        result="T"
    elif tax_s >= 4 and depth_s >=2 and tnf_s>=3:
        result="T"
    elif tax_s >= 3 and depth_s >=2 and tnf_s>=4:
        result="T"
    elif tax_s >= 2 and depth_s >=3 and tnf_s>=5:
        result="T"
    elif tax_s >= 1 and depth_s >=4 and tnf_s>=6:
        result="T"
    else:
        result="F"
    return result

cluster={}
cluster_scg={}
cluster_tnf={}
cluster_tax={}
cluster_gc={}
cluster_depth={}
cluster_len={}
tmp_ct1=0
tmp_ct2=0
score=0
counter=0

# get length rank
scaf_sorted = [i[0] for i in sorted(scaffold_fitered.items(), key=lambda d:d[1]["length"], reverse= True)]
for k in scaf_sorted:
    if counter == 0 : # set largest scaffold as a first cluster
        counter+=1
        clst_nm="Cluster"+str(counter)
        cluster[clst_nm]=[k]
        cluster_scg[clst_nm]=scg_dic[k]
        cluster_tnf[clst_nm]=scaffold[k]["TNF"]
        cluster_len[clst_nm]=scaffold[k]["length"]
        cluster_tax[clst_nm]=scaffold[k]["tax"]
        cluster_gc[clst_nm]=scaffold[k]["GC"]
        cluster_depth[clst_nm]=scaffold[k]["depth"]
        new_cls="F"
        continue
    for cls in cluster:
        try: ## for unclassified sequences
            tax_s = compare_tax(cluster_tax[cls], scaffold[k]["tax"])
        except KeyError: 
            tax_s = compare_tax(cluster_tax[cls], ["None", "None", "None", "None", "None", "None"])
            scaffold[k]["tax"] = ["None", "None", "None", "None", "None", "None"]
        tnf_s = compare_tnf(cluster_tnf[cls], scaffold[k]["TNF"])
        depth_s = compare_depth(cluster_depth[cls], float(scaffold[k]["depth"]))
        result = score_result(tax_s, tnf_s, depth_s)
        gc_r = compare_gc(scaffold[k]["GC"],cluster_gc[cls])
        tmp_sgc = plus_scg(cluster_scg[cls], scg_dic[k])
        mcg_count = summerize_scg(tmp_sgc)
        if result =="T" and gc_r == "T" and mcg_count <= 6:
            cluster_scg[cls] = tmp_sgc
            cluster_tnf[cls] = cluster_tnf_summarize(cluster_tnf[cls], cluster_len[cls], scaffold[k]["TNF"], scaffold[k]["length"])
            cluster_gc[cls] = (float(cluster_gc[cls])*float(cluster_len[cls])+float(scaffold[k]["GC"])*float(scaffold[k]["length"]))/(float(cluster_len[cls])+float(scaffold[k]["length"]))
            cluster_depth[cls] = (float(cluster_depth[cls])*float(cluster_len[cls])+float(scaffold[k]["depth"])*float(scaffold[k]["length"]))/(float(cluster_len[cls])+float(scaffold[k]["length"]))
            cluster[cls].append(k)
            cluster_len[cls] = float(cluster_len[cls]) + float(scaffold[k]["length"])
            new_cls = "F"
            break
        else : 
            new_cls = "T"
    if new_cls=="T":
        counter+=1
        clst_nm="Cluster"+str(counter)
        cluster[clst_nm]=[k]
        cluster_scg[clst_nm]=scg_dic[k]
        cluster_tnf[clst_nm]=scaffold[k]["TNF"]
        cluster_len[clst_nm]=scaffold[k]["length"]
        cluster_tax[clst_nm]=scaffold[k]["tax"]
        cluster_gc[clst_nm]=scaffold[k]["GC"]
        cluster_depth[clst_nm]=scaffold[k]["depth"]

print ("There are %s clusters" %len(cluster.keys()))

#write result 
print("-------------------------------writing bins ----------------------------------------")
out_subclass = open(outdir + "_cluster_subclass", 'w')
out_cls = open(outdir + "_cluster_length", 'w')


ct_bins=0
for clst in cluster:
    seq_num = 0
    clst_len=0
    for scfd in cluster[clst]:
        seq_num += 1
        clst_len+=float(scaffold[scfd]["length"])
        out_subclass.write(clst + "\t" + scfd + "\n")
    out_cls.write(clst + "\t" + str(seq_num) + "\t" + str(clst_len) + "\n")
    if clst_len >= 1000000:
        ct_bins += 1
    
print("There are %d bins large than 1 Mbp." %ct_bins)

out_subclass.close()
out_cls.close()
