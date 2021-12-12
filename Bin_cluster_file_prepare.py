import os
import glob
from pathlib import Path
import snakemake
import gzip

workdir: config['workdir']
parse_dia_out_script = "/mnt1/links/diamond_report_modified_contig.pl"
SGC_hmm_path = "/mnt1/database/SGC_db/bacteria_139_CSCG.hmm"


sample_list = {}
with open(config['file_names_txt'],'r') as f:
    for line in f:
        items = line.strip().split("\t")
        sample_list[items[0]] = items[1]

sample_id = list(sample_list.keys())
print(sample_list)
##
rule all:
    input:
        # expand("02_contig_taxonomy/{sample}_con.tax", sample = sample_id),
        "01_protein_all.gatk_con.tax.splt.done",
        expand("02_contig_5mer_freq/{sample}.5mer", sample = sample_id),
        expand("02_contig_SCG/{sample}.scg.bt.out", sample = sample_id)

rule run_prodigal:
    input:
        lambda wildcards: sample_list[wildcards.sample]
    output:
        "01_protein/{sample}.faa",
        "01_protein/{sample}.gff"
    threads:5
    shell:
        """
        prodigal -i {input} -f gff -a {output[0]} -o {output[1]} -q 
        """
rule run_cat_prot:
    input:
        expand("01_protein/{sample}.faa", sample = sample_id)
    output:
        temp("01_protein_all.faa")
    shell:
        """
        cat {input} > {output}
        """

rule run_get_faa_names:
    input:
        "01_protein_all.faa"
    output:
        "01_protein_all.faa.names"
    shell:
        """
        grep '>' {input} | perl -lane '@F[0]=~/>(.*)_(\d+)/; print "$1_$2\t$1" ' > {output}
        """

rule run_diamond:
    input:
        "01_protein_all.faa"
    output:
        "01_protein_all.gatk.dia.tsv"
    threads:100
    shell:
        """
        diamond blastp --quiet -k 1 --threads {threads} --db /mnt1/database/GTDB/protein_db/GTDB_rep_species.faa.dmnd --query {input} --outfmt 6 qseqid sseqid stitle pident qlen slen length mismatch gapopen qstart qend sstart send evalue bitscore --out {output}
        """

rule run_tax_parse:
    input:
        "01_protein_all.gatk.dia.tsv",
        "01_protein_all.faa.names"
    output:
        "01_protein_all.gatk.con.tax"
    params:
        parse_script = parse_dia_out_script
    shell:
        """
        {params.parse_script} {input} 01_protein_all.gatk
        """

rule run_get_scaf2sample:
    input:
        lambda wildcards: sample_list[wildcards.sample]
    output:
        temp("00_scaf2sampe.file/{sample}.fa2sample")
    shell:
        """
        Fasta2sample.sh -i {input} -n {wildcards.sample} > {output}
        """

rule run_merge_scaf2sample_files:
    input:
        expand("00_scaf2sampe.file/{sample}.fa2sample", sample = sample_id)
    output:
        "00_scaf2sampe.list"
    shell:
        """
        cat {input} > {output}
        """

checkpoint run_tax_split:
    input:
        rules.run_tax_parse.output,
        "00_scaf2sampe.list"
    output:
        touch("01_protein_all.gatk_con.tax.splt.done"),
        directory("02_contig_taxonomy")
    params:
        outdir = "02_contig_taxonomy"
    shell:
        """
        Bin_taxonomy_file_split.py -i {input[0]} -f {input[1]} -o {params.outdir}
        """

rule run_5mer_freq:
    input:
        lambda wildcards: sample_list[wildcards.sample]
    output:
        "02_contig_5mer_freq/{sample}.5mer"
    threads:3
    shell:
        """
        create_kmer_freq_modified.py {input} 5 > {output}
        """

rule run_SCG:
    input:
        "01_protein/{sample}.faa"
    output:
        temp("02_contig_SCG/{sample}.scg.hmm.out"), 
        "02_contig_SCG/{sample}.scg.bt.out"
    threads:2
    params:
        SGC_db = SGC_hmm_path
    shell:
        """
        hmmscan -o {output[0]} --tblout {output[1]} --cpu {threads} {params.SGC_db} {input}
        """