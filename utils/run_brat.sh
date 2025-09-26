#!/bin/bash

# script to run brat (bulk rna-seq analysis tool)

## the python script needs to be verified as an executable! :
##chmod +x brat.py
##mv brat.py brat

## also need to add path to folder where script lives to my global path like so:
##echo 'export PATH=/scratch/alpine/${USER}/test_snake_w_slurm/:$PATH' >> ~/.bashrc
## source ~/.bashrc

## actual command
brat \
    --raw_seq_dir 'test_seqs' \
    --metadata_file 'test_samp_metadata.tsv' \
    --out_dir_name 'test' \
    --genome_fasta '/scratch/alpine/mapgar@xsede.org/outside_rna_seq_refs/human_comprehensive/GRCh38.p14.genome.fa' \
    --genome_name 'GRCh38.p14' \
    --gtf_file '/scratch/alpine/mapgar@xsede.org/outside_rna_seq_refs/human_comprehensive/gencode.v48.annotation.gtf' \
    --annot_col_name 'gene_name' \
    --read_length 0 \
    --refflat '/scratch/alpine/mapgar@xsede.org/outside_rna_seq_refs/human_comprehensive/gencode.v48.annotation.refflat' \
    --ribosomal_int_list '/scratch/alpine/mapgar@xsede.org/outside_rna_seq_refs/human_comprehensive/GRCh38.p14_ribosomal_interval_list_with_header.txt' \
    --run_rsem \
    --extra_cutadapt_params ' ' \
    --extra_star_params '--sjdbGTFtagExonParentGene gene_name --outSAMunmapped Within KeepPairs --outReadsUnmapped Fastx --outFilterScoreMinOverLread 0.20 --outFilterMatchNminOverLread 0.20' \
    --extra_rsem_params ' ' \
    --latency_wait 60 \
    2>logs/brat.log
    ##--dry_run




