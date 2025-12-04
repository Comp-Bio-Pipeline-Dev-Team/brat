#!/bin/bash

# script to run brat (bulk rna-seq analysis tool)

## the python script needs to be verified as an executable!:
## changing script name to just "brat" for ease of use (and bc its fun)
## also need to add path to folder where script lives to my global path like so:
##chmod +x brat/brat.py && mv brat/brat.py brat/brat && echo 'export PATH=/scratch/alpine/${USER}/brat/:$PATH' >> ~/.bashrc && source ~/.bashrc

## actual command
brat \
    --raw_seq_dir 'rp_bulkRNAseq' \
    --metadata_file 'rp_bulkRNAseq/rp_mini_metadata.tsv' \
    --out_dir_name 'rp_mini_samples' \
    --profile 'files_for_users/slurm_profile.yml' \
    --fastq_screen_genomes 'files_for_users/custom_fastq_screen_genomes.csv' \
    --genome_fasta '/projects/mapgar@xsede.org/outside_rna_seq_refs/human_comprehensive/GRCh38.p14.genome.fa' \
    --genome_name 'GRCh38.p14' \
    --gtf_file '/projects/mapgar@xsede.org/outside_rna_seq_refs/human_comprehensive/gencode.v48.annotation.gtf' \
    --annot_col_name 'gene_name' \
    --read_length 0 \
    --refflat '/projects/mapgar@xsede.org/outside_rna_seq_refs/human_comprehensive/gencode.v48.annotation.refflat' \
    --ribosomal_int_list '/projects/mapgar@xsede.org/outside_rna_seq_refs/human_comprehensive/GRCh38.p14_ribosomal_interval_list_with_header.txt' \
    --run_rsem \
    --extra_cutadapt_params ' ' \
    --extra_star_params '--sjdbGTFtagExonParentGene gene_name --outSAMunmapped Within KeepPairs --outReadsUnmapped Fastx --outFilterScoreMinOverLread 0.20 --outFilterMatchNminOverLread 0.20' \
    --extra_rsem_params ' ' \
    --latency_wait 60 \
    --use_singularity \
    2>rp_mini_fastq_screen.log
    ##--dry_run




