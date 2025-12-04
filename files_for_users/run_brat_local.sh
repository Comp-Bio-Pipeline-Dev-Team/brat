#!/bin/bash

# script to run brat (bulk rna-seq analysis tool)

## the python script needs to be verified as an executable!:
## also need to add path to folder where script lives to your global path like so:
##chmod +x ~/brat/brat && echo 'export PATH=/$HOME/brat/:$PATH' >> ~/.bashrc && source ~/.bashrc

## activate brat conda environment
conda activate brat_v1.0.0

## run brat!
brat \
    --raw_seq_dir '/home/my_bulk_seqs' \
    --metadata_file 'metadata.csv' \
    --out_dir_name 'brat_out' \ ## this is the default name, can be changed 
    --profile 'mylocal_profile.yml' \ ## local execution profile
    --fastq_screen_genomes 'fastq_screen_genomes.csv' \ ## fastq screen contaminants 
    --genome_fasta 'GRCh38.p14.genome.fa' \ ## the genome FASTA for your samples
    --genome_name 'GRCh38.p14' \ ## the name of this genome
    --gtf_file 'gencode.v48.annotation.gtf' \ ## the annotation file for your genome
    --annot_col_name 'gene_name' \ ## default is gene_id, this will work for most annotation files except for gencode ones
    --read_length 0 \ ## default, we will calculate the read length for you unless specified
    --refflat 'gencode.v48.annotation.refflat' \ ## picard refFlat file
    --ribosomal_int_list 'GRCh38.p14_ribosomal_interval_list_with_header.txt' \ ## picard rRNA interval list
## optional parameters!! ##
    --run_rsem \ ## if included, will run RSEM as an optional quantification method 
    --extra_cutadapt_params ' ' \
    --extra_star_params ' ' \
    --extra_rsem_params ' ' \
    --latency_wait 60 \ ## default is 60s, can be increased/decreased
    --use_singularity \ ## if included, will run with singularity/apptainer
    --dry_run ## if included, will pretend to run brat to check for file path errors before you actually run it 