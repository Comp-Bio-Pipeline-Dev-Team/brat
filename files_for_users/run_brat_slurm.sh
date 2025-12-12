#!/bin/bash

#SBATCH --nodes=1 # use one node
#SBATCH --time=10:00:00 # 10h, make sure this is long enough to run the entire analysis
#SBATCH --account=amc-general # normal, amc, long, mem (use mem when using the amem partition) - THESE ARE AMC-ALPINE SPECIFIC
#SBATCH --partition=amilan # amilian, ami100, aa100, amem, amc - THESE ARE AMC-ALPINE SPECIFIC
#SBATCH --qos=normal 
#SBATCH --ntasks=1 # total processes/threads, shouldn't need more than this
#SBATCH --job-name=my_brat_run_DATE
#SBATCH --output=my_brat_run_DATE_%J.log ## consult this file if your brat run fails!
#SBATCH --mem=5G # suffix K,M,G,T can be used with default to M, brat uses <5G 
#SBATCH --mail-user=YOUR_EMAIL@SCHOOL.EDU
#SBATCH --mail-type=ALL ## will send you an email when this job starts, stops, or fails
#SBATCH --error=my_brat_run_DATE_%J.err ## consult this file if your brat run fails!

## THIS IS AMC-ALPINE SPECIFIC
module load miniforge/24.11.3-0
module load singularity/3.6.4

## install and activate brat conda environment
conda env create -f 'brat.yml' -n brat_v1.0.0
conda activate brat_v1.0.0

## if you're using relative paths: you may want to add this command to the sbatch script
## because if not, you have to submit job from within the directory your input files live in
## cd /home/mydir/

## THIS IS AMC-ALPINE SPECIFIC
## to use with apptainer/referencing singularity containers
## can add to your .bashrc so you dont have to include here
export SINGULARITY_CACHEDIR=/scratch/alpine/$USER 
export SINGULARITY_CACHDIR=/scratch/alpine/$USER
export SINGULARITY_TMPDIR=/scratch/alpine/$USER
export TMP=/scratch/alpine/$USER
export TMPDIR=/scratch/alpine/$USER
export TEMP=/scratch/alpine/$USER
export TEMPDIR=/scratch/alpine/$USER

## symlinking apptainer_cache so images aren't saved to home directory (not enough storage)
if [ -L "/scratch/alpine/$USER/apptainer_cache" ]; 
then
    echo "Singularity symlink already exists!"
else
    echo "Singularity symlink does not exist. Generating symlink..."
    ln -s /scratch/alpine/$USER/apptainer_cache ~/.singularity
fi

## adding brat directory to global path so you can execute brat from anywhere in your computer
## this can also be added to your .bashrc so you dont have to run it again like so:
## echo 'export PATH=/scratch/alpine/${USER}/brat/:$PATH' >> ~/.bashrc && source ~/.bashrc
export PATH=/scratch/alpine/${USER}/brat/:$PATH

## if you're getting this error: brat: command not found, run the command below
##chmod +x /scratch/alpine/$USER/brat/brat

## actually running brat
brat \
    --raw_seq_dir 'my_bulk_seqs' \
    --metadata_file 'metadata.csv' \
    --out_dir_name 'brat_out' \ ## this is the default name, can be changed 
    --profile 'slurm_profile.yml' \ ## slurm execution profile
    --fastq_screen_genomes 'fastq_screen_genomes.csv' \ ## fastq screen contaminants
    --fastq_screen_subset_num 0 \ ## number of reads to use for fastq screen  
    --genome_fasta 'GRCh38.p14.genome.fa' \ ## the genome FASTA for your samples
    --genome_name 'GRCh38.p14' \ ## the name of this genome
    --gtf_file 'gencode.v48.annotation.gtf' \ ## the annotation file for your genome
    --annot_col_name 'gene_name' \ ## default is gene_id, this will work for most annotation files except for gencode ones
    --transcript_type 'exon' \ ## default is exon, this will work for most annotation files except for bacterial/viral ones
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