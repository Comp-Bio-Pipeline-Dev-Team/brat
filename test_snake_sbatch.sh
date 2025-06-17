#!/bin/bash

#SBATCH --nodes=1 # use one node
#SBATCH --time=08:00:00 # 2h, madi edit: had to up the time so the job wasn't canceled
#SBATCH --account=amc-general # normal, amc, long, mem (use mem when using the amem partition)
#SBATCH --partition=amilan # amilian, ami100, aa100, amem, amc
#SBATCH --qos=normal 
#SBATCH --ntasks=1 # total processes/threads, madi edit: need n tasks to be the same as --threads in fastqc command
#SBATCH --job-name=test_star_genome_index_gen
#SBATCH --output=logs/06162025/test_star_index_%J.log
#SBATCH --mem=5G # suffix K,M,G,T can be used with default to M, madi edit: uses <5G 
#SBATCH --mail-user=MADISON.APGAR@CUANSCHUTZ.EDU
#SBATCH --mail-type=FAIL,END
#SBATCH --error=logs/06162025/test_star_index_%J.err

module load miniforge/24.11.3-0
module load singularity/3.6.4
conda activate snake

## may want to add this command to the sbatch script bc if not you have to submit job from within the test_snake_w_slurm dir
## cd /scratch/alpine/mapgar@xsede.org/test_snake_w_slurm

## to use with apptainer/referencing singularity containers
## can add to .bashrc
##export SINGULARITY_CACHEDIR=/scratch/alpine/$USER 
##export SINGULARITY_CACHDIR=/scratch/alpine/$USER
##export SINGULARITY_TMPDIR=/scratch/alpine/$USER
##export TMP=/scratch/alpine/$USER
##export TMPDIR=/scratch/alpine/$USER
##export TEMP=/scratch/alpine/$USER
##export TEMPDIR=/scratch/alpine/$USER

##/scratch/alpine/$USER

snakemake \
    -s workflow/snakefile \
    --configfile workflow/config_files/test_config.yml \
    --workflow-profile workflow/profiles/default \
    --jobs unlimited \
    --software-deployment-method apptainer

## add --rerun-incomplete flag eventually!
## may need to adjust the jobs parameter!