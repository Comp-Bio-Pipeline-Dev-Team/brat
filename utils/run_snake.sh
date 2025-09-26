#!/bin/bash

## put symlink ahead of this? - snakemake param where you want singularity images stored 
## if [ -L "/scratch/alpine/$USER/apptainer_cache" ]; 
##then
##    echo "Singularity symlink already exists!"
##else
##    echo "Singularity symlink does not exist. Generating symlink..."
##    ln -s /scratch/alpine/$USER/apptainer_cache ~/.singularity
##fi

snakeLog="logs/toySnake_rsem_snakeLog_09252025.log"

snakemake \
    -s workflow/toy_snakefile.smk \
    --configfile workflow/config_files/test_sample_config.yml \
    --jobs unlimited \
    --workflow-profile workflow/profiles/debug \
    --software-deployment-method apptainer \
    --rerun-incomplete \
    --dry-run \
    2>${snakeLog}
    ##--dry-run
    ##--unlock


    ##--singularity-args "-B /scratch/alpine/mapgar@xsede.org:/scratch/alpine/mapgar@xsede.org"
