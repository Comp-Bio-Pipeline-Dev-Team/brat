#!/bin/bash

## put symlink ahead of this? - snakemake param where you want singularity images stored 
## if [ -L "/scratch/alpine/$USER/apptainer_cache" ]; 
##then
##    echo "Singularity symlink already exists!"
##else
##    echo "Singularity symlink does not exist. Generating symlink..."
##    ln -s /scratch/alpine/$USER/apptainer_cache ~/.singularity
##fi

snakeLog="logs/2roberta_compGenome_snakeLog_07232025.log"

snakemake \
    -s workflow/snakefile \
    --configfile workflow/config_files/roberta_compGenome_config.yml \
    --jobs unlimited \
    --workflow-profile workflow/profiles/default \
    --software-deployment-method apptainer \
    --rerun-incomplete \
    2>${snakeLog}
    ##--dry-run
    ##--unlock


    ##--singularity-args "-B /scratch/alpine/mapgar@xsede.org:/scratch/alpine/mapgar@xsede.org"
