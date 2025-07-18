#!/bin/bash

## put symlink ahead of this? - snakemake param where you want singularity images stored 
## if [ -L "/scratch/alpine/$USER/apptainer_cache" ]; 
##then
##    echo "Singularity symlink already exists!"
##else
##    echo "Singularity symlink does not exist. Generating symlink..."
##    ln -s /scratch/alpine/$USER/apptainer_cache ~/.singularity
##fi

snakeLog="logs/07182025/snake_log_07182025.log"

snakemake \
    -s workflow/snakefile \
    --configfile workflow/config_files/test_config.yml \
    --jobs unlimited \
    --workflow-profile workflow/profiles/debug \
    --software-deployment-method apptainer \
    --rerun-incomplete \
    --dry-run
    ##2>${snakeLog}
    ##--unlock


    ##--singularity-args "-B /scratch/alpine/mapgar@xsede.org:/scratch/alpine/mapgar@xsede.org"
