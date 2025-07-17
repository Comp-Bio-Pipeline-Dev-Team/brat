#!/bin/bash

## put symlink ahead of this? - snakemake param where you want singularity images stored 

snakeLog="logs/07172025/snake_log_07172025.log"

snakemake \
    -s workflow/snakefile \
    --configfile workflow/config_files/test_config.yml \
    --jobs unlimited \
    --workflow-profile workflow/profiles/default \
    --software-deployment-method apptainer \
    --rerun-incomplete \
    --dry-run
    ##2>${snakeLog}
    ##--unlock
    ##--dry-run


    ##--unlock
    ##--singularity-args "-B /scratch/alpine/mapgar@xsede.org:/scratch/alpine/mapgar@xsede.org"
    ##--slurm-keep-successful-logs
