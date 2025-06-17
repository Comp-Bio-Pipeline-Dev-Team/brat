#!/bin/bash

## put symlink ahead of this? - snakemake param where you want singularity images stored 

snakemake \
    -s workflow/snakefile \
    --configfile workflow/config_files/test_config.yml \
    --jobs unlimited \
    --workflow-profile workflow/profiles/default \
    --software-deployment-method apptainer
    ##--dry-run


    ##--unlock
    ##--singularity-args "-B /scratch/alpine/mapgar@xsede.org:/scratch/alpine/mapgar@xsede.org"
    ##--slurm-keep-successful-logs
