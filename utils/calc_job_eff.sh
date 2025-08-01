#!/bin/bash

echo "gathering job efficiency stats...."

## relative to where you're calling the script from!! not where the script actually lives
slurmLogs="logs/slurm/"
jobStatsFile="logs/roberta_compGenome_jobs_07182025.csv"

python3 utils/madis_job_efficiency_v2.py --directory ${slurmLogs} \
				      	 --out_file ${jobStatsFile}
				      	 ##--keep_failed

echo "saving job efficiency stats to ${jobStatsFile}!"
