#!/bin/bash

echo "gathering job efficiency stats...."

## relative to where you're calling the script from!! not where the script actually lives
slurmLogs="logs/slurm/"
jobStatsFile="hyperthreaded_bratJobs_01152026.csv"

python3 madis_job_efficiency_v2.py --directory ${slurmLogs} \
				      	 		   --out_file ${jobStatsFile}
				      	 		   ##--keep_failed

echo "saving job efficiency stats to ${jobStatsFile}!"
