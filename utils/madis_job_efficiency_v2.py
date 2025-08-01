"""
This script parses the log files output by a head snakemake job and finds job efficiencies

Assumptions:
	1. This assumes that the log files have a certain structure relating to the "log" and submitting information.
	   This structure will be made by default using the slurm executor in snakemake version 8. I'm not sure about
	   the structure in other instances
	2. This assumes that the seff command is available on your system and has the format:
		```
		Job ID: 5915204
		Cluster: alpine
		User/Group: kwellswrasman@xsede.org/kwellswrasmanpgrp@xsede.org
		State: COMPLETED (exit code 0)
		Nodes: 1
		Cores: 11
		CPU Utilized: 05:18:00
		CPU Efficiency: 82.36% of 06:26:06 core-walltime
		Job Wall-clock time: 00:35:06
		Memory Utilized: 5.40 GB
		Memory Efficiency: 13.51% of 40.00 GB
		```

Still to do
	1. I want to clean up the output to show the percents more cleanly. This will mean better parsing the effiencies
	2. Eventually it might be helpful to plot these in some way. Maybe a bar plot where the color is the job name?
"""

import re
import sys
import argparse
import glob
import os
from os.path import basename
from os.path import dirname
import warnings
import subprocess
import pandas as pd


def main():
	options = setup()
	job_dict = get_job_information(options.in_dir)
	get_efficiencies(job_dict)
	make_data_frame(job_dict, options.out_file, options.keep_failed)

#########
# Setup #
#########

def setup():
    """
    Gets command line arguments and returns a Namespace object
    """

    # House keeping to read in arguments from the command line
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--directory", dest = "in_dir",
        help = "Path to the directory containing the logs from the master job",
        default = "none",
        action = "store",
        metavar = "\b")

    parser.add_argument("-o", "--out_file", dest = "out_file",
        help = "Path to the output file to save stats to",
        default = "none",
        action = "store",
        metavar = "\b")


    parser.add_argument("--keep_failed", action="store_true",
		help = "Keep failed jobs in final output")

    args = parser.parse_args()

    return(args)

#######################
# Get job information #
#######################

def get_job_information(in_dir):
	"""
	Parse subdirectories in the logs directory and pull out all job ids, job names, and run dates 

	This assumes that the log files are in subdirectories based on the date (MMDDYYYY) and that sbatch .log/.err file names are 
	formatted like so:
	- run_name_jobID.log
	- run_name_jobID.err

	A dictionary is made that has the job id as the key and a second dictionary with the job name as the
	value (and job_name as the key).

	This structure is so that we can easily add more information about the job and make it into a nice pandas data frame
	"""
	job_dict = {}
	all_files = [f for f in glob.glob(glob.escape(in_dir) + "/*/*.err")]

	for file_path in all_files:
		## pull file name from file path
		file_name = basename(file_path)

		## pull base directory of the file from file path (the job name)
		job_name = dirname(file_path).split("/")[-1]

		## pull the name of what was run from the file name
		## are there always 8 digits in job ids created by slurm??
		##job_name = re.sub("_[0-9]{8}.err", "", file_name)

		##\[(.*?)\]

		## pull name of wildcard from log file
		with open(file_path, "r") as input_file:
			for line in input_file:
				## regex to match any letter, number, or special character that isn't a front slash, newline, or return character
				match_sample_obj = re.search(r"wildcards: sample=([^/\n\r]*)", line)
				match_date_obj = re.search(r"^\[(.*?)\]", line)
				if match_sample_obj is not None:
					pre_sample_name = match_sample_obj.group(1)
					## regex to start from the beginning of the line and pull everything before the first comma
					sample_name = re.match(r"^([^,]+)", pre_sample_name).group(1)
					##print(sample_name)
				if match_date_obj is not None:
					run_date = match_date_obj.group(1)
					##print(run_date)

		## pull the job id out of the file name 
		##job_id = file_name.split("_", 3)[-1]
		job_id = re.sub(".err", "", file_name)

		job_dict[job_id] = {"job_name": job_name,
					  		"sample_id": sample_name,
							"date": run_date}
	
	return(job_dict)

def get_efficiencies(job_dict):
	"""
	This uses subprocess to run the `seff` command and captures the output. The output is then parsed
	for all of the helpful efficiency measures.

	This is added to the existing dictionary using the job id. The nested dictionaries are expanded with
	all of the additional efficiency measures
	"""
	for job_id in job_dict:
		command = ["seff", job_id]
		# Execute the command and capture its output
		command_output = subprocess.run(command, capture_output = True)

		# Grab the standard out, decode, and split
		output_lines = command_output.stdout.decode('utf-8').splitlines()

		# Iterate through list to pull out what we want
		for line in output_lines:
			if re.search("State", line):
				job_dict[job_id]["state"] = line.split(":")[1]
			elif re.search("Cores", line):
				job_dict[job_id]["cores"] = line.split(":")[1]	
			elif re.search("CPU Utilized", line):
				job_dict[job_id]["cpu_utilized"] = line.split(": ")[1]
			elif re.search("CPU Efficiency", line):
				job_dict[job_id]["cpu_efficiency"] = line.split(": ")[1]
			elif re.search("Memory Utilized", line):
				job_dict[job_id]["memory_utilized"] = line.split(":")[1]
			elif re.search("Memory Efficiency", line):
				job_dict[job_id]["memory_efficiency"] = line.split(":")[1]
			elif re.search("Job Wall-clock", line):
				job_dict[job_id]["run_time"] = line.split(": ")[1]


def make_data_frame(job_dict, output_file, keep_failed):
	"""
	This makes a pandas data frame using the layered defaut dictionary made previously. This is then
	transformed such that the job ids are the rows and the stats are the columns.

	You can chose (via command line options) to either keep or throw out the failed jobs. Default is
	to throw them out.

	The final data frame is saved as a csv file.
	"""
	df = pd.DataFrame(job_dict).T

	if(keep_failed):
		df.to_csv(output_file)
	else:
		filtered_df = df[~df["state"].str.contains("FAILED")]
		filtered_df.to_csv(output_file)
	

if __name__ == "__main__":
	main()