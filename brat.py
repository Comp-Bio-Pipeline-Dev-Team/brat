#!/usr/bin/env python3

import subprocess
import argparse
import os
from os.path import join as pj
from os.path import dirname
from ruamel.yaml import YAML
from ruamel.yaml.scalarstring import SingleQuotedScalarString
from importlib.resources import files

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw_seq_dir",
                        help="The filepath to the directory containing raw sequencing files.")
    parser.add_argument("--metadata_file",
                        help="The filepath to the metadata file as a .tsv. \
                              Your metadata file must contain at least three columns: \
                              'sampleid', 'forward', and 'reverse'.")
    parser.add_argument("--out_dir_name",
                        help="The name of the output directory, default is 'brat_out'.",
                        default="brat_out")
    parser.add_argument("--genome_fasta",
                        help="The filepath to the genome fasta file.")
    parser.add_argument("--genome_name",
                        help="The name of the genome, e.g. hg38, mm10, etc.")
    parser.add_argument("--gtf_file",
                        help="The filepath to the gtf annotation file.")
    parser.add_argument("--annot_col_name",
                        help="The column name in the gtf file that contains the gene names, default is gene_id.",
                        default="gene_id")
    parser.add_argument("--read_length",
                        help="The length of the reads, if not specified, the pipeline will infer read length.",
                        default=0)
    parser.add_argument("--refflat",
                        help="The filepath to the picard refflat file.")
    parser.add_argument("--ribosomal_int_list",
                        help="The filepath to the picard ribsomal interval list file.")
    parser.add_argument("--run_rsem",
                        help="Boolean, True or False, whether or not to run rsem quantification, default is False.",
                        action='store_true')
    parser.add_argument("--extra_cutadapt_params",
                        help="Extra parameters to pass to cutadapt, written on a single line and enclosed in single quotes.",
                        default=" ")
    parser.add_argument("--extra_star_params",
                        help="Extra parameters to pass to STAR, written on a single line and enclosed in single quotes.",
                        default=" ")
    parser.add_argument("--extra_rsem_params",
                        help="Extra parameters to pass to RSEM, written on a single line and enclosed in single quotes.",
                        default=" ")
    parser.add_argument("--latency_wait",
                        help="The amount of time (in seconds) to wait for files to appear",
                        default=60)
    parser.add_argument("--dry_run",
                        help="If this parameter is specified, you can practice running the workflow without \
                              actually starting it",
                        action='store_true')
    return parser.parse_args()

## functions
## dirname(__file__) gets the directory where this script (brat.py) is located
def get_snake_path():
    return pj(dirname(__file__), "workflow/toy_snakefile.smk")

def get_config_path():
    return pj(dirname(__file__), "workflow/config_files/test_brat_config.yml")

def get_profile_path():
    return pj(dirname(__file__), "workflow/profiles/default")


def create_config_file(config_path,
                       args):
    yaml = YAML()
    yaml.preserve_quotes = True
    yaml.default_flow_style = False

    config_params = {"raw_seq_in": args.raw_seq_dir,
                     "metadata_file": args.metadata_file,
                     "out_dir": args.out_dir_name,
                     "align_to_fasta": args.genome_fasta,
                     "align_to_name": args.genome_name,
                     "align_to_gtf": args.gtf_file,
                     "gtf_annot_col": args.annot_col_name,
                     "raw_seq_read_length": int(args.read_length),
                     "picard_refFlat": args.refflat,
                     "picard_rrna_interval_list": args.ribosomal_int_list,
                     "cutadapt_params": SingleQuotedScalarString(args.extra_cutadapt_params),
                     "star_params": SingleQuotedScalarString(args.extra_star_params),
                     "rsem_params": SingleQuotedScalarString(args.extra_rsem_params)}
    
    if args.run_rsem is True:
        config_params["run_rsem"] = True
    else:
        config_params["run_rsem"] = False
    
    with open(config_path, 'w') as outfile:
        yaml.dump(config_params, outfile)
    

def assemble_snake_command(snake_path,
                           config_path,
                           profile_path,
                           args):

    snake_command = ["snakemake",
                    "-s", snake_path,
                    "--configfile", config_path,
                    "--jobs", "unlimited",
                    "--workflow-profile", profile_path,
                    "--software-deployment-method", "apptainer",
                    "--rerun-incomplete",
                    "--latency-wait", str(args.latency_wait)]
    
    if args.dry_run is True:
        snake_command.append("--dry-run")
    
    return snake_command


def main():
    args = get_args()

    relative_config_path = "workflow/config_files/test_brat_config.yml"
    
    create_config_file(get_config_path(),
                       args)
    
    command = assemble_snake_command(get_snake_path(),
                                     get_config_path(),
                                     get_profile_path(),
                                     args)
    
    complete = subprocess.run(command)


if __name__=="__main__":
    main()
