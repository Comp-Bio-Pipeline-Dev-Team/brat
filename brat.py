#!/usr/bin/env python3

import subprocess
import argparse
import os
from os.path import join as pj
from os.path import dirname
from os.path import exists
from pathlib import Path
from ruamel.yaml import YAML
from ruamel.yaml.scalarstring import SingleQuotedScalarString

yaml = YAML()

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
    ## need to figure out the default profile for this!! (local execution?)
    parser.add_argument("--profile",
                        help="The filepath to the snakemake profile .yaml file you want to use.")
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
    parser.add_argument("--use_singularity",
                        help="If this parameter is specified, the workflow will run using singularity containers \
                              instead of conda environments (default). NOTE: apptainer MUST be installed to run this \
                              pipeline with singularity!",
                        action='store_true')
    parser.add_argument("--dry_run",
                        help="If this parameter is specified, you can practice running the workflow without \
                              actually starting it",
                        action='store_true')
    return parser.parse_args()

## functions
## dirname(__file__) gets the directory where this script (brat.py) is located
def get_snake_path():
    return pj(dirname(__file__), "workflow/snakefile")

def get_config_path():
    return pj(dirname(__file__), "workflow/config_files/brat_config.yml")

def get_multiqc_config_path():
    return pj(dirname(__file__), "workflow/config_files/multiqc_config.yaml")

def get_profile_path():
    return pj(dirname(__file__), "workflow/profiles/default")

## make sure profile path exists and is a file not a directory
## come up with multiple profiles and let user choose the best one for their use case
def move_workflow_profile(profile_path,
                          args):
    if exists(args.profile) and not os.path.isdir(args.profile):
        print("Setting up user profile...")
        fix_profile_config_name = ["cp", args.profile, pj(profile_path, "config.yaml")]
        subprocess.run(fix_profile_config_name) 
        print("Setup successful!")
    else:
        print("The specified profile path does not exist or is a directory, please check that your input is a .yaml file and try again")
        exit()

## attempting to move multiqc config file to tmp.brat of current working directory so singularity can find it
def move_multiqc_config(multiqc_config_path):
    multiqc_config_src = multiqc_config_path
    multiqc_config_dest = pj(os.getcwd(), "tmp.brat/multiqc_config.yaml")
    if not exists(pj(os.getcwd(), "tmp.brat")):
        os.makedirs(pj(os.getcwd(), "tmp.brat"))
    copy_multiqc_config = ["cp", multiqc_config_src, multiqc_config_dest]
    subprocess.run(copy_multiqc_config)
    print(f"multiqc config file copied to {multiqc_config_dest}")

## NOTE!!: we don't need to symlink raw seqs directory anymore bc singularity can find anything relative to where the script is run from
## symlinking raw seqs directory to working directory of pipeline so bind points dont break 
##def symlink_raw_seqs(args):
##    if exists(args.raw_seq_dir) and os.path.isdir(args.raw_seq_dir):
##        outside_dir = args.raw_seq_dir
##        #wanted_symlink = pj(dirname(__file__), Path(outside_dir).name)
##        wanted_symlink = pj(os.getcwd(), Path(outside_dir).name)
##        if not exists(wanted_symlink):
##            os.symlink(outside_dir, wanted_symlink,
##                       target_is_directory=True)
##            print(f"symlink created at {wanted_symlink}")
##    else:
##        print("--raw_seq_dir does not exist or is not a directory, please check your input and try again")
##        exit()


def create_config_file(config_path,
                       args):
    yaml.preserve_quotes = True
    yaml.default_flow_style = False

    ## trying to set multiqc config path to tmp.brat so singularity can find it 
    ## and so snakemake doesnt lose its mind about the absolute path 
    multiqc_config_dest = "tmp.brat/multiqc_config.yaml"

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
                     "rsem_params": SingleQuotedScalarString(args.extra_rsem_params),
                     "run_rsem": True if args.run_rsem else False,
                     "deployment_method": "singularity" if args.use_singularity else "conda",
                     "multiqc_config_file": multiqc_config_dest}
    
    with open(config_path, 'w') as outfile:
        yaml.dump(config_params, outfile)


## print run info to console and append to multiqc config file
def get_run_info(args):
    config_text = f"""Raw Sequence Directory:             {args.raw_seq_dir}
    Metadata File:                      {args.metadata_file}
    Snakemake Profile:                  {args.profile}
    Genome Used:                        {args.genome_name}
    Genome FASTA:                       {args.genome_fasta}
    Annotation File (GTF):              {args.gtf_file}
    Annotation Gene Name Column:        {args.annot_col_name}
    Raw Sequence Read Length:           {args.read_length}
    Picard RefFlat File:                {args.refflat}
    Picard Ribosomal Intervals File:    {args.ribosomal_int_list}
    MultiQC Config File:                tmp.brat/multiqc_config.yaml
    RSEM Quantification:                {'True' if args.run_rsem else 'False'}
    Additional Cutadapt Parameters:     {args.extra_cutadapt_params}
    Additional STAR Parameters:         {args.extra_star_params}
    Additional RSEM Parameters:         {args.extra_rsem_params}
    Latency Wait (secs):                {args.latency_wait}
    Deployment Method:                  {'singularity' if args.use_singularity else 'conda'}
    Dry Run:                            {'True' if args.dry_run else 'False'}"""

    ## print run config to console
    print(f"""Running with config:\n    {config_text}""")

    ## create multiqc config file with run config info appended to it
    config_dict = {"report_comment": f"you ran this workflow with:\n\n <pre><code>\n    {config_text}\n</code></pre>\n"}
    
    yaml.preserve_quotes = True
    yaml.default_flow_style = False
    yaml.width = 100000 ## make sure long lines dont get wrapped
    tmp_yml_path = "tmp.brat/tmp.yml"
    with open(tmp_yml_path, "w") as f:
        yaml.dump(config_dict, f)
    
    ## concatenate tmp yml file to multiqc config file in tmp.brat
    command = ["cat", tmp_yml_path, ">>", "tmp.brat/multiqc_config.yaml"]
    subprocess.run(" ".join(command), shell=True)



def assemble_snake_command(snake_path,
                           config_path,
                           profile_path,
                           args):

    snake_command = ["snakemake",
                    "-s", snake_path,
                    "--configfile", config_path,
                    "--workflow-profile", profile_path,
                    "--rerun-incomplete",
                    "--latency-wait", str(args.latency_wait)]
    
    if args.use_singularity is True:
        snake_command.extend(["--software-deployment-method", "apptainer"])
    else:
        snake_command.extend(["--software-deployment-method", "conda"])
    
    if args.dry_run is True:
        snake_command.append("--dry-run")
    
    return snake_command


def main():
    args = get_args()

    move_multiqc_config(get_multiqc_config_path())
    
    create_config_file(get_config_path(),
                       args)
    
    move_workflow_profile(get_profile_path(),
                          args)
    
    command = assemble_snake_command(get_snake_path(),
                                     get_config_path(),
                                     get_profile_path(),
                                     args)
    
    ## print snakemake command and run info to console upon runtime 
    print("Running brat with the following command:")
    print(" ".join(command))
    get_run_info(args)
    
    complete = subprocess.run(command)


if __name__=="__main__":
    main()
