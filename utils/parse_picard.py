import pandas as pd
from os.path import basename
from os.path import dirname
from os.path import exists
from os.path import join as pj
import re
import glob
import subprocess
import ruamel.yaml
from ruamel.yaml import YAML

## finds line number for the line containing a specified string and returns only the line number 
##def get_lines_to_skip(file_path,
##                      line_pattern):
##    with open(file_path, 'r') as raw_file:
##            lines_w_nums = [(line_num, line.strip()) for line_num,line in enumerate(raw_file, start=1)]
##            skip_line = list(filter(lambda item: line_pattern in item[1], lines_w_nums))[0][0]
##    return(skip_line)

## copilot's optimized version of the function:
## i do like this one better 
def get_lines_to_skip(file_path, 
                      line_pattern):
    with open(file_path, 'r') as raw_file:
        for line_num, line in enumerate(raw_file, start=1):
            if re.search(line_pattern, line):
                return line_num
    raise ValueError(f"Pattern '{line_pattern}' not found in file '{file_path}'")
    

## tell me which line ## METRICS is on and then skip that many lines 
def wrangle_picard_output(file_pattern,
                          out_file):
    picard_collectRNAseq_files = [f for f in glob.glob(file_pattern)]
    search_line_pattern = r"^## METRICS"

    picard_out_list = []

    for file in picard_collectRNAseq_files:

        skip_line_num = get_lines_to_skip(file_path=file,
                                          line_pattern=search_line_pattern)

        picard_file = pd.read_csv(file, 
                                 sep="\t",
                                 skiprows=skip_line_num,
                                 nrows=1)

        file_name = basename(file)
        sample = re.sub(".picard.metrics.txt", "", file_name)
        picard_file["sampleid"] = sample

        if (picard_file["PCT_R1_TRANSCRIPT_STRAND_READS"] >= 0.70).all():
            picard_file["strandedness"] = "R1"
            picard_file["rsem_strand_key"] = 1
        elif (picard_file["PCT_R2_TRANSCRIPT_STRAND_READS"] >= 0.70).all():
            picard_file["strandedness"] = "R2"
            picard_file["rsem_strand_key"] = 0
        else:
            picard_file["strandedness"] = "unstranded"
            picard_file["rsem_strand_key"] = 0.5

        picard_out_list.append(picard_file)

    all_picard_res = pd.concat(picard_out_list, ignore_index=True)
    all_picard_res.to_csv(out_file, sep="\t")

picard_fp = "/Users/apgarm/projects/pi_projects/rpelanda/bulkRNAseq_humanized_mice_autoreactive_vs_nonAutoreactive_bcells_01092025/comp_genome_redo/picard/*.picard.metrics.txt"
#wrangle_picard_output(file_pattern=picard_fp,
                      #out_file="test_picard_res.tsv")



def wrangle_picard_insert_size(file_pattern,
                               out_file):
    picard_collectInsertSize_files = [f for f in glob.glob(file_pattern)]
    search_line_pattern = r"^## METRICS"

    picard_out_list = []

    for file in picard_collectInsertSize_files:

        skip_line_num = get_lines_to_skip(file_path=file,
                                          line_pattern=search_line_pattern)
        picard_file = pd.read_csv(file, 
                                  sep="\t",
                                  skiprows=skip_line_num,
                                  nrows=1)

        file_name = basename(file)
        sample = re.sub(".picard.insertSize.txt", "", file_name)
        picard_file["sampleid"] = sample

        picard_out_list.append(picard_file)

    all_picard_res = pd.concat(picard_out_list, ignore_index=True)
    all_picard_res.to_csv(out_file, sep="\t")


##wrangle_picard_insert_size(file_pattern="/Users/apgarm/projects/pi_projects/rpelanda/bulkRNAseq_humanized_mice_autoreactive_vs_nonAutoreactive_bcells_01092025/comp_genome_redo/picard/*.picard.insertSize.txt",
                           ##out_file="test_picard_insertSize.tsv")



def wrangle_star_log(file_pattern,
                     out_file):
    star_log_files = [f for f in glob.glob(file_pattern)]

    star_out_list = []

    for file in star_log_files:
        star_file = pd.read_csv(file, 
                        sep="\t",
                        names=["key", "value"])
        
        ## remove the pipe symbol from the first column 
        star_file["key"] = star_file["key"].str.replace(" |", "")
        ## pull excess white space off of the values in the first column
        star_file["key"] = star_file["key"].str.strip().replace(r"[\s/]", "_", regex=True).replace(r"[,:-]", "", regex=True)
        ##star_file["key"] = star_file["key"].str.replace(r"[\s/]", "_", regex=True)
        ##star_file["key"] = star_file["key"].str.replace(r"[,:-]", "", regex=True)
        star_file["key"] = star_file["key"].str.replace("%", "pct")
        ## remove percent signs from the values in the second column 
        star_file["value"] = star_file["value"].str.replace("%", "")
        ## drop na values from the second column 
        star_file.dropna(subset=["value"], inplace=True)

        ## flip the two columns to rows 
        star_file = star_file.T.reset_index(drop=True)
        ## make the first row the column names
        star_file.rename(columns=star_file.iloc[0], inplace=True)
        ## drop the duplicate column name row 
        star_file.drop(star_file.index[0], inplace=True)

        star_file.columns = star_file.columns.str.lower()

        ## pull the sampleid from the file name and add it as a new column  
        file_name = basename(file)
        sample = re.sub(".Log.final.out", "", file_name)
        star_file["sampleid"] = sample

        star_out_list.append(star_file)

    all_star_res = pd.concat(star_out_list, ignore_index=True)
    all_star_res.to_csv(out_file, sep="\t") 


test_star_fp = "/Users/apgarm/projects/pi_projects/rpelanda/bulkRNAseq_humanized_mice_autoreactive_vs_nonAutoreactive_bcells_01092025/comp_genome_redo/star_alignment/*.Log.final.out"


##wrangle_star_log(file_pattern=test_star_fp,
                 ##out_file="test_star_file.tsv")

def comb_filepaths(filepath1,
                   filepath2):
    return(pj(filepath1, filepath2))


def make_fp_dict(metadata_df,
                 dataset_dir):
    raw_sample_dict_out = {}

    ## pulling the forward/reverse files for each sampleid an putting them in a dictionary
    for sample in metadata_df['sampleid']:
        sample_filt_meta = metadata_df.loc[metadata_df['sampleid'] == sample]
        ## need to have .iloc[0] to just pull character from column 
        raw_file_list = [sample_filt_meta['forward'].iloc[0], sample_filt_meta['reverse'].iloc[0]]
        ## raw files with directory!
        ## adding the raw_seq_in directory to the filepaths here to make my life easier
        proc_rawFile_list = [comb_filepaths(dataset_dir, filepath) for filepath in raw_file_list]
        raw_sample_dict_out.update({sample: proc_rawFile_list})

    return(raw_sample_dict_out)


def specified_strandedness(metadata_df,
                           sampleid):
    if 'strandedness' in metadata_df.columns:
        sample_filt_meta = metadata_df.loc[metadata_df['sampleid'] == sampleid]
        rsem_strandedness = sample_filt_meta['strandedness']
    else:
        rsem_strandedness=""
    return(rsem_strandedness)


## create test sample metadata file
##test_dict = {'sampleid': ['test'],
##             'forward': ['test_R1.fq.gz'],
##             'reverse': ['test_R2.fq.gz']}
##test_df = pd.DataFrame(test_dict)
##
##test_df.to_csv("test_samp_metadata.tsv", sep="\t")

fake_profile_path = "/Users/apgarm/fake_profile.yml"
fake_profile_dir = "/Users/apgarm/fake_profile_dir"

def move_workflow_profile(user_profile_path):
    ## make sure profile path exists and is a file not a directory
    if exists(user_profile_path) and not os.path.isdir(user_profile_path):
        print("Setting up user profile...")
        ## come up with multiple profiles and let user choose the best one for their use case 
        fix_profile_config_name = ["cp", user_profile_path, pj("workflow/profiles", "default/config.yaml")]
        subprocess.run(fix_profile_config_name) 
        print("Setup successful!")
    else:
        print("The specified profile path does not exist or is a directory, please check that your input is a .yaml file and try again.")
        exit()

##move_workflow_profile(user_profile_path=fake_profile_path)


## testing out creating run config and adding to multiqc report comment
def assemble_report_config():
    config_text = """Raw Sequence Directory:             /path/to/raw_seqs
    Metadata File:                      /path/to/metadata.tsv
    Snakemake Profile:                  /path/to/profile/config.yaml
    Genome FASTA:                       /path/to/genome.fasta
    Annotation File (GTF):              /path/to/annotations.gtf
    Annotation Gene Name Column:        gene_name
    Raw Sequence Read Length:           0
    Picard RefFlat File:                /path/to/refflat.txt
    Picard Ribosomal Intervals File:    /path/to/ribo_intervals.list
    MultiQC Config File:                tmp.brat/multiqc_config.yaml
    RSEM Quantification:                True
    Additional Cutadapt Parameters:     --nextseq-trim=20
    Additional STAR Parameters:         --outFilterMultimapNmax 20
    Additional RSEM Parameters:         --seed 12345
    Latency Wait (secs):                10
    Deployment Method:                  singularity
    Dry Run:                            True"""

    print(f"""Running with Config:\n    {config_text}""")

    config_dict = {"report_comment": f"you ran this workflow with:\n\n <pre><code>\n{config_text}\n</code></pre>\n"}
    yaml = YAML()
    yaml.preserve_quotes = True
    yaml.default_flow_style = False
    yaml.width = 10000
    with open("test_multiqc.yaml", "w") as f:
        yaml.dump(config_dict, f)
    command = ["cat", "test_multiqc.yaml", ">>", "final_test.yml"]
    subprocess.run(" ".join(command), shell=True)
    

assemble_report_config()