## attempting to put together overall functions to use in my pipeline
from os.path import join as pj
import pandas as pd
from os.path import basename
import re
import glob

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

## searchs through the lines of a file looking for a specific pattern and returns the line number of the first occurrence
## of that pattern
## to be used in the picard output file wrangling functions (below)!! 
def get_lines_to_skip(file_path, 
                      line_pattern):
    with open(file_path, 'r') as raw_file:
        for line_num, line in enumerate(raw_file, start=1):
            if re.search(line_pattern, line):
                return line_num
    raise ValueError(f"Pattern '{line_pattern}' not found in file '{file_path}'")


def picard_calculate_strandedness(file_pattern,
                                  out_file):
    picard_out_dict = {}

    picard_collectRNAseq_files = [f for f in glob.glob(file_pattern)]
    search_line_pattern = r"^## METRICS"

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

        ## print out the calculation here too!
        if (picard_file["PCT_R1_TRANSCRIPT_STRAND_READS"] >= 0.70).all():
            picard_file["strandedness"] = "R1"
            picard_file["rsem_strand_key"] = 1
        elif (picard_file["PCT_R2_TRANSCRIPT_STRAND_READS"] >= 0.70).all():
            picard_file["strandedness"] = "R2"
            picard_file["rsem_strand_key"] = 0
        else:
            picard_file["strandedness"] = "unstranded"
            picard_file["rsem_strand_key"] = 0.5

        picard_out_dict.update({f"{sample}": picard_file})

    all_picard_res = pd.concat(picard_out_dict, ignore_index=True)
    all_picard_res.to_csv(out_file, sep="\t")


def concat_picard_insert_size(file_pattern,
                              out_file):
    picard_out_dict = {}

    picard_collectInsertSize_files = [f for f in glob.glob(file_pattern)]
    search_line_pattern = r"^## METRICS"

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

        picard_out_dict.update({f"{sample}": picard_file})

    all_picard_res = pd.concat(picard_out_dict, ignore_index=True)
    all_picard_res.to_csv(out_file, sep="\t")


def concat_star_log(file_pattern,
                    out_file):
    star_log_files = [f for f in glob.glob(file_pattern)]
    star_out_dict = {}

    for file in star_log_files:
        star_file = pd.read_csv(file, 
                        sep="\t",
                        names=["key", "value"])
        ## remove the pipe symbol from the first column 
        star_file["key"] = star_file["key"].str.replace(" |", "")
        ## pull excess white space off of the values in the first column
        star_file["key"] = star_file["key"].str.strip()
        ## replace spaces and / with underscores
        star_file["key"] = star_file["key"].str.replace(r"[\s/]", "_", regex=True)
        ## replace commas, colons, and hyphens with nothing
        star_file["key"] = star_file["key"].str.replace(r"[,:-]", "", regex=True)
        ## replace % symbol with word pct
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

        star_out_dict.update({f"{sample}": star_file})

    all_star_res = pd.concat(star_out_dict, ignore_index=True)
    all_star_res.to_csv(out_file, sep="\t")


def specified_strandedness(metadata_df,
                           sampleid):
    if 'strandedness' in metadata_df.columns:
        sample_filt_meta = metadata_df.loc[metadata_df['sampleid'] == sampleid]
        rsem_strandedness = sample_filt_meta['strandedness']
    else:
        rsem_strandedness=""
    return(rsem_strandedness)






    

