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


def picard_calculate_strandedness(file_pattern,
                                  out_file):
    picard_out_dict = {}

    picard_collectRNAseq_files = [f for f in glob.glob(file_pattern)]

    for file in picard_collectRNAseq_files:
        picard_file = pd.read_csv(file, 
                                sep="\t",
                                skiprows=6,
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

        picard_out_dict.update({f"{sample}": picard_file})

    all_picard_res = pd.concat(picard_out_dict, ignore_index=True)
    all_picard_res.to_csv(out_file)







    

