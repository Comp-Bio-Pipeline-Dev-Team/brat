## attempting to put together overall functions to use in my pipeline
from os.path import join as pj
import pandas as pd

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







    

