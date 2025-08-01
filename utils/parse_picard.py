import pandas as pd
from os.path import basename
import re
import glob


def wrangle_picard_output(file_pattern,
                          out_file):
    picard_collectRNAseq_files = [f for f in glob.glob(file_pattern)]

    picard_out_dict = {}

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

    all_picard_res.to_csv(out_file, sep="\t")


##wrangle_picard_output(file_pattern="bulk_RNAseq_out/picard/*/*.picard.metrics.txt",
                      ##out_file="test_picard_res.csv")


def wrangle_picard_insert_size(file_pattern,
                               out_file):
    picard_collectInsertSize_files = [f for f in glob.glob(file_pattern)]

    picard_out_dict = {}

    for file in picard_collectInsertSize_files:
        picard_file = pd.read_csv(file, 
                                  sep="\t",
                                  skiprows=6,
                                  nrows=1)

        file_name = basename(file)
        sample = re.sub(".picard.insertSize.txt", "", file_name)
        picard_file["sampleid"] = sample

        picard_out_dict.update({f"{sample}": picard_file})


    all_picard_res = pd.concat(picard_out_dict, ignore_index=True)

    all_picard_res.to_csv(out_file, sep="\t")


##wrangle_picard_insert_size(file_pattern="/Users/apgarm/projects/pi_projects/rpelanda/bulkRNAseq_humanized_mice_autoreactive_vs_nonAutoreactive_bcells_01092025/comp_genome_redo/picard/*.picard.insertSize.txt",
                           ##out_file="test_picard_insertSize.tsv")



def wrangle_star_log(file_pattern,
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
        star_file["key"] = star_file["key"].str.replace(r"[\s/]", "_", regex=True)
        star_file["key"] = star_file["key"].str.replace(r"[,:-]", "", regex=True)
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


test_star_fp = "/Users/apgarm/projects/pi_projects/rpelanda/bulkRNAseq_humanized_mice_autoreactive_vs_nonAutoreactive_bcells_01092025/comp_genome_redo/star_alignment/*.Log.final.out"


wrangle_star_log(file_pattern=test_star_fp,
                 out_file="test_star_file.tsv")


##star_file = pd.read_csv(test_star_fp, 
##                        sep="\t",
##                        names=["key", "value"])
#### remove the pipe symbol from the first column 
##star_file["key"] = star_file["key"].str.replace(" |", "")
#### pull excess white space off of the values in the first column
##star_file["key"] = star_file["key"].str.strip()
#### remove percent signs from the values in the second column 
##star_file["value"] = star_file["value"].str.replace("%", "")
#### drop na values from the second column 
##star_file.dropna(subset=["value"], inplace=True)
#### flip the two columns to rows 
##final_star_file = star_file.T.reset_index(drop=True)
#### make the first row the column names
##final_star_file.rename(columns=final_star_file.iloc[0], inplace=True)
#### drop the duplicate column name row 
##final_star_file.drop(final_star_file.index[0], inplace=True)
##
#### pull the sampleid from the file name and add it as a new column 
##file_name = basename(test_star_fp)
##sample = re.sub(".Log.final.out", "", file_name)
##final_star_file["sampleid"] = sample
##
##print(final_star_file)

