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

    all_picard_res.to_csv(out_file)




##wrangle_picard_output(file_pattern="bulk_RNAseq_out/picard/*/*.picard.metrics.txt",
                      ##out_file="test_picard_res.csv")


