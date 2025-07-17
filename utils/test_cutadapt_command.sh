#!/bin/bash

## had to add set +o pipefail; to beginning of command so that snakemake doesn't detect non-zero error codes and fail
## (snakemake was detecing 141 error codes for first two steps in pipe bc pipe is so long they complete before later steps do)
## added '<' to zcat command so its portable to macos, macos prefers gzcat or zcat < but not plain zcat (its dumb)
## can use to check the error codes of all components of a pipe: echo ${PIPESTATUS[@]}
#totalScores=$( set +o pipefail; zcat < rp_bulkRNAseq/NA242plus_NA_S11_L004_R1_001.fastq.gz | awk 'NR%4==0' | sed 100q | grep -o . | sort | uniq | wc -l )

#if [ ${totalScores} -eq 4 ];
#then
    # NovaSeq now bins qc values to 2, 12, 23, 37 (if using NovaSeq), also use --nextseq-trim= instead of -q since NovaSeq is 2 color chemistry
    #echo "Quality scores are binned to 4 values (your score = ${totalScores}), therefore, using 2-color chemistry; running --nextseq-trim instead of -q"
    ## snakemake can insert {params.whatever} inside the double quotes and its evaluated as the value (30)
    #chemistryFlag="--nextseq-trim=30"
#else
    #echo "Quality scores are not binned (your score = ${totalScores}), therefore use -q"
    #chemistryFlag="-q 30"
#fi

## this works!! bash variable needs to be referenced exactly as it is or else it won't be evaluated by cutadapt!
#cutadapt \
    #-a AGATCGGAAGAG \
    #-A AGATCGGAAGAG \
    #"${chemistryFlag}"\
    #--minimum-length=10 \
    #--pair-filter=any \
    #-o bulk_RNAseq_out/cutadapt/NA242plus_NA_S11_L004/NA242plus_NA_S11_L004_R1_trimmed.fastq.gz \
    #-p bulk_RNAseq_out/cutadapt/NA242plus_NA_S11_L004/NA242plus_NA_S11_L004_R2_trimmed.fastq.gz \
    #rp_bulkRNAseq/NA242plus_NA_S11_L004_R1_001.fastq.gz rp_bulkRNAseq/NA242plus_NA_S11_L004_R2_001.fastq.gz>bulk_RNAseq_out/cutadapt/NA242plus_NA_S11_L004/NA242plus_NA_S11_L004_cutadapt.log \
    #2>bulk_RNAseq_out/cutadapt/NA242plus_NA_S11_L004/NA242plus_NA_S11_L004_cutadapt.err


## count read lengths in raw seq files for STAR
##check if read length set to zero or do want to infer? can take user input parameter (greater than 0)
readLength=$( awk '{print $5}' bulk_RNAseq_out/pretrimming_multiqc_report_data/multiqc_general_stats.txt | tail -n+2 | head -n 1 )
overheadLength=$( echo "$((${readLength}-1))" )
echo "Sequence read length is: ${overheadLength}"

## calculation of genomeLength - remove headers, remove whitespace and new lines, then count number of characters (bases)
genomeLength=$(cat test_refs/hg38.p14.fa | grep -v "^>" | tr -d [:space:] | wc -m)
echo "The length of the provided genome is: ${genomeLength} bases"

## calculation of genomeSAindex: min(14, log2(GenomeLength)/2 - 1)
genome_indexSize=$(printf %.0f $(echo "((l(${genomeLength})/l(2))/2)-1" | bc -l))
if [[ ${genome_indexSize} -ge 14 ]]; then
    genome_indexSize=14
fi
echo "Using genome index size of: ${genome_indexSize}"