#### NOTE: ALL GLOBAL VARIABLES USED IN RULES ARE DEFINED IN THE SNAKEFILE, NOT IN THIS RULE FILE, SO THEY CAN BE REFERENCED THROUGHOUT THE WORKFLOW ####
#### CHECK SNAKEFILE FOR DEFINITIONS OF ALL GLOBAL VARIABLES USED IN THIS RULE FILE (AND THE WORKFLOW IN GENERAL)                                    ####

#### OPTIONAL: subworkflow to generate gene and isoform level expression estimates with RSEM if user desires ####

#### GENERATE RSEM INDEX ####
## global variables within this rule: ALIGN_INDEX_NAME, REF_DIR, ANNOT_COL_NAME
## wont want to run this rule unless its needed - user param option to run rsem if gene counts from star are not desired (default will be False)
rule generate_rsem_index:
    input:
        fastaFile = NEW_FASTA_PATH,
        annotFile = NEW_ANNOT_PATH
    output:
        rsem_index_dir = directory("tmp.brat/reference_indices/rsem/")
    singularity:
        RSEM_SING
    conda:
        RSEM_CONDA
    shell:
        """
        if [ {ANNOT_COL_NAME} == "gene_id" ];
        then
            annotFile={input.annotFile}
            echo "Using gene_id column in .gtf file for sequence annotations!"
        else
            cp {input.annotFile} {REF_DIR}/tmp.gtf
            sed -i -e 's/gene_id/tmp_gene_id/g' -e 's/{ANNOT_COL_NAME}/gene_id/g' {REF_DIR}/tmp.gtf
            annotFile={REF_DIR}/tmp.gtf
            echo "Using {ANNOT_COL_NAME} column in .gtf file for sequence annotations!"
        fi

        mkdir -p {output.rsem_index_dir}

        rsem-prepare-reference --gtf ${{annotFile}} {input.fastaFile} {output.rsem_index_dir}/{ALIGN_INDEX_NAME}
        """


#### RUN RSEM QUANTIFICATION ####
## global variables within this rule: PROC_CMD, CORES, ALIGN_INDEX_NAME, EXTRA_RSEM_PARAMS
## need to include the genome name prefix for the rsem index files or else it loses its mind
## have to also provide the file prefix with the output directory path (I forgot this)
rule run_rsem_quantification:
    input:
        aligned_transcriptBam = pj(OUT_DIR_NAME, "star_alignment/{sample}/{sample}.Aligned.toTranscriptome.out.bam"),
        rsem_index_dir = "tmp.brat/reference_indices/rsem/",
        picard_metrics_file = pj(OUT_DIR_NAME, "picard_collectRnaSeq/allSample_picard_metrics.tsv")
    output:
        rsem_out_dir = directory(pj(OUT_DIR_NAME, "rsem_quantification/{sample}/"))
    singularity:
        RSEM_SING
    conda:
        RSEM_CONDA
    log:
        software_log = pj(SOFTWARE_LOG_DIR, "{sample}.rsem.log")
    params:
        current_sample = lambda wc: wc.get("sample"),
        user_strandedness = specified_strandedness(metadata_df=METADATA, sampleid=lambda wc: wc.get("sample"))
    shell:
        """ 
        isEmpty={params.user_strandedness}

        if [[ -z "${{isEmpty}}" ]];
        then
            ## use wildcard to pull the line for each sampleid and the value from the last column (rsem strandedness will always be in the last column)
            strandedness=$( grep "{params.current_sample}" {input.picard_metrics_file} | awk -F'\t' '{{print $NF}}' )

            ## grab second to last column value to print the strandedness for user understanding
            readableStrandedness=$( grep "{params.current_sample}" {input.picard_metrics_file} | awk -F'\t' '{{print $(NF-1)}}' )
            echo "Calculated strandedness of {params.current_sample} is: ${{readableStrandedness}}"
        else
            strandedness={params.user_strandedness}
            echo "User specified strandedness for {params.current_sample} is: ${{strandedness}}"
        fi
        
        ## getting rsem version for multiqc report
        ( echo -n "rsem: "; printf "\"%s\"\n" "$(rsem-calculate-expression --version 2>&1 | sed 's/[^0-9.]//g')" ) > {log.software_log} 2>&1

        ## calculating number of threads to use
        if [ {CORES} == "None" ]; then nThreads=$( echo "{PROC_CMD} * 2" | bc -l ); else nThreads={CORES}; fi
        echo "Requested {CORES} cores, using ${{nThreads}} threads for RSEM quantification"

        mkdir -p {output.rsem_out_dir}

        rsem-calculate-expression --paired-end \
                                  --bam \
                                  -p ${{nThreads}} \
                                  {EXTRA_RSEM_PARAMS} \
                                  --forward-prob ${{strandedness}} \
                                  --time {input.aligned_transcriptBam} {input.rsem_index_dir}/{ALIGN_INDEX_NAME} {output.rsem_out_dir}/{params.current_sample}
        
        """