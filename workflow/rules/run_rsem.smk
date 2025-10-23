## wont want to run this rule unless its needed - config option/looking at star log files for %uniquely_mapped v %percent_multimapped v %unmapped thresholds for when its run
## rsem genome index generation (probs will be moved to a subworkflow but its here for now) 
## finally works using outside fasta/gtf files with above rule (yay!!)
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
    params:
        rsem_index_name = ALIGN_INDEX_NAME,
        ref_directory = REF_DIR,
        genome_annot_col = ANNOT_COL_NAME ## gene_id should be the default!! only changed based on user input!
    shell:
        """
        if [ {params.genome_annot_col} == "gene_id" ];
        then
            annotFile={input.annotFile}
            echo "Using gene_id column in .gtf file for sequence annotations!"
        else
            cp {input.annotFile} {params.ref_directory}/tmp.gtf
            sed -i -e 's/gene_id/tmp_gene_id/g' -e 's/{params.genome_annot_col}/gene_id/g' {params.ref_directory}/tmp.gtf
            annotFile={params.ref_directory}/tmp.gtf
            echo "Using {params.genome_annot_col} column in .gtf file for sequence annotations!"
        fi

        mkdir -p {output.rsem_index_dir}

        rsem-prepare-reference --gtf ${{annotFile}} {input.fastaFile} {output.rsem_index_dir}/{params.rsem_index_name}
        """


## rsem 
## need to include the genome name prefix for the rsem index files or else it loses its mind
## have to also provide the file prefix with the output directory path (I forgot this)
## need to test that the user specified strandedness actually works as desired!!!
## look for rsem flag that can designate which column to pull 
rule run_rsem_quantification:
    input:
        aligned_transcriptBam = pj(OUT_DIR_NAME, "star_alignment/{sample}/{sample}.Aligned.toTranscriptome.out.bam"),
        rsem_index_dir = "tmp.brat/reference_indices/rsem/",
        picard_metrics_file = pj(OUT_DIR_NAME, "picard/allSample_picard_metrics.tsv")
    output:
        rsem_out_dir = directory(pj(OUT_DIR_NAME, "rsem_quantification/{sample}/"))
    singularity:
        RSEM_SING
    conda:
        RSEM_CONDA
    log:
        software_log = pj(SOFTWARE_LOG_DIR, "{sample}.rsem.log")
    params:
        nthreads = 6, ## nthreads = cpus_per_task*2 (for now)
        current_sample = lambda wc: wc.get("sample"),
        user_strandedness = specified_strandedness(metadata_df=METADATA, sampleid=lambda wc: wc.get("sample")), ## need to check that this will pull the wildcard appropriately
        rsem_index_name = ALIGN_INDEX_NAME,
        user_added_rsemParams = EXTRA_RSEM_PARAMS
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
        ( echo -n "rsem: "; printf "\"%s\"\n" "$(rsem-calculate-expression --version)" ) > {log.software_log}

        mkdir -p {output.rsem_out_dir}

        rsem-calculate-expression --paired-end \
                                  --bam \
                                  -p {params.nthreads} \
                                  {params.user_added_rsemParams} \
                                  --forward-prob ${{strandedness}} \
                                  --time {input.aligned_transcriptBam} {input.rsem_index_dir}/{params.rsem_index_name} {output.rsem_out_dir}/{params.current_sample}
        
        """