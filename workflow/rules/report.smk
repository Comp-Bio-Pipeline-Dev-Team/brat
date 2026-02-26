#### NOTE: ALL GLOBAL VARIABLES USED IN RULES ARE DEFINED IN THE SNAKEFILE, NOT IN THIS RULE FILE, SO THEY CAN BE REFERENCED THROUGHOUT THE WORKFLOW ####
#### CHECK SNAKEFILE FOR DEFINITIONS OF ALL GLOBAL VARIABLES USED IN THIS RULE FILE (AND THE WORKFLOW IN GENERAL)                                    ####

#### REQUIRED: sub-workflow to generate the overall multiqc report for all analysis run ####

## global variables within this rule: SOFTWARE_LOG_DIR, OUT_DIR_NAME, MULTIQC_CONFIG_PATH
## take previous multiqc outputs out of rule since they mess up the report (multiqc gets confused)
rule create_report:
    input:
        moduleDirs = expand(pj(OUT_DIR_NAME, "{module}/{sample}/"),
                            module=MODULE_LIST,
                            sample=SAMPLE_LIST),
        fastqScreenDirs = expand(pj(OUT_DIR_NAME, "fastq_screen/{sample}/{read}/"),
                                 sample=SAMPLE_LIST,
                                 read=TRIMMED_READS)
    output: 
        outFile = pj(OUT_DIR_NAME, "brat_run_report.html")
    singularity:
        MULTIQC_SING
    conda:
        MULTIQC_CONDA
    params:
        multiqcFilename = "brat_run_report.html"
    shell:
        """
        ## new way (a lot less code) - create a yaml file with software versions for multiqc to read directly
        cat {SOFTWARE_LOG_DIR}/*.log | sort | uniq > {OUT_DIR_NAME}/brat_mqc_versions.yml

        multiqc {input.moduleDirs} {input.fastqScreenDirs} {OUT_DIR_NAME}/brat_mqc_versions.yml --config {MULTIQC_CONFIG_PATH} -o {OUT_DIR_NAME} --filename {params.multiqcFilename}
        """