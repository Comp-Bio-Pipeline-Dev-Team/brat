## sub-workflow to generate the overall multiqc report for all analysis run 

## take previous multiqc outputs out of rule since they may mess up the report??
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
        softwareLogs = SOFTWARE_LOG_DIR,
        outDir = OUT_DIR_NAME,
        multiqcFilename = "brat_run_report.html",
        multiqc_config = MULTIQC_CONFIG_PATH
    shell:
        """
        ## new way (a lot less code) - create a yaml file with software versions for multiqc to read directly
        cat {params.softwareLogs}/*.log | sort | uniq > {params.outDir}/brat_mqc_versions.yml

        multiqc {input.moduleDirs} {input.fastqScreenDirs} {params.outDir}/brat_mqc_versions.yml --config {params.multiqc_config} -o {params.outDir} --filename {params.multiqcFilename}
        """