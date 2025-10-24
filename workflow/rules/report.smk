## sub-workflow to generate the overall multiqc report for all analysis run 

rule create_report:
    input:
        moduleDirs = expand(pj(OUT_DIR_NAME, "{module}/{sample}/"),
                            module=MODULE_LIST,
                            sample=SAMPLE_LIST),
        multiqcDirs = [pj(OUT_DIR_NAME, "pretrimming_multiqc_report_data"),
                       pj(OUT_DIR_NAME, "posttrimming_multiqc_report_data")],
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
        multiqc_config = MULTIQC_CONFIG_PATH ## workflow/config_files/multiqc_config.yaml
    shell:
        """
        ## pull all unique entries from log files to create software used log for multiqc
        cat {params.softwareLogs}/*.log | sort | uniq > {params.outDir}/brat_software_used.log

        ## insert software used into multiqc config file via yq
        yq eval -i -y '.software_versions = loadstr("{params.outDir}/brat_software_used.log")' {params.multiqc_config}

        multiqc {input.moduleDirs} {input.multiqcDirs} --config {params.multiqc_config} -o {params.outDir} --filename {params.multiqcFilename}
        """