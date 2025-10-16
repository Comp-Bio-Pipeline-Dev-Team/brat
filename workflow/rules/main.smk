##### SUBWORKFLOW ONE!! #####
## set up to run for all paired samples - in .csv file that tells you what pairs on, read in as a dict of the file names per sample
## output directory per samples
## can keep threads=2 and ntasks=1
##### madis notes: #####
## both examples under "input:" work to pull both forward/reverse reads in for each sample 
## could also do this: inFiles = expand(pj(RAW_SEQ_IN, "{{sample}}_{read}.fastq.gz"), read=READS) 
rule run_pretrimming_fastqc:
    input:
        inFiles = pull_rawSeq_fps
    output:
        outDir = directory(pj(OUT_DIR_NAME, "pretrimming_fastqc/{sample}/"))
    singularity:
        FASTQC_SING
    conda:
        FASTQC_CONDA
    params:
        threads = 2,
        software_log = SOFTWARE_LOG
    shell:
        """
        echo "running fastqc"

        ## getting fastqc version for multiqc report
        ( echo -n "fastqc: "; printf "\"%s\"\n" "$(fastqc --version)" ) >> {params.software_log}

        mkdir {output.outDir}

        fastqc {input.inFiles} -o {output.outDir} --threads {params.threads}
        """


## took '.' at end of multiqc command out bc i think its causing problems with the report generation - it was
## added "--force" to make multiqc overwrite directories that already exist bc snakemake likes to create the output directory before
## multiqc does so multiqc was naming it something different and snakemake was failing out bc it couldn't find the file (yay)
rule run_pretrimming_multiqc:
    input:
        ## tells snakemake to wait for all fastqc outputs to be made before starting this rule since all of them are required for this rule
        inDirs = expand(pj(OUT_DIR_NAME, "pretrimming_fastqc/{sample}/"),
                         sample=SAMPLE_LIST)
    output:
        outFile = pj(OUT_DIR_NAME, "pretrimming_multiqc_report.html"),
        pretrimming_multiqc_stats = pj(OUT_DIR_NAME, "pretrimming_multiqc_report_data/multiqc_general_stats.txt")
    singularity:
        MULTIQC_SING
    conda:
        MULTIQC_CONDA
    params:
        outDir = OUT_DIR_NAME, ## this might need the forward slash at the end again 
        multiqcFilename = "pretrimming_multiqc_report.html",
        software_log = SOFTWARE_LOG
    shell:
        """
        echo "running multiqc"

        ## getting multiqc version for multiqc report
        ( echo -n "multiqc: "; printf "\"%s\"\n" "$(multiqc --version)" ) >> {params.software_log}

        multiqc {input.inDirs} --force -o {params.outDir} --filename {params.multiqcFilename}
        """

## {params.inDir}*/


## had to add set +o pipefail; to beginning of totalScores command so that snakemake doesn't detect non-zero error codes and fail
## (snakemake was detecing 141 error codes for first two steps in pipe bc pipe is so long they complete before later steps do)
## added '<' to zcat command so its portable to macos, macos prefers gzcat or zcat < but not plain zcat (its dumb)
## can use to check the error codes of all components of a pipe: echo ${PIPESTATUS[@]}
## snakemake can insert {params.whatever} inside the double quotes of chemistryFlag variable and its evaluated as the value (30)
## add {params.user_added_cutadaptParams} to the command!
rule run_cutadapt:
    input:
        inFiles = pull_rawSeq_fps
    output:
        forwardTrimmed = pj(OUT_DIR_NAME, "cutadapt/{sample}/{sample}_R1_trimmed.fastq.gz"),
        reverseTrimmed = pj(OUT_DIR_NAME, "cutadapt/{sample}/{sample}_R2_trimmed.fastq.gz")
    singularity:
        CUTADAPT_SING
    conda:
        CUTADAPT_CONDA
    log:
        log = pj(OUT_DIR_NAME, "cutadapt/{sample}/{sample}_cutadapt.log"),
        error = pj(OUT_DIR_NAME, "cutadapt/{sample}/{sample}_cutadapt.err")
    params:
        adapter_threePrime = "AGATCGGAAGAG", ## do i need to put this in the config file?
        qualTrim = 30,
        minReadLen = 10,
        n_threads = 4, ## n_threads = cpus_per_task*2
        user_added_cutadaptParams = EXTRA_CUTADAPT_PARAMS,
        software_log = SOFTWARE_LOG
    shell:
        """
        totalScores=$( set +o pipefail; zcat < {input.inFiles[0]} | awk 'NR%4==0' | sed 100q | grep -o . | sort | uniq | wc -l )

        if [ ${{totalScores}} -eq 4 ];
        then
            # NovaSeq now bins qc values to 2, 12, 23, 37 (if using NovaSeq), also use --nextseq-trim= instead of -q since NovaSeq is 2 color chemistry
            echo "Quality scores are binned to 4 values (your score = ${{totalScores}}) therefore, using 2-color chemistry; running --nextseq-trim instead of -q"
            chemistryFlag="--nextseq-trim={params.qualTrim}"
        else
            echo "Quality scores are not binned (your score = ${{totalScores}}), therefore use -q"
            chemistryFlag="-q {params.qualTrim}"
        fi

        ## getting cutadapt version for multiqc report
        ( echo -n "cutadapt: "; printf "\"%s\"\n" "$(cutadapt --version)" ) >> {params.software_log}

        cutadapt -a {params.adapter_threePrime} \
                 -A {params.adapter_threePrime} \
                 "${{chemistryFlag}}" \
                 --minimum-length={params.minReadLen} \
                 -j {params.n_threads} \
                 --pair-filter=any \
                 -o {output.forwardTrimmed} \
                 -p {output.reverseTrimmed} {params.user_added_cutadaptParams} \
                 {input.inFiles[0]} {input.inFiles[1]}>{log.log} \
                 2>{log.error} 
        """

## can reuse previous rules this way!
## can alter the input, output, params, etc directives, anything not defined here will be inherited from the original rule
## CANNOT alter the execution step (i.e. shell)
## running fastqc on samples posttrimming via cutadapt
use rule run_pretrimming_fastqc as run_posttrimming_fastqc with:
    input:
        inFiles = [pj(OUT_DIR_NAME, "cutadapt/{sample}/{sample}_R1_trimmed.fastq.gz"), pj(OUT_DIR_NAME, "cutadapt/{sample}/{sample}_R2_trimmed.fastq.gz")]
    output:
        outDir = directory(pj(OUT_DIR_NAME, "posttrimming_fastqc/{sample}/"))


## running multiqc on samples posttrimming via cutadapt 
## if you include a directive to overwrite here (i.e. params), you'll need to specify all arguments underneath it,
## even if they don't change
use rule run_pretrimming_multiqc as run_posttrimming_multiqc with:
    input:
        inDirs = expand(pj(OUT_DIR_NAME, "posttrimming_fastqc/{sample}/"),
                        sample=SAMPLE_LIST)
    output:
        outFile = pj(OUT_DIR_NAME, "posttrimming_multiqc_report.html")
    params:
        outDir = OUT_DIR_NAME,
        multiqcFilename = "posttrimming_multiqc_report.html"



## copy input fasta/gtf files into the test_snake_w_slurm (or whatever working) directory so bind points dont break
## symlinks for individual files didnt work bc symlink names are copied over but not the actual files (i.e. containers break the symlinks)
rule copy_reference_files:
    input:
        fastaFile = ALIGN_FASTA,
        annotFile = ALIGN_ANNOT,
        refFlat = PICARD_REFFLAT,
        riboIntList = PICARD_RRNA_INTERVAL
    output:
        new_loc_fasta = NEW_FASTA_PATH,
        new_loc_annot = NEW_ANNOT_PATH,
        new_refFlat = NEW_PICARD_REFFLAT,
        new_riboIntList = NEW_PICARD_RRNA_INTERVAL
    singularity:
        UBUNTU_SING
    params:
        ref_directory = REF_DIR
    shell:
        """
        mkdir -p {params.ref_directory}

        refFileList=({input.fastaFile}
                     {input.annotFile}
                     {input.refFlat}
                     {input.riboIntList}
                    )

        for file in ${{refFileList[*]}};
        do
            fileName=$( basename ${{file}} )

            ## trying rsync with this instead of cp to see if it works better
            rsync -av --progress ${{file}} {params.ref_directory}/${{fileName}}

            echo "Copying ${{fileName}} to {params.ref_directory}/!"
        done

        """


## got errors about bc not being installed in container (it wasnt so thats fine) and out of memory error which was weird 
## increased memory to check if it does okay - needed 50GBs of memory 
## check if read length set to zero or do want to infer? can take user input parameter (greater than 0) - i cant remember why we were doing this
rule generate_star_index:
    input:
        pretrimming_multiqc_report = pj(OUT_DIR_NAME, "pretrimming_multiqc_report_data/multiqc_general_stats.txt"),
        fastaFile = NEW_FASTA_PATH,
        annotFile = NEW_ANNOT_PATH
    output:
        star_index_dir = directory("reference_indices/star/")
    singularity:
        STAR_SING
    conda:
        STAR_CONDA
    params:
        input_read_length = READ_LENGTH,
        n_threads = 22, ## n_threads = cpus_per_task*2
        feature_type = "exon"
    shell:
        """
        if [ {params.input_read_length} -eq 0 ];
        then
            ## count read lengths in raw seq files for STAR
            readLength=$( awk '{{print $5}}' {input.pretrimming_multiqc_report} | tail -n+2 | head -n 1 )
            overheadLength=$( echo "$((${{readLength}}-1))" )
            echo "Inferred sequence read length is: ${{overheadLength}}"
        else
            overheadLength={params.input_read_length}
            echo "User input sequence read length is: ${{overheadLength}}"
        fi

        ## calculation of genomeLength - remove headers, remove whitespace and new lines, then count number of characters (bases)
        genomeLength=$(cat {input.fastaFile} | grep -v "^>" | tr -d [:space:] | wc -m)
        echo "The length of the provided genome is: ${{genomeLength}} bases"

        ## calculation of genomeSAindex: min(14, log2(GenomeLength)/2 - 1)
        genome_indexSize=$(printf %.0f $(echo "((l(${{genomeLength}})/l(2))/2)-1" | bc -l))
        if [[ ${{genome_indexSize}} -ge 14 ]]; then
            genome_indexSize=14
        fi
        echo "Using genome index size of: ${{genome_indexSize}}"

        ## making output directory
        mkdir -p {output.star_index_dir}

        STAR --runThreadN {params.n_threads} \
             --runMode genomeGenerate \
             --genomeDir {output.star_index_dir} \
             --genomeFastaFiles {input.fastaFile} \
             --sjdbGTFfile {input.annotFile} \
             --sjdbGTFfeatureExon {params.feature_type} \
             --sjdbOverhang ${{overheadLength}} \
             --genomeSAindexNbases ${{genome_indexSize}}

        """

## had to put user added star params on a separate line so that snakemake would parse the commands correctly - was excluding the first command before
rule run_star_alignment:
    input:
        pretrimming_multiqc_report = pj(OUT_DIR_NAME, "pretrimming_multiqc_report_data/multiqc_general_stats.txt"),
        star_genome_index_dir = "reference_indices/star/",
        forwardTrimmed = pj(OUT_DIR_NAME, "cutadapt/{sample}/{sample}_R1_trimmed.fastq.gz"),
        reverseTrimmed = pj(OUT_DIR_NAME, "cutadapt/{sample}/{sample}_R2_trimmed.fastq.gz"),
        annotFile = NEW_ANNOT_PATH
    output:
        star_out_dir = directory(pj(OUT_DIR_NAME, "star_alignment/{sample}/")),
        star_out_coordFile = pj(OUT_DIR_NAME, "star_alignment/{sample}/{sample}.Aligned.sortedByCoord.out.bam"),
        star_out_transcriptFile = pj(OUT_DIR_NAME, "star_alignment/{sample}/{sample}.Aligned.toTranscriptome.out.bam"),
        star_out_finalLog = pj(OUT_DIR_NAME, "star_alignment/{sample}/{sample}.Log.final.out") 
    singularity:
        STAR_SING
    conda:
        STAR_CONDA
    params:
        input_read_length = READ_LENGTH,
        n_threads = 18, ## n_threads = cpus_per_task*1.5, madi note: doubling the nthreads per cpu caused out of memory errors within the first minute of the job running
        star_sample_prefix = pj(OUT_DIR_NAME, "star_alignment/{sample}/{sample}."), ## may need another decoy here bc snakemake takes the front slash off of the output fp above (now the file names look funky)
        user_added_starParams = EXTRA_STAR_PARAMS,
        software_log = SOFTWARE_LOG
    shell:
        """
        ## have to redo read length calculation, maybe not the best? but idk what else to do 
        if [ {params.input_read_length} -eq 0 ];
        then
            ## count read lengths in raw seq files for STAR
            readLength=$( awk '{{print $5}}' {input.pretrimming_multiqc_report} | tail -n+2 | head -n 1 )
            overheadLength=$( echo "$((${{readLength}}-1))" )
            echo "Inferred sequence read length is: ${{overheadLength}}"
        else
            overheadLength={params.input_read_length}
            echo "User input sequence read length is: ${{overheadLength}}"
        fi 

        ## getting star version for multiqc report
        ( echo -n "star: "; printf "\"%s\"\n" "$(STAR --version)" ) >> {params.software_log}

        mkdir -p {output.star_out_dir}

        STAR --runMode alignReads \
             --runThreadN {params.n_threads} \
             --genomeDir {input.star_genome_index_dir} \
             --readFilesCommand zcat \
             --sjdbGTFfile {input.annotFile} \
             --readFilesIn {input.forwardTrimmed} {input.reverseTrimmed} \
             --outFileNamePrefix {params.star_sample_prefix} \
             --outSAMtype BAM SortedByCoordinate \
             --quantMode TranscriptomeSAM GeneCounts \
             --twopassMode Basic \
             --sjdbOverhang ${{overheadLength}} \
             {params.user_added_starParams}
        """


## rule for picard collect rna seq metrics!
## this works
## picard.jar file needs to be local, NOT in the container, probs will download it in the move_reference_files rule or have it prewrapped in the pipeline
rule run_picard_collect_rna_seq:
    input:
        aligned_coordBam = pj(OUT_DIR_NAME, "star_alignment/{sample}/{sample}.Aligned.sortedByCoord.out.bam"),
        refFlat_file = NEW_PICARD_REFFLAT,
        riboIntList_file = NEW_PICARD_RRNA_INTERVAL
    output:
        collect_rnaSeq_file = pj(OUT_DIR_NAME, "picard/{sample}/{sample}.picard.metrics.txt")
    singularity:
        PICARD_SING
    conda:
        PICARD_CONDA
    params:
        sample_out_dir = pj(OUT_DIR_NAME, "picard/{sample}/"),
        command = PICARD_CMD,
        jar_file = "/opt/picard-2.27.5/picard.jar", ## tell it to look in the container for this, change permissions of jar file in container (exec java -jar picard.jar in the container)
        strandedness = "NONE",
        software_log = SOFTWARE_LOG
    shell:
        """
        ## get picard version for multiqc report
        ( echo -n "picard: "; printf "\"%s\"\n" "$({params.command} CollectRnaSeqMetrics --version)" ) >> {params.software_log}

        mkdir -p {params.sample_out_dir}

        {params.command} CollectRnaSeqMetrics \
             I={input.aligned_coordBam} \
             O={output.collect_rnaSeq_file} \
             REF_FLAT={input.refFlat_file} \
             STRAND={params.strandedness} \
             RIBOSOMAL_INTERVALS={input.riboIntList_file}
        """



## rule for picard collect insert size metrics !
## this works
rule run_picard_collect_insert_size:
    input:
        aligned_coordBam = pj(OUT_DIR_NAME, "star_alignment/{sample}/{sample}.Aligned.sortedByCoord.out.bam")
    output:
        insertSize_file = pj(OUT_DIR_NAME, "picard/{sample}/{sample}.picard.insertSize.txt"),
        insertSize_histogram = pj(OUT_DIR_NAME, "picard/{sample}/{sample}.picard.insertSize_histogram.pdf")
    singularity:
        PICARD_SING
    conda:
        PICARD_CONDA
    params:
        sample_out_dir = pj(OUT_DIR_NAME, "picard/{sample}/"),
        command = PICARD_CMD,
        jar_file = "/opt/picard-2.27.5/picard.jar" ## this is in the container!!
    shell:
        """
        # check if outDir exists; if it doesn't it will make the output directory
        ## idk if i need to do this or not 
        if [ -d {params.sample_out_dir} ];
        then
            echo "Directory {params.sample_out_dir} exists. Proceeding with analysis."
        else
            echo "Directory {params.sample_out_dir} does not exist. Generating directory..."
            mkdir -p {params.sample_out_dir}
        fi

        {params.command} CollectInsertSizeMetrics \
             I={input.aligned_coordBam} \
             O={output.insertSize_file} \
             H={output.insertSize_histogram}
        """


## need a script here to calculate strandedness from picard collect rna seq metrics output BEFORE rsem!!
## write python function and put it in snake_functions.py to reference here 
## add automated strand interpretation column in picard collect rna seq metrics table (sanity checking)
## all picard for all samples should complete before we put them all into rule to wrangle picard results into same df and save as .csv
## use run: instead of shell: 
## do I need to import needed libraries in the run directive or will it be okay? - it wil be okay 
## i should probably name this rule something different but i cant think of anything right now
rule calculate_strandedness:
    input:
        collect_rnaSeq_files = expand(pj(OUT_DIR_NAME, "picard/{sample}/{sample}.picard.metrics.txt"),
                                      sample=SAMPLE_LIST),
        collect_insertSize_files = expand(pj(OUT_DIR_NAME, "picard/{sample}/{sample}.picard.insertSize.txt"),
                                          sample=SAMPLE_LIST),
        star_log_files = expand(pj(OUT_DIR_NAME, "star_alignment/{sample}/{sample}.Log.final.out"),
                                sample=SAMPLE_LIST)
    output:
        allSample_picard_strandedness = pj(OUT_DIR_NAME, "picard/allSample_picard_metrics.tsv"),
        allSample_picard_insertSize = pj(OUT_DIR_NAME, "picard/allSample_picard_insertSize.tsv"),
        allSample_star_logs = pj(OUT_DIR_NAME, "star_alignment/allSample_star_logs.tsv")
    params:
        picard_metrics_pattern = pj(OUT_DIR_NAME, "picard/*/*.picard.metrics.txt"),
        picard_insertSize_pattern = pj(OUT_DIR_NAME, "picard/*/*.picard.insertSize.txt"),
        star_logs_pattern = pj(OUT_DIR_NAME, "star_alignment/*/*.Log.final.out")
    run:
        picard_calculate_strandedness(file_pattern=params.picard_metrics_pattern,
                                      out_file=output.allSample_picard_strandedness)
        concat_picard_insert_size(file_pattern=params.picard_insertSize_pattern,
                                  out_file=output.allSample_picard_insertSize)
        concat_star_log(file_pattern=params.star_logs_pattern,
                        out_file=output.allSample_star_logs)