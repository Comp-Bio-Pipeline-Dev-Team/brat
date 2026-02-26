#### NOTE: ALL GLOBAL VARIABLES USED IN RULES ARE DEFINED IN THE SNAKEFILE, NOT IN THIS RULE FILE, SO THEY CAN BE REFERENCED THROUGHOUT THE WORKFLOW ####
#### CHECK SNAKEFILE FOR DEFINITIONS OF ALL GLOBAL VARIABLES USED IN THIS RULE FILE (AND THE WORKFLOW IN GENERAL)                                    ####

#### REQUIRED: main subworkflow to for sequence quality control (fastqc/multiqc), trimming (cutadapt), possible contamination screening (fastq screen), ####
#### alignment (star), and gene count quantification (star/picard)                                                                                      ####


#### PRETRIMMING FASTQC ####
## global variables used in this rule: CORES
## could also do this instead of input function: inFiles = expand(pj(RAW_SEQ_IN, "{{sample}}_{read}.fastq.gz"), read=READS) 
rule run_pretrimming_fastqc:
    input:
        inFiles = pull_rawSeq_fps
    output:
        outDir = directory(pj(OUT_DIR_NAME, "pretrimming_fastqc/{sample}/"))
    singularity:
        FASTQC_SING
    conda:
        FASTQC_CONDA
    log:
        software_log = pj(SOFTWARE_LOG_DIR, "{sample}.fastqc1.log")
    shell:
        """
        echo "running fastqc"

        ## getting fastqc version for multiqc report
        ( echo -n "fastqc: "; printf "\"%s\"\n" "$(fastqc --version 2>&1 | sed 's/[^0-9.]//g')" ) > {log.software_log} 2>&1

        ## calculating number of threads to use
        if [ {CORES} == "None" ]; then nThreads=$( echo "{PROC_CMD} * 2" | bc -l ); else nThreads={CORES}; fi 
        echo "Requested {CORES} cores, using ${{nThreads}} threads for fastqc"

        mkdir {output.outDir}

        fastqc {input.inFiles} -o {output.outDir} --threads ${{nThreads}}
        """

#### PRETRIMMING MULTIQC ####
## global variables within this rule: OUT_DIR_NAME
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
        pretrimming_multiqc_stats = pj(OUT_DIR_NAME, "pretrimming_multiqc_report_data/multiqc_general_stats.txt"),
        outDir = directory(pj(OUT_DIR_NAME, "pretrimming_multiqc_report_data/"))
    singularity:
        MULTIQC_SING
    conda:
        MULTIQC_CONDA
    log:
        software_log = pj(SOFTWARE_LOG_DIR, "multiqc1.log")
    params:
        multiqcFilename = "pretrimming_multiqc_report.html"
    shell:
        """
        echo "running multiqc"

        ## getting multiqc version for multiqc report
        ( echo -n "multiqc: "; printf "\"%s\"\n" "$(multiqc --version 2>&1 | sed 's/[^0-9.]//g')" ) > {log.software_log} 2>&1

        multiqc {input.inDirs} --force -o {OUT_DIR_NAME} --filename {params.multiqcFilename}
        """

#### SEQUENCE TRIMMING (CUTADAPT) ####
## global variables within this rule: CORES, TRIM_ADAPTER, TRIM_QUAL_SCORE, MIN_READ_LEN, PROC_CMD, EXTRA_CUTADAPT_PARAMS
## had to add set +o pipefail; to beginning of totalScores command so that snakemake doesn't detect non-zero error codes and fail
## (snakemake was detecing 141 error codes for first two steps in pipe bc pipe is so long they complete before later steps do)
## added '<' to zcat command so its portable to macos, macos prefers gzcat or zcat < but not plain zcat (its dumb)
## can use to check the error codes of all components of a pipe: echo ${PIPESTATUS[@]}
## snakemake can insert {params.whatever} inside the double quotes of chemistryFlag variable and its evaluated as the value (30)
rule run_cutadapt:
    input:
        inFiles = pull_rawSeq_fps
    output:
        forwardTrimmed = pj(OUT_DIR_NAME, "cutadapt/{sample}/{sample}_R1_trimmed.fastq.gz"),
        reverseTrimmed = pj(OUT_DIR_NAME, "cutadapt/{sample}/{sample}_R2_trimmed.fastq.gz"),
        cutadapt_dir = directory(pj(OUT_DIR_NAME, "cutadapt/{sample}/"))
    singularity:
        CUTADAPT_SING
    conda:
        CUTADAPT_CONDA
    log:
        log = pj(OUT_DIR_NAME, "cutadapt/{sample}/{sample}_cutadapt.log"),
        error = pj(OUT_DIR_NAME, "cutadapt/{sample}/{sample}_cutadapt.err"),
        software_log = pj(SOFTWARE_LOG_DIR, "{sample}.cutadapt.log")
    shell:
        """
        totalScores=$( set +o pipefail; zcat < {input.inFiles[0]} | awk 'NR%4==0' | sed 100q | grep -o . | sort | uniq | wc -l )

        if [ ${{totalScores}} -eq 4 ];
        then
            # NovaSeq now bins qc values to 2, 12, 23, 37 (if using NovaSeq), also use --nextseq-trim= instead of -q since NovaSeq is 2 color chemistry
            echo "Quality scores are binned to 4 values (your score = ${{totalScores}}) therefore, using 2-color chemistry; running --nextseq-trim instead of -q"
            chemistryFlag="--nextseq-trim={TRIM_QUAL_SCORE}"
        else
            echo "Quality scores are not binned (your score = ${{totalScores}}), therefore use -q"
            chemistryFlag="-q {TRIM_QUAL_SCORE}"
        fi

        ## getting cutadapt version for multiqc report
        ( echo -n "cutadapt: "; printf "\"%s\"\n" "$(cutadapt --version 2>&1 | sed 's/[^0-9.]//g')" ) > {log.software_log} 2>&1

        ## calculating number of threads to use
        if [ {CORES} == "None" ]; then nThreads=$( echo "{PROC_CMD} * 2" | bc -l ); else nThreads={CORES}; fi
        echo "Requested {CORES} cores, using ${{nThreads}} threads for cutadapt"

        cutadapt -a {TRIM_ADAPTER} \
                 -A {TRIM_ADAPTER} \
                 "${{chemistryFlag}}" \
                 --minimum-length={MIN_READ_LEN} \
                 -j ${{nThreads}} \
                 --pair-filter=any \
                 -o {output.forwardTrimmed} \
                 -p {output.reverseTrimmed} {EXTRA_CUTADAPT_PARAMS} \
                 {input.inFiles[0]} {input.inFiles[1]}>{log.log} \
                 2>{log.error} 
        """

#### POSTTRIMMING FASTQC ####
## can reuse previous rules this way!
## can alter the input, output, params, etc directives, anything not defined here will be inherited from the original rule
## CANNOT alter the execution step (i.e. shell)
use rule run_pretrimming_fastqc as run_posttrimming_fastqc with:
    input:
        inFiles = [pj(OUT_DIR_NAME, "cutadapt/{sample}/{sample}_R1_trimmed.fastq.gz"), pj(OUT_DIR_NAME, "cutadapt/{sample}/{sample}_R2_trimmed.fastq.gz")]
    output:
        outDir = directory(pj(OUT_DIR_NAME, "posttrimming_fastqc/{sample}/"))
    log:
        software_log = "tmp.brat/software_used/{sample}.fastqc2.log"


#### POSTTRIMMING MULTIQC ####
## if you include a directive to overwrite here (i.e. params), you'll need to specify all arguments underneath it,
## even if they don't change
use rule run_pretrimming_multiqc as run_posttrimming_multiqc with:
    input:
        inDirs = expand(pj(OUT_DIR_NAME, "posttrimming_fastqc/{sample}/"),
                        sample=SAMPLE_LIST)
    output:
        outFile = pj(OUT_DIR_NAME, "posttrimming_multiqc_report.html"),
        outDir = directory(pj(OUT_DIR_NAME, "posttrimming_multiqc_report_data/"))
    log:
        software_log = pj(SOFTWARE_LOG_DIR, "multiqc2.log")
    params:
        multiqcFilename = "posttrimming_multiqc_report.html"


#### GENERATE FASTQ SCREEN CONFIG FILE ####
## global variables within this rule: FQSCREEN_GENOMES
## update: snakemake was mad ab list input so added list to params and took out input
## how to pull everything from the gene names column on the command line (i probs wont need this): awk -F , '{{print $1}}' fastq_screen_genomes.csv | sed '1d'
rule create_fastq_screen_config: 
    output:
        fqScreen_config = "tmp.brat/fastq_screen.conf.mod"
    singularity:
        UBUNTU_SING
    params:
        index_loc = "tmp.brat/reference_indices/bowtie2"
    shell:
        """
        ## adding call for bowtie2
        printf "BOWTIE2\t bowtie2\n" > {output.fqScreen_config}

        bashGenomes=({FQSCREEN_GENOMES})

        for genome in ${{bashGenomes[*]}};
        do
            printf "DATABASE\t ${{genome}}\t {params.index_loc}/${{genome}}/${{genome}}\n" >> {output.fqScreen_config}
        done
        """


#### BOWTIE2 INDICES FOR FASTQ SCREEN ####
## this is based on the genome names and links provided in the csv file specified by the user (FQSCREEN_FILE)
## had to remove carriage return characters from links
## bowtie2 version isnt pulled correctly when pipeline is run with conda in alpine (due to perl being mad about the conda prefix including "@xsede.org")
rule generate_fastq_screen_index:
    input:
        linkFile = FQSCREEN_FILE
    output:
        index = directory("tmp.brat/reference_indices/bowtie2/{fq_screen_genome}/")
    singularity:
        FQSCREEN_SING
    conda:
        FQSCREEN_CONDA
    params:
        genome_name = lambda wc: wc.get("fq_screen_genome"),
        genome_fasta_dir = "tmp.brat/reference_indices/fastq_screen_fnas/",
        index_name = "tmp.brat/reference_indices/bowtie2/{fq_screen_genome}/{fq_screen_genome}"
    log:
        software_log = pj(SOFTWARE_LOG_DIR, "{fq_screen_genome}.bowtie2.log")
    shell:
        """
        ## get appropriate link 
        link=$( grep {params.genome_name} {input.linkFile} | sed 's/^[^,]*,//' | tr -d '\r' )

        ## and download it to specified directory via wget 
        wget -P {params.genome_fasta_dir} ${{link}}

        ## pulling off name of file 
        fileName=$( basename ${{link}} )

        ## pulling bowtie2 version
        ## im not using sed for this one bc it was a bit more complicated to pull the actual version off the first line
        ( echo -n "bowtie2: "; printf "\"%s\"\n" "$(bowtie2 --version 2>&1 | grep "version" | head -n 1 | awk -F ' ' '{{print $NF}}')" ) > {log.software_log} 2>&1

        ## make index output directory or else bowtie2 will fail 
        mkdir -p {output.index}

        ## bowtie2 command to build index from downloaded fasta
        bowtie2-build -f {params.genome_fasta_dir}${{fileName}} \
                         {params.index_name}
        """


#### CONTAMINANT IDENTIFICATION (FASTQ SCREEN) ####
## global variables within this rule: FQSCREEN_SUBSET_NUM, CORES, PROC_CMD
## every single trimmed fastq file needs to go into this rule separately 
rule run_fastq_screen:
    input:
        builtIndices = expand("tmp.brat/reference_indices/bowtie2/{fq_screen_genome}/",
                               fq_screen_genome=FQSCREEN_GENOMES),
        fqScreen_config = "tmp.brat/fastq_screen.conf.mod",
        trimmedReads = pj(OUT_DIR_NAME, "cutadapt/{sample}/{sample}_{read}_trimmed.fastq.gz")
    output:
        outDir = directory(pj(OUT_DIR_NAME, "fastq_screen/{sample}/{read}/"))
    singularity:
        FQSCREEN_SING
    conda:
        FQSCREEN_CONDA
    log:
        software_log = pj(SOFTWARE_LOG_DIR, "{sample}.{read}.fastq_screen.log")
    shell:
        """
        ## pulling fastq screen version for report
        ( echo -n "fastq screen: "; printf "\"%s\"\n" "$(fastq_screen --version 2>&1 | sed 's/[^0-9.]//g')" ) > {log.software_log} 2>&1

        ## calculating number of threads to use
        if [ {CORES} == "None" ]; then nThreads=$( echo "{PROC_CMD} * 2" | bc -l ); else nThreads={CORES}; fi
        echo "Requested {CORES} cores, using ${{nThreads}} threads for fastq screen"

        fastq_screen --outdir {output.outDir} \
                     --subset {FQSCREEN_SUBSET_NUM} \
                     --conf {input.fqScreen_config} \
                     --threads ${{nThreads}} \
                     --nohits \
                     {input.trimmedReads}
        """


#### COPY USER INPUT REFERENCE FILES TO WORKING DIRECTORY ####
## global variables within this rule: REF_DIR, NEW_FASTA_PATH, NEW_ANNOT_PATH, NEW_PICARD_REFFLAT, NEW_PICARD_RRNA_INTERVAL
## copy input fasta/gtf files into the working directory so bind points dont break
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
    shell:
        """
        mkdir -p {REF_DIR}

        refFileList=({input.fastaFile}
                     {input.annotFile}
                     {input.refFlat}
                     {input.riboIntList}
                    )

        for file in ${{refFileList[*]}};
        do
            fileName=$( basename ${{file}} )

            ## trying rsync with this instead of cp to see if it works better
            rsync -av --progress ${{file}} {REF_DIR}/${{fileName}}

            echo "Copying ${{fileName}} to {REF_DIR}/!"
        done

        """


#### GENERATE STAR INDEX FOR INPUT GENOME ####
## global variables within this rule: READ_LENGTH, CORES, PROC_CMD, STAR_FEATURE_TYPE
## increased memory to check if it does okay - needed 50GBs of memory 
## check if read length set to zero or do want to infer? can take user input parameter (greater than 0)
rule generate_star_index:
    input:
        pretrimming_multiqc_report = pj(OUT_DIR_NAME, "pretrimming_multiqc_report_data/multiqc_general_stats.txt"),
        fastaFile = NEW_FASTA_PATH,
        annotFile = NEW_ANNOT_PATH
    output:
        star_index_dir = directory("tmp.brat/reference_indices/star/")
    singularity:
        STAR_SING
    conda:
        STAR_CONDA
    shell:
        """
        if [ {READ_LENGTH} -eq 0 ];
        then
            ## count read lengths in raw seq files for STAR
            readLength=$( awk '{{print $5}}' {input.pretrimming_multiqc_report} | tail -n+2 | head -n 1 )
            overheadLength=$( echo "$((${{readLength}}-1))" )
            echo "Inferred sequence read length is: ${{overheadLength}}"
        else
            overheadLength={READ_LENGTH}
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

        ## calculating number of threads to use
        if [ {CORES} == "None" ]; then nThreads=$( echo "{PROC_CMD} * 1.5" | bc -l ); else nThreads={CORES}; fi
        echo "Requested {CORES} cores, using ${{nThreads}} threads for STAR genome generation"

        ## making output directory
        mkdir -p {output.star_index_dir}

        STAR --runThreadN ${{nThreads}} \
             --runMode genomeGenerate \
             --genomeDir {output.star_index_dir} \
             --genomeFastaFiles {input.fastaFile} \
             --sjdbGTFfile {input.annotFile} \
             --sjdbGTFfeatureExon {STAR_FEATURE_TYPE} \
             --sjdbOverhang ${{overheadLength}} \
             --genomeSAindexNbases ${{genome_indexSize}}

        """

#### SEQUENCE ALIGNMENT (STAR) ####
## global variables within this rule: READ_LENGTH, CORES, PROC_CMD, ANNOT_COL_NAME, EXTRA_STAR_PARAMS
## had to put user added star params on a separate line so that snakemake would parse the commands correctly - was excluding the first command before
## doubling the nthreads per cpu caused out of memory errors within the first minute of the job running (1.5 is fine, wouldnt go above that)
rule run_star_alignment:
    input:
        pretrimming_multiqc_report = pj(OUT_DIR_NAME, "pretrimming_multiqc_report_data/multiqc_general_stats.txt"),
        star_genome_index_dir = "tmp.brat/reference_indices/star/",
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
    log:
        software_log = pj(SOFTWARE_LOG_DIR, "{sample}.star.log")
    params:
        star_sample_prefix = pj(OUT_DIR_NAME, "star_alignment/{sample}/{sample}.")
    shell:
        """
        ## have to redo read length calculation, maybe not the best? but idk what else to do 
        if [ {READ_LENGTH} -eq 0 ];
        then
            ## count read lengths in raw seq files for STAR
            readLength=$( awk '{{print $5}}' {input.pretrimming_multiqc_report} | tail -n+2 | head -n 1 )
            overheadLength=$( echo "$((${{readLength}}-1))" )
            echo "Inferred sequence read length is: ${{overheadLength}}"
        else
            overheadLength={READ_LENGTH}
            echo "User input sequence read length is: ${{overheadLength}}"
        fi 

        ## getting star version for multiqc report
        ( echo -n "star: "; printf "\"%s\"\n" "$(STAR --version 2>&1 | sed 's/[^0-9.]//g')" ) > {log.software_log} 2>&1

        ## calculating number of threads to use
        if [ {CORES} == "None" ]; then nThreads=$( echo "{PROC_CMD} * 1.5" | bc -l ); else nThreads={CORES}; fi
        echo "Requested {CORES} cores, using ${{nThreads}} threads for STAR alignment"

        mkdir -p {output.star_out_dir}

        STAR --runMode alignReads \
             --runThreadN ${{nThreads}} \
             --genomeDir {input.star_genome_index_dir} \
             --readFilesCommand zcat \
             --sjdbGTFfile {input.annotFile} \
             --sjdbGTFtagExonParentGene {ANNOT_COL_NAME} \
             --readFilesIn {input.forwardTrimmed} {input.reverseTrimmed} \
             --outFileNamePrefix {params.star_sample_prefix} \
             --outSAMtype BAM SortedByCoordinate \
             --quantMode TranscriptomeSAM GeneCounts \
             --twopassMode Basic \
             --sjdbOverhang ${{overheadLength}} \
             {EXTRA_STAR_PARAMS}
        """


#### GET RNA SEQ METRICS (PICARD) ####
## global variables within this rule: PICARD_CMD
## picard jar file is included in singularity image but not needed for conda version of picard so global variable allows for flexibility between deployment methods
rule run_picard_collect_rna_seq:
    input:
        aligned_coordBam = pj(OUT_DIR_NAME, "star_alignment/{sample}/{sample}.Aligned.sortedByCoord.out.bam"),
        refFlat_file = NEW_PICARD_REFFLAT,
        riboIntList_file = NEW_PICARD_RRNA_INTERVAL
    output:
        collect_rnaSeq_file = pj(OUT_DIR_NAME, "picard_collectRnaSeq/{sample}/{sample}.picard.metrics.txt"),
        collect_rnaSeq_dir = directory(pj(OUT_DIR_NAME, "picard_collectRnaSeq/{sample}/"))
    singularity:
        PICARD_SING
    conda:
        PICARD_CONDA
    log:
        software_log = pj(SOFTWARE_LOG_DIR, "{sample}.picard.log")
    params:
        sample_out_dir = pj(OUT_DIR_NAME, "picard_collectRnaSeq/{sample}/"),
        strandedness = "NONE"
    shell:
        """
        ## get picard version for multiqc report
        ( echo -n "picard: "; printf "\"%s\"\n" "$({params.command} CollectRnaSeqMetrics --version 2>&1 | sed 's/[^0-9.]//g')" ) > {log.software_log} 2>&1

        mkdir -p {output.collect_rnaSeq_dir}

        {PICARD_CMD} CollectRnaSeqMetrics \
             I={input.aligned_coordBam} \
             O={output.collect_rnaSeq_file} \
             REF_FLAT={input.refFlat_file} \
             STRAND={params.strandedness} \
             RIBOSOMAL_INTERVALS={input.riboIntList_file}
        """


#### GET INSERT SIZE METRICS (PICARD) ####
## global variables within this rule: PICARD_CMD
rule run_picard_collect_insert_size:
    input:
        aligned_coordBam = pj(OUT_DIR_NAME, "star_alignment/{sample}/{sample}.Aligned.sortedByCoord.out.bam")
    output:
        insertSize_file = pj(OUT_DIR_NAME, "picard_collectInsertSize/{sample}/{sample}.picard.insertSize.txt"),
        insertSize_histogram = pj(OUT_DIR_NAME, "picard_collectInsertSize/{sample}/{sample}.picard.insertSize_histogram.pdf"),
        insertSize_dir = directory(pj(OUT_DIR_NAME, "picard_collectInsertSize/{sample}/"))
    singularity:
        PICARD_SING
    conda:
        PICARD_CONDA
    params:
        sample_out_dir = pj(OUT_DIR_NAME, "picard_collectInsertSize/{sample}/")
    shell:
        """
        mkdir -p {output.insertSize_dir}

        {PICARD_CMD} CollectInsertSizeMetrics \
             I={input.aligned_coordBam} \
             O={output.insertSize_file} \
             H={output.insertSize_histogram}
        """


#### CALCULATE STRANDEDNESS FROM PICARD METRICS AND WRANGLE STAR OUTPUTS INTO ONE FILE ####
## calculate strandedness from picard collect rna seq metrics output BEFORE rsem!!
## add automated strand interpretation column in picard collect rna seq metrics table (sanity checking)
## all picard for all samples should complete before we put them all into rule to wrangle picard results into same df and save as .csv
## dont need to import needed libraries in the run directive (for now)
rule calculate_strandedness:
    input:
        collect_rnaSeq_files = expand(pj(OUT_DIR_NAME, "picard_collectRnaSeq/{sample}/{sample}.picard.metrics.txt"),
                                      sample=SAMPLE_LIST),
        collect_insertSize_files = expand(pj(OUT_DIR_NAME, "picard_collectInsertSize/{sample}/{sample}.picard.insertSize.txt"),
                                          sample=SAMPLE_LIST),
        star_log_files = expand(pj(OUT_DIR_NAME, "star_alignment/{sample}/{sample}.Log.final.out"),
                                sample=SAMPLE_LIST)
    output:
        allSample_picard_strandedness = pj(OUT_DIR_NAME, "picard_collectRnaSeq/allSample_picard_metrics.tsv"),
        allSample_picard_insertSize = pj(OUT_DIR_NAME, "picard_collectInsertSize/allSample_picard_insertSize.tsv"),
        allSample_star_logs = pj(OUT_DIR_NAME, "star_alignment/allSample_star_logs.tsv")
    params:
        picard_metrics_pattern = pj(OUT_DIR_NAME, "picard_collectRnaSeq/*/*.picard.metrics.txt"),
        picard_insertSize_pattern = pj(OUT_DIR_NAME, "picard_collectInsertSize/*/*.picard.insertSize.txt"),
        star_logs_pattern = pj(OUT_DIR_NAME, "star_alignment/*/*.Log.final.out")
    run:
        picard_calculate_strandedness(file_pattern=params.picard_metrics_pattern,
                                      out_file=output.allSample_picard_strandedness)
        concat_picard_insert_size(file_pattern=params.picard_insertSize_pattern,
                                  out_file=output.allSample_picard_insertSize)
        concat_star_log(file_pattern=params.star_logs_pattern,
                        out_file=output.allSample_star_logs)