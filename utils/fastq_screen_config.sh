#!/bin/bash

## how to pull url for a specific genome from fastq_screen_genomes.csv
## genome will probably be a wildcard in snakefile
genome="acholeplasma_laidlawii_PG_8A"
link=$( grep ${genome} ../files_for_users/fastq_screen_genomes.csv | sed 's/^[^,]*,//' )

cleanedLink=$( echo ${link} | tr -d '\r' )

echo ${cleanedLink}

## and download it to specified directory via wget 
wget -P tmp.brat/reference_indices/fastq_screen_fnas/ ${link}

fileName=$( basename ${link} )

## pulling bowtie2 version
## im not using sed for this one bc it was a bit more complicated to pull the actual version off the first line
( echo -n "bowtie2: "; printf "\"%s\"\n" "$(bowtie2 --version 2>&1 | head -n 1 | awk -F ' ' '{{print $NF}}')" ) > {log.software_log} 2>&1

## bowtie2 command to build index from downloaded fasta
bowtie2-build -f tmp.brat/reference_indices/fastq_screen_fnas/${fileName} \
              tmp.brat/reference_indices/bowtie2/${genome}/${genome}



## assemble fastq screen config file
## how to pull everything from the gene names column on the command line (i probs wont need this)
## awk -F , '{print $1}' fastq_screen_genomes.csv | sed '1d'
genomeList=("acholeplasma_laidlawii_PG_8A"
              "mesomycoplasma_hyorhinis_MCLD"
              "metamycoplasma_hominis_ATCC_23114"
              "mycoplasmopsis_fermentans_M64"
             )

## adding call for bowtie2 (hopefully this works)
printf "BOWTIE2\t bowtie2\n" > brat.fastq_screen.conf.mod

for genome in ${genomeList[*]};
do
    printf "DATABASE\t ${genome}\t tmp.brat/reference_indices/bowtie2/${genome}/${genome}\n" >> brat.fastq_screen.conf.mod
done


## actual command to run fastq screen (in separate rule)
## i think each trimmed sequence goes in individually regardless of being R1/R2?
## might have to create separate subdirectories for R1/R2 so snakemake doesnt lose its mind
## may need to make all the bowtie2 indices an input to this rule so that snakemake waits for those to be built first 

## pulling fastq screen version for report
( echo -n "fastq screen: "; printf "\"%s\"\n" "$(fastq_screen --version 2>&1 | sed 's/[^0-9.]//g')" ) > {log.software_log} 2>&1

fastq_screen --outdir "brat_out/fastq_screen/{sample}/" \
             --subset 0 \
             --conf "tmp.brat/brat.fastq_screen.conf.mod" \
             --threads 10 \
             --nohits \
             "brat_out/cutadapt/{sample}/{sample}_{read}_trimmed.fastq.gz