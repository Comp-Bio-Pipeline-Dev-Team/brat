#!/bin/bash

testRefFileList=("/scratch/alpine/mapgar@xsede.org/outside_rna_seq_refs/user_reference_files/hg38.refGene.gtf"
                 "/scratch/alpine/mapgar@xsede.org/outside_rna_seq_refs/user_reference_files/hg38.p14.fa"
                )
symlinkTo="test_reference_files"

mkdir -p ${symlinkTo}

for file in ${testRefFileList[*]};
do
    fileName=$( basename ${file} )

    if [ -L ${symlinkTo}/${fileName} ]; 
    then
        echo "Symlink already exists for ${fileName}!"
    else
        echo "Symlink for ${fileName} does not exist. Generating symlink..."
        ln -s ${file} ${symlinkTo}/${fileName}
    fi
done