#!/bin/bash

################################################
##   SETUP OF PIPELINE FOLDER AND YML CONDA   ##
################################################

#create folder for Pipeline
mkdir ./Pipeline_Virus_Assembly_V0
cd ./Pipeline_Virus_Assembly_V0

#create .yml with all tools

echo ""
echo "Install environment.."
echo ""

echo "name: Pipeline_V0
channels:
- defaults
- bioconda
- conda-forge
dependencies:
- samtools=1.6
- bedtools=2.25.0
- fastqc=0.11.7
- seqtk=1.3
- mafft=7.471
- bbmap=37.99
- megahit=1.2.9
- multiqc=1.9
- qualimap=2.2.2d
- blast=2.7.1
- git=2.23.0
- cutadapt=1.18
" > ./Pipeline_Virus_Assembly_V0.yml

################################################
##   MAKE ENVIRONMENT                         ##
################################################

# make new environment

conda env create -f ./Pipeline_Virus_Assembly_V0.yml

rm Pipeline_Virus_Assembly_V0.yml
