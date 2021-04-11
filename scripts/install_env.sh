#!/bin/bash

################################################
##   SETUP OF PIPELINE FOLDER AND YML CONDA   ##
################################################

#create .yml with all tools

echo ""
echo "Install environment.."
echo ""

echo "name: virus_assembly
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
" > ./virus_assembly.yml

################################################
##   MAKE ENVIRONMENT                         ##
################################################

# make new environment

conda env create -f ./virus_assembly.yml

rm ./virus_assembly.yml

# Download github 

git clone \
  https://raw.githubusercontent.com/jsede/virus_assembly.git \
  $CONDA_PREFIX/virus_assembly


  # Create links to pipeline
  
find \
  $CONDA_PREFIX/virus_assembly/ \
  -name "*.sh" \
  -exec chmod +x {} \;
  
ln -s \
  $CONDA_PREFIX/virus_assembly/virus_assembly.sh \
  $CONDA_PREFIX/bin/virus_assembly