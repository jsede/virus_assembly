# virus_assembly

A bash pipeline for de novo assembly of viral genomes generated via Illumina NGS. Currently handles the following viruses: HIV-1, RRV, HMPV.

Current version: V1

**Table of contents**
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)

## Requirements
A conda package manager like [Miniconda3](https://docs.conda.io/en/latest/miniconda.html).

## Installation

1.  Download the initial environment installation file 
   ```
   wget https://raw.githubusercontent.com/jsede/virus_assembly/main/scripts/install_env.sh
   ```
2. Run the script in the terminal 
  ```
   bash ./install_env.sh
   ```
3. Check if installation worked
  ```
   conda activate virus_assembly
  ```
 ## Usage
   ```
 Pipeline: NGS pipeline for viral assembly.
usage: virus_assembly [-h -v -p -q] (-i dir -m value -t value )
(-s string) 
with:
    -h  Show help text
    -v  Version of the pipeline
    -n  Name of RUN.
    -i  Input directory
    -s  Viral species
    -c  Perform clipping of primers
    -q  Perform quality check using fastQC
    -m  Memory
    -t  Number of threads
   ```
