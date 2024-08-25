#!/bin/bash

cd ~/GitHub/hello-nextflow/day2

bash_scripts/00_index.sh
bash_scripts/01_fastqc.sh gut
bash_scripts/02_quant.sh gut
bash_scripts/03_multiqc.sh
