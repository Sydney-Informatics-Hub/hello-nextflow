#!/bin/bash

cd ~/GitHub/hello-nextflow/day2

# Samples to iterate through
samples=("gut" "lung" "liver")

bash_scripts/00_index.sh

for sample in "${samples[@]}"; do
	bash_scripts/01_fastqc.sh $sample
	bash_scripts/02_quant.sh $sample
done

bash_scripts/03_multiqc.sh
