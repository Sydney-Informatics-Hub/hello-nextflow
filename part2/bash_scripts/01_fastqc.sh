#!/usr/bin/bash

SAMPLE_ID=$1
READS_1="data/ggal/${SAMPLE_ID}_1.fq"
READS_2="data/ggal/${SAMPLE_ID}_2.fq"

mkdir -p "results/fastqc_${SAMPLE_ID}_logs"
fastqc \
	--outdir "results/fastqc_${SAMPLE_ID}_logs" \
	--format fastq ${READS_1} ${READS_2}
