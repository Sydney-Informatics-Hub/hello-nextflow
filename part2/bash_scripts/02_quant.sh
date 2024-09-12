#!/usr/bin/bash

SAMPLE_ID="gut"
READS_1="data/ggal/${SAMPLE_ID}_1.fq"
READS_2="data/ggal/${SAMPLE_ID}_2.fq"

salmon quant \
	--libType=U \
	-i results/salmon_index \
	-1 ${READS_1} \
	-2 ${READS_2} \
	-o results/${SAMPLE_ID}
