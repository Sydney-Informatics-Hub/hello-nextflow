#!/usr/bin/bash

mkdir "results"
salmon index --transcripts data/ggal/transcriptome.fa --index results/salmon_index
