#!/usr/bin/bash

mkdir -p "results"
salmon index -t data/ggal/transcriptome.fa -i results/salmon_index
