#!/usr/bin/bash

mkdir "results"
salmon index -t data/ggal/transcriptome.fa -i results/salmon_index
