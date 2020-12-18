#!/bin/bash

# This script needs a file with SRA accession
source activate full-pipe-main-env

cat $1 | xargs -I{} prefetch {}

find * -type d -print0 | xargs -I{} -0 fastq-dump -I --split-files {}


