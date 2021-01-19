#!/usr/bin/bash

source activate full-pipe-main-env

#Get the sra fastq
cat $1 | xargs -I{} sh -c "prefetch --max-size 80G -c {}"

#Get the normal fastq format
cat $1 | xargs -I{} sh -c "fasterq-dump {}/{}.sra"

#gzip it
cat $1 | xargs -I{} sh -c "gzip {}.sra_1.fastq"
cat $1 | xargs -I{} sh -c "gzip {}.sra_2.fastq"

#Get the sra associated metadata
cat $1 | xargs -I{} sh -c "wget -O ./{}_info.csv 'http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term= {}'"

accession=$(realpath $1)
fastq_dir=$(realpath .)
python3 ../scripts/create_sample_json_from_SRA.py --accession_file $accession --fastq_dir $fastq_dir



