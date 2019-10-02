#!/bin/bash

# Stop on error
set -e

#Install UCSC utilities

wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/bigWigToWig -O scripts/bigWigToWig
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/bigBedToBed -O scripts/bigBedToBed
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/bedToBigBed -O scripts/bedToBigBed
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/wigToBigWig -O scripts/wigToBigWig

chmod +x scripts/bigWigToWig scripts/bigBedToBed scripts/bedToBigBed scripts/wigToBigWig


# Install conda env

echo "Creating conda env for chromHMM"
conda env create -f envs/full-pipe-chromhmm.yml

echo "Creating conda env for phantompeakqualtools"
conda env create -f envs/full-pipe-spp.yml

echo "Creating conda env for most tools in the pipeline"
conda env create -f envs/full-pipe-main-env.yml

echo "Creating conda env for MACS1.4 and MACS2"
conda env create -f envs/full-pipe-macs.yml



#get hs TSS
#wget -q -O - "http://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz"  | gunzip -c | grep -v "#" | awk '($3=="gene")' | grep protein_coding | awk '{OFS="\t"};{if($7 == "+"){start = $4} else if($7 == "-"){start = $5}};{print $1, start-1, start}'
