#!/bin/bash


#Install UCSC utilities
echo -e "\e[33m#################################################\e[0m"
echo -e "\e[33mInstalling dependencies for pipeline full-chipseq\e[0m"
echo -e "\e[33m#################################################\n\e[0m"


wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/bigBedToBed -O scripts/bigBedToBed
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/bedToBigBed -O scripts/bedToBigBed

chmod +x scripts/bigBedToBed scripts/bedToBigBed


# Install conda env
echo -e "\e[33m#################################################\e[0m"
echo -e "\e[33mInstalling conda environments\e[0m"
echo -e "\e[33m#################################################\n\e[0m"

echo -e "Creating conda env for \e[32mbowtie2 samtools bedtools deepTools sambamba fastqc multiqc featureCounts chromHMM homer \e[0m in the pipeline"
mamba env create -f envs/full-pipe-main-env.yml --name full-pipe-main-env

echo -e "Creating conda env for \e[32mphantompeakqualtools\e[0m"
mamba env create -f envs/full-pipe-spp.yml --name full-pipe-spp

echo -e "Creating conda env for \e[32mMACS1.4\e[0m and \e[32mMACS2\e[0m"
mamba env create -f envs/full-pipe-macs.yml --name full-pipe-macs

echo -e "Creating conda env for \e[32mtss plots\e[0m"
mamba env create -f envs/full-pipe-tss.yml --name full-pipe-tss

echo -e "\e[33m#################################################\e[0m"
echo -e "\e[33mInstalling data\e[0m"
echo -e "\e[33m#################################################\n\e[0m"

# TSS
# For hs
echo -e "\e[95mFor homo sapiens GRCh37\e[0m"

echo -e "\e[95m Annotation file and TSS beds\e[0m"
wget -q -O - "http://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz" | gunzip -c > ./data/Homo_sapiens.GRCh37.87.gtf
grep -v "#" ./data/Homo_sapiens.GRCh37.87.gtf | awk '($3=="gene")' \
| grep protein_coding | awk '{OFS="\t"};{if($7 == "+"){start = $4} else if($7 == "-"){start = $5}};{print $1, start-1, start}' \
| grep -v "\." | grep -v "_" | grep -v "MT" | grep -v "X" | grep -v "Y" > ./data/GRCh37_TSS.bed



echo -e "\e[95m Blacklist regions from ENCODE DCC\e[0m"
curl -LJO https://github.com/Boyle-Lab/Blacklist/archive/v2.0.tar.gz
tar -xvf Blacklist-2.0.tar.gz -C ./data
gzip -d data/Blacklist-2.0/lists/*.bed.gz

