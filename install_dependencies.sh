#!/bin/bash

# Stop on error
set -e

#Install UCSC utilities

wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/bigWigToWig -O scripts/bigWigToWig
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/bigBedToBed -O scripts/bigBedToBed
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/bedToBigBed -O scripts/bedToBigBed
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/wigToBigWig -O scripts/wigToBigWig

chmod +x scripts/bigWigToWig scripts/bigBedToBed scripts/bedToBigBed scripts/wigToBigWig
