#!/usr/bin/env python3


import json
import os
import csv
import re
from os.path import join
import argparse
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument("--fastq_dir", help="Required. the FULL path to the fastq folder")
parser.add_argument("--meta", help="Required. the FULL path to the tab delimited meta file")
args = parser.parse_args()

assert args.fastq_dir is not None, "please provide the path to the fastq folder"
assert args.meta is not None, "please provide the path to the meta file"


## collect all the fastq.gz full path in to a list
fastq_paths = []
paired_dic = {}
for root, dirs, files in os.walk(args.fastq_dir):
    for file in files:
        if file.endswith("fastq.gz"):
            full_path = join(root, file)
            fastq_paths.append(full_path)


FILES = defaultdict(lambda: defaultdict(list))

with open(args.meta, "r") as f:
    reader = csv.reader(f, delimiter = "\t")
    # skip the header
    header = next(reader)
    for row in reader:
        sample_name = row[0].strip()
        fastq_name = row[1].strip()
        sample_type = row[2].strip()
        paired = row[3].strip()
        
        if paired == "no" :
    	## now just assume the file name in the metafile contained in the fastq file path
            paired_dic[sample_name] = 0
            fastq_full_path = [x for x in fastq_paths if fastq_name in x]
            if fastq_full_path:
                FILES[sample_name][sample_type].extend(fastq_full_path)
            else:
                print("sample {sample_name} missing {sample_type} {fastq_name} fastq files".format(sample_name = sample_name, sample_type = sample_type, fastq_name = fastq_name))
        elif paired == "yes":
            paired_dic[sample_name] = 1
            fastq_full_path = [x for x in fastq_paths if fastq_name in x]
            fastq_R1_full_path = [x for x in fastq_full_path if "R1" in x]
            fastq_R2_full_path = [x for x in fastq_full_path if "R2" in x]
            if fastq_full_path:
                FILES[sample_name][sample_type] = {"R1" : fastq_R1_full_path ,"R2" : fastq_R2_full_path}
            else:
                print("sample {sample_name} missing {sample_type} {fastq_name} fastq files".format(sample_name = sample_name, sample_type = sample_type, fastq_name = fastq_name))
            
print()
sample_num = len(FILES.keys())

print ("In total, {} unique samples will be processed".format(sample_num))
print ("------------------------------------------")

for sample in FILES.keys():
	print ("Sample {sample} has {n} mark or TF".format(sample = sample, n = len(FILES[sample])))
print ("------------------------------------------")
for sample_name in sorted(FILES.keys()):
	for sample_type in FILES[sample_name]:
            if paired_dic[sample_name]:
                R1 = "".join(FILES[sample_name][sample_type]["R1"])
                R2 = "".join(FILES[sample_name][sample_type]["R2"])
                print("Paired sample {sample_name}'s {sample_type} fastq path is {R1} for forward strand and {R2} for reverse".format(sample_name = sample_name, sample_type = sample_type, R1 = R1, R2 = R2))
            else:
                fastq_file = " ".join(FILES[sample_name][sample_type])
                print("Sample {sample_name}'s {sample_type} fastq path is {fastq_file}".format(sample_name = sample_name, sample_type = sample_type, fastq_file = fastq_file))
print ("------------------------------------------")

print("Output : samples_from_{meta}.json".format(meta = os.path.splitext(args.meta)[0]))
print ("------------------------------------------")

print("Please check the samples_from_{meta}.json file for any errors before launching the pipeline".format(meta = os.path.splitext(args.meta)[0]))
print()

js = json.dumps(FILES, indent = 4, sort_keys=True)
open('samples_from_{meta}.json'.format(meta = os.path.splitext(args.meta)[0]),'w').writelines(js)




