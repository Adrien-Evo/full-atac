
import json
import os
import re
from os.path import join
import argparse
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument("--accession_file", help="Required. the FULL path to the accession file with all the SRA ids")
parser.add_argument("--fastq_dir", help="Required. the FULL path to the root folder where fastq can be found. It will dive into all subfolder matching fastqs with their accession name. This si only for paired end.")
args = parser.parse_args()

assert args.accession_file is not None, "please provide the path to the SRA accession file"
assert args.fastq_dir is not None, "please provide the path to the fastq folder"	


## collect all the fastq.gz full path in to a list
fastq_paths = []
FILES = defaultdict(lambda: defaultdict(list))

with open(args.accession_file, "r") as f:
    lines = f.readlines()
    for root, dirs, files in os.walk(args.fastq_dir):
        for file in files:
            if file.endswith("fastq.gz"):
                full_path = join(root, file)
                fastq_paths.append(full_path)
		
    for accession in lines:
        sample_name = accession.strip()
        ## now just assume the file name in the accession file  iscontained in the fastq file path
        fastq_full_path = [x for x in fastq_paths if sample_name in x]
        if fastq_full_path:
            R1 = [x for x in fastq_full_path if "_1" in x]
            R2 = [x for x in fastq_full_path if "_2" in x]
            R1.sort()
            R2.sort()
            FILES[sample_name]["R1"]= R1
            FILES[sample_name]["R2"]= R2
        else:
            print("sample {sample_name} missing its fastq files".format(sample_name = sample_name))

sample_num = len(FILES.keys())
output_file = os.path.join(os.path.dirname(args.accession_file),"samples_test.json")

print ("A total {} unique samples will be processed".format(sample_num))
print( "Output file in json format available here:\n{output_file}".format(output_file = output_file))
print ("------------------------------------------")
js = json.dumps(FILES, indent = 4, sort_keys=True)
open('samples_test.json', 'w').writelines(js)
