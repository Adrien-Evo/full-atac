import csv
import os
import json
import yaml
import numpy as np
from snakemake.logging import logger
import re

#  Safe execution of scripts  #
shell.prefix("set -eo pipefail; echo BEGIN at $(date); ")

# ======================================================== #
# ================== Config file loading ================= #
# ======================================================== #

# Loading config file items
FILES = json.load(open(config['SAMPLES_JSON']))
WORKDIR = os.path.abspath(config["OUTPUT_DIR"])
PROJECT_NAME = config['PROJECT_NAME']
BAM_INPUT = config['bam']
CUTRUN_SAMPLES = config['cutrun']
NARROW_BROAD = yaml.safe_load(open(config['narrow_broad']))
# -- genome related config - #
GENOME_FASTA = config['genome_fasta']
GENOME_GTF = config['genome_gtf']
GENOME_SIZE = config['genome_size']
GENOME_TSS = config['genome_tss']
BLACKLIST_DCC = config['blacklist']


###########################################################################
#################### Defining samples, cases, controls ####################
###########################################################################

# ======================================================== #
# ==== Defining for the first time CASES samples    ====== #
# ======================================================== #

SAMPLE_MARK = []
SAMPLES_NAMES = sorted(FILES.keys())

# Create sample_Marks list for all samples
# e.g. Mousekidney01_H3K27, Mousekidney01_H3K27me3, Mouseliver04_H3K27, Mouseliver04_H3K27me3
for sample in SAMPLES_NAMES:
    for mark_or_tf in FILES[sample].keys():
        SAMPLE_MARK.append(sample + "_" + mark_or_tf)


# Which sample_type is used as control for calling peaks: e.g. Input, IgG...
CONTROL_NAME = config["control"]

#  defining CASES samples : can be CUTRUN or CHIPSEQ #

CASES = [sample for sample in SAMPLE_MARK if CONTROL_NAME not in sample]


# ======================================================== #
# ============== Defining CUTRUN samples ================= #
# ======================================================== #

CUTRUN_MARK = []
#for sample in CUTRUN_SAMPLES:
#    for mark_or_tf in FILES[sample].keys():
#        CUTRUN_MARK.append(sample + "_" + mark_or_tf)


# ======================================================== #
# =========== Defining normal CHipSeqsamples ============= #
# ======================================================== #
CHIPSEQ_SAMPLES = [sample for sample in SAMPLES_NAMES if sample not in CUTRUN_SAMPLES]

CHIPSEQ_MARK = [sample for sample in SAMPLE_MARK if sample not in CUTRUN_MARK]

CHIPSEQ_CASES = [sample for sample in CHIPSEQ_MARK if CONTROL_NAME not in sample]

# ======================================================== #
# =========  Defining Control/Input samples   ============ #
# ======================================================== #


CONTROL_SAMP = [sample for sample in SAMPLE_MARK if CONTROL_NAME in sample]

# Create a dictionary linking each sample with their control fastq or bam e.g. { Mousekidney01 : controlIG16_TCGCTAGA_L001_R1_002.fastq.gz 16_TCGCTAGA_L001_R1_001.fastq.gz}  #
# Here we want to check if input/control are used multiple times for different sample. We then create aliases for those input/control to ensure not having to align them twice. We use the fastq path as a dictionary key

controlFile = dict()
for samp in CONTROL_SAMP:
    sample = "".join(samp.split("_")[0:-1])
    #Here mark should contain the Input name
    mark = samp.split("_")[-1]
    # If Input is paired end 
    if('R1' in FILES[sample][mark] or 'R2' in FILES[sample][mark]):
        controlFile[sample] = " ".join(FILES[sample][mark]['R1'])
    else:
        controlFile[sample] = " ".join(FILES[sample][mark])


#Checking if control has been found:
if not bool(controlFile):
    logger.warning("Can't file any controls/input named " + CONTROL_NAME + ". Will work without controls.")


# Finding duplicate values from controlFile by flipping the dictionary
controlFileFlipped = {} 
  
for key, value in controlFile.items():
    if value not in controlFileFlipped: 
        controlFileFlipped[value] = [key] 
    else: 
        controlFileFlipped[value].append(key) 
# controlFileFlipped dict now is of length the number of unique controls, with the samples using those controls as values

# mergedInputDict allows to create a generic name for the Inputs ( Input1, Input2 etc etc)
mergedInputDict = controlFileFlipped
i = 1
for key, value in controlFileFlipped.items():
    inputname = "Input" + str(i)
    mergedInputDict[key] = inputname
    i = i + 1

# Now creating CONTROL_SAMPLE_DICT,  linking sample with their unique Input using the generic input name e.g. {Mousekidney: Input1}
CONTROL_SAMPLE_DICT = {}
for key, value in controlFile.items():
    CONTROL_SAMPLE_DICT[key] = mergedInputDict[value]


CONTROL_SAMPLE_MARK_DICT = {}
for key, value in CONTROL_SAMPLE_DICT.items():
    for samp in FILES[key].keys():
        if CONTROL_NAME not in samp: 
            CONTROL_SAMPLE_MARK_DICT[key + "_" + samp ]  = value


# Create a dict with Input name to the fastq file
CONTROL_MERGED_FILES = {}
for samp in CONTROL_SAMP:
    sample = "".join(samp.split("_")[0:-1])
    #Here mark should contain the Input name
    mark = samp.split("_")[-1]
    CONTROL_MERGED_FILES[CONTROL_SAMPLE_DICT[sample]] = FILES[sample][mark]



# Creating the dictionnary for all sample name linked to their input file
#CASES_SAMPLE_FILES becomes CHIPSEQ_SAMPLE_FILES
CHIPSEQ_SAMPLE_FILES = {}
for case in CHIPSEQ_CASES:
    sample = "_".join(case.split("_")[0:-1])
    mark = case.split("_")[-1]
    CHIPSEQ_SAMPLE_FILES[case] = FILES[sample][mark]

CUTRUN_SAMPLE_FILES = {}
#for case in CUTRUN_MARK:
#    sample = "_".join(case.split("_")[0:-1])
#    mark = case.split("_")[-1]
#    CUTRUN_SAMPLE_FILES[case] = FILES[sample][mark]



ALL_SAMPLE_FILES = {**CHIPSEQ_SAMPLE_FILES, **CONTROL_MERGED_FILES, **CUTRUN_SAMPLE_FILES }
ALL_CHIPSEQ_FILES = {**CHIPSEQ_SAMPLE_FILES, **CONTROL_MERGED_FILES}
CONTROLS = list(CONTROL_MERGED_FILES.keys())

# ~~~~~~~~~~~~~~ All samples ~~~~~~~~~~~~~~ #
ALL_SAMPLES = CASES + CONTROLS
ALL_CHIPSEQ = CHIPSEQ_CASES + CONTROLS

# ======================================================== #
# =========== Defining paired/single samples ============= #
# ======================================================== #

PAIRED_SAMPLES = []
SINGLE_SAMPLES = []

for samples, files in ALL_SAMPLE_FILES.items():
    if 'R1' in files and 'R2' in files:
        PAIRED_SAMPLES.append(samples)
    else:
        SINGLE_SAMPLES.append(samples)

# ======================================================== #
# ============= Creating helper dictionaries ============= #
# ======================================================== # 

# ~~~~~~~~~~~~~~ Samples dict ~~~~~~~~~~~~~ #

# Regroup Marks or TF by sample
# e.g. Mousekidney01: [H3K27, H3K27me3, Input1], Mouseliver04: [H3K27, H3K27me3, Input2]
SAMPLES = dict()
for sample in sorted(FILES.keys()):
    for mark in FILES[sample].keys():
        #if(mark not in CONTROL_NAME):
        SAMPLES.setdefault(sample, []).append(mark)

# Adding the proper merged input to the marks of each samples
for sample in CHIPSEQ_SAMPLES:
    if sample in CONTROL_SAMPLE_DICT.keys():
        SAMPLES.setdefault(sample,[]).append(CONTROL_SAMPLE_DICT[sample])

# ~~~~~~~~~~~ Samples_name dict ~~~~~~~~~~~ #

# Regroup Marks or TF by sample using the full sample name SAMPLE_MARK
# e.g. Mousekidney01: [Mousekidney01_H3K27, Mousekidney01_H3K27me3], Mouseliver04: [_Mouseliver04_H3K27, Mouseliver04_H3K27me3]
SAMPLES_COMPLETE_NAME = dict()
for sample in sorted(FILES.keys()):
    for mark in FILES[sample].keys():
        if(mark not in CONTROL_NAME):
            SAMPLES_COMPLETE_NAME.setdefault(sample, []).append(sample + "_" + mark)

# Adding the proper merged input to the marks of each samples
for sample in CHIPSEQ_SAMPLES:
    if sample in CONTROL_SAMPLE_DICT.keys():
        SAMPLES_COMPLETE_NAME.setdefault(sample,[]).append(CONTROL_SAMPLE_DICT[sample])

# ~~~~~~~~~~~~~~ Marks dicts ~~~~~~~~~~~~~~ #

# Regroup samples per marks or TF
# e.g. H3K27: [Mousekidney01, Mouseliver04], H3K27me3: [Mousekidney01, Mouseliver04]
MARKS = dict()
for sample in sorted(FILES.keys()):
    for mark in FILES[sample].keys():
        if(mark not in CONTROL_NAME):
            MARKS.setdefault(mark, []).append(sample)

# Here create a list with all marks without Input/Control mentionned, before adding the controls to the dict
MARKS_NO_CONTROL = list(MARKS.keys())

# Adding the key for the merged input
for key, value in CONTROL_SAMPLE_DICT.items():
    MARKS.setdefault(value, []).append(key)

# ~~~~~~~~~~~~~~ Marks_name dicts ~~~~~~~~~~~~~~ #

# Regroup samples per marks or TF
# e.g. H3K27: [Mousekidney01_H3K27, Mouseliver04_H3K27], H3K27me3: [Mousekidney01_H3K27me3, Mouseliver04_H3K27me3]
MARKS_COMPLETE_NAME = dict()
for sample in sorted(FILES.keys()):
    for mark in FILES[sample].keys():
        if(mark not in CONTROL_NAME):
            MARKS_COMPLETE_NAME.setdefault(mark, []).append(sample + "_"  + mark)

# Adding the key for the merged input
for key, value in CONTROL_SAMPLE_DICT.items():
    MARKS_COMPLETE_NAME.setdefault(value, []).append(key)


# Adding complete file name for SINGLE and PAIRED end samples

PAIRED_SAMPLES_DICT_FASTQ_R1 = {}
PAIRED_SAMPLES_DICT_FASTQ_R2 = {}
SINGLE_SAMPLES_DICT_FASTQ = {}

for sample in PAIRED_SAMPLES :
    PAIRED_SAMPLES_DICT_FASTQ_R1[sample] = ALL_SAMPLE_FILES[sample]['R1']
    PAIRED_SAMPLES_DICT_FASTQ_R2[sample] = ALL_SAMPLE_FILES[sample]['R2']

for sample in SINGLE_SAMPLES :
    SINGLE_SAMPLES_DICT_FASTQ[sample] = ALL_SAMPLE_FILES[sample]

FASTQ_NAME_WITH_SUFFIX = expand("{sample}", sample = SINGLE_SAMPLES)
FASTQ_NAME_WITH_SUFFIX.extend(expand("{sample}_R1", sample = PAIRED_SAMPLES))
FASTQ_NAME_WITH_SUFFIX.extend(expand("{sample}_R2", sample = PAIRED_SAMPLES))


# Checking dict
#print("PAIRED      ", PAIRED_SAMPLES)
#print("SAMPLES with full fastq path     ",ALL_SAMPLE_FILES)
#print("ALL_CONTROL     ", CONTROL_MERGED_FILES)
#print("MARKS     ", MARKS)
#print("MARKS_NO_CONTROL     ", MARKS_NO_CONTROL)
#print("MARKS_COMPLETE_NAME     ", MARKS_COMPLETE_NAME)
#print("CONTROL_SAMPLE_DICT     ",CONTROL_SAMPLE_DICT)
#print("CONTROL_SAMPLE_MARK_DICT     ",CONTROL_SAMPLE_MARK_DICT)
#print("R1 Paired SAMPLES with full fastq path     ",PAIRED_SAMPLES_DICT_FASTQ_R1)
#print("R2 Paired SAMPLES with full fastq path     ",PAIRED_SAMPLES_DICT_FASTQ_R2)
#print("Single end SAMPLES with full fastq path     ",SINGLE_SAMPLES_DICT_FASTQ)
###########################################################################
############################# Helper functions ############################
###########################################################################

# --- only canonical chr --- #
def get_canonical_chromSize(genomeSize):
    with open(genomeSize, 'r') as fi:
        size = [line for line in fi if not re.search('\.|_', line)]

    output_file_name = "data/canonical_genome_size.txt"
    with open(output_file_name, 'w') as fo:
        for item in size:
            fo.write("%s" % item)
    return(output_file_name)

# --- Getting all chr in a genome size file in a list. Usefull for chromHMM to predict how many files will be produced --- #
def get_chr(chromSize):
    with open(chromSize, 'r') as fs:
        chr = [line.rstrip().split('\t')[0] for line in fs]
    return(chr)

CANONICAL_CHR = get_chr(get_canonical_chromSize(GENOME_SIZE))
print(CANONICAL_CHR)
###########################################################################
########################### Listing OUTPUT FILES ##########################
###########################################################################

# Not all of those output files will be use in the snakemake rules but it's good to have them all at the same place #


# ~~~~~~ files with case and control ~~~~~~ #
ALL_PEAKS = []
ALL_BIGBED = []
ALL_FEATURECOUNTS = []
ALL_BROADPEAKCOUNTS = []
ALL_NARROWPEAKCOUNTS = []
ALL_RELAXEDPEAKCOUNTS = []
ALL_STRINGENTPEAKCOUNTS = []
ALL_ANNOTATED_PEAKS = []

# going through all cases samples (sample_mark) #
for sample_mark in CASES:
    ALL_PEAKS.append(os.path.join(WORKDIR, "peak_calling/macs2_broad/{}-macs2_peaks.broadPeak").format(sample_mark))
    #ALL_PEAKS.append(os.path.join(WORKDIR, "peak_calling/macs2_narrow/{}-macs2_peaks.narrowPeak").format(sample_mark))
    ALL_BIGBED.append(os.path.join(WORKDIR, "visualisation/bigbeds/{}-macs2_peaks.bb").format(sample_mark))
    ALL_FEATURECOUNTS.append(os.path.join(WORKDIR, "QC/{}.FRiP.summary").format(sample_mark))
    ALL_BROADPEAKCOUNTS.append(os.path.join(WORKDIR, "QC/{}-broadpeak-count_mqc.json").format(sample_mark))
    ALL_NARROWPEAKCOUNTS.append(os.path.join(WORKDIR, "QC/{}-narrowpeak-count_mqc.json").format(sample_mark))
    ALL_RELAXEDPEAKCOUNTS.append(os.path.join(WORKDIR, "QC/{}-seacr-relaxed-count_mqc.json").format(sample_mark))
    ALL_STRINGENTPEAKCOUNTS.append(os.path.join(WORKDIR, "QC/{}-seacr-stringent-count_mqc.json").format(sample_mark))
    
    ALL_ANNOTATED_PEAKS.append(os.path.join(WORKDIR, "annotation/{}-peaks_annotated.txt").format(sample_mark))
   

ALL_PEAKCOUNTS =  ALL_BROADPEAKCOUNTS + ALL_NARROWPEAKCOUNTS + ALL_RELAXEDPEAKCOUNTS + ALL_STRINGENTPEAKCOUNTS 

# ~~~~~~~~~~~~~~~ Bam files ~~~~~~~~~~~~~~~ #
CONTROL_BAM = expand(os.path.join(WORKDIR, "alignment/bams/{sample}.sorted.bam"), sample = CONTROL_MERGED_FILES)
CASE_BAM = expand(os.path.join(WORKDIR, "alignment/bams/{sample}.sorted.bam"), sample = CASES)
ALL_BAM     = CONTROL_BAM + CASE_BAM


# ~~ All samples files (cases + control) ~~ #
ALL_DOWNSAMPLE_BAM = expand(os.path.join(WORKDIR, "alignment/downsampling/{sample}-downsample.sorted.bam"), sample = ALL_SAMPLES)
ALL_FASTQ  = expand(os.path.join(WORKDIR, "alignment/fastq/{fname}.fastq.gz"), fname = FASTQ_NAME_WITH_SUFFIX)
ALL_FASTQC  = expand(os.path.join(WORKDIR, "QC/fastqc/{fname}_fastqc.zip"), fname = FASTQ_NAME_WITH_SUFFIX)
ALL_BOWTIE_LOG = expand(os.path.join(WORKDIR, "logs/{sample}.align"), sample = ALL_SAMPLES)
ALL_INDEX = expand(os.path.join(WORKDIR, "alignment/bams/{sample}.sorted.bam.bai"), sample = ALL_SAMPLES)
ALL_DOWNSAMPLE_INDEX = expand(os.path.join(WORKDIR, "downsampling/{sample}-downsample.sorted.bam.bai"), sample = ALL_SAMPLES)
ALL_SAMTOOLSSTATS = expand(os.path.join(WORKDIR, "QC/samtoolsstats/{sample}.sorted.bam.stats"), sample = ALL_SAMPLES)
ALL_SAMTOOLSSTATS.extend(expand(os.path.join(WORKDIR, "QC/samtoolsstats/{sample}.raw.bam.stats"), sample = ALL_SAMPLES))
ALL_PHANTOM = expand(os.path.join(WORKDIR, "QC/phantompeakqualtools/{sample}.spp.out"), sample = ALL_SAMPLES)
ALL_BIGWIG = expand(os.path.join(WORKDIR, "visualisation/bigwigs/{sample}.bw"), sample = ALL_SAMPLES)

ALL_SAMBLASTER_LOG = expand(os.path.join(WORKDIR, "logs/{sample}.samblaster.log"), sample = ALL_SAMPLES)

# ~~~~~~~~~~~ ENCODE QC specific ~~~~~~~~~~ #
ALL_TSS = expand(os.path.join(WORKDIR, "QC/encodeQC/{sample}.tss.txt"), sample = ALL_SAMPLES)
ALL_TSS.extend(expand(os.path.join(WORKDIR, "QC/encodeQC/{sample}.max.txt"), sample = ALL_SAMPLES))
ALL_ENCODE = expand(os.path.join(WORKDIR, "QC/encodeQC/{sample}.encodeQC.txt"), sample = ALL_SAMPLES)

# ~~ Just sample with control ~~ #

ALL_BIGWIG_INPUT = expand(os.path.join(WORKDIR, "visualisation/bigwigs_with_control/{casewithcontrol}.bw"), casewithcontrol = CONTROL_SAMPLE_MARK_DICT)

# ~~~~~~~~~~~ Deeptools specific ~~~~~~~~~~ #
# ---- Grouped by marks ---- #
ALL_COMPUTEMATRIX = expand(os.path.join(WORKDIR, "QC/{mark}.computeMatrix.gz"), mark = MARKS)
ALL_DPQC_PLOT = expand(os.path.join(WORKDIR, "QC/plots/tss/{mark}.plotHeatmap.png"), mark = MARKS)
ALL_DPQC_PLOT.extend(expand(os.path.join(WORKDIR, "QC/plots/profile/{mark}.plotProfile.png"), mark = MARKS))
ALL_DPQC = expand(os.path.join(WORKDIR, "QC/{mark}.plotProfileOutFileNameData.txt"), mark = MARKS)


#Temp commenting of allCorrelation since it takes ages with a lot of samples
# --- Grouped by samples --- #
ALL_DPQC_PLOT.extend(expand(os.path.join(WORKDIR, "QC/plots/fingerprint/{sample}.fingerprint.png"), sample = ALL_SAMPLES))
#ALL_DPQC_PLOT.extend([os.path.join(WORKDIR, "QC/plots/correlation/plotCorrelation.png")])

ALL_DPQC.extend(expand(os.path.join(WORKDIR, "QC/{sample}.plotFingerprintOutRawCounts.txt"), sample = ALL_SAMPLES))
ALL_DPQC.extend(expand(os.path.join(WORKDIR, "QC/{sample}.plotFingerprintOutQualityMetrics.txt"), sample = ALL_SAMPLES))
#ALL_DPQC.extend([os.path.join(WORKDIR, "QC/outFileCorMatrix.txt")])


###TEMPORARY SEACR peaks

ALL_SEACR_peaks = expand(os.path.join(WORKDIR, "peak_calling/seacr/{case}.stringent.bed"), case = CASES)

ALL_SEACR_peaks.extend(expand(os.path.join(WORKDIR, "peak_calling/seacr/{case}.relaxed.bed"), case = CASES))


# ~~~~~~~~~~~ ChromHMM specific ~~~~~~~~~~~ #
if config["chromHMM"]:

    
    CANONICAL_GENOME_SIZE = get_canonical_chromSize(GENOME_SIZE)
    CHR_FOR_CHROMHMM = get_chr(CANONICAL_GENOME_SIZE)
    # Read histone
    HISTONE_INCLUDED = config["histone_for_chromHMM"]
    HISTONE_FOR_CHROMHMM = [histone for histone in MARKS if histone in HISTONE_INCLUDED ]
    #sample_mark for chromHMM in case not all samples have all the marks
    SAMPLE_MARK_FOR_CHROMHMM  = [MARKS_COMPLETE_NAME[histone] for histone in MARKS_COMPLETE_NAME if histone in HISTONE_INCLUDED]
    #Flattening the list
    SAMPLE_MARK_FOR_CHROMHMM = [sample for sublist in SAMPLE_MARK_FOR_CHROMHMM for sample in sublist]
    #sample for chromHMM in case not all samples have all the marks
    SAMPLE_FOR_CHROMHMM  = [MARKS[histone] for histone in MARKS if histone in HISTONE_INCLUDED]
    #Flattening the list and removing duplicates
    SAMPLE_FOR_CHROMHMM = list(set([sample for sublist in SAMPLE_FOR_CHROMHMM for sample in sublist]))

    BAM_TO_BED = expand(os.path.join(WORKDIR, "bamtobed/{sample}.bed"), sample = SAMPLE_MARK_FOR_CHROMHMM + CONTROLS)
    CHROMHMM = expand(os.path.join(WORKDIR, "chromHMM/learn_{nb_state}_states/{sample}_{nb_state}_segments.bed"), sample = SAMPLE_FOR_CHROMHMM, nb_state = config["state"])
    CHROMHMM_TABLE = [os.path.join(WORKDIR, "chromHMM/cellmarkfiletable.txt")]

# ~~~~~~~~~~~~~~~~~~ Misc ~~~~~~~~~~~~~~~~~ #
ALL_MULTIQC = [os.path.join(WORKDIR, "multiQC/multiqc_report.html")]

ALL_CONFIG= [os.path.join(WORKDIR, "alignment/bams.json")]

HUB_FOLDER = os.path.join(WORKDIR, "visualisation/{}_UCSC_hub").format(PROJECT_NAME)

ALL_HUB = [os.path.join(HUB_FOLDER,"{}.hub.txt").format(PROJECT_NAME)]

# ======================================================== #
# ==================== MULTIQC inputs ==================== #
# ======================================================== #

# Depends on bam or fastq as input #
ALL_MULTIQC_INPUT = ALL_PHANTOM + ALL_DPQC + ALL_FEATURECOUNTS +ALL_PEAKCOUNTS + ALL_ENCODE + ALL_TSS + ALL_SAMTOOLSSTATS
if not BAM_INPUT:
    ALL_MULTIQC_INPUT.extend(ALL_ENCODE + ALL_FASTQC + ALL_BOWTIE_LOG + ALL_SAMBLASTER_LOG)

###########################################################################
########################### Targets for rule all ##########################
###########################################################################
TARGETS = []
TARGETS.extend(ALL_ANNOTATED_PEAKS)
TARGETS.extend(ALL_MULTIQC)
TARGETS.extend(ALL_HUB)
#print(ALL_MULTIQC_INPUT)
# Since output from bam input are not used as input, needs to be put in the rule all for execution #
if not BAM_INPUT:
    TARGETS.extend(ALL_CONFIG)

#Here some of the ouputs of the rules are not used by multiQC and need adding to the rule all
#TARGETS.extend(ALL_DPQC_PLOT)
TARGETS.extend(ALL_PEAKS)

# ~~~~~~~~~~~~~~~~ ChromHMM ~~~~~~~~~~~~~~~ #
if config["chromHMM"]:
    TARGETS.extend(CHROMHMM)

############################################################################
############################ Wildcards functions ###########################
############################################################################

# ~ Aggregation of bigwigs by Marks or TF ~ #
def get_big_wig_with_mark_or_tf(wildcards):
    samples = MARKS_COMPLETE_NAME[wildcards.mark]
    bigwigs = list()
    if wildcards.mark in CONTROLS:
            bigwigs.append(os.path.join(WORKDIR, "visualisation/bigwigs/" + wildcards.mark + ".bw"))
    else:
        for s in samples:
            bigwigs.append(os.path.join(WORKDIR, "visualisation/bigwigs/" + s + ".bw"))
    return bigwigs

# ~ Aggregation of bigwigs by sample for each of its marks or TF ~ #
def get_big_wig_per_sample(wildcards):
    samples = SAMPLES_COMPLETE_NAME[wildcards.samp]
    bigwigs = list()
    if wildcards.samp in CONTROLS:
            bigwigs.append(os.path.join(WORKDIR, "visualisation/bigwigs/" + wildcards.samp + ".bw"))
    else:
        for s in samples:
            bigwigs.append(os.path.join(WORKDIR, "visualisation/bigwigs/" + s + ".bw"))
    return bigwigs

# ~~~~~ Aggregation of bams per sample ~~~~ #
def get_bams_per_sample(wildcards):
    marks = SAMPLES_COMPLETE_NAME[wildcards.samp]
    bams = list()
    for s in marks:
        bams.append(os.path.join(WORKDIR, "alignment/bams/" + s + ".sorted.bam"))
    return bams

# ~~~ Aggregation of bam idx per sample ~~~ #
def get_bam_index_per_sample(wildcards):
    marks = SAMPLES_COMPLETE_NAME[wildcards.samp]
    bams = list()
    for s in marks:
        bams.append(os.path.join(WORKDIR, "alignment/bams/" + s + ".sorted.bam.bai"))
    return bams

# ~~~~~~~~ marks or tf per samples ~~~~~~~~ #
def get_all_marks_per_sample(wildcards):
    return SAMPLES_COMPLETE_NAME[wildcards.samp]

# ~~~~~~~~ get fastq for single or paired end ~~~~~~~~ #
def get_fastq(wildcards):
    output = []
    if wildcards.sample in SINGLE_SAMPLES:
        output = os.path.join(WORKDIR, "alignment/fastq/" + wildcards.sample + ".fastq.gz")
    elif wildcards.sample in PAIRED_SAMPLES:
        output.append(os.path.join(WORKDIR, "alignment/fastq/" + wildcards.sample + "_R1.fastq.gz"))
        output.append(os.path.join(WORKDIR, "alignment/fastq/" + wildcards.sample + "_R2.fastq.gz"))
    return output


# ~~~~~~~ bam files for sample_mark ~~~~~~~ #
def get_bams(wildcards):
    return ALL_SAMPLE_FILES[wildcards.sample]

# ~~~~~~~ peak sets per marks or tf ~~~~~~~ #
def get_peaks(wildcards):
    for key, value in MARKS_COMPLETE_NAME.items():
        if any(wildcards.case in sample_mark for sample_mark in MARKS_COMPLETE_NAME[key]):
            if key in NARROW_BROAD:
                if(NARROW_BROAD[key] == 'narrow'):
                    return os.path.join(WORKDIR, "peak_calling/macs2_narrow/" + wildcards.case + "-macs2_peaks.narrowPeak")
                elif(NARROW_BROAD[key] == 'broad'):
                    return os.path.join(WORKDIR, "peak_calling/macs2_broad/" + wildcards.case + "-macs2_peaks.broadPeak")
            else:
                logger.warning("Marks or TF not in the {} for {}. Will work with narrow peaks".format(config['narrow_broad'],wildcards.case))
                return os.path.join(WORKDIR, "peak_calling/macs2_narrow/" + wildcards.case + "-macs2_peaks.narrowPeak")




def get_control_downsampled_bais(wildcards):
    return os.path.join(WORKDIR, "alignment/downsampling/"+ CONTROL_SAMPLE_MARK_DICT[wildcards.case] + "-downsample.sorted.bam.bai"),

def get_control_downsampled_bams_with_input(wildcards):
    return os.path.join(WORKDIR, "alignment/downsampling/"+ CONTROL_SAMPLE_MARK_DICT[wildcards.casewithcontrol] + "-downsample.sorted.bam"),

def get_control_downsampled_bais_with_input(wildcards):
    return os.path.join(WORKDIR, "alignment/downsampling/"+ CONTROL_SAMPLE_MARK_DICT[wildcards.casewithcontrol] + "-downsample.sorted.bam.bai"),



###########################################################################
################################# rule all ################################
###########################################################################

#TODO see if this is usefull
localrules: all


rule all:
    input: TARGETS
#rule all:
#    input: ALL_SEACR_peaks

# Ordering rule for ambiguous cases for paired and single end management
ruleorder: merge_fastqs_paired > merge_fastqs_single > fastqc

###########################################################################
######################### Alignement using BOWTIE2 ########################
###########################################################################

# Those rules are only executed if inputs   are fastq #
if BAM_INPUT == False:

    #Now only for single-end ChIPseq
    rule merge_fastqs_single:
        input: lambda wildcards: SINGLE_SAMPLES_DICT_FASTQ[wildcards.sample]
        output: temp(os.path.join(WORKDIR, "alignment/fastq/{sample}.fastq.gz"))
        message: "Merging fastqs."
        shell: 
            """
            cat {input} > {output}
            """

    rule merge_fastqs_paired:
        input:
            R1 = lambda wildcards: PAIRED_SAMPLES_DICT_FASTQ_R1[wildcards.sample],
            R2 = lambda wildcards: PAIRED_SAMPLES_DICT_FASTQ_R2[wildcards.sample]
        output:
            R1 = temp(os.path.join(WORKDIR,"alignment/fastq/{sample}_R1.fastq.gz")),
            R2 = temp(os.path.join(WORKDIR,"alignment/fastq/{sample}_R2.fastq.gz"))
        message: "Merging fastqs."
        shell:
            """
            cat {input.R1} > {output.R1}
            cat {input.R2} > {output.R2}
            """

    # Simple alignment with bowtie 2 followed by sorting #
    rule align:
        input:  get_fastq
        output: temp(os.path.join(WORKDIR, "alignment/{sample}.raw.bam"))
        threads: 16
        log:    
            bowtie = os.path.join(WORKDIR, "logs/{sample}.align"),
            samblaster = os.path.join(WORKDIR, "logs/{sample}.samblaster.log")
        params:
            input= (
                lambda wildcards, input: ["-U", input]
                if wildcards.sample in SINGLE_SAMPLES
                else ["-1", input[0], "-2", input[1], "-I", "0", "-X", "2000", "--local", "--dovetail"]
                )
        shell:
            """
            source activate full-pipe-main-env
            bowtie2 --threads {threads} -x {config[idx_bt2]} {params.input} 2> {log.bowtie} \
            | samblaster --ignoreUnmated 2> {log.samblaster}\
            | samtools view -bu `{CANONICAL_CHR}`\
            | samtools sort -m 8G -@ 4 -T {output}.tmp -o {output}
            """


    # Get the duplicates marked sorted bam, remove unmapped reads by samtools view -F 1804 #
    rule filter_alignment:
        input:  os.path.join(WORKDIR, "alignment/{sample}.raw.bam")
        output: os.path.join(WORKDIR, "alignment/bams/{sample}.sorted.bam")
        log:    os.path.join(WORKDIR, "logs/{sample}.filter")
        shell:
            """
            source activate full-pipe-main-env
            bedtools intersect -v -abam {input} -b {BLACKLIST_DCC} | samtools view -bu -F 1804 -q {config[MQ]} \
            | samtools sort -m 8G -@ 4 -T {output}.tmp -o {output} 2> {log}
            """
    
    # This rule is not followed by other rules, so its output has to be added to the rule all conditionally on BAM_INPUT, if Bam are used as input or not #
    rule create_bam_json:
        input: expand(os.path.join(WORKDIR, "alignment/bams/{sample}.sorted.bam"), sample = ALL_SAMPLES),
        output: os.path.join(WORKDIR, "alignment/bams.json")
        params: SAMPLES
        run:
            #sample dictionary that will be dumped as a json
            dict_for_json = {}
            #Going through all samples usign the sample dictionary -> {sample1: [mark1, mark2],sample2: [mark1, mark2])
            for samp in SAMPLES:
                #Mini dictionary containing marks and associated sorted bam filepath
                mini_dict = {}
                for mark in SAMPLES[samp]:
                    
                    filepath  =  os.path.join(WORKDIR, "alignment/bams/" + samp + "_" + mark + ".sorted.bam")
                    mini_dict.update({mark : filepath})

                #Adding the mini dictionary to the sample dictionary
                dict_for_json[samp] = mini_dict

            with open(output[0],'w') as outFile:
                outFile.write(json.dumps(dict_for_json, indent = 4))


###########################################################################
################################ BAM INPUT ################################
###########################################################################

# Pipeline can use bams as input #
if BAM_INPUT:

    rule symlink_bam:
        input: get_bams
        output: os.path.join(WORKDIR, "alignment/bams/{sample}.sorted.bam")
        shell:
            """
            ln -s {input} {output}
            """

# ~~~~~~~~~~~~~ Indexing bams ~~~~~~~~~~~~~ #
rule index_bam:
    input:  os.path.join(WORKDIR, "alignment/bams/{sample}.sorted.bam")
    output: os.path.join(WORKDIR, "alignment/bams/{sample}.sorted.bam.bai")
    log:    os.path.join(WORKDIR, "logs/{sample}.bam.index")
    threads: 1
    shell:
        """
        source activate full-pipe-main-env
        samtools index {input} 2> {log}
        """

###########################################################################
############################### Downsampling ##############################
###########################################################################
#  Using user provided parameters, bam will be downsampled. Samtools stats is used for read counting and fed to sambamba  #

rule samtoolsstats_raw_bam:
    input:  os.path.join(WORKDIR, "alignment/{sample}.raw.bam")
    output: os.path.join(WORKDIR, "QC/samtoolsstats/{sample}.raw.bam.stats")
    log:    os.path.join(WORKDIR, "logs/{sample}.bam.stats")
    threads: 1
    params: jobname = "{sample}"
    shell:
        """
        source activate full-pipe-main-env
        samtools stats {input} > {output} 2> {log}
        """
rule samtoolsstats_sorted_bam:
    input:  os.path.join(WORKDIR, "alignment/bams/{sample}.sorted.bam")
    output: os.path.join(WORKDIR, "QC/samtoolsstats/{sample}.sorted.bam.stats")
    log:    os.path.join(WORKDIR, "logs/{sample}.bam.stats")
    threads: 1
    params: jobname = "{sample}"
    shell:
        """
        source activate full-pipe-main-env
        samtools stats {input} > {output} 2> {log}
        """


#downsampling
rule down_sample:
    input: 
        bam = os.path.join(WORKDIR, "alignment/bams/{sample}.sorted.bam"), 
        bai = os.path.join(WORKDIR, "alignment/bams/{sample}.sorted.bam.bai"), 
        stats = os.path.join(WORKDIR, "QC/stats/{sample}.sorted.bam.stats")
    output: 
        bam = os.path.join(WORKDIR, "alignment/downsampling/{sample}-downsample.sorted.bam"), 
        bai = os.path.join(WORKDIR, "alignment/downsampling/{sample}-downsample.sorted.bam.bai")
    log: os.path.join(WORKDIR, "logs/{sample}.downsample.log")
    threads: 4
    params: 
    log: os.path.join(WORKDIR, "logs/{sample}.phantompeakqual.log")
    message: "downsampling for {sample}"
    shell:
        """
        source activate full-pipe-main-env
        sambamba view -f bam -t {threads} --subsampling-seed=3 -s `grep "reads mapped:" {input.stats} | cut -f3 | awk '{{ratio = {config[target_reads]}/$0}};{{if(ratio < 1 )print ratio; else print 1}}'` {input.bam} | samtools sort -m 2G -@ 5 -T {output.bam}.tmp > {output.bam} 2> {log}
        samtools index {output.bam}
        """

###########################################################################
#################################### QC ###################################
###########################################################################

#Library complexity
if BAM_INPUT == False:
    rule encode_complexity:
        input: 
            bam = os.path.join(WORKDIR, "alignment/{sample}.raw.bam") 
        output: 
            os.path.join(WORKDIR, "QC/encodeQC/{sample}.encodeQC.txt")
        shell:
            """
            source activate full-pipe-main-env
            bedtools bamtobed -i {input} | awk 'BEGIN{{OFS="\t"}}{{print $1,$2,$3,$6}}' | grep -v 'chrM' | sort | uniq -c \
            | awk 'BEGIN{{mt=0;m0=0;m1=0;m2=0;OFS="\t"}} ($1==1){{m1=m1+1}} ($1==2){{m2=m2+1}} {{m0=m0+1}} {{mt=mt+$1}} \
            END{{m1_m2=-1.0; if(m2>0) m1_m2=m1/m2; print "Sample Name","NRF","PBC1","PBC2"; print "{wildcards.sample}",m0/mt,m1/m0,m1_m2}}' > {output}
            """

    rule fastqc:
        input:  os.path.join(WORKDIR, "alignment/fastq/{fname}.fastq.gz")
        output: os.path.join(WORKDIR, "QC/fastqc/{fname}_fastqc.zip"), os.path.join(WORKDIR, "QC/fastqc/{fname}_fastqc.html")
        log:    os.path.join(WORKDIR, "logs/{fname}.fastqc")
        params:
            output_dir = os.path.join(WORKDIR, "QC/fastqc")
        shell:
            """
            source activate full-pipe-main-env
            fastqc -o {params.output_dir} -f fastq --noextract {input} 2> {log}
            """

# Phantompeakqualtools computes a robust fragment length using the cross correlation (xcor) metrics.
rule phantom_peak_qual:
    input: 
        bam = os.path.join(WORKDIR, "alignment/bams/{sample}.sorted.bam"), 
        bai = os.path.join(WORKDIR, "alignment/bams/{sample}.sorted.bam.bai")
    output: os.path.join(WORKDIR, "QC/phantompeakqualtools/{sample}.spp.out")
    log: os.path.join(WORKDIR, "logs/{sample}.phantompeakqual.log")
    threads: 4
    params: os.path.join(WORKDIR, "QC/phantompeakqualtools/")
    message: "phantompeakqual for {input}"
    shell:
        """
        source activate full-pipe-spp
        run_spp.R -c={input.bam} -savp -rf -p={threads} -odir={params}  -out={output} -tmpdir={params} 2> {log}
        """
rule get_tss_info:
    input:
        bam = os.path.join(WORKDIR, "alignment/bams/{sample}.sorted.bam"),
	bai = os.path.join(WORKDIR, "alignment/bams/{sample}.sorted.bam.bai"),
    output:
        tss_count = os.path.join(WORKDIR, "QC/encodeQC/{sample}.tss.txt"),
        tss_stat = os.path.join(WORKDIR, "QC/encodeQC/{sample}.max.txt")
    params:
        tss = GENOME_TSS,
        genome_size = GENOME_SIZE,
        plot_folder = os.path.join(WORKDIR, "QC/encodeQC/"),
        sample_name = "{sample}",
        fragment_length = 200
    threads: 8
    shell:
        """
        source activate full-pipe-tss
        python2.7 scripts/tss_plots.py --bam {input.bam} --tss {params.tss} --sampleName {params.sample_name} --chromSizes {params.genome_size} \
        --plotFolder {params.plot_folder} --outputmultiQC_plot {output.tss_count} --outputmultiQC_stat {output.tss_stat} --fragLength {params.fragment_length} \
        --bins 400 --window 2000 --threads {threads}
        """

# Deeptools QC
rule computeMatrix_QC:
    input : get_big_wig_with_mark_or_tf 
    output : os.path.join(WORKDIR, "QC/{mark}.computeMatrix.gz")
    params : GENOME_TSS
    threads: 4
    shell:
        """
        source activate full-pipe-main-env
	    computeMatrix reference-point -S {input} -R {params} -a 2000 -b 2000 -out {output} --numberOfProcessors {threads}
        """

# Deeptools QC
rule plotHeatmap:
    input :  os.path.join(WORKDIR, "QC/{mark}.computeMatrix.gz")
    output : os.path.join(WORKDIR, "QC/plots/tss/{mark}.plotHeatmap.png")
    shell:
        """
        source activate full-pipe-main-env
        plotHeatmap -m {input} -out {output} --colorMap jet
        """

# Deeptools QC
rule plotProfile:
    input :  os.path.join(WORKDIR, "QC/{mark}.computeMatrix.gz")
    output : 
        plot = os.path.join(WORKDIR, "QC/plots/profile/{mark}.plotProfile.png"),
        outFileNameData = os.path.join(WORKDIR, "QC/{mark}.plotProfileOutFileNameData.txt")
    shell:
        """
        source activate full-pipe-main-env
        plotProfile -m {input} -out {output.plot} --outFileNameData {output.outFileNameData} --refPointLabel TSS
        """


# Deeptools correlation of all samples grouped together for multiQC
rule all_multiBigwigSummary:
    input : ALL_BIGWIG
    output : os.path.join(WORKDIR, "QC/multiBigwigSummary.npz")
    threads: 4
    shell:
        """
        source activate full-pipe-main-env
	    multiBigwigSummary bins --bwfiles {input} -out {output} --region chr1 --numberOfProcessors {threads}
        """

rule all_plotCorrelation:
    input :  os.path.join(WORKDIR, "QC/multiBigwigSummary.npz")
    output : 
        plot = os.path.join(WORKDIR, "QC/plots/correlation/plotCorrelation.png"),
        outFileNameData = os.path.join(WORKDIR, "QC/outFileCorMatrix.txt")
    shell:
        """
        source activate full-pipe-main-env
        plotCorrelation --corData {input} --plotFile {output.plot} --outFileCorMatrix {output.outFileNameData} --corMethod pearson --whatToPlot heatmap --plotNumbers --skipZeros
        """

# ChipSeq QCs plots from deeptools. Plotfingerprints are really usefull to see focal enrichment of your Chip-Seq enrichment
rule plotFingerPrint:
    input:
        bam = os.path.join(WORKDIR, "alignment/bams/{sample}.sorted.bam"), 
        bai = os.path.join(WORKDIR, "alignment/bams/{sample}.sorted.bam.bai")
    output:
        plot = os.path.join(WORKDIR, "QC/plots/fingerprint/{sample}.fingerprint.png"), 
        rawCounts = os.path.join(WORKDIR, "QC/{sample}.plotFingerprintOutRawCounts.txt"),
        qualityMetrics = os.path.join(WORKDIR, "QC/{sample}.plotFingerprintOutQualityMetrics.txt")
    shell:
        """
        source activate full-pipe-main-env
        plotFingerprint -b {input.bam} --plotFile {output.plot} --region chr1 --skipZeros --numberOfSamples 100000 --minMappingQuality {config[MQ]} --plotTitle {wildcards.sample} --outRawCounts {output.rawCounts} --outQualityMetrics {output.qualityMetrics}
        """

# Here featureCounts is used to compute the Fraction of reads in peaks. In case there is no peak called, a fake peak is added in the saf file to avoid featureCounts throwing an error and stopping the snakemake
rule get_FRiP_for_multiqc:
    input:
        peaks = get_peaks,
        bam = os.path.join(WORKDIR, "alignment/bams/{case}.sorted.bam"), 
    output:
        os.path.join(WORKDIR, "QC/{case}.FRiP.summary")
    params:
        saf = os.path.join(WORKDIR, "QC/{case}.saf"),
        outputName = os.path.join(WORKDIR, "QC/{case}.FRiP")
    shell:
        """
        source activate full-pipe-main-env
        awk 'BEGIN{{OFS="\t";print "GeneID", "Chr","Start","End","Strand";print "FUBAR", "1","1000","1001","."}}{{print $4,$1,$2,$3,$6}}' {input.peaks} > {params.saf}
        featureCounts -a {params.saf} -F SAF -o {params.outputName} {input.bam}
        """

rule get_broad_peak_counts_for_multiqc:
    input:
        peaks = os.path.join(WORKDIR, "peak_calling/macs2_broad/{case}-macs2_peaks.broadPeak"),
    output:
        os.path.join(WORKDIR, "QC/{case}-broadpeak-count_mqc.json")
    params:
        peakType = "MACS2 broadPeak"
    shell:
        """
        source activate full-pipe-main-env
        python3 scripts/count_peaks.py --peak_type {params.peakType} --peaks {input.peaks} --sample_name {wildcards.case} > {output}
        """

rule get_narrow_peak_counts_for_multiqc:
    input:
        peaks = os.path.join(WORKDIR, "peak_calling/macs2_narrow/{case}-macs2_peaks.narrowPeak"),
    output:
        os.path.join(WORKDIR, "QC/{case}-narrowpeak-count_mqc.json")
    params:
        peakType = "MACS2 narrowPeak"
    shell:
        """
        source activate full-pipe-main-env
        python3 scripts/count_peaks.py --peak_type {params.peakType} --peaks {input.peaks} --sample_name {wildcards.case} > {output}
        """

#SEACR relaxed
rule get_relaxed_peak_counts_for_multiqc:
    input:
        peaks = os.path.join(WORKDIR, "peak_calling/seacr/{case}.relaxed.bed")
    output:
        os.path.join(WORKDIR, "QC/{case}-seacr-relaxed-count_mqc.json")
    params:
        peakType = "SEACR relaxed"
    shell:
        """
        source activate full-pipe-main-env
        python3 scripts/count_peaks.py --peak_type {params.peakType} --peaks {input.peaks} --sample_name {wildcards.case} > {output}
        """

#SEACR stringent
rule get_stringent_peak_counts_for_multiqc:
    input:
        peaks = os.path.join(WORKDIR, "peak_calling/seacr/{case}.stringent.bed")
    output:
        os.path.join(WORKDIR, "QC/{case}-seacr-stringent-count_mqc.json")
    params:
        peakType = "SEACR stringent"
    shell:
        """
        source activate full-pipe-main-env
        python3 scripts/count_peaks.py --peak_type {params.peakType} --peaks {input.peaks} --sample_name {wildcards.case} > {output}
        """


###########################################################################
############################### PEAK CALLING ##############################
###########################################################################



def get_bam_peak_calling_input(wildcards):
    
    input_for_peak_calling = []
    #Checking if paired or not
    input_for_peak_calling.append(os.path.join(WORKDIR, "alignment/downsampling/" + wildcards.case + "-downsample.sorted.bam"))
        #Checking if sample has a control
    if wildcards.case in CONTROL_SAMPLE_MARK_DICT.keys() :
        input_for_peak_calling.append(os.path.join(WORKDIR, "alignment/downsampling/"+ CONTROL_SAMPLE_MARK_DICT[wildcards.case] + "-downsample.sorted.bam"))
    
    return input_for_peak_calling


# Peak calling using MACS 2
rule call_broad_peaks_macs2:
    input: get_bam_peak_calling_input
    output:
        broad = os.path.join(WORKDIR, "peak_calling/macs2_broad/{case}-macs2_peaks.broadPeak")
    log: os.path.join(WORKDIR, "logs/{case}-call-peaks_macs2.log")
    params:
        name = "{case}-macs2", 
        outdir = os.path.join(WORKDIR, "peak_calling/macs2_broad"),
        control = (
            lambda wildcards, input: ["-c", input[1]]
            if wildcards.case in CONTROL_SAMPLE_MARK_DICT
            else [""]
        ),
        format = (
            lambda wildcards: ["BAM"]
            if wildcards.case in SINGLE_SAMPLES
            else ["BAMPE"]
        )
    message: "Calling broadpeaks with macs2."
    shell:
        """
        source activate full-pipe-macs
        macs2 callpeak -t {input[0]} \
            {params.control} --keep-dup all -f {params.format} -g {config[macs2_g]} \
            --outdir {params.outdir} -n {params.name} -p {config[macs2_pvalue]} --broad --broad-cutoff {config[macs2_pvalue_broad_cutoff]} &> {log}
        """


# Peak calling using MACS 2
rule call_narrow_peaks_macs2:
    input: 
        get_bam_peak_calling_input
    output:
        narrow = os.path.join(WORKDIR, "peak_calling/macs2_narrow/{case}-macs2_peaks.narrowPeak")
    log: os.path.join(WORKDIR, "logs/{case}-call-narrowpeaks_macs2.log")
    params:
        name = "{case}-macs2", 
        jobname = "{case}", 
        outdir = os.path.join(WORKDIR, "peak_calling/macs2_narrow"),
        control = (
            lambda wildcards, input: ["-c", input[1]]
            if wildcards.case in CONTROL_SAMPLE_MARK_DICT
            else [""]
        ),
        format = (
            lambda wildcards: ["BAM"]
            if wildcards.case in SINGLE_SAMPLES
            else ["BAMPE"]
        )
    message: "Calling broadpeaks with macs2."
    shell:
        """
        source activate full-pipe-macs
        macs2 callpeak -t {input[0]} \
            {params.control} --nomodel --keep-dup all -f {params.format} -g {config[macs2_g]} \
            --outdir {params.outdir} -n {params.name} &> {log}
        """


# Peak calling using seacr for cutrun samples 
rule get_input_for_seacr:
    input: 
        bam = os.path.join(WORKDIR, "alignment/bams/{sample}.sorted.bam") 
    output: 
        bedgraph = os.path.join(WORKDIR, "peak_calling/seacr/{sample}.fragments.bedgraph")
    params:
        prefix = os.path.join(WORKDIR, "peak_calling/seacr/{sample}"),
        genome_size = GENOME_SIZE
    shell:
        """
        source activate full-pipe-main-env
        bedtools bamtobed -bedpe -i {input.bam} > {params.prefix}.bed
        awk '$1==$4 && $6-$2 < 1000 {{print $0}}' {params.prefix}.bed > {params.prefix}.clean.bed
        cut -f 1,2,6 {params.prefix}.clean.bed | sort -k1,1 -k2,2n -k3,3n > {params.prefix}.fragments.bed
        bedtools genomecov -bg -i {params.prefix}.fragments.bed -g {params.genome_size} > {output.bedgraph}
        """

def get_bam_peak_calling_input_seacr(wildcards):
    
    input_for_peak_calling = []
    #Checking if paired or not
    input_for_peak_calling.append(os.path.join(WORKDIR, "peak_calling/seacr/" + wildcards.case + ".fragments.bedgraph"))
    #Checking if sample has a control
    if wildcards.case in CONTROL_SAMPLE_MARK_DICT.keys() :
        input_for_peak_calling.append(os.path.join(WORKDIR, "peak_calling/seacr/"+ CONTROL_SAMPLE_MARK_DICT[wildcards.case] + ".fragments.bedgraph"))
    
    return input_for_peak_calling

#Here if the case sample has a control, then its going to be used. Else it's the numeric threshold
rule call_stringent_peaks_seacr:
    input: 
        get_bam_peak_calling_input_seacr 
    output: 
        stringent = os.path.join(WORKDIR, "peak_calling/seacr/{case}.stringent.bed")
    params:
        prefix = os.path.join(WORKDIR, "peak_calling/seacr/{case}"),
        control = (
            lambda wildcards, input: [input[1]]
            if wildcards.case in CONTROL_SAMPLE_MARK_DICT
            else ["0.01"]
        )
    shell:
        """
        source activate full-pipe-main-env
        SEACR_1.3.sh {input[0]} {params.control} norm stringent {params.prefix}
        """

rule call_relaxed_peaks_seacr:
    input: 
        get_bam_peak_calling_input_seacr 
    output: 
        stringent = os.path.join(WORKDIR, "peak_calling/seacr/{case}.relaxed.bed")
    params:
        prefix = os.path.join(WORKDIR, "peak_calling/seacr/{case}"),
        control = (
            lambda wildcards, input: [input[1]]
            if wildcards.case in CONTROL_SAMPLE_MARK_DICT
            else ["0.01"]
        )
    shell:
        """
        source activate full-pipe-main-env
        SEACR_1.3.sh {input[0]} {params.control} norm relaxed {params.prefix}
        """

###########################################################################
####### VISUALIZATION bigWig and bigBed generation and HUB creation #######
###########################################################################
# rule get_bigwigs_using_inputs:
#     input : 
#         case =  os.path.join(WORKDIR, "alignment/downsampling/{casewithcontrol}-downsample.sorted.bam"),
#         bai_case = os.path.join(WORKDIR, "alignment/downsampling/{casewithcontrol}-downsample.sorted.bam.bai"),
#         control = get_control_downsampled_bams_with_input, 
#         bai_control = get_control_downsampled_bais_with_input,
#         spp = os.path.join(WORKDIR, "QC/phantompeakqualtools/{casewithcontrol}.spp.out")
#     output:  os.path.join(WORKDIR, "visualisation/bigwigs_with_control/{casewithcontrol}.bw")
#     log: os.path.join(WORKDIR, "logs/{casewithcontrol}.makebw")
#     threads: 4
#     params: jobname = "{casewithcontrol}"
#     message: "Making bigwig of {casewithcontrol} log2 fold change versus its control"
#     shell:
#         """
#         source activate full-pipe-main-env
#         bamCompare --bamfile1 {input.case} --bamfile2 {input.control} \
#         --normalizeUsing RPKM  --operation log2 --operation first --scaleFactorsMethod None --binSize 30 --smoothLength 150 --numberOfProcessors {threads} \
#         --extendReads `cut -f3 {input.spp} | awk 'BEGIN{{FS=","}}{{print $1}}'` -o {output} 2> {log}
#         """

rule get_bigwigs:
    input : 
        bam = os.path.join(WORKDIR, "alignment/downsampling/{sample}-downsample.sorted.bam"),
        bai = os.path.join(WORKDIR, "alignment/downsampling/{sample}-downsample.sorted.bam.bai"),
        spp = os.path.join(WORKDIR, "QC/phantompeakqualtools/{sample}.spp.out")
    output: os.path.join(WORKDIR, "visualisation/bigwigs/{sample}.bw")
    log: os.path.join(WORKDIR, "logs/{sample}.makebw")
    threads: 4
    params: jobname = "{sample}"
    message: "Making bigwig of {sample}"
    shell:
        """
        source activate full-pipe-main-env
        bamCoverage -b {input.bam} --normalizeUsing RPKM --binSize 30 --smoothLength 150 --numberOfProcessors {threads} \
        --extendReads `cut -f3 {input.spp} | awk 'BEGIN{{FS=","}}{{print $1}}'` -o {output} 2> {log}
        """

# Cleaning peaks by taking only columns 1, 2, 3 because narrowpeaks and broadpeaks are different.
rule get_bigbeds:
    input: get_peaks
    output: os.path.join(WORKDIR, "visualisation/bigbeds/{case}-macs2_peaks.bb")
    params : 
        bed1 = temp(os.path.join(WORKDIR, "visualisation/bigbeds/{case}-macs2_peaks_temp.bed"))
    shell:
        """
        source activate full-pipe-main-env
        cut -f1,2,3 {input} | LC_COLLATE=C sort -k1,1 -k2,2n > {params.bed1}
        scripts/bedToBigBed {params.bed1} {GENOME_SIZE} {output} -type=bed3
        """

# Creating a Hub for UCSC
rule get_UCSC_hub:
    input:  
        bed = ALL_BIGBED, 
        bigwig = ALL_BIGWIG
    output:
        ALL_HUB
    log: os.path.join(WORKDIR, "logs/log.trackhub")
    params:
        output_dir = HUB_FOLDER, 
        sample_name = list(SAMPLES.keys()), 
        categories = MARKS_NO_CONTROL
    shell:
        """
        source activate full-pipe-main-env
        python3 scripts/makeUCSCtrackHub.py --hub_name {PROJECT_NAME} --sample_name {params.sample_name} --categories {params.categories} \
        --output_dir {params.output_dir} --peaks {input.bed} --bw {input.bigwig} 2> {log}
        """

########################################################################### 
################################# MULTIQC #################################
###########################################################################

# MultiQC: Moving the config file in the multiqc dir + adding some custom info in the header
rule multiQC_config:
    input : "config/multiqc_config.yaml"
    output: os.path.join(WORKDIR, "multiQC/multiqc_config.yaml")
    message: "Moving multiqc config"
    shell:
        """
        sed "s/DATE/$(date)/g" {input} | sed "s/PROJECTNAME/{PROJECT_NAME}/g" > {output}
        """

# MultiQC: Takes fastqc, phantompeakqualtools, deeptools, bowtie, fastQC and feature counts for FRiP and custom MACS2 peak counts
rule multiQC:
    input : 
        multiqc_files = ALL_MULTIQC_INPUT,
        multiqc_config = os.path.join(WORKDIR, "multiQC/multiqc_config.yaml")
    output: ALL_MULTIQC
    params: os.path.join(WORKDIR, "multiQC/")
    log: os.path.join(WORKDIR, "logs/multiqc.log")
    message: "multiqc for all logs"
    shell:
        """
        source activate full-pipe-main-env
        multiqc {input.multiqc_files} -o {params} --config {input.multiqc_config} -v -f 2> {log}
        """

###########################################################################
########################### Fonctional analysis ###########################
###########################################################################

# Annotation with homer
rule homer_annotate:
    input: 
        peaks = get_peaks
    output:
        annotated_peaks = os.path.join(WORKDIR, "annotation/{case}-peaks_annotated.txt")
    params:
        genome = GENOME_FASTA, 
        gtf = GENOME_GTF
    log: os.path.join(WORKDIR,"logs/{case}-homer-annotated.log")
    message: "Annotating peak files with HOMER"
    shell:
        """
        source activate full-pipe-main-env
        annotatePeaks.pl {input.peaks} {params.genome} -gtf {params.gtf} > {output.annotated_peaks} 2> {log}
        """

###########################################################################
################################# ChromHMM ################################
###########################################################################

# This allow for all necessary steps for ChromHMM execution with the number of states declared in the config file 
if config["chromHMM"]:

    rule bam2bed:
        input :
            os.path.join(WORKDIR, "alignment/downsampling/{sample}-downsample.sorted.bam")
        output:
            os.path.join(WORKDIR, "bamtobed/{sample}.bed")
        params: jobname = "{sample}"
        log: os.path.join(WORKDIR, "logs/{sample}-bam2bed.log")
        message: "converting bam to bed for {input}"
        shell:
            """
            source activate full-pipe-main-env
            bedtools bamtobed -i {input} > {output}
            """
    #This rule does not need an input, but it's nice to have it wrapped up like that anyway.
    rule make_table:
        input : 
            expand(os.path.join(WORKDIR, "bamtobed/{sample}.bed"), sample = SAMPLE_MARK_FOR_CHROMHMM + CONTROLS)
        output : 
            os.path.join(WORKDIR, "chromHMM/cellmarkfiletable.txt")
        log: os.path.join(WORKDIR, "logs/make_table_chromHMM.log")
        message: "Making the cellmark table for chromHMM"
        run:
            import os
            from os.path import join
            with open (output[0], "w") as f:
                for histone in HISTONE_FOR_CHROMHMM:
                    for sample in MARKS[histone]:
                        control = CONTROL_SAMPLE_DICT[sample]
                        case_bed = sample + "_"+ histone + ".bed"
                        if os.path.exists(join(os.path.join(WORKDIR, "bamtobed"), case_bed)):
                            f.write(sample + "\t" +  histone + "\t" + case_bed + "\t" + control + ".bed" + "\n")
    
    rule chromHMM_binarize:
        input :
            cellmarkfiletable = os.path.join(WORKDIR, "chromHMM/cellmarkfiletable.txt"), 
            beds = expand(os.path.join(WORKDIR, "bamtobed/{sample}.bed"), sample = SAMPLE_MARK_FOR_CHROMHMM + CONTROLS)
        output:
            expand(os.path.join(WORKDIR, "chromHMM/binarizedData/{sample}_{chr}_binary.txt"), sample = SAMPLE_FOR_CHROMHMM, chr = CHR_FOR_CHROMHMM)
        log:
            os.path.join(WORKDIR, "logs/chromhmm_bin.log")
        params:
            folder = os.path.join(WORKDIR, "chromHMM/binarizedData/"), 
            bamtobed_folder = os.path.join(WORKDIR, "bamtobed/"), 
            memory = "32G"
        shell:
            """
            source activate full-pipe-main-env
            ChromHMM.sh -Xmx{params.memory} BinarizeBed -b {config[binsize]} {CANONICAL_GENOME_SIZE} {params.bamtobed_folder} {input.cellmarkfiletable} {params.folder} 2> {log}
            """

    rule chromHMM_learn:
        input:
            expand(os.path.join(WORKDIR, "chromHMM/binarizedData/{sample}_{chr}_binary.txt"), sample = SAMPLE_FOR_CHROMHMM, chr = CHR_FOR_CHROMHMM)
        output: 
            expand(os.path.join(WORKDIR, "chromHMM/learn_{nb_state}_states/{sample}_{nb_state}_segments.bed"), sample = SAMPLE_FOR_CHROMHMM, nb_state = config["state"])
        log: os.path.join(WORKDIR, "logs/chromhmm_learn.log")
        params:
            input_folder = os.path.join(WORKDIR, "chromHMM/binarizedData/"), 
            output_folder = expand(os.path.join(WORKDIR, "chromHMM/learn_{nb_state}_states/"), nb_state = config["state"]), 
            memory = "32G"
        shell:
            """
            source activate full-pipe-main-env
            unset DISPLAY && ChromHMM.sh -Xmx{params.memory} LearnModel -p 0 -b {config[binsize]} \
            {params.input_folder} {params.output_folder} {config[state]} {CANONICAL_GENOME_SIZE} 2> {log}
            """
