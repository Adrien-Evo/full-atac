import csv
import os
import json
import yaml
import numpy as np
from snakemake.logging import logger

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
NARROW_BROAD = yaml.load(open(config['narrow_broad']))
# -- genome related config - #
GENOME_FASTA = config['genome_fasta']
GENOME_GTF = config['genome_gtf']
GENOME_SIZE = config['genome_size']
GENOME_TSS = config['genome_tss']


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
    for sample_type in FILES[sample].keys():
        SAMPLE_MARK.append(sample + "_" + sample_type)


# Which sample_type is used as control for calling peaks: e.g. Input, IgG...
CONTROL_NAME = config["control"]

#  defining CASES samples  #
CASES = [sample for sample in SAMPLE_MARK if CONTROL_NAME not in sample]


# ======================================================== #
# =========  Defining Control/Input samples   ============ #
# ======================================================== #


#  Create a dictionary linking each sample with their control fastq e.g. { Mousekidney01 : controlIG16_TCGCTAGA_L001_R1_002.fastq.gz 16_TCGCTAGA_L001_R1_001.fastq.gz}  #
# Joining the list to allow for usage as a dictionary key
controlFastq = dict()
for samp in SAMPLE_MARK:
    if CONTROL_NAME in samp:
        sample = "".join(samp.split("_")[0:-1])
        mark = samp.split("_")[-1]
        controlFastq[sample] = " ".join(FILES[sample][mark])
 
# Finding duplicate values from controlFastq by flipping it 
controlFastqFlipped = {} 
  
for key, value in controlFastq.items():
    if value not in controlFastqFlipped: 
        controlFastqFlipped[value] = [key] 
    else: 
        controlFastqFlipped[value].append(key) 
  
# controlFastqFlipped dict now is of length the number of unique controls, with the samples using those controls as values

# mergedInputDit allows to create a generic name for the Inputs ( Input1, Input2 etc etc)
mergedInputDict = controlFastqFlipped
i = 1
for key, value in controlFastqFlipped.items():
    inputname = "Input" + str(i)
    mergedInputDict[key] = inputname
    i = i + 1

# Now creating CONTROL_SAMPLE_DICT,  linking sample with their unique Input using the generic input name e.g. {Mousekidney: Input1}
CONTROL_SAMPLE_DICT = {}
for key, value in controlFastq.items():
    CONTROL_SAMPLE_DICT[key] = mergedInputDict[value]


# Flipping the flipped dictionary to have a link between generic input name and their corresponding fastq/Bam files. Splitting back the file name
CONTROL_MERGED_FILES = {}
for key, value in controlFastqFlipped.items():
    CONTROL_MERGED_FILES[value] = key.split(" ")

# Creating the dictionnary for all sample name linked to their input file
CASES_SAMPLE_FILES = {}
for case in CASES:
    sample = "_".join(case.split("_")[0:-1])
    mark = case.split("_")[-1]
    CASES_SAMPLE_FILES[case] = FILES[sample][mark]

ALL_SAMPLE_FILES = {**CASES_SAMPLE_FILES, **CONTROL_MERGED_FILES}

CONTROLS = list(CONTROL_MERGED_FILES.keys())

# ~~~~~~~~~~~~~~ All samples ~~~~~~~~~~~~~~ #
ALL_SAMPLES = CASES + CONTROLS

# ======================================================== #
# ============= Creating helper dictionaries ============= #
# ======================================================== # 

# ~~~~~~~~~~~~~~ Samples dict ~~~~~~~~~~~~~ #

# Regroup Marks or TF by sample
# e.g. Mousekidney01: [H3K27, H3K27me3], Mouseliver04: [H3K27, H3K27me3]
SAMPLES = dict()
for sample in sorted(FILES.keys()):
    for mark in FILES[sample].keys():
        if(mark not in CONTROL_NAME):
            SAMPLES.setdefault(sample, []).append(mark)

# Adding the proper merged input to the marks of each samples
for key in SAMPLES.keys():
    SAMPLES.setdefault(key,[]).append(CONTROL_SAMPLE_DICT[key])

# ~~~~~~~~~~~ Samples_name dict ~~~~~~~~~~~ #

# Regroup Marks or TF by sample using the full sample name SAMPLE_MARK
# e.g. Mousekidney01: [Mousekidney01_H3K27, Mousekidney01_H3K27me3], Mouseliver04: [_Mouseliver04_H3K27, Mouseliver04_H3K27me3]
SAMPLES_COMPLETE_NAME = dict()
for sample in sorted(FILES.keys()):
    for mark in FILES[sample].keys():
        if(mark not in CONTROL_NAME):
            SAMPLES_COMPLETE_NAME.setdefault(sample, []).append(sample + "_" + mark)

# Adding the proper merged input to the marks of each samples
for key in SAMPLES_COMPLETE_NAME.keys():
    SAMPLES_COMPLETE_NAME.setdefault(key,[]).append(CONTROL_SAMPLE_DICT[key])

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

# Here create a list with all marks without Input/Control mentionned, before adding the controls to the dict
MARKS_NO_CONTROL_COMPLETE_NAME = list(MARKS_COMPLETE_NAME.keys())

# Adding the key for the merged input
for key, value in CONTROL_SAMPLE_DICT.items():
    MARKS_COMPLETE_NAME.setdefault(value, []).append(key)


###########################################################################
########################### Listing OUTPUT FILES ##########################
###########################################################################

# Not all of those output files will be use in the snakemake rules but it's good to have them all at the same place #


# ~~~~~~ files with case and control ~~~~~~ #
ALL_PEAKS = []
ALL_BIGWIG_INPUT = []
ALL_BIGBED = []
ALL_FEATURECOUNTS = []
ALL_BROADPEAKCOUNTS = []
ALL_NARROWPEAKCOUNTS = []
ALL_ANNOTATED_PEAKS = []

# going through all cases samples (sample_mark) #
for case in CASES:
    sample = "_".join(case.split("_")[0:-1])
    control = CONTROL_SAMPLE_DICT[sample]
    if control in CONTROLS:
        ALL_PEAKS.append(os.path.join(WORKDIR, "08peak_macs1/{}-vs-{}-macs1-narrow_peaks.bed").format(case, control))
        ALL_PEAKS.append(os.path.join(WORKDIR, "09peak_macs2/{}-vs-{}-macs2_peaks.broadPeak").format(case, control))
        ALL_BIGWIG_INPUT.append(os.path.join(WORKDIR, "06bigwig_input/{}-vs-{}.bw").format(case, control))
        ALL_BIGBED.append(os.path.join(WORKDIR, "12UCSC/{}-vs-{}-macs2_peaks.bb").format(case, control))
        ALL_FEATURECOUNTS.append(os.path.join(WORKDIR, "DPQC/{}-vs-{}.FRiP.summary").format(case,control))
        ALL_BROADPEAKCOUNTS.append(os.path.join(WORKDIR, "DPQC/{}-vs-{}-broadpeak-count_mqc.json").format(case,control))
        ALL_NARROWPEAKCOUNTS.append(os.path.join(WORKDIR, "DPQC/{}-vs-{}-narrowpeak-count_mqc.json").format(case,control))
        ALL_ANNOTATED_PEAKS.append(os.path.join(WORKDIR, "annotate/{}-vs-{}-peaks_annotated.txt").format(case,control))
# ~~~~~~~~~~~~~~~ Bam files ~~~~~~~~~~~~~~~ #
CONTROL_BAM = expand(os.path.join(WORKDIR, "03aln/{sample}.sorted.bam"), sample = CONTROL_MERGED_FILES)
CASE_BAM = expand(os.path.join(WORKDIR, "03aln/{sample}.sorted.bam"), sample = CASES)
ALL_BAM     = CONTROL_BAM + CASE_BAM


# ~~ All samples files (cases + control) ~~ #
ALL_DOWNSAMPLE_BAM = expand(os.path.join(WORKDIR, "04aln_downsample/{sample}-downsample.sorted.bam"), sample = ALL_SAMPLES)
ALL_FASTQ   = expand(os.path.join(WORKDIR, "01seq/{sample}.fastq"), sample = ALL_SAMPLES)
ALL_FASTQC  = expand(os.path.join(WORKDIR, "02fqc/{sample}_fastqc.zip"), sample = ALL_SAMPLES)
ALL_BOWTIE_LOG = expand(os.path.join(WORKDIR, "00log/{sample}.align"), sample = ALL_SAMPLES)
ALL_INDEX = expand(os.path.join(WORKDIR, "03aln/{sample}.sorted.bam.bai"), sample = ALL_SAMPLES)
ALL_DOWNSAMPLE_INDEX = expand(os.path.join(WORKDIR, "04aln_downsample/{sample}-downsample.sorted.bam.bai"), sample = ALL_SAMPLES)
ALL_FLAGSTAT = expand(os.path.join(WORKDIR, "03aln/{sample}.sorted.bam.flagstat"), sample = ALL_SAMPLES)
ALL_PHANTOM = expand(os.path.join(WORKDIR, "05phantompeakqual/{sample}.spp.out"), sample = ALL_SAMPLES)
ALL_BIGWIG = expand(os.path.join(WORKDIR, "07bigwig/{sample}.bw"), sample = ALL_SAMPLES)
ALL_ENCODE = expand(os.path.join(WORKDIR, "DPQC/{sample}.encodeQC.txt"), sample = ALL_SAMPLES)


# ~~~~~~~~~~~ Deeptools specific ~~~~~~~~~~ #
# ---- Grouped by marks ---- #
ALL_COMPUTEMATRIX = expand(os.path.join(WORKDIR, "DPQC/{mark}.computeMatrix.gz"), mark = MARKS)
ALL_DPQC_PLOT = expand(os.path.join(WORKDIR, "DPQC/{mark}.plotHeatmap.png"), mark = MARKS)

# --- Grouped by samples --- #
ALL_DPQC_PLOT.extend(expand(os.path.join(WORKDIR, "DPQC/{samp}.fingerprint.png"), samp = SAMPLES))
ALL_DPQC = expand(os.path.join(WORKDIR, "DPQC/{samp}.plotFingerprintOutRawCounts.txt"), samp = SAMPLES)
ALL_DPQC.extend(expand(os.path.join(WORKDIR, "DPQC/{samp}.plotFingerprintOutQualityMetrics.txt"), samp = SAMPLES))

# ~~~~~~~~~~~ ChromHMM specific ~~~~~~~~~~~ #
if config["chromHMM"]:

    def get_chr(chromSize):
        with open(chromSize, 'r') as fs:
            chr = [line.rstrip().split('\t')[0] for line in fs]
            return(chr)
    # Read histone
    HISTONE_INCLUDED = config["histone_for_chromHMM"].split(" ")
    HISTONE_CASES = [sample for sample in CASES if sample.split("_")[-1] in HISTONE_INCLUDED ]
    HISTONE_SAMPLE = list(set([sample.split("_")[0] for sample in CASES if sample.split("_")[-1] in HISTONE_INCLUDED ]))
    BAM_TO_BED = expand(os.path.join(WORKDIR, "bamtobed/{sample}.bed"), sample = HISTONE_CASES + CONTROLS)
    CHROMHMM = expand(os.path.join(WORKDIR, "chromHMM/learn_{nb_state}_states/{sample}_{nb_state}_segments.bed"), sample = HISTONE_SAMPLE, nb_state = config["state"])
    CHRHMM = get_chr(config['chromHmm_g'])
    CHROMHMM_TABLE = [os.path.join(WORKDIR, "chromHMM/cellmarkfiletable.txt")]

# ~~~~~~~~~~~~~~~~~~ Misc ~~~~~~~~~~~~~~~~~ #
ALL_MULTIQC = [os.path.join(WORKDIR, "10multiQC/multiqc_report.html")]

ALL_CONFIG= [os.path.join(WORKDIR, "03aln/bams.json")]

HUB_FOLDER = os.path.join(WORKDIR, "UCSC_HUB")
ALL_HUB = [os.path.join(HUB_FOLDER,"{}.hub.txt").format(PROJECT_NAME)]

# ======================================================== #
# ==================== MULTIQC inputs ==================== #
# ======================================================== #

# Depends on bam or fastq as input #
if not BAM_INPUT:
    ALL_MULTIQC_INPUT = ALL_PHANTOM + ALL_DPQC + ALL_FEATURECOUNTS + ALL_BROADPEAKCOUNTS + ALL_NARROWPEAKCOUNTS + ALL_ENCODE + ALL_FASTQC + ALL_BOWTIE_LOG 
else:
    ALL_MULTIQC_INPUT = ALL_PHANTOM + ALL_DPQC + ALL_FEATURECOUNTS + ALL_ENCODE + ALL_BROADPEAKCOUNTS + ALL_NARROWPEAKCOUNTS

###########################################################################
########################### Targets for rule all ##########################
###########################################################################
TARGETS = []
TARGETS.extend(ALL_ANNOTATED_PEAKS)
TARGETS.extend(ALL_MULTIQC)
TARGETS.extend(ALL_HUB)

# Since output from bam input are not used as input, needs to be put in the rule all for execution #
if not BAM_INPUT:
    TARGETS.extend(ALL_CONFIG)

#TEMP
TARGETS.extend(ALL_DPQC_PLOT)


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
            bigwigs.append(os.path.join(WORKDIR, "07bigwig/" + wildcards.mark + ".bw"))
    else:
        for s in samples:
            bigwigs.append(os.path.join(WORKDIR, "07bigwig/" + s + ".bw"))
    return bigwigs

# ~~~~~ Aggregation of bams per sample ~~~~ #
def get_bams_per_sample(wildcards):
    marks = SAMPLES_COMPLETE_NAME[wildcards.samp]
    bams = list()
    for s in marks:
        bams.append(os.path.join(WORKDIR, "03aln/" + s + ".sorted.bam"))
    return bams

# ~~~ Aggregation of bam idx per sample ~~~ #
def get_bam_index_per_sample(wildcards):
    marks = SAMPLES_COMPLETE_NAME[wildcards.samp]
    bams = list()
    for s in marks:
        bams.append(os.path.join(WORKDIR, "03aln/" + s + ".sorted.bam.bai"))
    return bams

# ~~~~~~~~ marks or tf per samples ~~~~~~~~ #
def get_all_marks_per_sample(wildcards):
    return SAMPLES_COMPLETE_NAME[wildcards.samp]

# ~~~~~~ fastq files for sample_mark ~~~~~~ #
def get_fastq(wildcards):
    return ALL_SAMPLE_FILES[wildcards.sample]

# ~~~~~~~ bam files for sample_mark ~~~~~~~ #
def get_bams(wildcards):
    return ALL_SAMPLE_FILES[wildcards.sample]

# ~~~~~~~ peak sets per marks or tf ~~~~~~~ #
def get_peaks(wildcards):
    for key, value in MARKS_COMPLETE_NAME.items():
        if any(wildcards.case in sample_mark for sample_mark in MARKS_COMPLETE_NAME[key]):
            if key in NARROW_BROAD:
                if(NARROW_BROAD[key] == 'narrow'):
                    return os.path.join(WORKDIR, "08peak_macs1/" + wildcards.case + "-vs-" + wildcards.control + "-macs1-narrow_peaks.bed")
                elif(NARROW_BROAD[key] == 'broad'):
                    return os.path.join(WORKDIR, "09peak_macs2/" + wildcards.case + "-vs-" + wildcards.control + "-macs2_peaks.broadPeak")
            else:
                logger.warning("Marks or TF not in the {} for {}. Will work with narrow peaks".format(config['narrow_broad'],wildcards.case))
                return os.path.join(WORKDIR, "08peak_macs1/" + wildcards.case + "-vs-" + wildcards.control + "-macs1-narrow_peaks.bed")

#TODO see if this is usefull
localrules: all

rule all:
    input: TARGETS

###########################################################################
######################### Alignement using BOWTIE2 ########################
###########################################################################

# Those rules are only executed if inputs are fastq #
if BAM_INPUT == False:

    #Now only for single-end ChIPseq
    rule merge_fastqs:
        input: get_fastq
        output: os.path.join(WORKDIR, "01seq/{sample}.fastq")
        log: os.path.join(WORKDIR, "00log/{sample}.unzip")
        params: jobname = "{sample}"
        shell: "gunzip -c {input} > {output} 2> {log}"

    rule fastqc:
        input:  os.path.join(WORKDIR, "01seq/{sample}.fastq")
        output: os.path.join(WORKDIR, "02fqc/{sample}_fastqc.zip"), os.path.join(WORKDIR, "02fqc/{sample}_fastqc.html")
        log:    os.path.join(WORKDIR, "00log/{sample}.fastqc")
        params:
            output_dir = os.path.join(WORKDIR, "02fqc")
        shell:
            """
            source activate full-pipe-main-env
            fastqc -o {params.output_dir} -f fastq --noextract {input} 2> {log}
            """

    # Simple alignment with bowtie 2 followed by sorting #
    rule align:
        input:  os.path.join(WORKDIR, "01seq/{sample}.fastq")
        output: os.path.join(WORKDIR, "03aln/{sample}.temp.bam")
        log:    os.path.join(WORKDIR, "00log/{sample}.align")
        shell:
            """
            source activate full-pipe-main-env
            bowtie2 -p 8 -x {config[idx_bt2]} -q {input} 2> {log} \
            | samblaster \
            | samtools view -bu - \
            | samtools sort -m 8G -@ 4 -T {output}.tmp -o {output}
            """

    # Get the duplicates marked sorted bam, remove unmapped reads by samtools view -F 1804 #
    # Samblaster should run before samtools sort #
    rule filter_alignment:
        input:  os.path.join(WORKDIR, "03aln/{sample}.temp.bam")
        output: os.path.join(WORKDIR, "03aln/{sample}.sorted.bam")
        log:    os.path.join(WORKDIR, "00log/{sample}.filter")
        shell:
            """
            source activate full-pipe-main-env
            samtools view -bu -F 1804 {input} -q 30 \
            | samtools sort -m 8G -@ 4 -T {output}.tmp -o {output} 2> {log}
            """
    
    # This rule is not followed by other rules, so its output has to be added to the rule all conditionally on BAM_INPUT, if Bam are used as input or not #
    rule create_bam_json:
        input: expand(os.path.join(WORKDIR, "03aln/{sample}.sorted.bam"), sample = ALL_SAMPLES),
        output: os.path.join(WORKDIR, "03aln/bams.json")
        params: SAMPLES
        run:
            #sample dictionary that will be dumped as a json
            dict_for_json = {}
            #Going through all samples usign the sample dictionary -> {sample1: [mark1, mark2],sample2: [mark1, mark2])
            for samp in SAMPLES:
                #Mini dictionary containing marks and associated sorted bam filepath
                mini_dict = {}
                for mark in SAMPLES[samp]:

                    filepath  =  os.path.join(WORKDIR, "03aln/" + samp + "_" + mark + ".sorted.bam")
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
        output: os.path.join(WORKDIR, "03aln/{sample}.sorted.bam")
        shell:
            """
            ln -s {input} {output}
            """

# ~~~~~~~~~~~~~ Indexing bams ~~~~~~~~~~~~~ #
rule index_bam:
    input:  os.path.join(WORKDIR, "03aln/{sample}.sorted.bam")
    output: os.path.join(WORKDIR, "03aln/{sample}.sorted.bam.bai")
    log:    os.path.join(WORKDIR, "00log/{sample}.bam.index")
    threads: 1
    shell:
        """
        source activate full-pipe-main-env
        samtools index {input} 2> {log}
        """

###########################################################################
############################### Downsampling ##############################
###########################################################################
#  Using user provided parameters, bam will be downsampled. Flagstat is used for read counting and fed to sambamba  #

# flagstat
rule flagstat_bam:
    input:  os.path.join(WORKDIR, "03aln/{sample}.sorted.bam")
    output: os.path.join(WORKDIR, "03aln/{sample}.sorted.bam.flagstat")
    log:    os.path.join(WORKDIR, "00log/{sample}.bam.flagstat")
    threads: 1
    params: jobname = "{sample}"
    message: "flagstat_bam {input}: {threads} threads"
    shell:
        """
        source activate full-pipe-main-env
        samtools flagstat {input} > {output} 2> {log}
        """

#downsampling
rule down_sample:
    input: 
        bam = os.path.join(WORKDIR, "03aln/{sample}.sorted.bam"), 
        bai = os.path.join(WORKDIR, "03aln/{sample}.sorted.bam.bai"), 
        flagstat = os.path.join(WORKDIR, "03aln/{sample}.sorted.bam.flagstat")
    output: 
        bam = os.path.join(WORKDIR, "04aln_downsample/{sample}-downsample.sorted.bam"), 
        bai = os.path.join(WORKDIR, "04aln_downsample/{sample}-downsample.sorted.bam.bai")
    log: os.path.join(WORKDIR, "00log/{sample}.downsample.log")
    threads: 5
    params: 
    log: os.path.join(WORKDIR, "00log/{sample}.phantompeakqual.log")
    message: "downsampling for {input}"
    shell:
        """
        source activate full-pipe-main-env
        sambamba view -f bam -t 5 --subsampling-seed=3 -s `sed '5q;d' {input.flagstat} | cut -d" " -f1 | awk '{{ratio = {config[target_reads]}/$0}};{{if(ratio < 1 )print ratio; else print 1}}'` {input.bam} | samtools sort -m 2G -@ 5 -T {output.bam}.tmp > {output.bam} 2> {log}
        samtools index {output.bam}
        """

###########################################################################
#################################### QC ###################################
###########################################################################

#Library complexity
rule encode_complexity:
    input: 
        bam = os.path.join(WORKDIR, "03aln/{sample}.temp.bam") 
    output: 
        os.path.join(WORKDIR, "DPQC/{sample}.encodeQC.txt")
    threads: 4
    shell:
        """
        source activate full-pipe-main-env
        bedtools bamtobed -i {input} | awk 'BEGIN{{OFS="\t"}}{{print $1,$2,$3,$6}}' | grep -v 'chrM' | sort | uniq -c \
        | awk 'BEGIN{{mt=0;m0=0;m1=0;m2=0;OFS="\t"}} ($1==1){{m1=m1+1}} ($1==2){{m2=m2+1}} {{m0=m0+1}} {{mt=mt+$1}} \
        END{{m1_m2=-1.0; if(m2>0) m1_m2=m1/m2; print "Sample Name","NRF","PBC1","PBC2"; print "{wildcards.sample}",m0/mt,m1/m0,m1_m2}}' > {output}
        """

# Phantompeakqualtools computes a robust fragment length using the cross correlation (xcor) metrics.
rule phantom_peak_qual:
    input: 
        bam = os.path.join(WORKDIR, "03aln/{sample}.sorted.bam"), 
        bai = os.path.join(WORKDIR, "03aln/{sample}.sorted.bam.bai")
    output: os.path.join(WORKDIR, "05phantompeakqual/{sample}.spp.out")
    log: os.path.join(WORKDIR, "00log/{sample}.phantompeakqual.log")
    threads: 4
    params: os.path.join(WORKDIR, "05phantompeakqual/")
    message: "phantompeakqual for {input}"
    shell:
        """
        source activate full-pipe-spp
        run_spp -c={input.bam} -savp -rf -p=4 -odir={params}  -out={output} -tmpdir={params} 2> {log}
        """

# Deeptools QC
rule computeMatrix_QC:
    input : get_big_wig_with_mark_or_tf 
    output : os.path.join(WORKDIR, "DPQC/{mark}.computeMatrix.gz")
    params : GENOME_TSS
    shell:
        """
        source activate full-pipe-main-env
	    computeMatrix reference-point -S {input} -R {params} -a 3000 -b 3000 -out {output} --numberOfProcessors max/2
        """

# Deeptools QC
rule plotHeatmap:
    input :  os.path.join(WORKDIR, "DPQC/{mark}.computeMatrix.gz")
    output : os.path.join(WORKDIR, "DPQC/{mark}.plotHeatmap.png")
    shell:
        """
        source activate full-pipe-main-env
        plotHeatmap -m {input} -out {output} --colorMap jet
        """

# ChipSeq QCs plots from deeptools. Plotfingerprints are really usefull to see focal enrichment of your CHip-Seq enrichment
rule plotFingerPrint:
    input:
        bam = get_bams_per_sample, 
        bai = get_bam_index_per_sample
    output:
        plot = os.path.join(WORKDIR, "DPQC/{samp}.fingerprint.png"), 
        rawCounts = os.path.join(WORKDIR, "DPQC/{samp}.plotFingerprintOutRawCounts.txt"),
        qualityMetrics = os.path.join(WORKDIR, "DPQC/{samp}.plotFingerprintOutQualityMetrics.txt")
    params: 
        labels = get_all_marks_per_sample
    shell:
        """
        source activate full-pipe-main-env
        plotFingerprint -b {input.bam} --plotFile {output.plot} --labels {params.labels} --region chr1 --skipZeros --numberOfSamples 100000 --minMappingQuality 30 --plotTitle {wildcards.samp} --outRawCounts {output.rawCounts} --outQualityMetrics {output.qualityMetrics}
        """

rule get_FRiP_for_multiqc:
    input:
        peaks = get_peaks,
        bam = os.path.join(WORKDIR, "03aln/{case}.sorted.bam"), 
    output:
        os.path.join(WORKDIR, "DPQC/{case}-vs-{control}.FRiP.summary")
    params:
        saf = os.path.join(WORKDIR, "DPQC/{case}.saf"),
        outputName = os.path.join(WORKDIR, "DPQC/{case}-vs-{control}.FRiP")
    shell:
        """
        source activate full-pipe-main-env
        awk 'BEGIN{{OFS="\t";print "GeneID", "Chr","Start","End","Strand"}}{{print $4,$1,$2,$3,$6}}' {input.peaks} > {params.saf}
        featureCounts -a {params.saf} -F SAF -o {params.outputName} {input.bam}
        """

rule get_broad_peak_counts_for_multiqc:
    input:
        peaks = os.path.join(WORKDIR, "09peak_macs2/{case}-vs-{control}-macs2_peaks.broadPeak"),
    output:
        os.path.join(WORKDIR, "DPQC/{case}-vs-{control}-broadpeak-count_mqc.json")
    params:
        peakType = "broadPeak"
    shell:
        """
        source activate full-pipe-main-env
        python3 scripts/count_peaks.py --peak_type {params.peakType} --peaks {input.peaks} --sample_name {wildcards.case} > {output}
        """

rule get_narrow_peak_counts_for_multiqc:
    input:
        peaks = os.path.join(WORKDIR, "08peak_macs1/{case}-vs-{control}-macs1-narrow_peaks.bed"),
    output:
        os.path.join(WORKDIR, "DPQC/{case}-vs-{control}-narrowpeak-count_mqc.json")
    params:
        peakType = "narrowPeak"
    shell:
        """
        source activate full-pipe-main-env
        python3 scripts/count_peaks.py --peak_type {params.peakType} --peaks {input.peaks} --sample_name {wildcards.case} > {output}
        """


###########################################################################
############################### PEAK CALLING ##############################
###########################################################################

# Peak calling using MACS
rule call_narrow_peaks_macs1:
    input: 
        control = os.path.join(WORKDIR, "04aln_downsample/{control}-downsample.sorted.bam"), 
        case = os.path.join(WORKDIR, "04aln_downsample/{case}-downsample.sorted.bam"),
        spp = os.path.join(WORKDIR, "05phantompeakqual/{case}.spp.out")
    output:
        os.path.join(WORKDIR, "08peak_macs1/{case}-vs-{control}-macs1-narrow_peaks.bed")
    log:
        macs1_nomodel = os.path.join(WORKDIR, "00log/{case}-vs-{control}-call-peaks-macs1-nomodel.log")
    params:
        name = "{case}-vs-{control}-macs1-narrow",
        jobname = "{case}", 
        outdir = os.path.join(WORKDIR, "08peak_macs1/")
    message: "Calling narrow peaks with macs14."
    shell:
        """
        source activate full-pipe-macs
        # nomodel with shiftsize half of the estimated fragment length from phantompeakqual.
        macs -t {input.case} \
            -c {input.control} --keep-dup all -f BAM -g {config[macs_g]} \
            --outdir {params.outdir} -n {params.name} --shiftsize `cut -f3 {input.spp} | awk 'BEGIN{{FS=","}}{{printf "%.0f",($1+1)/2}}'` --nomodel -p {config[macs_pvalue]} &> {log.macs1_nomodel}
        """
        
# Peak calling using MACS 2
rule call_broad_peaks_macs2:
    input: 
        control = os.path.join(WORKDIR, "04aln_downsample/{control}-downsample.sorted.bam"), 
        case = os.path.join(WORKDIR, "04aln_downsample/{case}-downsample.sorted.bam"),
        spp = os.path.join(WORKDIR, "05phantompeakqual/{case}.spp.out")
    output:
        broad = os.path.join(WORKDIR, "09peak_macs2/{case}-vs-{control}-macs2_peaks.broadPeak")
    log: os.path.join(WORKDIR, "00log/{case}-vs-{control}-call-peaks_macs2.log")
    params:
        name = "{case}-vs-{control}-macs2", 
        jobname = "{case}", 
        outdir = os.path.join(WORKDIR, "09peak_macs2")
    message: "Calling broadpeaks with macs2."
    shell:
        """
        source activate full-pipe-macs
        ## for macs2, when nomodel is set, --extsize is default to 200bp, this is the same as 2 * shift-size in macs14.
        macs2 callpeak -t {input.case} \
            -c {input.control} --keep-dup all -f BAM -g {config[macs2_g]} \
            --outdir {params.outdir} -n {params.name} --extsize `cut -f3 {input.spp} | awk 'BEGIN{{FS=","}}{{print $1}}'` -p {config[macs2_pvalue]} --broad --broad-cutoff {config[macs2_pvalue_broad]} --nomodel &> {log}
        """


###########################################################################
####### VISUALIZATION bigWig and bigBed generation and HUB creation #######
###########################################################################

rule get_bigwigs_using_inputs:
    input : 
        case =  os.path.join(WORKDIR, "04aln_downsample/{case}-downsample.sorted.bam"),
        bai_case = os.path.join(WORKDIR, "04aln_downsample/{case}-downsample.sorted.bam.bai"),
        control = os.path.join(WORKDIR, "04aln_downsample/{control}-downsample.sorted.bam"), 
        bai_control = os.path.join(WORKDIR, "04aln_downsample/{control}-downsample.sorted.bam.bai"),
        spp = os.path.join(WORKDIR, "05phantompeakqual/{case}.spp.out")
    output:  os.path.join(WORKDIR, "06bigwig_input/{case}-vs-{control}.bw")
    log: os.path.join(WORKDIR, "00log/{case}-vs-{control}.makebw")
    threads: 5
    params: jobname = "{case}"
    message: "Making bigwig of {case} log2 fold change versus {control}"
    shell:
        """
        source activate full-pipe-main-env
        bamCompare --bamfile1 {input.case} --bamfile2 {input.control} \
        --normalizeUsing RPKM  --operation log2 --operation first --scaleFactorsMethod None --binSize 10 --smoothLength 30 --numberOfProcessors max/2 \
        --extendReads `cut -f3 {input.spp} | awk 'BEGIN{{FS=","}}{{print $1}}'` -o {output} 2> {log}
        """

rule get_bigwigs:
    input : 
        bam = os.path.join(WORKDIR, "04aln_downsample/{sample}-downsample.sorted.bam"),
        bai = os.path.join(WORKDIR, "04aln_downsample/{sample}-downsample.sorted.bam.bai"),
        spp = os.path.join(WORKDIR, "05phantompeakqual/{sample}.spp.out")
    output: os.path.join(WORKDIR, "07bigwig/{sample}.bw")
    log: os.path.join(WORKDIR, "00log/{sample}.makebw")
    threads: 5
    params: jobname = "{sample}"
    message: "Making bigwig of {sample}"
    shell:
        """
        source activate full-pipe-main-env
        bamCoverage -b {input.bam} --normalizeUsing RPKM --binSize 10 --smoothLength 30 -p 5 --numberOfProcessors max/2 \
        --extendReads `cut -f3 {input.spp} | awk 'BEGIN{{FS=","}}{{print $1}}'` -o {output} 2> {log}
        """

# Cleaning peaks by taking only columns 1, 2, 3 because narrowpeaks and broadpeaks are different.
rule get_bigbeds:
    input: get_peaks
    output: os.path.join(WORKDIR, "12UCSC/{case}-vs-{control}-macs2_peaks.bb")
    params : 
        bed1 = temp(os.path.join(WORKDIR, "12UCSC_broad/{case}-vs-{control}-macs2_peaks.bed"))
    shell:
        """
        cut -f1,2,3 {input} | LC_COLLATE=C sort -k1,1 -k2,2n > {params.bed1}
        scripts/bedToBigBed {params.bed1} {GENOME_SIZE} {output} -type=bed3
        """

# Creating a Hub for UCSC
rule get_UCSC_hub:
    input:  
        bed = ALL_BIGBED, 
        bigwig = ALL_BIGWIG_INPUT
    output:
        ALL_HUB
    log: os.path.join(WORKDIR, "00log/log.trackhub")
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
    input : "multiqc_config.yaml"
    output: os.path.join(WORKDIR, "10multiQC/multiqc_config.yaml")
    message: "Moving multiqc config"
    shell:
        """
        sed "s/DATE/$(date)/g" {input} | sed "s/PROJECTNAME/{PROJECT_NAME}/g"  > {output}
        """

# MultiQC: Takes fastqc, phantompeakqualtools, deeptools, bowtie, fastQC and feature counts for FRiP and custom MACS2 peak counts
rule multiQC:
    input : 
        multiqc_files = ALL_MULTIQC_INPUT,
        multiqc_config = os.path.join(WORKDIR, "10multiQC/multiqc_config.yaml")
    output: ALL_MULTIQC
    params: os.path.join(WORKDIR, "10multiQC/")
    log: os.path.join(WORKDIR, "00log/multiqc.log")
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
        annotated_peaks = os.path.join(WORKDIR, "annotate/{case}-vs-{control}-peaks_annotated.txt")
    params:
        genome = GENOME_FASTA, 
        gtf = GENOME_GTF
    log: os.path.join(WORKDIR,"00log{case}-vs-{control}-homer-annotated.log"
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
            os.path.join(WORKDIR, "04aln_downsample/{sample}-downsample.sorted.bam")
        output:
            os.path.join(WORKDIR, "bamtobed/{sample}.bed")
        params: jobname = "{sample}"
        log: os.path.join(WORKDIR, "00log/{sample}-bam2bed.log")
        message: "converting bam to bed for {input}"
        shell:
            """
            source activate full-pipe-chromhmm
            bedtools bamtobed -i {input} > {output}
            """

    rule make_table:
        input : 
            expand(os.path.join(WORKDIR, "bamtobed/{sample}.bed"), sample = HISTONE_CASES + CONTROLS)
        output : 
            os.path.join(WORKDIR, "chromHMM/cellmarkfiletable.txt")
        log: os.path.join(WORKDIR, "00log/make_table_chromHMM.log")
        message: "making a table for chromHMM"
        run:
            import os
            from os.path import join
            with open (output[0], "w") as f:
                for case in HISTONE_CASES:
                    sample = "_".join(case.split("_")[0:-1])
                    mark = case.split("_")[-1]
                    control = sample + "_" + CONTROL_NAME
                    case_bed = case + ".bed"
                    if os.path.exists(join(os.path.join(WORKDIR, "bamtobed"), case_bed)):
                        f.write(sample + "\t" +  mark + "\t" + case + ".bed" + "\t" + control + ".bed" + "\n")

    rule chromHMM_binarize:
        input :
            cellmarkfiletable = os.path.join(WORKDIR, "chromHMM/cellmarkfiletable.txt"), 
            beds = expand(os.path.join(WORKDIR, "bamtobed/{sample}.bed"), sample = HISTONE_CASES + CONTROLS)
        output:
            expand(os.path.join(WORKDIR, "chromHMM/binarizedData/{sample}-{chr}-binary.txt"), sample = HISTONE_SAMPLE, chr = CHRHMM)
        log:
            os.path.join(WORKDIR, "00log/chromhmm_bin.log")
        params:
            folder = os.path.join(WORKDIR, "chromHMM/binarizedData/"), 
            bamtobed_folder = os.path.join(WORKDIR, "bamtobed/"), 
            memory = "32G"
        shell:
            """
            source activate full-pipe-chromhmm
            ChromHMM.sh -Xmx{params.memory} BinarizeBed -b {config[binsize]} {config[chromHmm_g]} {params.bamtobed_folder} {input.cellmarkfiletable} {params.folder} 2> {log}
            """

    rule chromHMM_learn:
        input:
            expand(os.path.join(WORKDIR, "chromHMM/binarizedData/{sample}-{chr}-binary.txt"), sample = HISTONE_SAMPLE, chr = CHRHMM)
        output: 
            expand(os.path.join(WORKDIR, "chromHMM/learn_{nb_state}_states/{sample}_{nb_state}_segments.bed"), sample = HISTONE_SAMPLE, nb_state = config["state"])
        log: os.path.join(WORKDIR, "00log/chromhmm_learn.log")
        params:
            input_folder = os.path.join(WORKDIR, "chromHMM/binarizedData/"), 
            output_folder = expand(os.path.join(WORKDIR, "chromHMM/learn_{nb_state}_states/"), nb_state = config["state"]), 
            memory = "32G"
        shell:
            """
            source activate full-pipe-chromhmm
            unset DISPLAY && ChromHMM.sh -Xmx{params.memory} LearnModel -p 0 -b {config[binsize]} {params.input_folder} {params.output_folder} {config[state]} {config[chromHmm_g]} 2> {log}
            """
