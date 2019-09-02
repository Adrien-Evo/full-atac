import csv
import os
import json
import numpy as np


shell.prefix("set -eo pipefail; echo BEGIN at $(date); ")

configfile: "config_bam.yaml"
localrules: all
# localrules will let the rule run locally rather than submitting to cluster
# computing nodes, this is for very small jobs

# load cluster config file
CLUSTER = json.load(open(config['CLUSTER_JSON']))
FILES = json.load(open(config['SAMPLES_JSON']))
TSS_BED = config['tss_bed']
WORKDIR = os.path.abspath(config["OUTPUT_DIR"])
PROJECT_NAME = config['PROJECT_NAME']
BAM_INPUT = config['bam']

# ======================================================== #
#     Defining for the first time CASES and CONTROLS        #
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
# ====== Checking for number of Input/control sample ===== #
# ======================================================== #


#  Create a dictionary linking each sample with their control fastq e.g. { Mousekidney01 : controlIG16_TCGCTAGA_L001_R1_002.fastq.gz 16_TCGCTAGA_L001_R1_001.fastq.gz}  #
# joining the list to allow for usage asa dictionary key
controlFastq = dict()
for samp in SAMPLE_MARK:
    if CONTROL_NAME in samp:
        sample = "".join(samp.split("_")[0:-1])
        mark = samp.split("_")[-1]
        controlFastq[sample] = " ".join(FILES[sample][mark])
 
# finding duplicate values from controlFastq by flipping it 
controlFastqFlipped = {} 
  
for key, value in controlFastq.items():
    if value not in controlFastqFlipped: 
        controlFastqFlipped[value] = [key] 
    else: 
        controlFastqFlipped[value].append(key) 
  
# controlFastqFlipped dict now is of length the number of unique controls, with the samples using those controls as values

#mergedInputDit allows to create a generic name for the Inputs
mergedInputDict = controlFastqFlipped
i = 1
for key, value in controlFastqFlipped.items():
    inputname = "Input" + str(i)
    mergedInputDict[key] = inputname
    i = i + 1

#Now creating CONTROL_SAMPLE_DICT,  linking sample with their unique Input using the generic input name e.g. {Mousekidney: Input1}
CONTROL_SAMPLE_DICT = {}
for key, value in controlFastq.items():
    CONTROL_SAMPLE_DICT[key] = mergedInputDict[value]


#Flipping the flipped dictionary to have a link between generic input name and their corresponding file. Splitting back the file name
CONTROL_MERGED_FILES = {}
for key, value in controlFastqFlipped.items():
    CONTROL_MERGED_FILES[value] = key.split(" ")

#Creating the dictionnary for all sample name linked to their input file
CASES_SAMPLE_FILES = {}
for case in CASES:
    sample = "_".join(case.split("_")[0:-1])
    mark = case.split("_")[-1]
    CASES_SAMPLE_FILES[case] = FILES[sample][mark]

ALL_SAMPLE_FILES = {**CASES_SAMPLE_FILES, **CONTROL_MERGED_FILES}

CONTROLS = list(CONTROL_MERGED_FILES.keys())

# ======================================================== #
#  Creating dictionaries to associate Marks with their samples and samples with their marks  #
# ======================================================== #

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

SAMPLES_COMPLETE_NAME = dict()
for sample in sorted(FILES.keys()):
    for mark in FILES[sample].keys():
        if(mark not in CONTROL_NAME):
            SAMPLES_COMPLETE_NAME.setdefault(sample, []).append(sample + "_"+ mark)

# Adding the proper merged input to the marks of each samples
for key in SAMPLES_COMPLETE_NAME.keys():
    SAMPLES_COMPLETE_NAME.setdefault(key,[]).append(CONTROL_SAMPLE_DICT[key])


# Regroup samples per marks or TF
# e.g. H3K27: [Mousekidney01, Mouseliver04], H3K27me3: [Mousekidney01, Mouseliver04]
MARKS = dict()
for sample in sorted(FILES.keys()):
    for mark in FILES[sample].keys():
        if(mark not in CONTROL_NAME):
            MARKS.setdefault(mark, []).append(sample)

#Here create a list with all makrs without the control
MARKS_NO_CONTROL = list(MARKS.keys())

#Adding the key for the merged input
for key, value in CONTROL_SAMPLE_DICT.items():
    MARKS.setdefault(value, []).append(key)



#  list BAM files
CONTROL_BAM = expand(os.path.join(WORKDIR, "03aln/{sample}.sorted.bam"), sample = CONTROL_MERGED_FILES)
CASE_BAM = expand(os.path.join(WORKDIR, "03aln/{sample}.sorted.bam"), sample = CASES)

#  peaks and bigwigs
ALL_PEAKS = []
ALL_inputSubtract_BIGWIG = []
ALL_BROADPEAK = []
ALL_BIGWIGUCSC = []
ALL_FEATURECOUNTS = []

for case in CASES:
    sample = "_".join(case.split("_")[0:-1])
    control = CONTROL_SAMPLE_DICT[sample]
    if control in CONTROLS:
        ALL_PEAKS.append(os.path.join(WORKDIR, "08peak_macs1/{}-vs-{}-macs1_peaks.bed").format(case, control))
        ALL_PEAKS.append(os.path.join(WORKDIR, "08peak_macs1/{}-vs-{}-macs1-nomodel_peaks.bed").format(case, control))
        ALL_PEAKS.append(os.path.join(WORKDIR, "09peak_macs2/{}-vs-{}-macs2_peaks.xls").format(case, control))
        ALL_PEAKS.append(os.path.join(WORKDIR, "09peak_macs2/{}-vs-{}-macs2_peaks.broadPeak").format(case, control))
        ALL_inputSubtract_BIGWIG.append(os.path.join(WORKDIR, "06bigwig_inputSubtract/{}-subtract-{}.bw").format(case, control))
        ALL_BROADPEAK.append(os.path.join(WORKDIR, "12UCSC_broad/{}-vs-{}-macs2_peaks.broadPeak").format(case, control))
        ALL_BIGWIGUCSC.append(os.path.join(WORKDIR, "UCSC_compatible_bigWig/{}-subtract-{}.bw").format(case, control))
        ALL_FEATURECOUNTS.append(os.path.join(WORKDIR, "DPQC/{}-vs-{}.FRiP.summary").format(case,control))


ALL_SAMPLES = CASES + CONTROLS
ALL_BAM     = CONTROL_BAM + CASE_BAM
ALL_DOWNSAMPLE_BAM = expand(os.path.join(WORKDIR, "04aln_downsample/{sample}-downsample.sorted.bam"), sample = ALL_SAMPLES)
ALL_FASTQ   = expand(os.path.join(WORKDIR, "01seq/{sample}.fastq"), sample = ALL_SAMPLES)
ALL_FASTQC  = expand(os.path.join(WORKDIR, "02fqc/{sample}_fastqc.zip"), sample = ALL_SAMPLES)
ALL_BOWTIE_LOG = expand(os.path.join(WORKDIR, "00log/{sample}.align"), sample = ALL_SAMPLES)
ALL_INDEX = expand(os.path.join(WORKDIR, "03aln/{sample}.sorted.bam.bai"), sample = ALL_SAMPLES)
ALL_DOWNSAMPLE_INDEX = expand(os.path.join(WORKDIR, "04aln_downsample/{sample}-downsample.sorted.bam.bai"), sample = ALL_SAMPLES)
ALL_FLAGSTAT = expand(os.path.join(WORKDIR, "03aln/{sample}.sorted.bam.flagstat"), sample = ALL_SAMPLES)
ALL_PHANTOM = expand(os.path.join(WORKDIR, "05phantompeakqual/{sample}.spp.out"), sample = ALL_SAMPLES)
ALL_BIGWIG = expand(os.path.join(WORKDIR, "07bigwig/{sample}.bw"), sample = ALL_SAMPLES)

ALL_COMPUTEMATRIX = expand(os.path.join(WORKDIR, "DPQC/{mark}.computeMatrix.gz"), mark = MARKS)
ALL_DPQC_PLOT = expand(os.path.join(WORKDIR, "DPQC/{mark}.plotHeatmap.png"), mark = MARKS)
ALL_DPQC_PLOT.extend(expand(os.path.join(WORKDIR, "DPQC/{samp}.fingerprint.png"), samp = SAMPLES))
ALL_DPQC = expand(os.path.join(WORKDIR, "DPQC/{samp}.plotFingerprintOutRawCounts.txt"), samp = SAMPLES)
ALL_DPQC.extend(expand(os.path.join(WORKDIR, "DPQC/{samp}.plotFingerprintOutQualityMetrics.txt"), samp = SAMPLES))


ALL_QC = [os.path.join(WORKDIR, "10multiQC/multiQC_log.html")]

ALL_CONFIG= [os.path.join(WORKDIR, "03aln/bams.json")]
HUB_FOLDER = os.path.join(WORKDIR, "UCSC_HUB")
ALL_HUB = [os.path.join(HUB_FOLDER,"{}.hub.txt").format(PROJECT_NAME)]




if not BAM_INPUT:
    ALL_MULTIQC_INPUT = ALL_FLAGSTAT + ALL_PHANTOM + ALL_DPQC + ALL_FEATURECOUNTS + ALL_FASTQC + ALL_BOWTIE_LOG 
else:
    ALL_MULTIQC_INPUT = ALL_FLAGSTAT + ALL_PHANTOM + ALL_DPQC + ALL_FEATURECOUNTS



TARGETS = []
TARGETS.extend(ALL_PEAKS)
TARGETS.extend(ALL_QC)
TARGETS.extend(ALL_HUB)

if not BAM_INPUT:
    TARGETS.extend(ALL_CONFIG)


# Output files for ChromHMM
if config["chromHMM"]:

    def get_chr(chromSize):
        with open(chromSize, 'r') as fs:
            chr = [line.rstrip().split('\t')[0] for line in fs]
            return(chr)
    #Read histone
    HISTONE_INCLUDED = config["histone_for_chromHMM"].split(" ")
    HISTONE_CASES = [sample for sample in CASES if sample.split("_")[-1] in HISTONE_INCLUDED ]
    HISTONE_SAMPLE = list(set([sample.split("_")[0] for sample in CASES if sample.split("_")[-1] in HISTONE_INCLUDED ]))
    BAM_TO_BED = expand(os.path.join(WORKDIR, "bamtobed/{sample}.bed"), sample = HISTONE_CASES + CONTROLS)
    CHROMHMM = expand(os.path.join(WORKDIR, "chromHMM/learn_{nb_state}_states/{sample}_{nb_state}_segments.bed"), sample = HISTONE_SAMPLE, nb_state = config["state"])
    CHRHMM = get_chr(config['chromHmm_g'])
    CHROMHMM_TABLE = [os.path.join(WORKDIR, "chromHMM/cellmarkfiletable.txt")]
    #Adding the ChromHMM outputs to rule all
    TARGETS.extend(CHROMHMM)

# Aggregation of bigwigs by Marks or TF
def get_big_wig_with_mark_or_tf(wildcards):
    samples = MARKS[wildcards.mark]
    bigwigs = list()
    for s in samples:
        bigwigs.append(os.path.join(WORKDIR, "07bigwig/"+s+"_"+wildcards.mark+".bw"))
    return bigwigs

# Aggregation of bams per sample
def get_bams_per_sample(wildcards):
    marks = SAMPLES_COMPLETE_NAME[wildcards.samp]
    bams = list()
    for s in marks:
        bams.append(os.path.join(WORKDIR, "03aln/"+s+".sorted.bam"))
    return bams

# Aggregation of bam index per sample, necessary to avoid launching without bai
def get_bam_index_per_sample(wildcards):
    marks = SAMPLES_COMPLETE_NAME[wildcards.samp]
    bams = list()
    for s in marks:
        bams.append(os.path.join(WORKDIR, "03aln/"+s+".sorted.bam.bai"))
    return bams

# Just return all marks or tf for a single sample using wildcards
def get_all_marks_per_sample(wildcards):
    return SAMPLES[wildcards.samp]


#  Get a list of fastq.gz files for the same mark, same sample
def get_fastq(wildcards):
    return ALL_SAMPLE_FILES[wildcards.sample]

def get_bams(wildcards):
    return ALL_SAMPLE_FILES[wildcards.sample]



localrules: all

rule all:
    input: TARGETS

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
        conda:
            "envs/full-atac-main-env.yml"
        params :
            output_dir = os.path.join(WORKDIR, "02fqc")
        shell:
            """
            fastqc -o {params.output_dir} -f fastq --noextract {input} 2> {log}
            """

    # Get the duplicates marked sorted bam, remove unmapped reads by samtools view -F 4 and duplicated reads by samblaster -r
    # Samblaster should run before samtools sort
    rule align:
        input:  os.path.join(WORKDIR, "01seq/{sample}.fastq")
        output: os.path.join(WORKDIR, "03aln/{sample}.sorted.bam"), os.path.join(WORKDIR, "00log/{sample}.align")
        threads: CLUSTER["align"]["cpu"]
        conda:
            "envs/full-atac-main-env.yml"
        params:
            bowtie = " --chunkmbs 320 -m 1 --best -p 5 ", 
            jobname = "{sample}"
        message: "aligning {input}: 16 threads"
        log:
            bowtie = os.path.join(WORKDIR, "00log/{sample}.align"), 
            markdup = os.path.join(WORKDIR, "00log/{sample}.markdup")
        shell:
            """
            bowtie2 -p 4 -x {config[idx_bt1]} -q {input} 2> {log.bowtie} \
            | samblaster --removeDups \
        | samtools view -Sb -F 4 - \
        | samtools sort -m 8G -@ 4 -T {output[0]}.tmp -o {output[0]} 2> {log.markdup}
        """

    
    # This rule is not followed by other rules, so its output has to be added to the rule all conditionally on BAM_INPUT, if Bam are used as input or not
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



if BAM_INPUT:

    rule symlink_bam:
        input: get_bams
        output: os.path.join(WORKDIR, "03aln/{sample}.sorted.bam")
        shell:
            """
            ln -s {input} {output}
            """


# Simply indexing bams
rule index_bam:
    input:  os.path.join(WORKDIR, "03aln/{sample}.sorted.bam")
    output: os.path.join(WORKDIR, "03aln/{sample}.sorted.bam.bai")
    log:    os.path.join(WORKDIR, "00log/{sample}.bam.index")
    threads: 1
    conda:
        "envs/full-atac-main-env.yml"
    params: jobname = "{sample}"
    message: "index_bam {input}: {threads} threads"
    shell:
        """
        samtools index {input} 2> {log}
        """

# Check number of reads mapped by samtools flagstat, the output will be used for downsampling
rule flagstat_bam:
    input:  os.path.join(WORKDIR, "03aln/{sample}.sorted.bam")
    output: os.path.join(WORKDIR, "03aln/{sample}.sorted.bam.flagstat")
    log:    os.path.join(WORKDIR, "00log/{sample}.bam.flagstat")
    threads: 1
    conda:
        "envs/full-atac-main-env.yml"
    params: jobname = "{sample}"
    message: "flagstat_bam {input}: {threads} threads"
    shell:
        """
        samtools flagstat {input} > {output} 2> {log}
        """

# Phantompeakqualtools computes a robust fragment length using the cross correlation (xcor) metrics.
rule phantom_peak_qual:
    input: 
        bam = os.path.join(WORKDIR, "03aln/{sample}.sorted.bam"), 
        bai = os.path.join(WORKDIR, "03aln/{sample}.sorted.bam.bai")
    output: os.path.join(WORKDIR, "05phantompeakqual/{sample}.spp.out")
    log: os.path.join(WORKDIR, "00log/{sample}.phantompeakqual.log")
    threads: 4
    conda:
        "envs/spp.yml"
    params: os.path.join(WORKDIR, "05phantompeakqual/")
    message: "phantompeakqual for {input}"
    shell:
        """
        run_spp -c={input.bam} -savp -rf -p=4 -odir={params}  -out={output} -tmpdir={params} 2> {log}
        """
# Downsampling for relative comparisons of samples
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
    conda:
        "envs/full-atac-main-env.yml"
    params: 
    log: os.path.join(WORKDIR, "00log/{sample}.phantompeakqual.log")
    message: "downsampling for {input}"
    shell:
        """
        sambamba view -f bam -t 5 --subsampling-seed=3 -s `sed '5q;d' {input.flagstat} | cut -d" " -f1 | awk '{{ratio = {config[target_reads]}/$0}};{{if(ratio < 1 )print ratio; else print 1}}'` {input.bam} | samtools sort -m 2G -@ 5 -T {output.bam}.tmp > {output.bam} 2> {log}
        samtools index {output.bam}
        """
rule make_inputSubtract_bigwigs:
    input : 
        control = os.path.join(WORKDIR, "04aln_downsample/{control}-downsample.sorted.bam"), 
        case =  os.path.join(WORKDIR, "04aln_downsample/{case}-downsample.sorted.bam")
    output:  os.path.join(WORKDIR, "06bigwig_inputSubtract/{case}-subtract-{control}.bw")
    log: os.path.join(WORKDIR, "00log/{case}-vs-{control}inputSubtract.makebw")
    threads: 5
    conda:
        "envs/full-atac-main-env.yml"
    params: jobname = "{case}"
    message: "making input subtracted bigwig for {input}"
    shell:
        """
        bamCompare --bamfile1 {input.case} --bamfile2 {input.control} --normalizeUsing RPKM  --operation log2 --operation first --scaleFactorsMethod None --binSize 30 --smoothLength 300 -p 5  --extendReads 200 -o {output} 2> {log}
        """

rule make_bigwigs:
    input : os.path.join(WORKDIR, "04aln_downsample/{sample}-downsample.sorted.bam"), os.path.join(WORKDIR, "04aln_downsample/{sample}-downsample.sorted.bam.bai")
    output: os.path.join(WORKDIR, "07bigwig/{sample}.bw")
    log: os.path.join(WORKDIR, "00log/{sample}.makebw")
    threads: 5
    conda:
        "envs/full-atac-main-env.yml"
    params: jobname = "{sample}"
    message: "making bigwig for {input}"
    shell:
        """
        bamCoverage -b {input[0]} --normalizeUsing RPKM --binSize 30 --smoothLength 300 -p 5 --extendReads 200 -o {output} 2> {log}
        """

# Deeptools QC
rule computeMatrix_QC:
    input : get_big_wig_with_mark_or_tf 
    output : os.path.join(WORKDIR, "DPQC/{mark}.computeMatrix.gz")
    conda:
        "envs/full-atac-main-env.yml"
    params : TSS_BED
    shell:
        """
        computeMatrix reference-point -S {input} -R {TSS_BED} -a 3000 -b 3000 -out {output}
        """

# Deeptools QC
rule plotHeatmap:
    input :  os.path.join(WORKDIR, "DPQC/{mark}.computeMatrix.gz")
    output : os.path.join(WORKDIR, "DPQC/{mark}.plotHeatmap.png")
    conda:
        "envs/full-atac-main-env.yml"
    shell:
        """
        plotHeatmap -m {input} -out {output} --colorMap jet
        """

# # ChipSeq QCs plots from deeptools. Plotfingerprints are really usefull to see focal enrichment of your CHip-Seq enrichment
rule plotFingerPrint:
    input:
        bam = get_bams_per_sample, 
        bai = get_bam_index_per_sample
    output:
        plot = os.path.join(WORKDIR, "DPQC/{samp}.fingerprint.png"), 
        rawCounts = os.path.join(WORKDIR, "DPQC/{samp}.plotFingerprintOutRawCounts.txt"), 
        qualityMetrics = os.path.join(WORKDIR, "DPQC/{samp}.plotFingerprintOutQualityMetrics.txt")
    conda:
        "envs/full-atac-main-env.yml"
    params: 
        labels = get_all_marks_per_sample
    shell:
        """
        plotFingerprint -b {input.bam} --plotFile {output.plot} --labels {params.labels} --region chr1 --skipZeros --numberOfSamples 100000 --minMappingQuality 30 --plotTitle {wildcards.samp} --outRawCounts {output.rawCounts} --outQualityMetrics {output.qualityMetrics}
        """

rule get_UCSC_bigwig:
    input : os.path.join(WORKDIR, "06bigwig_inputSubtract/{case}-subtract-{control}.bw")
    output : os.path.join(WORKDIR, "UCSC_compatible_bigWig/{case}-subtract-{control}.bw")
    params : 
        wig1 = os.path.join(WORKDIR, "06bigwig_inputSubtract/{case}-subtract-{control}.temp.wig"), 
        wig2 = os.path.join(WORKDIR, "06bigwig_inputSubtract/{case}-subtract-{control}.temp2.wig")
    conda:
        "envs/full-atac-main-env.yml"
    shell:
        """
        bigWigToWig {input} {params.wig1}
        sed -r 's/^[0-9]|^X|^Y|^MT/chr&/g' {params.wig1} | LC_COLLATE=C sort -k1,1 -k2,2n > {params.wig2}
        wigToBigWig {params.wig2} ~/genome_size_UCSC_compatible_GRCh37.75.txt {output}
        rm {params.wig1} {params.wig2}
        """

# Peak calling using MACS
rule call_peaks_macs1:
    input: 
        control = os.path.join(WORKDIR, "04aln_downsample/{control}-downsample.sorted.bam"), 
        case = os.path.join(WORKDIR, "04aln_downsample/{case}-downsample.sorted.bam")
    output: os.path.join(WORKDIR, "08peak_macs1/{case}-vs-{control}-macs1_peaks.bed"), os.path.join(WORKDIR, "08peak_macs1/{case}-vs-{control}-macs1-nomodel_peaks.bed")
    log:
        macs1 = os.path.join(WORKDIR, "00log/{case}-vs-{control}-call-peaks_macs1.log"), 
        macs1_nomodel = os.path.join(WORKDIR, "00log/{case}-vs-{control}-call-peaks-macs1-nomodel.log")
    params:
        name1 = "{case}-vs-{control}-macs1", 
        name2 = "{case}-vs-{control}-macs1-nomodel", 
        jobname = "{case}", 
        outdir = os.path.join(WORKDIR, "08peak_macs1/")
    conda:
        "envs/macs.yml"
    message: "call_peaks macs14 {input}: {threads} threads"
    shell:
        """
        macs -t {input.case} \
            -c {input.control} --keep-dup all -f BAM -g {config[macs_g]} \
            --outdir {params.outdir} -n {params.name1} --single-profile -p {config[macs_pvalue]} &> {log.macs1}

        # nomodel for macs14, shift-size will be 100 bp (e.g. fragment length of 200bp)
        # can get fragment length from the phantompeakqual. Now set it to 200 bp for all.
        macs -t {input.case} \
            -c {input.control} --keep-dup all -f BAM -g {config[macs_g]} \
            --outdir {params.outdir} -n {params.name2} --nomodel -p {config[macs_pvalue]} &> {log.macs1_nomodel}
        """
        

# Peak calling using MACS 2

rule call_peaks_macs2:
    input: 
        control = os.path.join(WORKDIR, "04aln_downsample/{control}-downsample.sorted.bam"), 
        case = os.path.join(WORKDIR, "04aln_downsample/{case}-downsample.sorted.bam")
    output:
        bed = os.path.join(WORKDIR, "09peak_macs2/{case}-vs-{control}-macs2_peaks.xls"), 
        broad = os.path.join(WORKDIR, "09peak_macs2/{case}-vs-{control}-macs2_peaks.broadPeak")
    log: os.path.join(WORKDIR, "00log/{case}-vs-{control}-call-peaks_macs2.log")
    params:
        name = "{case}-vs-{control}-macs2", 
        jobname = "{case}", 
        outdir = os.path.join(WORKDIR, "09peak_macs2")
    conda:
        "envs/macs.yml"
    message: "call_peaks macs2 {input}: {threads} threads"
    shell:
        """
        ## for macs2, when nomodel is set, --extsize is default to 200bp, this is the same as 2 * shift-size in macs14.
        macs2 callpeak -t {input.case} \
            -c {input.control} --keep-dup all -f BAM -g {config[macs2_g]} \
            --outdir {params.outdir} -n {params.name} -p {config[macs2_pvalue]} --broad --broad-cutoff {config[macs2_pvalue_broad]} --nomodel &> {log}
        """

rule get_FRiP:
    input:
        peaks = os.path.join(WORKDIR, "09peak_macs2/{case}-vs-{control}-macs2_peaks.broadPeak"),
        bam = os.path.join(WORKDIR, "03aln/{case}.sorted.bam"), 
    output:
        os.path.join(WORKDIR, "DPQC/{case}-vs-{control}.FRiP.summary")
    params:
        saf = os.path.join(WORKDIR, "DPQC/{case}.saf"),
        outputName = os.path.join(WORKDIR, "DPQC/{case}-vs-{control}.FRiP")
    conda:
        "envs/full-atac-main-env.yml"
    shell:
        """
        awk 'BEGIN{{OFS="\t";print "GeneID", "Chr","Start","End","Strand"}}{{print $4,$1,$2,$3,$6}}' {input.peaks} > {params.saf}
        featureCounts -a {params.saf} -F SAF -o {params.outputName} {input.bam}
        """

# Cleaning broadGapped peak by adding "chr" on chr, sorting, setting score > 1000 to 1000 with awk then converting to bigbed
rule get_UCSC_bigBed:
    input: os.path.join(WORKDIR, "09peak_macs2/{case}-vs-{control}-macs2_peaks.broadPeak")
    output: os.path.join(WORKDIR, "12UCSC_broad/{case}-vs-{control}-macs2_peaks.broadPeak")
    params : 
        bed1 = temp(os.path.join(WORKDIR, "12UCSC_broad/{case}-vs-{control}-macs2_peaks.bed"))
    conda:
        "envs/full-atac-main-env.yml"
    shell:
        """
        sed -r 's/^[0-9]|^X|^Y|^MT/chr&/g' {input} | LC_COLLATE=C sort -k1,1 -k2,2n | awk '{{if($5 > 1000) $5=1000}}; {{print $0}}' > {params.bed1}
        bedToBigBed {params.bed1} ~/genome_size_UCSC_compatible_GRCh37.75.txt {output} -type=bed6+3
        """

# Creating a Hub for UCSC
rule get_UCSC_hub:
    input:  
        bed = ALL_BROADPEAK, 
        bigwig = ALL_BIGWIGUCSC
    output:
        ALL_HUB
    params:
        output_dir = HUB_FOLDER, 
        sample_name = list(SAMPLES.keys()), 
        categories = MARKS_NO_CONTROL
    conda:
        "envs/full-atac-main-env.yml"
    shell:
        """
        python3 scripts/makeUCSCtrackHub.py --hub_name {PROJECT_NAME} --sample_name {params.sample_name} --categories {params.categories} --output_dir {params.output_dir} --peaks {input.bed} --bw {input.bigwig} 2> makehub.err
        """

# MultiQC: Takes flagstats, fastqc, phantompeakqualtools and the deeptools outputs

rule multiQC:
    input : ALL_MULTIQC_INPUT
    output: ALL_QC
    params: os.path.join(WORKDIR, "10multiQC/")
    conda:
        "envs/full-atac-main-env.yml"
    log: os.path.join(WORKDIR, "00log/multiqc.log")
    message: "multiqc for all logs"
    shell:
        """
        multiqc {input} -o {params} -f -v -n multiQC_log 2> {log}
        """

# ChromHMM section.
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
        conda:
            "envs/chromhmm.yml"
        shell:
            """
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
        conda:
            "envs/chromhmm.yml"
        shell:
            """
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
        conda:
            "envs/chromhmm.yml"
        shell:
            """
            unset DISPLAY && ChromHMM.sh -Xmx{params.memory} LearnModel -p 0 -b {config[binsize]} {params.input_folder} {params.output_folder} {config[state]} {config[chromHmm_g]} 2> {log}
            """
