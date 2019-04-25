import csv
import os
import json

shell.prefix("set -eo pipefail; echo BEGIN at $(date); ")
shell.suffix("; exitstat=$?; echo END at $(date); echo exit status was $exitstat; exit $exitstat")

configfile: "config_all.yaml"
localrules: all
# localrules will let the rule run locally rather than submitting to cluster
# computing nodes, this is for very small jobs

# load cluster config file
CLUSTER = json.load(open(config['CLUSTER_JSON']))
FILES = json.load(open(config['SAMPLES_JSON']))
TSS_BED = config['tss_bed']
WORKDIR = os.path.abspath(config["OUTPUT_DIR"])
PROJECT_NAME = config['PROJECT_NAME']
## list all samples by sample_name and sample_type (sampleName_MARK is the basis of most of this pipeline atm)
MARK_SAMPLES = []
SAMPLES = sorted(FILES.keys())

##Create sample_Marks list for all samples
for sample in SAMPLES:
    for sample_type in FILES[sample].keys():
        MARK_SAMPLES.append(sample + "_" + sample_type)

##Regroup Marks or TF by sample
SAMPLES = dict()
for sample in sorted(FILES.keys()):
    for sample_type in FILES[sample].keys():
        SAMPLES.setdefault(sample,[]).append(sample_type)

##Regroup samples per marks or TF
MARKS = dict()
for sample in sorted(FILES.keys()):
    for sample_type in FILES[sample].keys():
            MARKS.setdefault(sample_type,[]).append(sample)




##Aggregation of bigwigs by Marks or TF
def getBigWigWithMarkorTF(wildcards):
    samples = MARKS[wildcards.mark]
    bigwigs = list()
    for s in samples:
        bigwigs.append(os.path.join(WORKDIR,"07bigwig/"+s+"_"+wildcards.mark+".bw"))
    return bigwigs

##Aggregation of bams per sample
def getBamsperSample(wildcards):

    marks = SAMPLES[wildcards.samp]
    bams = list()
    for s in marks:
        bams.append(os.path.join(WORKDIR,"03aln/"+wildcards.samp+"_"+s+".sorted.bam"))
    return bams

##Aggregation of bam index per sample, necessary to avoid lanching without bai
def getBaisperSample(wildcards):
    marks = SAMPLES[wildcards.samp]
    bais = list()
    for s in marks:
        bais.append(os.path.join(WORKDIR,"03aln/"+wildcards.samp+"_"+s+".sorted.bam.bai"))
    return bais

# which sample_type is used as control for calling peaks: e.g. Input, IgG...
CONTROL = config["control"]
MARKS_NO_CONTROL = list(MARKS.keys())
MARKS_NO_CONTROL.remove(CONTROL)
CONTROLS = [sample for sample in MARK_SAMPLES if CONTROL in sample]
CASES = [sample for sample in MARK_SAMPLES if CONTROL not in sample]

## multiple samples may use the same control input/IgG files
CONTROLS_UNIQUE = list(set(CONTROLS))


## list BAM files
CONTROL_BAM = expand(os.path.join(WORKDIR,"03aln/{sample}.sorted.bam"), sample=CONTROLS_UNIQUE)
CASE_BAM = expand(os.path.join(WORKDIR,"03aln/{sample}.sorted.bam"), sample=CASES)

## peaks and bigwigs
ALL_PEAKS = []
ALL_inputSubtract_BIGWIG = []
ALL_BROADPEAK = []
ALL_BIGWIGUCSC = []
ALL_COMPUTEMATRIX = []
ALL_PLOTS = []

for case in CASES:
    sample = "_".join(case.split("_")[0:-1])
    control = sample + "_" + CONTROL
    if control in CONTROLS:
        ALL_PEAKS.append(os.path.join(WORKDIR,"08peak_macs1/{}_vs_{}_macs1_peaks.bed").format(case, control))
        ALL_PEAKS.append(os.path.join(WORKDIR,"08peak_macs1/{}_vs_{}_macs1_nomodel_peaks.bed").format(case, control))
        ALL_PEAKS.append(os.path.join(WORKDIR,"09peak_macs2/{}_vs_{}_macs2_peaks.xls").format(case, control))
        ALL_PEAKS.append(os.path.join(WORKDIR,"09peak_macs2/{}_vs_{}_macs2_peaks.broadPeak").format(case, control))
        ALL_inputSubtract_BIGWIG.append(os.path.join(WORKDIR,"06bigwig_inputSubtract/{}_subtract_{}.bw").format(case, control))
        ALL_BROADPEAK.append(os.path.join(WORKDIR,"12UCSC_broad/{}_vs_{}_macs2_peaks.broadPeak").format(case, control))
        ALL_BIGWIGUCSC.append(os.path.join(WORKDIR,"UCSC_compatible_bigWig/{}_subtract_{}.bw").format(case,control))


ALL_SAMPLES = CASES + CONTROLS_UNIQUE
ALL_BAM     = CONTROL_BAM + CASE_BAM
ALL_DOWNSAMPLE_BAM = expand(os.path.join(WORKDIR,"04aln_downsample/{sample}-downsample.sorted.bam"), sample = ALL_SAMPLES)
ALL_FASTQ   = expand(os.path.join(WORKDIR,"01seq/{sample}.fastq"), sample = ALL_SAMPLES)
ALL_FASTQC  = expand(os.path.join(WORKDIR,"02fqc/{sample}_fastqc.zip"), sample = ALL_SAMPLES)
ALL_INDEX = expand(os.path.join(WORKDIR,"03aln/{sample}.sorted.bam.bai"), sample = ALL_SAMPLES)
ALL_DOWNSAMPLE_INDEX = expand(os.path.join(WORKDIR,"04aln_downsample/{sample}-downsample.sorted.bam.bai"), sample = ALL_SAMPLES)
ALL_FLAGSTAT = expand(os.path.join(WORKDIR,"03aln/{sample}.sorted.bam.flagstat"), sample = ALL_SAMPLES)
ALL_PHANTOM = expand(os.path.join(WORKDIR,"05phantompeakqual/{sample}.spp.out"), sample = ALL_SAMPLES)
ALL_BIGWIG = expand(os.path.join(WORKDIR,"07bigwig/{sample}.bw"), sample = ALL_SAMPLES)
ALL_COMPUTEMATRIX = expand(os.path.join(WORKDIR,"DPQC/{mark}.computeMatrix.gz"), mark = MARKS)
ALL_PLOTS = expand(os.path.join(WORKDIR,"DPQC/{mark}.plotHeatmap.png"), mark = MARKS)
ALL_PLOTS.extend(expand(os.path.join(WORKDIR,"DPQC/{samp}.fingerprint.png"),samp = SAMPLES))
ALL_DPQC = expand(os.path.join(WORKDIR,"DPQC/{samp}.plotFingerprintOutRawCounts.txt"),samp = SAMPLES)
ALL_DPQC.extend(expand(os.path.join(WORKDIR,"DPQC/{samp}.plotFingerprintOutQualityMetrics.txt"),samp = SAMPLES))
ALL_QC = [os.path.join(WORKDIR,"10multiQC/multiQC_log.html")]



TARGETS = []
TARGETS.extend(ALL_FASTQC)
TARGETS.extend(ALL_BAM)
TARGETS.extend(ALL_DOWNSAMPLE_BAM)
TARGETS.extend(ALL_INDEX)
TARGETS.extend(ALL_PHANTOM)
TARGETS.extend(ALL_PEAKS)
TARGETS.extend(ALL_BIGWIG)
TARGETS.extend(ALL_inputSubtract_BIGWIG)
TARGETS.extend(ALL_FASTQ)
TARGETS.extend(ALL_FLAGSTAT)
TARGETS.extend(ALL_QC)
TARGETS.extend(ALL_BROADPEAK)
TARGETS.extend(ALL_BIGWIGUCSC)
TARGETS.extend(ALL_COMPUTEMATRIX)
TARGETS.extend(ALL_PLOTS)

if config["chromHMM"]:

    def get_chr(chromSize):
        with open(chromSize,'r') as fs:
            chr = [line.rstrip().split('\t')[0] for line in fs]
            return(chr)

    HISTONE_INCLUDED = config["histone_for_chromHMM"].split(" ")
    HISTONE_CASES = [sample for sample in CASES if sample.split("_")[-1] in HISTONE_INCLUDED ]
    HISTONE_SAMPLE = list(set([sample.split("_")[0] for sample in CASES if sample.split("_")[-1] in HISTONE_INCLUDED ]))
    ALL_BED = expand(os.path.join(WORKDIR,"bamtobed/{sample}.bed"), sample = HISTONE_CASES)
    CHROMHMM = expand(os.path.join(WORKDIR,"chromHMM/learn_{nb_state}_states/{sample}_{nb_state}_segments.bed"),sample = HISTONE_SAMPLE,nb_state=config["state"])
    CHRHMM = get_chr(config['chromHmm_g'])
    CHROMHMM_TABLE = [os.path.join(WORKDIR,"chromHMM/cellmarkfiletable.txt")]
    TARGETS.extend(ALL_BED)
    TARGETS.extend(CHROMHMM)
    TARGETS.extend(CHROMHMM_TABLE)

localrules: all
rule all:
    input: TARGETS


## get a list of fastq.gz files for the same mark, same sample
def get_fastq(wildcards):
    sample = "_".join(wildcards.sample.split("_")[0:-1])
    mark = wildcards.sample.split("_")[-1]
    return FILES[sample][mark]

## now only for single-end ChIPseq,
rule merge_fastqs:
    input: get_fastq
    output: os.path.join(WORKDIR,"01seq/{sample}.fastq")
    log: os.path.join(WORKDIR,"00log/{sample}_unzip")
    threads: CLUSTER["merge_fastqs"]["cpu"]
    params: jobname = "{sample}"
    message: "merging fastqs gunzip -c {input} > {output}"
    shell: "gunzip -c {input} > {output} 2> {log}"

rule fastqc:
    input:  os.path.join(WORKDIR,"01seq/{sample}.fastq")
    output: os.path.join(WORKDIR,"02fqc/{sample}_fastqc.zip"), os.path.join(WORKDIR,"02fqc/{sample}_fastqc.html")
    log:    os.path.join(WORKDIR,"00log/{sample}_fastqc")
    threads: CLUSTER["fastqc"]["cpu"]
    conda:
        "envs/multifastqc.yml"
    params :
        output_dir = os.path.join(WORKDIR,"02fqc")
    message: "fastqc {input}: {threads}"
    shell:
        """
        fastqc -o {params.output_dir} -f fastq --noextract {input} 2> {log}
        """

# get the duplicates marked sorted bam, remove unmapped reads by samtools view -F 4 and dupliated reads by samblaster -r
# samblaster should run before samtools sort
rule align:
    input:  os.path.join(WORKDIR,"01seq/{sample}.fastq")
    output: os.path.join(WORKDIR,"03aln/{sample}.sorted.bam"), os.path.join(WORKDIR,"00log/{sample}.align")
    threads: CLUSTER["align"]["cpu"]
    conda:
        "envs/alignment.yml"
    params:
        bowtie = "--chunkmbs 320 -m 1 --best -p 5 ",
        jobname = "{sample}"
    message: "aligning {input}: 16 threads"
    log:
        bowtie = os.path.join(WORKDIR,"00log/{sample}.align"),
        markdup = os.path.join(WORKDIR,"00log/{sample}.markdup")
    shell:
        """
        bowtie2 -p 4 -x {config[idx_bt1]} -q {input} 2> {log.bowtie} \
        | samblaster --removeDups \
    | samtools view -Sb -F 4 - \
    | samtools sort -m 8G -@ 4 -T {output[0]}.tmp -o {output[0]} 2> {log.markdup}
    """

rule index_bam:
    input:  os.path.join(WORKDIR,"03aln/{sample}.sorted.bam")
    output: os.path.join(WORKDIR,"03aln/{sample}.sorted.bam.bai")
    log:    os.path.join(WORKDIR,"00log/{sample}.index_bam")
    threads: 1
    conda:
        "envs/alignment.yml"
    params: jobname = "{sample}"
    message: "index_bam {input}: {threads} threads"
    shell:
        """
        samtools index {input} 2> {log}
        """

# check number of reads mapped by samtools flagstat, the output will be used for downsampling
rule flagstat_bam:
    input:  os.path.join(WORKDIR,"03aln/{sample}.sorted.bam")
    output: os.path.join(WORKDIR,"03aln/{sample}.sorted.bam.flagstat")
    log:    os.path.join(WORKDIR,"00log/{sample}.flagstat_bam")
    threads: 1
    conda:
        "envs/alignment.yml"
    params: jobname = "{sample}"
    message: "flagstat_bam {input}: {threads} threads"
    shell:
        """
        samtools flagstat {input} > {output} 2> {log}
        """
   
rule plotFingerPrint:
    input: 
        bam = getBamsperSample,
        bai = getBaisperSample
    output: 
        plot = os.path.join(WORKDIR,"DPQC/{samp}.fingerprint.png"),
        rawCounts = os.path.join(WORKDIR,"DPQC/{samp}.plotFingerprintOutRawCounts.txt"),
        qualityMetrics = os.path.join(WORKDIR,"DPQC/{samp}.plotFingerprintOutQualityMetrics.txt")
    conda:
        "envs/deeptools.yml"
    params: 
        labels = lambda wildcards : [wildcards.samp + "_" + marks for marks in SAMPLES[wildcards.samp]]
    shell:
        """
        plotFingerprint -b {input.bam} --plotFile {output.plot} --labels {params.labels} --region chr1 --skipZeros --numberOfSamples 100000 --minMappingQuality 30 --plotTitle {wildcards.samp} --outRawCounts {output.rawCounts} --outQualityMetrics {output.qualityMetrics}
        """

rule phantom_peak_qual:
    input: 
        bam = os.path.join(WORKDIR,"03aln/{sample}.sorted.bam"), 
        bai = os.path.join(WORKDIR,"03aln/{sample}.sorted.bam.bai")
    output: os.path.join(WORKDIR,"05phantompeakqual/{sample}.spp.out")
    log: os.path.join(WORKDIR,"00log/{sample}_phantompeakqual.log")
    threads: 4
    conda:
        "envs/spp.yml"
    params: os.path.join(WORKDIR,"05phantompeakqual/")
    message: "phantompeakqual for {input}"
    shell:
        """
        run_spp -c={input.bam} -savp -rf -p=4 -odir={params}  -out={output} -tmpdir={params} 2> {log}
        """

rule down_sample:
    input: 
        bam = os.path.join(WORKDIR,"03aln/{sample}.sorted.bam"), 
        bai = os.path.join(WORKDIR,"03aln/{sample}.sorted.bam.bai"), 
        flagstat = os.path.join(WORKDIR,"03aln/{sample}.sorted.bam.flagstat")
    output: 
        bam = os.path.join(WORKDIR,"04aln_downsample/{sample}-downsample.sorted.bam"), 
        bai = os.path.join(WORKDIR,"04aln_downsample/{sample}-downsample.sorted.bam.bai")
    log: os.path.join(WORKDIR,"00log/{sample}_downsample.log")
    threads: 5
    conda:
        "envs/alignment.yml"
    params: 
    log: os.path.join(WORKDIR,"00log/{sample}.phantompeakqual.log")
    message: "downsampling for {input}"
    shell:
        """
        sambamba view -f bam -t 5 --subsampling-seed=3 -s `sed '5q;d' {input.flagstat} | cut -d" " -f1 | awk '{{ratio = {config[target_reads]}/$0}};{{if(ratio < 1 )print ratio; else print 1}}'` {input.bam} | samtools sort -m 2G -@ 5 -T {output.bam}.tmp > {output.bam} 2> {log}
        samtools index {output.bam}
        """

rule make_inputSubtract_bigwigs:
    input : 
        control = os.path.join(WORKDIR,"04aln_downsample/{control}-downsample.sorted.bam"),
        case =  os.path.join(WORKDIR,"04aln_downsample/{case}-downsample.sorted.bam")
    output:  os.path.join(WORKDIR,"06bigwig_inputSubtract/{case}_subtract_{control}.bw")
    log: os.path.join(WORKDIR,"00log/{case}_{control}inputSubtract.makebw")
    threads: 5
    conda:
        "envs/deeptools.yml"
    params: jobname = "{case}"
    message: "making input subtracted bigwig for {input}"
    shell:
        """
        bamCompare --bamfile1 {input.case} --bamfile2 {input.control} --normalizeUsing RPKM  --operation log2 --operation first --scaleFactorsMethod None --binSize 30 --smoothLength 300 -p 5  --extendReads 200 -o {output} 2> {log}
        """

rule make_bigwigs:
    input : os.path.join(WORKDIR,"04aln_downsample/{sample}-downsample.sorted.bam"), os.path.join(WORKDIR,"04aln_downsample/{sample}-downsample.sorted.bam.bai")
    output: os.path.join(WORKDIR,"07bigwig/{sample}.bw")
    log: os.path.join(WORKDIR,"00log/{sample}.makebw")
    threads: 5
    conda:
        "envs/deeptools.yml"
    params: jobname = "{sample}"
    message: "making bigwig for {input}"
    shell:
        """
        bamCoverage -b {input[0]} --normalizeUsing RPKM --binSize 30 --smoothLength 300 -p 5 --extendReads 200 -o {output} 2> {log}
        """

rule computeMatrix_QC:
    input : getBigWigWithMarkorTF 
    output : os.path.join(WORKDIR,"DPQC/{mark}.computeMatrix.gz")
    conda:
        "envs/deeptools.yml"
    params : TSS_BED
    shell:
        """
        computeMatrix reference-point -S {input} -R {TSS_BED} -a 3000 -b 3000 -out {output}
        """

rule plotHeatmap:
    input :  os.path.join(WORKDIR,"DPQC/{mark}.computeMatrix.gz")
    output : os.path.join(WORKDIR,"DPQC/{mark}.plotHeatmap.png")
    conda:
        "envs/deeptools.yml"
    shell:
        """
        plotHeatmap -m {input} -out {output} --colorMap jet
        """

rule get_UCSC_bigwig:
    input : os.path.join(WORKDIR,"06bigwig_inputSubtract/{case}_subtract_{control}.bw")
    output : os.path.join(WORKDIR,"UCSC_compatible_bigWig/{case}_subtract_{control}.bw")
    params : 
        wig1 = os.path.join(WORKDIR,"06bigwig_inputSubtract/{case}_subtract_{control}_temp.wig"),
        wig2 = os.path.join(WORKDIR,"06bigwig_inputSubtract/{case}_subtract_{control}_temp2.wig")
    conda:
        "envs/ucsc-utilities.yml"
    shell:
        """
        bigWigToWig {input} {params.wig1}
        sed -r 's/^[0-9]|^X|^Y|^MT/chr&/g' {params.wig1} | LC_COLLATE=C sort -k1,1 -k2,2n > {params.wig2}
        wigToBigWig {params.wig2} ~/genome_size_UCSC_compatible_GRCh37.75.txt {output}
        """

rule call_peaks_macs1:
    input: 
        control = os.path.join(WORKDIR,"04aln_downsample/{control}-downsample.sorted.bam"),
        case = os.path.join(WORKDIR,"04aln_downsample/{case}-downsample.sorted.bam")
    output: os.path.join(WORKDIR,"08peak_macs1/{case}_vs_{control}_macs1_peaks.bed"), os.path.join(WORKDIR,"08peak_macs1/{case}_vs_{control}_macs1_nomodel_peaks.bed")
    log:
        macs1 = os.path.join(WORKDIR,"00log/{case}_vs_{control}_call_peaks_macs1.log"),
        macs1_nomodel = os.path.join(WORKDIR,"00log/{case}_vs_{control}_call_peaks_macs1_nomodel.log")
    params:
        name1 = "{case}_vs_{control}_macs1",
        name2 = "{case}_vs_{control}_macs1_nomodel",
        jobname = "{case}",
        outdir = os.path.join(WORKDIR,"08peak_macs1/")
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
        

rule call_peaks_macs2:
    input: 
        control = os.path.join(WORKDIR,"04aln_downsample/{control}-downsample.sorted.bam"),
        case = os.path.join(WORKDIR,"04aln_downsample/{case}-downsample.sorted.bam")
    output:
        bed = os.path.join(WORKDIR,"09peak_macs2/{case}_vs_{control}_macs2_peaks.xls"),
        broad = os.path.join(WORKDIR,"09peak_macs2/{case}_vs_{control}_macs2_peaks.broadPeak")
    log: os.path.join(WORKDIR,"00log/{case}_vs_{control}_call_peaks_macs2.log")
    params:
        name = "{case}_vs_{control}_macs2",
        jobname = "{case}",
        outdir = os.path.join(WORKDIR,"09peak_macs2")
    conda:
        "envs/macs.yml"
    message: "call_peaks macs2 {input}: {threads} threads"
    shell:
        """
        ## for macs2, when nomodel is set, --extsize is default to 200bp, this is the same as 2 * shift-size in macs14.
        macs2 callpeak -t {input.case} \
            -c {input.control} --keep-dup all -f BAM --bdg -g {config[macs2_g]} \
            --outdir {params.outdir} -n {params.name} -p {config[macs2_pvalue]} --broad --broad-cutoff {config[macs2_pvalue_broad]} --nomodel &> {log}
        """

#Cleaning broadGapped peak by adding "chr" on chr, sorting, setting score > 1000 to 1000 with awk then converting to bigbed
rule get_UCSC_bigBed:
    input: os.path.join(WORKDIR,"09peak_macs2/{case}_vs_{control}_macs2_peaks.broadPeak")
    output: os.path.join(WORKDIR,"12UCSC_broad/{case}_vs_{control}_macs2_peaks.broadPeak")
    params : 
        bed1 = temp(os.path.join(WORKDIR,"12UCSC_broad/{case}_vs_{control}_macs2_peaks.bed"))
    conda:
        "envs/ucsc-utilities.yml"
    shell:
        """
        sed -r 's/^[0-9]|^X|^Y|^MT/chr&/g' {input} | LC_COLLATE=C sort -k1,1 -k2,2n | awk '{{if($5 > 1000) $5 = 1000}}; {{print $0}}' > {params.bed1}
        bedToBigBed {params.bed1} ~/genome_size_UCSC_compatible_GRCh37.75.txt {output} -type=bed6+3
        """
        
rule get_UCSC_hub:
    input:  
        bed = ALL_BROADPEAK,
        bigwig = ALL_BIGWIGUCSC
    params:
        output_dir = os.path.join(WORKDIR,"HUB/"),
        sample_name = list(SAMPLES.keys()),
        categories = MARKS_NO_CONTROL
    conda:
        "envs/trackhub.yml"
    shell:
        """
        python3 scripts/makeUCSCtrackHub.py --hub_name {PROJECT_NAME} --sample_name {params.sample_name} --categories {params.categories} --output_dir {params.output_dir} --peaks {input.bed} --bw {input.bigwig} 2> makehub.err
        """


rule multiQC:
    input :
        expand(os.path.join(WORKDIR,"00log/{sample}.align"), sample = ALL_SAMPLES),
        ALL_FLAGSTAT,
        ALL_FASTQC,
        ALL_PHANTOM,
        ALL_DPQC
    output: os.path.join(WORKDIR,"10multiQC/multiQC_log.html")
    params: os.path.join(WORKDIR,"10multiQC/")
    conda:
        "envs/multifastqc.yml"
    log: os.path.join(WORKDIR,"00log/multiqc.log")
    message: "multiqc for all logs"
    shell:
        """
        multiqc {input} -o {params} -d -f -v -n multiQC_log 2> {log}
        """


if config["chromHMM"]:
    rule bam2bed:
        input :
            os.path.join(WORKDIR,"04aln_downsample/{sample}-downsample.sorted.bam")
        output:
            os.path.join(WORKDIR,"bamtobed/{sample}.bed")
        params: jobname = "{sample}"
        log: os.path.join(WORKDIR,"00log/{sample}_bam2bed.log")
        message: "converting bam to bed for {input}"
        conda:
            "envs/chromhmm.yml"
        shell:
            """
            bedtools bamtobed -i {input} > {output}
            """

    rule make_table:
        input : 
            expand(os.path.join(WORKDIR,"bamtobed/{sample}.bed"), sample = HISTONE_CASES)
        output : 
            os.path.join(WORKDIR,"chromHMM/cellmarkfiletable.txt")
        log: os.path.join(WORKDIR,"00log/make_table_chromHMM.log")
        message: "making a table for chromHMM"
        run:
            import os
            from os.path import join
            with open (output[0], "w") as f:
                for case in HISTONE_CASES:
                    sample = "_".join(case.split("_")[0:-1])
                    mark = case.split("_")[-1]
                    control = sample + "_" + CONTROL
                    case_bed = case + ".bed"
                    if os.path.exists(join(os.path.join(WORKDIR,"bamtobed"), case_bed)):
                        f.write(sample + "\t" +  mark + "\t" + case + ".bed" + "\t" + control + ".bed" + "\n")

    rule chromHMM_binarize:
        input :
            cellmarkfiletable = os.path.join(WORKDIR,"chromHMM/cellmarkfiletable.txt"),
            beds = expand(os.path.join(WORKDIR,"bamtobed/{sample}.bed"), sample = HISTONE_CASES)
        output:
            expand(os.path.join(WORKDIR,"chromHMM/binarizedData/{sample}_{chr}_binary.txt"), sample = HISTONE_SAMPLE, chr=CHRHMM)
        log:
            os.path.join(WORKDIR,"00log/chromhmm_bin.log")
        params:
            folder = os.path.join(WORKDIR,"chromHMM/binarizedData/"),
            bamtobed_folder = os.path.join(WORKDIR,"bamtobed/"),
            memory = "32G"
        conda:
            "envs/chromhmm.yml"
        shell:
            """
            ChromHMM.sh -Xmx{params.memory} BinarizeBed -b {config[binsize]} {config[chromHmm_g]} {params.bamtobed_folder} {input.cellmarkfiletable} {params.folder} 2> {log}
            """

    rule chromHMM_learn:
        input:
            expand(os.path.join(WORKDIR,"chromHMM/binarizedData/{sample}_{chr}_binary.txt"), sample = HISTONE_SAMPLE, chr=CHRHMM)
        output: 
            expand(os.path.join(WORKDIR,"chromHMM/learn_{nb_state}_states/{sample}_{nb_state}_segments.bed"),sample = HISTONE_SAMPLE,nb_state=config["state"])
        log: os.path.join(WORKDIR,"00log/chromhmm_learn.log")
        params:
            input_folder = os.path.join(WORKDIR,"chromHMM/binarizedData/"),
            output_folder = expand(os.path.join(WORKDIR,"chromHMM/learn_{nb_state}_states/"),nb_state=config["state"]),
            memory = "32G"
        conda:
            "envs/chromhmm.yml"
        shell:
            """
            unset DISPLAY && ChromHMM.sh -Xmx{params.memory} LearnModel -p 0 -b {config[binsize]} {params.input_folder} {params.output_folder} {config[state]} {config[chromHmm_g]} 2> {log}
            """