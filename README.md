# full-chipseq


**A snakemake pipeline to process single end human ChIP-seq data on a linux machine.**

This pipeline will allow you to do the following:
* Quality control of your samples
* Alignment and peak calling
* Generate a hub to visualize your samples on the UCSC genome browser
* Fonctionnal annotation of your peaks
* Motif Finding TODO
* Super enhancer TODO


This pipeline is based on the [pyflow-ChIP-seq](https://github.com/crazyhottommy/pyflow-ChIPseq) from the talented [Ming Tang](https://github.com/crazyhottommy). I repackaged most of the steps, added some, deleted some, created new conda environments and packaged all QC outputs in multiQC.

This pipeline should be easy to install with snakemake and conda present on the system. At the moment it support running on a sge cluster.



# Installation
## Prerequisite
* Anaconda installed
* Snakemake installed
## Steps

On a linux terminal, please run the following commands :

```
git clone https://gitlab.univ-nantes.fr/foucal-a/full-chipseq.git
cd full-chipseq
./install_dependencies
```

Those commands will clone the git directory, move into the folder and install a few things :
* UCSC executables like *bedTbigBed*
* Conda environments
* Genomic data (for *homo sapiens* only at the moment, with hg19)


# Preparing to launch the pipeline

To launch the pipeline you will need :
* One json file that details where the samples can be found and what are the marks. An example can be found in "samples_from_meta.json"
This file can be created manually or you can use a meta file to be used with the sample2json script to create your json file this way :
python sample2json.py --meta meta.txt --fastq_dir directory

where directory contains all the fastq.gz needed

* One config files detailling all the different options and data needed to launch the pipe. One example is provided as config_example.yaml in the config directory



# Launching the pipeline


# Inspecting the output

# Delving into the pipeline : details of the rules

## Alignement and filtering steps


* <span style="color:green">rule </span> <span style="color:yellow"> merge_fastqs</span>: Merging of fastq files together into unzipped fastq. To save disk space, the fastq output is temporary.
* <span style="color:green">rule </span> <span style="color:yellow"> align</span>: Alignement using bowtie2, with bowtie index to be provided in the parameters. Then samblaster is used to mark duplicates, followed by removing of patches and alternate loci. To save disk space, the bam output is temporary. This steps has been added to have an unfiltered bam file to compute some QC metrics.
* <span style="color:green">rule </span> <span style="color:yellow"> filter_alignment</span>: Alignment filtering removing blacklist regions, reads below the mapping quality threshold specified in the parameter file and duplicate removal, followed by sorting. Blacklist regions are downloaded with the install_dependencies scripts from their [github repo](https://github.com/Boyle-Lab/). Thread on why it's still best practise to remove blacklist regions : 
https://www.biostars.org/p/184537/ and https://www.biostars.org/p/361297/
. Filtering is done with -F 1804 to get uniquely mapped reads, PCR duplicate removed and QC passed.

* <span style="color:green">rule </span> <span style="color:yellow"> create_bam_json</span>: simple rule to create a json file integrating the filtered bam. When starting the pipeline with the bam option = True and this json as a sample_json it avoids having to realign the fastqs.
*  <span style="color:green">rule </span> <span style="color:yellow"> index_bam</span>: Indexing filtered bams with samtools index
*  <span style="color:green">rule </span> <span style="color:yellow"> flagstat_bam</span>: Flagstat is used to get the read numbers for each samples
*   <span style="color:green">rule </span> <span style="color:yellow"> down_sample</span>: filtered bams are downsampled to the number of specified reads in the parameters files. You can put an arbitrary large number to avoid downsampling (i.e. 200000000). The rationale for downsampling is from the [ENCODE Roadmap](https://www.nature.com/articles/nature14248.pdf) paper. The blog from [Ming Tang](http://crazyhottommy.blogspot.com/2016/05/downsampling-for-bam-files-to-certain.html) also highlights some methods and excerpt from the paper. This is probably not usefull unless you have replicates since you prolly won't reach the recommanded 30 millions reads. The pipeline does not support merging replicates at the moment.
## QC steps
The metrics are described in more details in the QC part of the documentation
* <span style="color:green">rule </span> <span style="color:yellow"> encode_complexity</span>: This will compute the NRF, PCB1 and PBC2 using the beds transformed from the raw bams with bedtools bamtobed. It then uses an awk command from the [ENCODE_DCC ATAC SEQ](https://github.com/ENCODE-DCC/atac-seq-pipeline) pipeline to compute the metrics.
* <span style="color:green">rule </span> <span style="color:yellow"> fastqc</span>: fastqc tools on the merged fastqc files.
* <span style="color:green">rule </span> <span style="color:yellow"> phantom_peak_qual</span>: phantompeakqualtools is used to compute its metrics and the fragment length for MACS2 peak calling.
* <span style="color:green">rule </span> <span style="color:yellow"> computeMatrix_QC</span>: This comes from the [deepTools](https://deeptools.readthedocs.io/en/develop/content/list_of_tools.html) suite of tools. A necessary step to compute plotProfile and plotHeatmap. The mode used here is reference point, using the Transcription Start Site (TSS) from genes downloaded from the ensemble annotation. Youc an use anything you want here, just put in the the GENOME_TSS part of the config file.
* <span style="color:green">rule </span> <span style="color:yellow"> plotHeatmap</span>: plotHeatmap from deepTools, centered around the GENOME_TSS points, for each mark or TF. The profile on top of the heatmap has a scale problem at the moment. 
* <span style="color:green">rule </span> <span style="color:yellow"> plotProfile</span>: deepTools plotProfile. Plots are grouped by mark or TF, but you'll find everything together in the multiQC report.
* <span style="color:green">rule </span> <span style="color:yellow"> plotFingerPrint</span>: deepTools plotFingerPrint. This gives you all the metrics for the chromosome 1. There is one plot per sample with all its marks or TF, but you'll find everything together in the multiQC report.
* <span style="color:green">rule </span> <span style="color:yellow"> get_FRiP_for_multiqc</span>: Fraction of reads in peak using the featureCOunts module of subRead. This will use the appropriate peaks (narrow or broad) and the filtered bams
* <span style="color:green">rule </span> <span style="color:yellow"> get_broad_peak_counts_for_multiqc</span>: A simple extraction of broad peak counts with a python script to create custome data for MultiQC.
* <span style="color:green">rule </span> <span style="color:yellow"> get_narrow_peak_counts_for_multiqc</span>:  simple extraction of narrow peak counts with a python script to create custome data for MultiQC.

## Peak calling
Details of the overall peak calling strategy is int he Peak calling strategies

* <span style="color:green">rule </span> <span style="color:yellow"> call_narrow_peaks_macs1</span>: MACS1 narrow peak calling with shitsize from the phantompeaqualtools estimation, nomodel (no fragment length estimation) with the p-value from the config files macs_pvalue.
* <span style="color:green">rule </span> <span style="color:yellow"> call_narrow_peaks_macs2</span>: MACS2 narrow peak calling with extsize ( shiftsize in MACS1) from the phantompeaqualtools estimation, nomodel  (no fragment length estimation) and p-value from the config file macs2_pvalue.
* <span style="color:green">rule </span> <span style="color:yellow"> call_broad_peaks_macs2</span>: MACS2 broad peak calling with extsize ( shiftsize in MACS1) from the phantompeaqualtools estimation, nomodel  (no fragment length estimation) and p-value from the config file macs2_pvalue and the broad peak cutoff from the config file macs2_pvalue_broad_cutoff.

## Visualization
Those steps will ultimately create a HUB that can be used on the UCSC genome browser
* <span style="color:green">rule </span> <span style="color:yellow"> get_bigwigs_using_inputs</span>: Sample vs INput bigWig file using bamCompare. This will create a bigwig signal file in log2 fold change of sample vs input. RPKM is used as a scaling method. binSize is 10 and smoothLentgh is 30. SES scaling might be more appropriate.
* <span style="color:green">rule </span> <span style="color:yellow"> get_bigwigs</span>: Single sample bigWig signal file. RPKM is used as a scaling method. binSize is 10 and smoothLentgh is 30.
* <span style="color:green">rule </span> <span style="color:yellow"> get_bigbeds</span>: Uses MACS2 peak file to create bigbeds file with UCSC scripts bedToBigBed.
* <span style="color:green">rule </span> <span style="color:yellow"> get_UCSC_hub</span>: Create a HUB with all peaks, sample vs input bigwig. Uses a python package from [Ryan Dale](https://github.com/daler/trackhub). 


## MultiQC
* <span style="color:green">rule </span> <span style="color:yellow"> multiQC_config</span>: Copying the multiQC config file in the config output folder.
* <span style="color:green">rule </span> <span style="color:yellow"> multiQC</span>: Using phantompeakqualtools, deepTools, featureCOunts, custome peaks counting, ENCODE, fastqc and BOWTIE logs to create  the multiQC report. If BAM im=nput is true some rules wont be executed : FASTQC BOWTIE and ENCODE.

## Annotation
* <span style="color:green">rule </span> <span style="color:yellow"> homer_annotate</span>: Annotate the peak files with HOMER. It uses the annotation file in the config file genome_gtf and genome_fasta.

## ChromHMM
This will only be laucnhed if the ChromHMM flag in the config fil is set to True
* <span style="color:green">rule </span> <span style="color:yellow"> bam2bed</span>: Convert filtered bam into bed files using bedtools bamtobed.
* <span style="color:green">rule </span> <span style="color:yellow"> make_table</span>: Create the cellmarkfiletable.txt necessary for ChromHMM using some small python commands and the dictionaries of samples.
* <span style="color:green">rule </span> <span style="color:yellow"> chromHMM_binarize</span>: Binarize steps of ChromHMM. This will use the binsize in the config file, a fixed 32G of memory for Java and the list of chromosomes without patches and alternate loci.
* <span style="color:green">rule </span> <span style="color:yellow"> chromHMM_learn</span>: LearnModel steps of ChromHMM. This will use the binsize in the config file, a fixed 32G of memory for Java and the list of chromosomes without patches and alternate loci.


# Peak calling strategies

The tool of choice here is MACS. Parameters are controlled by the p-value and for broad calls, the broad cutoff threshold. MACS2 ouputs will give you excel ouputs that will allow to select on the q value if needed. For the moment, there is only individual calls

# QC

A slew of variables are produced by the pipeline to help you assess the quality of your samples. They are summarized in the multiQC report.


## ENCODE metrics

They are summarized [here](https://www.encodeproject.org/data-standards/terms/)
The ENCODE consortium published guidelines in the following [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3431496/). This is a must read for anyone processing ChIP-Seq data. Here I quickly recapitulate ENCODE metrics: 

* **Non-Redundant Fraction (NRF)** : 

	$`NRF=\frac{\text{Nunique}}{\text{Ntotal}}`$
    
    where:
     * Nunique: Number of distinct uniquely mapping reads (i.e. after removing duplicates)
     * Ntotal: Total number of reads
	
* **PCR Bottlenecking Coefficient 1 (PBC1)**: 
     
	 $`PBC1=\frac{M1}{MDISTINCT}`$
     
     where:
     * M1: the number of genomic locations where exactly one read maps uniquely
     * MDISTINCT: the number of distinct genomic locations to which some read maps uniquely

* **PCR Bottlenecking Coefficient 1 (PBC2)**:

	$`PBC2=\frac{M1}{M2}`$ 

	where:
    * M1: number of genomic locations where only one read maps uniquely
    * M2: number of genomic locations where two reads map uniquely

Those are all computed by an awk command on bam file transformed into a bed file by the UCSC utility bamtobed

## phantompeakqualtools metrics
From [ENCODE](https://genome.ucsc.edu/ENCODE/qualityMetrics.html#chipSeq) and the [phantompeakqualtool](https://github.com/kundajelab/phantompeakqualtools) software
* **Normalized Strand Cross-correlation coefficient (NSC)**:
A measure of enrichment derived without dependence on prior determination of enriched regions. Forward and reverse strand read coverage signal tracks are computed (number of unique mapping read starts at each base in the genome on the + and - strand counted separately). The forward and reverse tracks are shifted towards and away from each other by incremental distances and for each shift, the Pearson correlation coefficient is computed. In this way, a cross-correlation profile is computed, representing the correlation between forward and reverse strand coverage at different shifts. The highest cross-correlation value is obtained at a strand shift equal to the predominant fragment length in the dataset as a result of clustering/enrichment of relative fixed-size fragments around the binding sites of the target factor or feature.
The NSC is the ratio of the maximal cross-correlation value (which occurs at strand shift equal to fragment length) divided by the background cross-correlation (minimum cross-correlation value over all possible strand shifts). Higher values indicate more enrichment, values less than 1.1 are relatively low NSC scores, and the minimum possible value is 1 (no enrichment). This score is sensitive to technical effects; for example, high-quality antibodies such as H3K4me3 and CTCF score well for all cell types and ENCODE production groups, and variation in enrichment in particular IPs is detected as stochastic variation. This score is also sensitive to biological effects; narrow marks score higher than broad marks (H3K4me3 vs H3K36me3, H3K27me3) for all cell types and ENCODE production groups, and features present in some individual cells, but not others, in a population are expected to have lower scores.

* **Relative Strand Cross-correlation coefficient (RSC)**:
A measure of enrichment derived without dependence on prior determination of enriched regions. Forward and reverse strand read coverage signal tracks are computed (number of unique mapping read starts at each base in the genome on the + and - strand counted separately). The forward and reverse tracks are shifted towards and away from each other by incremental distances and for each shift, the Pearson correlation coefficient is computed. In this way, a cross-correlation profile is computed representing the correlation values between forward and reverse strand coverage at different shifts. The highest cross-correlation value is obtained at a strand shift equal to the predominant fragment length in the dataset as a result of clustering/enrichment of relative fixed-size fragments around the binding sites of the target factor. For short-read datasets (< 100 bp reads) and large genomes with a significant number of non-uniquely mappable positions (e.g., human and mouse), a cross-correlation phantom-peak is also observed at a strand-shift equal to the read length. This read-length peak is an effect of the variable and dispersed mappability of positions across the genome. For a significantly enriched dataset, the fragment length cross-correlation peak (representing clustering of fragments around target sites) should be larger than the mappability-based read-length peak.
The RSC is the ratio of the fragment-length cross-correlation value minus the background cross-correlation value, divided by the phantom-peak cross-correlation value minus the background cross-correlation value. The minimum possible value is 0 (no signal), highly enriched experiments have values greater than 1, and values much less than 1 may indicate low quality.


## deepTools
[deepTools](https://deeptools.readthedocs.io/en/develop/index.html) is a suite of tools to process biological data and is especially useful for Chip-Seq. Numerous plots and metrics are computed to give you an idea on the quality of your data.
* Transcription Start Site (TSS) Enrichment plots : 
Give you an indication of enrichment around TSS for your marks or TF. Produced by deeptools using TSS start sites that are specified in the config file. You can provide your own data from the table browsing from ucsc or use the data from the pipeline. The Profile plot is included in the multiQC report.
![Alt](http://deeptools.ie-freiburg.mpg.de/etc/galaxy/web/welcome_hm_GC.png)

* DeepTools fingerprint metrics : 
 This helps you evaluate the focal enrichment of your mark or TF. It's explained in more details in [plotFingerprint](https://deeptools.readthedocs.io/en/develop/content/tools/plotFingerprint.html). 
 
## Metric using peaks called by MACS
* **Fraction of reads in peaks (FRiP)** : Computed using subCount Feature Count on the filtered bams and the appropriate MACS peaks
* **Counts of broadpeaks from MACS2**
* **Counts of narrow peaks from MACS1**




