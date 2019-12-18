# full-chipseq

A snakemake pipeline to process ChIP-seq data. This pipeline is based on the [pyflow-ChIP-seq](https://github.com/crazyhottommy/pyflow-ChIPseq) from the talented [Ming Tang](https://github.com/crazyhottommy). I repackaged most of the steps into conda environments and removed/added some steps.

This pipeline should be easy to install with snakemake and conda present on the system. At the moment it support running on a sge cluster.


#TODO
Create a input meta file that is processed into a json. Can also integrate the narrow mars into the json too
Load UCSC HUB to check for it
Verify that all files are commited and excluded
Check for macs peak calling about the discovery of peaks with different p-value# full-chipseq

**A snakemake pipeline to process single end human ChIP-seq data on a linux machine.**

This pipeline is based on the [pyflow-ChIP-seq](https://github.com/crazyhottommy/pyflow-ChIPseq) from the talented [Ming Tang](https://github.com/crazyhottommy). I repackaged most of the steps, added some, deleted some, created new conda environments and packaged all QC outputs in multiQC.

This pipeline should be easy to install with snakemake and conda present on the system. At the moment it support running on a sge cluster.



# Installation
On a linux terminal, please run the following commands :
```
git clone https://gitlab.univ-nantes.fr/foucal-a/full-chipseq.git
cd full-chipseq
./install_dependencies
```
Those commands will clone the git directory, move into the folder and install a few things :
* UCSC executables like *bedTbigBed*
* Conda environments
* Genomic data (for *homo sapiens* only at the moment)


##TEMP blacklist filtering

Thread on why it's still best practise to remove blacklist regions
https://www.biostars.org/p/184537/
https://www.biostars.org/p/361297/

https://github.com/Boyle-Lab/Blacklist/releases

# QC variable

A slew of variables are produced by the pipeline to help you assess the quality of your samples. They are summarized in the multiQC report.


## ENCODE metrics

They are summarized [here](https://www.encodeproject.org/data-standards/terms/)
The ENCODE consortium published guidelines in the following [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3431496/). This is a must read for anyone processing ChIP-Seq data. 

* **Non-Redundant Fraction (NRF)** : 

	$\frac{\text{Number of distinct uniquely mapping reads (i.e. after removing duplicates)}}{\text{Total number of reads}}$
	
* **PCR Bottlenecking Coefficient 1 (PBC1)**: 
     
	 $PBC1=\frac{M1}{MDISTINCT}$
     
     where:
     * M1: the number of genomic locations where exactly one read maps uniquely
     * MDISTINCT: the number of distinct genomic locations to which some read maps uniquely

* **PCR Bottlenecking Coefficient 1 (PBC2)**:

	$PBC2=\frac{M1}{M2}$ 

	where:
    * M1: number of genomic locations where only one read maps uniquely
    * M2: number of genomic locations where two reads map uniquely

Computed by an awk command on bamtobed
## phantompeakqualtools metrics
From [ENCODE](https://genome.ucsc.edu/ENCODE/qualityMetrics.html#chipSeq) and the [phantompeakqualtool](https://github.com/kundajelab/phantompeakqualtools) software
* **Normalized Strand Cross-correlation coefficient (NSC)**:
A measure of enrichment derived without dependence on prior determination of enriched regions. Forward and reverse strand read coverage signal tracks are computed (number of unique mapping read starts at each base in the genome on the + and - strand counted separately). The forward and reverse tracks are shifted towards and away from each other by incremental distances and for each shift, the Pearson correlation coefficient is computed. In this way, a cross-correlation profile is computed, representing the correlation between forward and reverse strand coverage at different shifts. The highest cross-correlation value is obtained at a strand shift equal to the predominant fragment length in the dataset as a result of clustering/enrichment of relative fixed-size fragments around the binding sites of the target factor or feature.
The NSC is the ratio of the maximal cross-correlation value (which occurs at strand shift equal to fragment length) divided by the background cross-correlation (minimum cross-correlation value over all possible strand shifts). Higher values indicate more enrichment, values less than 1.1 are relatively low NSC scores, and the minimum possible value is 1 (no enrichment). This score is sensitive to technical effects; for example, high-quality antibodies such as H3K4me3 and CTCF score well for all cell types and ENCODE production groups, and variation in enrichment in particular IPs is detected as stochastic variation. This score is also sensitive to biological effects; narrow marks score higher than broad marks (H3K4me3 vs H3K36me3, H3K27me3) for all cell types and ENCODE production groups, and features present in some individual cells, but not others, in a population are expected to have lower scores.

* **Relative Strand Cross-correlation coefficient (RSC)**:
A measure of enrichment derived without dependence on prior determination of enriched regions. Forward and reverse strand read coverage signal tracks are computed (number of unique mapping read starts at each base in the genome on the + and - strand counted separately). The forward and reverse tracks are shifted towards and away from each other by incremental distances and for each shift, the Pearson correlation coefficient is computed. In this way, a cross-correlation profile is computed representing the correlation values between forward and reverse strand coverage at different shifts. The highest cross-correlation value is obtained at a strand shift equal to the predominant fragment length in the dataset as a result of clustering/enrichment of relative fixed-size fragments around the binding sites of the target factor. For short-read datasets (< 100 bp reads) and large genomes with a significant number of non-uniquely mappable positions (e.g., human and mouse), a cross-correlation phantom-peak is also observed at a strand-shift equal to the read length. This read-length peak is an effect of the variable and dispersed mappability of positions across the genome. For a significantly enriched dataset, the fragment length cross-correlation peak (representing clustering of fragments around target sites) should be larger than the mappability-based read-length peak.
The RSC is the ratio of the fragment-length cross-correlation value minus the background cross-correlation value, divided by the phantom-peak cross-correlation value minus the background cross-correlation value. The minimum possible value is 0 (no signal), highly enriched experiments have values greater than 1, and values much less than 1 may indicate low quality.

* **Fraction of reads in peaks (FRiP)**: Computed using subCount Feature Count on downsampled bams and MACS peaks

## deepTools
[deepTools](https://deeptools.readthedocs.io/en/develop/index.html) is a suite of tools to process biological data. It provides usefull tools for QC.
* Transcription Start Site (TSS) Enrichment plots : 
Give you an indication of enrichment around TSS for your marks or TF. Produced by deeptools using TSS start sites
![Alt](H3K4me3.plotHeatmap.png)

* DeepTools fingerprint metrics : 
 Computed by deeptools [plotFingerprint](https://deeptools.readthedocs.io/en/develop/content/tools/plotFingerprint.html). This helps you evaluate the focal enrichment of your mark or TF.
 
#TODO
Create a input meta file that is processed into a json. Can also integrate the narrow mars into the json too
Load UCSC HUB to check for it
Verify that all files are commited and excluded
Check for macs peak calling about the discovery of peaks with different p-value
What can this pipeline do :



#Limitations/TODO
the MultiQC modules for deeptools plotProfile can't handle reference-point mode for computematrix.

