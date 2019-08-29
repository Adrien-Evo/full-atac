# full-chipseq

A snakemake pipeline to process ChIP-seq data. This pipeline is based on the [pyflow-ChIP-seq](https://github.com/crazyhottommy/pyflow-ChIPseq) from the talented [Ming Tang](https://github.com/crazyhottommy). I repackaged most of the steps into conda environments and removed/added some steps.

This pipeline should be easy to install with snakemake and conda present on the system. At the moment it support running on a sge cluster.


#TODO
Add the QC options
Add the detection of identical input to save computing time and space

