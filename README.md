# full-chipseq

A snakemake pipeline to process ChIP-seq data. This pipeline is based on the [pyflow-ChIP-seq](https://github.com/crazyhottommy/pyflow-ChIPseq) from the talented [Ming Tang](https://github.com/crazyhottommy). I repackaged most of the steps into conda environments and removed/added some steps.

This pipeline should be easy to install with snakemake and conda present on the system. At the moment it support running on a sge cluster.


#TODO
Work on the specificites of each histone marks : Downsmalping requicrement from ENCODE
WOrk on the specific of peak calling. Narrow or broad for different kinds of histones. Conditional input is tough though
