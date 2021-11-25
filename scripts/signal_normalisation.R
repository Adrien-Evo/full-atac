library(DESeq2)

args = commandArgs(trailingOnly=TRUE)

# read in the raw counts
rawcounts = read.table(args[1],h=T,stringsAsFactor=F, check.names = FALSE)

# Create a deseq2 object
dds <- DESeqDataSetFromMatrix(rawcounts,data.frame(row.names=seq_len(ncol(rawcounts))), ~1)

# estimate size factor
dds <-estimateSizeFactors(dds)

# store them in a simple named list
sizefactors <- 1/sizeFactors(dds)
#sizefactor are meant to be used as a divider by DESEQ2 (counts/scaling factor) while DeepTools bamCoverage scaling uses multiplication.
# It nneds to be reversed
# get output folder
output_folder <- dirname(args[1])
# write the scaling into a sample_name_scalingfactor.txt file
for (i in 1:length(sizefactors)){
    write.table(sizefactors[i],file.path(output_folder,paste0(names(sizefactors)[i],"_scalingfactor.txt")),quote=F, col.names = F, row.names = F)
}

