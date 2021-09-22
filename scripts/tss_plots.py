#! /usr/bin/env python

__author__ = 'AF'

# This script is adapted from the the ATAqC https://github.com/kundajelab/ataqc

#Order of import are important to set up the backend for plotting
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib import mlab
import argparse
import os
import numpy as np
import pybedtools
import metaseq
import logging


#######################################################
parser = argparse.ArgumentParser()
parser.add_argument("--bam", help="Required. Bam file, sorted and indexed")
parser.add_argument("--tss", help = "Required. Bed files with TSS coordinate")
parser.add_argument("--sampleName", help="Required. Name of your samples")
parser.add_argument("--chromSizes", help = "Required. genome size")
parser.add_argument("--fragLength",type=int, help= "Required. Average fragment length")
parser.add_argument("--bins",type=int, help="Size of the bins. Default 400", default="400")
parser.add_argument("--window",type=int, help="Distance from TSS", default="2000")
parser.add_argument("--threads",type=int, help="Thread numbers. Default 8", default="8")
parser.add_argument("--plotFolder", help="Required. Output Folder.")
parser.add_argument("--outputmultiQC_plot", help="Required. Output for multiQC plotting.")
parser.add_argument("--outputmultiQC_stat", help="Required. Output for multiQC stat.")



args = parser.parse_args()


bam_file = args.bam
tss = args.tss
prefix = args.sampleName
chromsizes = args.chromSizes
read_len = args.fragLength
plotFolder = args.plotFolder
outputmultiQC_plot = args.outputmultiQC_plot
outputmultiQC_stat = args.outputmultiQC_stat
bins = 400
bp_edge = 2000
processes = args.threads
greenleaf_norm=True


'''
Take bootstraps, generate tss plots, and get a mean and
standard deviation on the plot. Produces 2 plots. One is the
aggregation plot alone, while the other also shows the signal
at each TSS ordered by strength.
'''

logging.info('Generating tss plot...')
tss_plot_file = os.path.join(plotFolder,prefix + "_tss-enrich.pdf")
tss_plot_large_file = os.path.join(plotFolder,prefix + "_large_tss-enrich.pdf")

# Load the TSS file
tss = pybedtools.BedTool(tss)
tss_ext = tss.slop(b=bp_edge, g=chromsizes)
# Load the bam file
bam = metaseq.genomic_signal(bam_file, 'bam') # Need to shift reads and just get ends, just load bed file?
bam_array = bam.array(tss_ext, bins=bins, shift_width = -read_len/2, # Shift to center the read on the cut site
                    processes=processes, stranded=True)

if greenleaf_norm:
    # Use enough bins to cover 100 bp on either end
    num_edge_bins = int(100/(2*bp_edge/bins))
    bin_means = bam_array.mean(axis=0)
    avg_noise = (sum(bin_means[:num_edge_bins]) +
                sum(bin_means[-num_edge_bins:]))/(2*num_edge_bins)
    bam_array /= avg_noise
else:
    bam_array /= bam.mapped_read_count() / 1e6

# Generate a line plot
fig = plt.figure()
ax = fig.add_subplot(111)
x = np.linspace(-bp_edge, bp_edge, bins)

ax.plot(x, bam_array.mean(axis=0), color='r', label='Mean')
ax.axvline(0, linestyle=':', color='k')

# Note the middle high point (TSS)
tss_point_val = max(bam_array.mean(axis=0))

ax.set_xlabel('Distance from TSS (bp)')
ax.set_ylabel('Average read coverage (per million mapped reads)')
ax.legend(loc='best')

fig.savefig(tss_plot_file)

# Print a more complicated plot with lots of info

# Find a safe upper percentile - we can't use X if the Xth percentile is 0


upper_prct=99
if mlab.prctile(bam_array.ravel(), upper_prct) == 0.0:
    upper_prct = 100.0

plt.rcParams['font.size'] = 8
fig = metaseq.plotutils.imshow(bam_array,
                                x=x,
                                figsize=(5, 10),
                                vmin=5, vmax=upper_prct, percentile=True,
                                line_kwargs=dict(color='k', label='All'),
                                fill_kwargs=dict(color='k', alpha=0.3),
                                sort_by=bam_array.mean(axis=1))

# And save the file
fig.savefig(tss_plot_large_file)

f1 = open(outputmultiQC_plot, "w")
j = -2000
for i in bam_array.mean(axis=0):
    f1.write("%s\t%s\n" % (j, i) )
    j += 10

f1.close()


f2 = open(outputmultiQC_stat, "w")
f2.write("Sample Name\tTSS\n")
f2.write("%s\t%s\n" % (prefix, tss_point_val) )

f2.close()

