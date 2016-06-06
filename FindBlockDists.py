from pysam import Samfile
from numpy import histogram, arange

dists = []
for read in Samfile('assigned_dmelR_wasp_dedup.bam'):
    if read.mpos - read.pos > 0:
        dists.append(read.mpos - read.reference_end)

print(histogram(dists, arange(-100, 500, 10)))


