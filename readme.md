######################################################
# Spatially varying cis-regulatory divergence in     #
# Drosophila embryos elucidates cis-regulatory logic #
######################################################
_Peter A. Combs and Hunter B. Fraser_

_Department of Biology, Stanford University_

This repository contains the analysis code for Combs and Fraser 2017 ([bioRxiv
preprint](http://www.biorxiv.org/content/early/2017/08/10/175059); currently
submitted for review). Raw and processed data files available from the [Gene
Expression Omnibus](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE102233).

Very briefly, the data analyses in these scripts do two things:

1. Process RNA-seq data from cryosliced D. melanogaster x simulans hybrid
   embryos. We are looking for genes with spatially varying allele-specific
    expression (svASE) is different in one part of the embryo compared to the
    other.

2. Perform modeling of the cis-regulatory input functions of genes using a
   modeling approach inspired in large part by Ilsley, et al (2013). Using
   genome alignments and motif searches, we can then make inferences about which
   cis-regulatory changes actually produced the spatially varying ASE that we
   observed in part 1.

Almost all of the code is written in either Python 3 or
[Snakemake](http://snakemake.readthedocs.io/en/stable/).  Known dependencies
include:

* SRA Tools
* Snakemake
* Bowtie2 (Reference gDNA alignment)
* STAR (RNA-seq alignment)
* Cufflinks (RNA-seq quantification)
* Bedtools
* Samtools
* ASEr (https://github.com/thefraserlab/aser)
* Hornet (https://github.com/thefraserlab/hornet)
* Various python modules
	- pysam
	- progressbar
	- numpy/scipy/pandas/matplotlib
	- svgwrite
	- pyemd
	- BioPython
	- statsmodels


####
# RNA-seq and finding svASE
###

You should be able to go from raw reads to summary data tables by doing: 

    snakemake

Though you probably want to run this on a pretty high-powered machine---or,
ideally, a compute cluster. The basic steps for a single sample are:

1. Download reads from SRA
2. Map reads to the D. melanogaster genome (the resulting file is called
   `assigned_dmelR.bam` for various historical reasons)
3. Calculate absolute abundances using Cufflinks
4. Filter out potentially ambiguous/mismapped reads using Hornet to implement
   the WASP pipeline.
5. Count allele-specific reads for each gene using ASEr

Then, all of the expresion and ASE data are combined into summary tables.

We found that fitting either a logistic or gaussian function to the data found
all of the genes whose allele-specific expression had spatial patterns.
