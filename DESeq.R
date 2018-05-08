suppressPackageStartupMessages(library("DESeq2"))

args <- commandArgs(trailingOnly=TRUE)
print(args)
counts_orig <- read.table('analysis_godot/ase_summary_refalt.tsv', header=TRUE)
samples_orig <- read.table('DEseqDesign.txt', header=TRUE, stringsAsFactors=FALSE)
counts <- counts_orig[, samples_orig$FullSample]

sampleCounts <- colSums(counts)
refCounts <- sampleCounts[seq(2, length(sampleCounts), 2)]
altCounts <- sampleCounts[seq(1, length(sampleCounts), 2)]
sampleCounts[seq(1, length(sampleCounts), 2)] <- altCounts + refCounts
sampleCounts[seq(2, length(sampleCounts), 2)] <- altCounts + refCounts

samples <- samples_orig[sampleCounts > 1e5, ]
counts <- counts[ apply(counts, 1, max) >20, sampleCounts > 1e5]

sel = seq(1, nrow(counts), 1)
mel = grepl("melXsim", samples$FullSample)
sim = !mel

mel_mom <- DESeqDataSetFromMatrix(countData=counts[sel, samples$FullSample[mel]],
                                  colData=samples[mel, ],
                                  design=~Sample + AlignsTo)
mel_mom <- estimateSizeFactors(mel_mom)
sf <- sizeFactors(mel_mom)
alt_sf <- as.numeric(sf[seq(1, length(sf), 2)])
ref_sf <- as.numeric(sf[seq(2, length(sf), 2)])
sf[seq(1, length(sf), 2)] <- ref_sf + alt_sf
sf[seq(2, length(sf), 2)] <- ref_sf + alt_sf
sizeFactors(mel_mom) <- sf
mel_mom <- DESeq(mel_mom, parallel=TRUE)
res_mel <- results(mel_mom)
res_mel$gene <- rownames(res_mel)
res_mel <- res_mel[ , c(7,1,2,3,4,5,6)]
write.table(res_mel, 'analysis/mel_deseq.tsv', sep='\t', quote=FALSE, row.names=FALSE)

sim_mom <- DESeqDataSetFromMatrix(countData=counts[sel, samples$FullSample[sim]],
                                  colData=samples[sim, ],
                                  design=~Sample + AlignsTo)
sim_mom <- estimateSizeFactors(sim_mom)
sf <- sizeFactors(sim_mom)
alt_sf <- sf[seq(1, length(sf), 2)]
ref_sf <- sf[seq(2, length(sf), 2)]
sf[seq(1, length(sf), 2)] <- ref_sf + alt_sf
sf[seq(2, length(sf), 2)] <- ref_sf + alt_sf
sizeFactors(sim_mom) <- sf
sim_mom <- DESeq(sim_mom, parallel=TRUE)
res_sim <- results(sim_mom)
res_sim$gene <- rownames(res_sim)
res_sim <- res_sim[ , c(7,1,2,3,4,5,6)]
write.table(res_sim, 'analysis/sim_deseq.tsv', sep='\t', quote=FALSE, row.names=FALSE)

#maternal <- DESeqDataSetFromMatrix(countData=counts[sel, samples$FullSample], 
#                                   colData=samples, 
#                                   design=~ Sample + IsMaternal)
#
#sizeFactors(maternal) <- rep(1, nrow(samples))
#
#maternal <- DESeq(maternal, parallel=TRUE, test="LRT", reduced=~Sample)
#res <- results(maternal)
#res$gene <- rownames(res)
#res <- res[, c(7,1,2,3,4,5,6)]
#write.table(res, 'analysis/maternal_deseq.tsv', sep='\t', quote=FALSE, row.names=FALSE)
#
#
#samples$at <- factor(samples$AlignsTo)
#melsim <- DESeqDataSetFromMatrix(countData=counts[sel, samples$FullSample], colData=samples, design=~Mother + at)
#melsim <- DESeq(melsim, parallel=TRUE)
#res_ms<- results(melsim)
#res_ms$gene <- rownames(res_ms)
#res_ms <- res_ms[,c(7,1,2,3,4,5,6)]
#write.table(res_ms, 'analysis/species_deseq.tsv', sep='\t', quote=FALSE, row.names=FALSE)
