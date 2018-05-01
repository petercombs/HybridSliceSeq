library("DESeq2")

args <- commandArgs(trailingOnly=TRUE)
print(args)
counts_orig <- read.table('analysis_godot/ase_summary_refalt.tsv', header=TRUE)
samples_orig <- read.table('DEseqDesign.txt', header=TRUE, stringsAsFactors=FALSE)
counts <- counts_orig[, samples_orig$Sample]

samples <- samples_orig[colSums(counts) > 1e5, ]
counts <- counts[ apply(counts, 1, max) >20, colSums(counts) > 1e5]

maternal <- DESeqDataSetFromMatrix(countData=counts[ , samples$Sample], colData=samples, design=~Mother + IsMaternal)
maternal <- DESeq(maternal)
res <- results(maternal, contrast=c('IsMaternal', TRUE, FALSE))
res$gene <- rownames(res)
res <- res[, c(7,1,2,3,4,5,6)]
write.table(res, 'analysis/maternal_deseq.tsv', sep='\t', quote=FALSE, row.names=FALSE)


samples$at <- factor(samples$AlignsTo)
melsim <- DESeqDataSetFromMatrix(countData=counts[ , samples$Sample], colData=samples, design=~Mother + at)
melsim <- DESeq(melsim)
res_ms<- results(melsim)
res_ms$gene <- rownames(res_ms)
res_ms <- res[,c(7,1,2,3,4,5,6)]
write.table(res_ms, 'analysis/species_deseq.tsv', sep='\t', quote=FALSE, row.names=FALSE)
