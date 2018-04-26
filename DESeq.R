library("DESeq2")

args <- commandArgs(trailingOnly=TRUE)
print(args)
counts <- read.table('analysis_godot/ase_summary_refalt.tsv', header=TRUE)
samples <- read.table('DEseqDesign.txt', header=TRUE)

maternal <- DESeqDataSetFromMatrix(countData=counts[ , samples$Sample], colData=samples, design=~Mother + IsMaternal)
maternal <- DESeq(maternal)
res <- results(maternal, contrast=c('IsMaternal', TRUE, FALSE))
res$gene <- rownames(res)
res <- res[, c(7,1,2,3,4,5,6)]
write.table(res, 'analysis/maternal_deseq.tsv', sep='\t', quote=FALSE, row.names=FALSE)

