
## Make DESeq2 object and store it
## The users must install DESeq2 beforehand
## BiocManager::install("DESeq2")

library(DESeq2)

argv <- commandArgs(trailingOnly = TRUE)

abundance <- argv[1]
meta <- argv[2]
condition <- argv[3]
outputDir <- argv[4]


abund <- read.table(abundance, sep="\t", row.names=1, header=1)
meta <- read.table(meta, sep="\t", row.names=1, header=1)
abund_col <- colnames(abund)

abund <- apply(abund, 1, function(x) as.integer(x))
row.names(abund) <- abund_col

inc_samples <- intersect(row.names(meta), colnames(abund))
meta <- data.frame(meta[inc_samples, condition])
colnames(meta) <- condition

## Controlling of the other parameters
meta[["condition"]] <- factor(meta[[condition]])

dds <- DESeqDataSetFromMatrix(countData = abund[, inc_samples],
                              colData = meta,
                              design = ~ condition)

res <- DESeq(dds)
save(res, file=paste0(outputDir, "/dds_", condition, ".rda"))