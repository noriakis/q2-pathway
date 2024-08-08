
## Run DESeq2 and return statistics
## The users must install DESeq2 beforehand
## BiocManager::install("DESeq2")

library(DESeq2)

argv <- commandArgs(trailingOnly = TRUE)

abundance <- argv[1]
meta <- argv[2]
condition <- argv[3]
outPath <- argv[4]
outImagePath <- argv[5]
lvl1 <- argv[6]
lvl2 <- argv[7]
outputDir <- argv[8]


# abund <- read.table(abundance, sep="\t", row.names=1, header=1)
# meta <- read.table(meta, sep="\t", row.names=1, header=1)
# abund_col <- colnames(abund)

# abund <- apply(abund, 1, function(x) as.integer(x))
# row.names(abund) <- abund_col

# inc_samples <- intersect(row.names(meta), colnames(abund))
# meta <- data.frame(meta[inc_samples, condition])
# colnames(meta) <- condition


# ## Controlling of the other parameters
# meta[["condition"]] <- factor(meta[[condition]])

# dds <- DESeqDataSetFromMatrix(countData = abund[, inc_samples],
#                               colData = meta,
#                               design = ~ condition)

# res <- DESeq(dds)

load(paste0(outputDir, "/dds_", condition, ".rda"))
x <- as.data.frame(results(res, contrast=c("condition", lvl1, lvl2)))
write.table(x, file=outPath, sep="\t", quote=FALSE)

## Image
png(outImagePath, width=5, height=5, res=96, units="in")
plotMA(res)
dev.off()

sess <- sessionInfo()
save(sess, file=paste0(outputDir, "/session_deseq2.rda"))