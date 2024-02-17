
## Run ALDEx2 and return statistics
## The users must install ALDEx2 beforehand
## BiocManager::install("ALDEx2")

library(ALDEx2)

argv <- commandArgs(trailingOnly = TRUE)

abundance <- argv[1]
meta <- argv[2]
condition <- argv[3]
outPath <- argv[4]

abund <- read.table(abundance, sep="\t", row.names=1, header=1)
meta <- read.table(meta, sep="\t", row.names=1, header=1)
abund_col <- colnames(abund)

abund <- apply(abund, 1, function(x) as.integer(x))
row.names(abund) <- abund_col

inc_samples <- intersect(row.names(meta), colnames(abund))
meta <- data.frame(meta[inc_samples, condition])
colnames(meta) <- condition

## Controlling of the other parameters
res <- aldex(abund[,inc_samples], as.character(meta[[condition]]))
res <- as.data.frame(res)
write.csv(res, file=outPath)