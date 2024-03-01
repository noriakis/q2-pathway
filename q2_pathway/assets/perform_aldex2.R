
## Run ALDEx2 and return statistics
## The users must install ALDEx2 beforehand
## BiocManager::install("ALDEx2")

library(ALDEx2)

argv <- commandArgs(trailingOnly = TRUE)

abundance <- argv[1]
meta <- argv[2]
condition <- argv[3]
outPath <- argv[4]
outImagePath <- argv[5]
mcSamples <- as.integer(argv[6])

abund <- read.table(abundance, sep="\t", row.names=1, header=1)
meta <- read.table(meta, sep="\t", row.names=1, header=1)
abund_col <- colnames(abund)

abund <- apply(abund, 1, function(x) as.integer(x))
row.names(abund) <- abund_col

inc_samples <- intersect(row.names(meta), colnames(abund))
meta <- data.frame(meta[inc_samples, condition])
colnames(meta) <- condition

## Controlling of the other parameters

res <- aldex(abund[,inc_samples], as.character(meta[[condition]]), mcSamples)
x <- as.data.frame(res)
write.table(x, file=outPath, sep="\t", quote=FALSE)

## Image
png(outImagePath, width=10, height=5, res=96, units="in")
par(mfrow=c(1,2))
aldex.plot(res, type="MA", test="welch", xlab="Log-ratio abundance",
    ylab="Difference", main='Bland-Altman plot')
aldex.plot(res, type="MW", test="welch", xlab="Dispersion",
    ylab="Difference", main='Effect plot')
dev.off()