library(piphir)
library(dplyr)


argv <- commandArgs(trailingOnly = TRUE)
outputDir <- argv[1]
cnTablePath <- argv[2]
cn16STablePath <- argv[3]
full <- argv[4]

if (full=="False") {
	full <- FALSE
} else {
	full <- TRUE
}

cat("Reading in KO copy number\n")
kodf <- read.table(cnTablePath, row.names=1, header=1)
cat("Reading in 16S rRNA gene copy number\n")
cn <- read.table(cn16STablePath, row.names=1, sep="\t", header=FALSE)
cn[,1] <- as.numeric(cn[,1])
cat("Reading in BLAST results\n")
blast <- read.table(paste0(outputDir, "/blast_out.txt"), sep="\t", header=FALSE)
blast$V2 <- blast$V2 %>% strsplit(":") %>% vapply("[", 1, FUN.VALUE="a")
cat("Example reference ID:", blast$V2[1], "\n")
seqtab <- read.table(paste0(outputDir, "/seqtab.txt"), sep="\t", row.names=1, header=1)
result <- profileMetagenome(seqtab, cn, kodf, blast, full=full)

write.table(result, paste0(outputDir, "/ko_table.tsv"), sep="\t")

## Example
# ex <- loadExample()
# result <- profileMetagenome(ex$seqtab, ex$cn16s, ex$cnko, ex$blast, full=full)
