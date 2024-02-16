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

ex <- loadExample()
result <- profileMetagenome(ex$seqtab, ex$cn16s, ex$cnko, ex$blast, full=full)
write.table(result, paste0(outputDir, "/ko_table.tsv"), sep="\t")