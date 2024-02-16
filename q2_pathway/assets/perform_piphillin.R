library(piphir)
library(dplyr)


argv <- commandArgs(trailingOnly = TRUE)
outputDir <- argv[1]

ex <- loadExample()
result <- profileMetagenome(ex$seqtab, ex$cn16s, ex$cnko, ex$blast)
write.table(result, paste0(outputDir, "/ko_table.tsv"), sep="\t")