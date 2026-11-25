library(piphir)
library(dplyr)


argv <- commandArgs(trailingOnly = TRUE)
outputDir <- argv[1]
cnTablePath <- argv[2]
cn16STablePath <- argv[3]
full <- argv[4]
fullID <- argv[5]
fullKOID <- argv[6]

if (full=="False") {
	full <- FALSE
} else {
	full <- TRUE
	if (fullID=="None") {
		fullID <- NULL
	}
	if (fullKOID=="None") {
		fullKOID <- NULL
	}
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


dummySample <- NULL
if (ncol(seqtab) == 1) {
	dummySample <- "q2pathway_dummy_sample"
	while (dummySample %in% colnames(seqtab)) {
		dummySample <- paste0(dummySample, "_")
	}
	seqtab[[dummySample]] <- 0
}

if (full && !is.null(fullKOID)) {
	if (!(fullKOID %in% colnames(kodf))) {
		stop(paste0("Requested full_ko_id is absent from the KO copy-number table: ", fullKOID))
	}
	kodf <- kodf[, fullKOID, drop=FALSE]
	kodf[["q2pathway_dummy_KO"]] <- 0
}

if (full) {
	result <- profileMetagenome(
		seqtab, cn, kodf, blast,
		full = TRUE,
		fullOutput = "inmemory"
	)
} else {
	result <- profileMetagenome(seqtab, cn, kodf, blast, full = FALSE)
}
if (!is.null(dummySample) && dummySample %in% colnames(result)) {
	result <- result[, colnames(result) != dummySample, drop = FALSE]
}
if (!full) {
	write.table(result, paste0(outputDir, "/ko_table.tsv"), sep="\t")
} else {
	if (!is.null(fullID)) {
		result <- result[result$ID == fullID, , drop = FALSE]
	}
	if (!is.null(fullKOID)) {
		result <- result[result$KO == fullKOID, , drop = FALSE]
	}
	sample_cols <- setdiff(colnames(result), c("KO", "ID"))
	result_long <- tidyr::pivot_longer(
		result,
		cols = dplyr::all_of(sample_cols),
		names_to = "sample",
		values_to = "value"
	)
	write.table(
		result_long,
		paste0(outputDir, "/ko_table.txt"),
		sep = "\t",
		row.names = FALSE,
		quote = FALSE
	)
}

## Example
# ex <- loadExample()
# result <- profileMetagenome(ex$seqtab, ex$cn16s, ex$cnko, ex$blast, full=full)
	
