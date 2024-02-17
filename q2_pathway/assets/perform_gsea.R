library(fgsea)

# pathway_name_url = "https://rest.kegg.jp/list/pathway"
# pnm <- data.table::fread(pathway_name_url, header=FALSE)
# pnm$V1 <- gsub("map", "ko", pnm$V1)
# namec <- pnm$V2
# names(namec) <- pnm$V1

argv <- commandArgs(trailingOnly = TRUE)
outputDir <- argv[1]
pref <- argv[2]

cat("Beginning GSEA by fgsea\n")
cat(pref, "\n")

## Should be done in python
pathmap <- read.table(paste0(outputDir, "/pathway_map.tsv"), sep="\t", header=FALSE)
paths <- unique(pathmap$V1)

nl <- lapply(paths, function(p) {
	tmp <- subset(pathmap, pathmap$V1 == p)
	tmp$V2
})

names(nl) <- paths

ranks <- read.table(paste0(outputDir, "/values_",pref,".tsv"), sep="\t", header=FALSE)

## Remove NA
ranks <- ranks[!is.na(ranks$V2),]

vals <- ranks$V2
names(vals) <- ranks$V1
vals <- sort(vals, decreasing=TRUE)
res <- fgsea::fgsea(nl, vals)
ores <- data.frame(res)

## Changing list to char
ores[[8]] <- unlist(lapply(ores[[8]], function(x) paste0(x, collapse="/")))
ores <- ores[order(ores$padj), ]
# ores[["pathway"]] <- namec[ores[[1]]]

write.table(ores, paste0(outputDir, "/gsea_res_", pref, ".tsv"), sep="\t")


topPathwaysUp <- res[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- res[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

## output image
png(paste0(outputDir, "/gsea_res_", pref, ".png"),
	width=8, height=8, res=96, units="in")
plotGseaTable(nl[topPathways], vals, res, gseaParam=0.5)
dev.off()