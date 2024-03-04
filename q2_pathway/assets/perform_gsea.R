library(fgsea)

argv <- commandArgs(trailingOnly = TRUE)
outputDir <- argv[1]
pref <- argv[2]
kon <- as.integer(argv[3])
bg <- argv[4]

cat("Beginning GSEA by fgsea\n")
cat(pref, "\n")


ranks <- read.table(paste0(outputDir, "/values_",pref,".tsv"), sep="\t", header=FALSE)

## Remove NA
ranks <- ranks[!is.na(ranks$V2),]

vals <- ranks$V2
names(vals) <- ranks$V1
vals <- sort(vals, decreasing=TRUE)


## Should be done in python
pathmap <- read.table(paste0(outputDir, "/pathway_map.tsv"), sep="\t", header=FALSE)


if (bg!="all") {
    bgt <- read.table(paste0(outputDir, "/", pref, "_all_KO.txt"), sep="\t", header=FALSE)
    pathmap <- pathmap[pathmap$V2 %in% bgt$V1, ]
}

cat(dim(pathmap)[1], "\n")

if (kon==1) {
    change <- read.table(paste0(outputDir, "/pathway_names.tsv"), sep="\t", header=FALSE)
    desc <- change[,2]
    names(desc) <- change[,1]
    pathmap$V1 <- desc[pathmap$V1]
    
    rev <- change[,1]
    names(rev) <- change[,2]
}

paths <- unique(pathmap$V1)

nl <- lapply(paths, function(p) {
	tmp <- subset(pathmap, pathmap$V1 == p)
	tmp$V2
})

names(nl) <- paths



res <- fgsea::fgsea(nl, vals)
ores <- data.frame(res)

## Changing list to char
ores[[8]] <- unlist(lapply(ores[[8]], function(x) paste0(x, collapse="/")))
ores <- ores[order(ores$padj), ]
if (kon==1) {
    ores[["description"]] <- ores[["pathway"]]
    ores[["pathway"]] <- rev[ores[["pathway"]]]
}

write.table(ores, paste0(outputDir, "/gsea_res_", pref, ".tsv"), sep="\t")


topPathwaysUp <- res[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- res[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

## output image
png(paste0(outputDir, "/gsea_res_", pref, ".png"),
	width=14, height=8, res=96, units="in")
plotGseaTable(nl[topPathways], vals, res, gseaParam=0.5)
dev.off()