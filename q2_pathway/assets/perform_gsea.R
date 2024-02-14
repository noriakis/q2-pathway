
library(fgsea)
ranks <- read.table("ranks.txt", sep="\t")

res <- fgsea::fgsea(fgsea::examplePathways, fgsea::exampleRanks)
write.table(paste0("gsea_res_", condition, ".txt"), sep="\t")