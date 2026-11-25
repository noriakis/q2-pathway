suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: summarize_ko.R <out_dir> <label>")
}
out_dir <- args[[1]]
label <- args[[2]]

per_sample <- fread(file.path(out_dir, "abundance_agreement_per_sample.tsv"))
per_feature <- fread(file.path(out_dir, "abundance_agreement_per_feature.tsv"))

per_sample[, comparison := paste(method1, method2, sep = " vs ")]
per_feature[, comparison := paste(method1, method2, sep = " vs ")]

p_sample <- ggplot(per_sample, aes(x = comparison, y = rho)) +
  geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.3) +
  geom_boxplot(outlier.size = 0.4) +
  coord_flip() +
  theme_bw(base_size = 10) +
  labs(
    x = NULL,
    y = "Spearman rho (across shared KOs)",
    title = sprintf("KO abundance agreement per sample [%s]", label)
  )

p_feature <- ggplot(per_feature, aes(x = comparison, y = rho)) +
  geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.3) +
  geom_boxplot(outlier.size = 0.4) +
  coord_flip() +
  theme_bw(base_size = 10) +
  labs(
    x = NULL,
    y = "Spearman rho (across shared samples)",
    title = sprintf("KO abundance agreement per KO [%s]", label)
  )

ggsave(file.path(out_dir, "KO_agreement_per_sample.pdf"), p_sample,
       width = 8, height = 5)
ggsave(file.path(out_dir, "KO_agreement_per_feature.pdf"), p_feature,
       width = 8, height = 5)
