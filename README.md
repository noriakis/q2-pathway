# q2-pathway

QIIME2 plugin for visualizing and analyzing pathway information based on gene abundances

## Installation

After activating the QIIME2 environment, ALDEx2 shuold be installed if you are to use the statistics from the package.

```r
install.packages("BiocManager")
BiocManager::install("ALDEx2")
```

```shell
mamba install -c noriakisato q2-pathway
```

## q2-pathway

This plugin is used to analyze the functional analysis results from metagenomic dataset. 

`infer` module can perform an inferrence based on Piphillin or Tax4Fun2. For Tax4Fun2, the users should install the R package in the QIIME2 environment following [this tutorial](https://github.com/songweizhi/Tax4Fun2_short_tutorial), and download the reference database.
The database path should be set to `--p-reference-database`.

```shell
qiime pathway infer \
    --i-sequences rep-seqs.qza \
    --i-seq-table table-dada2.qza \
    --o-table infer_piphillin.qza \
    --p-threads 0
```

Also, `gsea` module is implemented for performing GSEA using `fgsea`, based on the KEGG PATHWAY mapping. The users should perform with `--verbose` for inspecting the GSEA output (like the existing of the ties). Althoug there is already a plugin (`q2-aldex2`), the function can rank the genes based on the statistics from ALDEx2. One should install ALDEx2 (`BiocManager::install("ALDEx2")`) beforehand.


```shell
qiime pathway gsea \
    --i-ko-table ko_metagenome.qza \
    --m-metadata-file metadata.tsv \
    --o-visualization gsea_output
```

The `summarize` module reports and compares the KO abundance table between the tables produced by multiple inference method. Using [`q2-sapienns`](https://github.com/gregcaporaso/q2-sapienns), the results from the shotgun metagenomics data can also be compared. The correlation metrics can be chosen from `spearman`, `pearson`, `kendall` by `--p-method`.

```shell
qiime pathway summarize \
    --i-ko-table ko_metagenome.qza \
    --m-metadata-file metadata.tsv \
    --o-visualization vis_output
```

`kegg` module is implemented for visualization of KEGG PATHWAY images colored by the statistics calculated from comparing the categorical variables in the metadata.

```shell
qiime pathway kegg \
    --i-ko-table ko_metagenome.qza \
    --m-metadata-file metadata.tsv \
    --p-pathway-id ko00240 \
    --o-visualization pathway_output
```


