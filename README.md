# q2-pathway

QIIME2 plugin for visualizing and analyzing pathway information based on gene abundances

## Installation

```shell
mamba install -c noriakisato q2-pathway
```

## q2-pathway

This plugin is used to analyze the functional analysis results from metagenomic dataset. Currently, `kegg` module is implemented for visualization of KEGG PATHWAY images colored by the statistics calculated from comparing the categorical variables in the metadata. Also, `gsea` module is implemented for performing GSEA using `fgsea`, based on the KEGG PATHWAY mapping. The users should perform with `--verbose` for inspecting the GSEA output (like the existing of the ties). The `summarize` module reports and compares the KO abundance table between the tables produced by multiple inference method. Using [`q2-sapienns`](https://github.com/gregcaporaso/q2-sapienns), the results from the shotgun metagenomics data can also be compared (WIP).


```shell
qiime pathway kegg \
    --i-ko-table ko_metagenome.qza \
    --m-metadata-file metadata.tsv \
    --p-pathway-id ko00240 \
    --o-visualization pathway_output
```

```shell
qiime pathway gsea \
    --i-ko-table ko_metagenome.qza \
    --m-metadata-file metadata.tsv \
    --o-visualization gsea_output
```

```shell
qiime pathway summarize \
    --i-ko-table ko_metagenome.qza \
    --m-metadata-file metadata.tsv \
    --o-visualization vis_output
```