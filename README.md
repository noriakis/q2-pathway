# q2-pathway

QIIME2 plugin for visualizing and analyzing pathway information based on gene abundances

## Installation

```shell
mamba install -c noriakisato q2-pathway
```

## q2-pathway

This plugin is used to analyze the functional analysis results from metagenomic dataset. Currently, `kegg` module is implemented for visualization of KEGG PATHWAY images colored by the statistics calculated from comparing the categorical variables in the metadata. Also, `gsea` module is implemented for performing GSEA using `fgsea`, based on the KEGG PATHWAY mapping. The users should perform with `--verbose` for inspecting the GSEA output (like the existing of the ties).

The example usage:

```shell
qiime pathway kegg \
    --i-ko-table ko_metagenome.qza \
    --m-metadata-file metadata.tsv \
    --p-pathway-id ko00240 \
    --o-visualization pathway_output
```

The example output:

<p align="center">
<img src="https://github.com/noriakis/software/blob/main/images/q2_pathway_ex2.png?raw=true" width="800px">
</p>

The example usage:

```shell
qiime pathway gsea \
    --i-ko-table ko_metagenome.qza \
    --m-metadata-file metadata.tsv \
    --o-visualization gsea_output
```

The example output:

<p align="center">
<img src="https://github.com/noriakis/software/blob/main/images/q2_pathway_ex1.png?raw=true" width="800px">
</p>