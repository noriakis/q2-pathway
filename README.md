# q2-pathway

QIIME 2 plugin for analyzing gene family abundances and biological pathway information obtained from 16S rRNA gene sequencing data.


<p align="center">
<img src="https://github.com/noriakis/software/blob/main/images/q2-pathway-diagram.png?raw=true" width="500px">
</p>

## Installation

After activating the QIIME 2 environment:

```shell
mamba install -c noriakisato -c bioconda q2-pathway
```

`ALDEx2` and `DESeq2` shuold be installed in the QIIME 2 environment if you are to use the statistics from the packages:

```r
install.packages("BiocManager")
BiocManager::install("ALDEx2")
BiocManager::install("DESeq2")
```

If `kegg` module is to be used, `pykegg` must be installed:

```shell
pip install pykegg
```

## q2-pathway

This plugin is used to analyze the functional prediction results from 16S rRNA gene sequencing dataset and optionally the profile from shotgun metagenomes.

`infer` module can perform an inferrence based on Piphillin or Tax4Fun2 algorithm. For Tax4Fun2, the users should install the R package in the QIIME 2 environment following [this tutorial](https://github.com/songweizhi/Tax4Fun2_short_tutorial), and download the artifact containing reference database [here](URLTBU). For `Piphillin`, the prebuilt database is bundled with conda installation.

Suppose that `rep-seqs.qza` is an artifact containing the representative sequence (of ASV or OTU) and `table-dada2.qza` contains the count table of the sequences:

```shell
## This runs the Piphillin algorithm
qiime pathway infer \
    --i-sequences rep-seqs.qza \
    --i-seq-table table-dada2.qza \
    --o-table infer_piphillin.qza \
    --p-threads 0
```

```shell
## This runs the Tax4Fun2 algorithm
qiime pathway infer-t4f2 \
    --i-sequences rep-seqs.qza \
    --i-seq-table table-dada2.qza \
    --i-database Tax4Fun2_ReferenceData_v2_artifact.qza \
    --o-table infer_t4f2.qza
```

Also, `gsea` module is implemented for performing GSEA using `fgsea`, based on the KEGG PATHWAY mapping. The users should perform with `--verbose` for inspecting the GSEA output (like the existing of the ties). Althoug there is already a plugin (`q2-aldex2`), the function can rank the genes based on the statistics from ALDEx2. One should install ALDEx2 (`BiocManager::install("ALDEx2")`) beforehand.

```shell
## This will run GSEA based on the log2 fold changes calculated from DESeq2
qiime pathway gsea \
    --i-ko-table ko_metagenome.qza \
    --m-metadata-file metadata.tsv \
    --o-visualization gsea_output \
    --p-method deseq2
```

The `summarize` module reports and compares the gene family abundance table between the tables produced by multiple inference methods including [`q2-picrust2`](https://github.com/gavinmdouglas/q2-picrust2). Using [`q2-sapienns`](https://github.com/gregcaporaso/q2-sapienns), the results from the shotgun metagenomics data can also be compared. In that case, the converting table for shotgun accession to 16S accession (`Metadata` containing `converted` column) must be provided to `--p-convert-table`. The correlation metrics can be chosen from `spearman`, `pearson`, `kendall` by `--p-method`. Also, the correlation based on the p-values proposed in Sun et al. 2020. can be calculated by specifying `--p-use-p`.

```shell
qiime pathway summarize \
    --i-ko-table ko_metagenome.qza \
    --m-metadata-file metadata.tsv \
    --o-visualization vis_output
```

Finally, the `kegg` module is implemented for visualization of KEGG PATHWAY images colored by the statistics calculated from comparing the categorical variables in the metadata. The plugin uses `pykegg` to produce the images.

```shell
qiime pathway kegg \
    --i-ko-table ko_metagenome.qza \
    --m-metadata-file metadata.tsv \
    --p-pathway-id ko00240 \
    --o-visualization pathway_output
```

## Bugs and errors

If you find bugs, suggestions, or errors, please kindly report them to Issues, or make a pull request, or report it directly to [e-mail](nori@hgc.jp).