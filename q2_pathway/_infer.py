from os import path
import pandas as pd
import qiime2
import q2templates
import numpy as np
import os
import subprocess
from q2_types.feature_data import DNAFASTAFormat
import pkg_resources
TEMPLATES = pkg_resources.resource_filename('q2_pathway', 'assets')

## Option to use q2-gcn-norm
def infer(sequences: pd.Series, reference_sequences: str = path.join(TEMPLATES, "ref_seqs.fna.gz"),
        cn_table: str = path.join(TEMPLATES, "cn_ko.tsv.gz"),
        cn_16s_table: str = path.join(TEMPLATES, "cn_16s.tsv.gz"), threads: int = 1) -> pd.DataFrame:
    with TemporaryDirectory() as temp_dir:
        cmd = ["vsearch", "--usearch_global"]
        
        cmd = ["Rscript", path.join(TEMPLATES, "perform_piphillin.R"), temp_dir]
        ko = pd.read_csv(path.join(temp_dir, "ko_table.tsv"), sep="\t", index_col=0, header=0)
        return(ko)