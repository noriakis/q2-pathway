from os import path
import pandas as pd
import qiime2
import q2templates
import numpy as np
import os
import subprocess
from q2_types.feature_data import DNAFASTAFormat

## Option to use q2-gcn-norm
def infer(sequences: DNAFASTAFormat, reference_sequences: DNAFASTAFormat,
        cn_table: pd.DataFrame, cn_16s_table: pd.DataFrame, threads: int = 1):
    cmd = ["vsearch", "--usearch_global"]