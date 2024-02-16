from os import path
import pandas as pd
import qiime2
import q2templates
import numpy as np
import os
import subprocess
from tempfile import TemporaryDirectory
import pkg_resources
TEMPLATES = pkg_resources.resource_filename('q2_pathway', 'assets')

## Option to use q2-gcn-norm for normalizing by rrnDB
def infer(sequences: pd.Series, reference_sequences: str = path.join(TEMPLATES, "ref_seqs.fna.gz"),
        cn_table: str = path.join(TEMPLATES, "cn_ko.tsv.gz"),
        cn_16s_table: str = path.join(TEMPLATES, "cn_16s.tsv.gz"), threads: int = 1, full: bool = False,
        pct_id: float = 0.99) -> pd.DataFrame:
    with TemporaryDirectory() as temp_dir:
    	repseq = path.join(temp_dir, "rep_seqs.fna")
        cmd = ["vsearch", "--usearch_global", repseq, "--db", reference_sequences, "--id", pct_id,
            '--top_hits_only',
            '--maxaccepts', '0',
            '--maxrejects', '0',
            '--uc_allhits',
            '--blast6out',
            path.join(temp_dir, "blast_out.txt"),
            "--threads", threads]
        # try:
        #     res = subprocess.run(cmd, check=True)
        # except subprocess.CalledProcessError as e:
        #     raise ValueError("Error running vsearch.")
                    
        cmd = ["Rscript", path.join(TEMPLATES, "perform_piphillin.R"), temp_dir]
        try:
            res = subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as e:
            raise ValueError("Error running piphillin algorithm.")
        ko = pd.read_csv(path.join(temp_dir, "ko_table.tsv"), sep="\t", index_col=0, header=0)
        return(ko.T)