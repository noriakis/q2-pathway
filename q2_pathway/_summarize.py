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

def summarize(output_dir: str, ko_table: pd.DataFrame,
    metadata: qiime2.Metadata) -> None:
    ## Assuming stratified output for ko_table
    ## Filter columns
    metadata = metadata.filter_ids(ko_table.index)
    metadata = metadata.filter_columns(column_type='categorical')
    metadata = metadata.filter_columns(
        drop_all_unique=True, drop_zero_variance=True, drop_all_missing=True)
    filenames = []

    ## save out metadata for download in viz
    metadata.save(os.path.join(output_dir, 'metadata.tsv'))

    ## Use in script
    metadata_df = metadata.to_dataframe()
    all_cols = ko_table.columns.values
    
    ## Loose checking (tax1_K00001)
    if "_" not in all_cols[0]:
        raise ValueError("The input ko_table seems to be not stratified.")
        
    all_taxs = list(set([i.split("_")[0] for i in all_cols]))
    all_kos = list(set([i.split("_")[1] for i in all_cols]))
    
    for column in metadata_df.columns:
        metadata_df_filt = metadata_df[metadata_df[column].notna()]
        
        data = pd.concat([ko_table, metadata_df_filt],
                         axis=1, join='inner')
        kos = []
        names = []
        groups = []
        
        for ko in all_kos:
            candidate_columns = [i for i in all_cols if ko in i]
            tmp_sp = {}
            for name, group in data.groupby(column):
                kos.append(ko)
                names.append('%s (n=%d)' % (name, len(group)))
                for i in candidate_columns:
                    tmp_sp[i.split("_")[0]] = list(group[i])
                groups.append(tmp_sp)
        ## We have KO list, along with group length and 
        ## Per-KO stratified abundance per species
