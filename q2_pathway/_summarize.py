from os import path
import pandas as pd
import qiime2
import shutil
import q2templates
import numpy as np
import os
import json
import base64
import subprocess
from tempfile import TemporaryDirectory
import seaborn as sns
from PIL import Image
from urllib.parse import quote
import io
import pkg_resources
TEMPLATES = pkg_resources.resource_filename('q2_pathway', 'assets')

def summarize(output_dir: str, ko_table: pd.DataFrame,
    metadata: qiime2.Metadata) -> None:

    strat = False
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
    filenames = []
    
    ## Loose checking (tax1_K00001)
    if "_" in all_cols[0]:
        strat = True
    if strat:
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
    else:
        all_kos = list(set(all_cols))
        for column in metadata_df.columns:
            metadata_df_filt = metadata_df[metadata_df[column].notna()]
            
            data = pd.concat([ko_table, metadata_df_filt],
                             axis=1, join='inner')
            for ko in all_kos[1:10]:
                prefix = column+"_"+ko
                bp = sns.boxplot(data=data, x=column, y=ko)
                fig = bp.get_figure()
                fig.savefig(path.join(output_dir, prefix + ".png"))
                
                img = Image.open(path.join(output_dir, column+"_"+ko+".png"))
                img_byte_arr = io.BytesIO()
                img.save(img_byte_arr, format='PNG')
                img_byte_arr = base64.b64encode(img_byte_arr.getvalue()).decode("utf-8")

                csv_path = os.path.join(output_dir, prefix + ".csv")
                data.loc[:,ko].to_csv(csv_path)
            
                ## Save the image
                jsonp = prefix + ".jsonp"
                with open(os.path.join(output_dir, jsonp), 'w') as fh:
                    fh.write('load_data("%s",' % prefix)
                    json.dump(img_byte_arr, fh)
                    fh.write(",'")
                    table = q2templates.df_to_html(pd.DataFrame(data.loc[:,ko]), escape=False)
                    fh.write(table.replace('\n', '').replace("'", "\\'"))
                    fh.write("','")
                    fh.write(quote(prefix+".csv"))
                    fh.write("');")
                filenames.append(jsonp)
    index = os.path.join(TEMPLATES, 'index.html')
    q2templates.render(index, output_dir, context={
        'columns': [quote(fn) for fn in filenames]
    })
    shutil.copytree(
        os.path.join(TEMPLATES, 'dist'),
        os.path.join(output_dir, 'dist'))
