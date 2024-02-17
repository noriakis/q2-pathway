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
import matplotlib.pyplot as plt
import pkg_resources
TEMPLATES = pkg_resources.resource_filename('q2_pathway', 'assets')
ko_url = "https://www.genome.jp/dbget-bin/www_bget?"


def summarize(output_dir: str,
    ko_table: pd.DataFrame,
    metadata: qiime2.Metadata,
    first: int = 100,
    tss: bool = False,
    ko_table2: pd.DataFrame = None,
    method: str = "pearson") -> None:

    if tss:
        ko_table = ko_table.apply(lambda x: x / sum(x), axis=1)

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
    mean_val = ko_table.apply(lambda x: np.nanmean(x), axis=0).sort_values(ascending=False)
    all_cols = mean_val.head(first).index.values

    # all_cols = ko_table.columns.values
    filenames = []
    
        
    
    ## If table2 is given, all the intersection will be used
    if ko_table2 is not None:
        
        if tss:
            ko_table2 = ko_table2.apply(lambda x: x / sum(x), axis=1)

        all_cols = list(set(ko_table2.columns.values) & set(ko_table.columns.values))
        
        ## First table will be used for subset
        mean_val = ko_table.loc[:, all_cols].apply(lambda x: np.nanmean(x), axis=0).sort_values(ascending=False)
        all_cols = mean_val.head(first).index.values
        
        all_samples = list(set(ko_table2.index.values) & set(ko_table.index.values))

    
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
            levels = metadata_df_filt[column].unique()
            data = pd.concat([ko_table, metadata_df_filt],
                             axis=1, join='inner')
            if ko_table2 is not None:
                data2 = pd.concat([ko_table2, metadata_df_filt],
                    axis=1, join="inner")
            
            for ko in all_kos:
                prefix = column+"_"+ko
                
                ## Clear plot
                
                if ko_table2 is None:
                    plt.figure()
                    bp = sns.boxplot(data=data, x=column, y=ko)
                    fig = bp.get_figure()
                    fig.savefig(path.join(output_dir, prefix + ".png"))
                    plt.close()
                    
                    img = Image.open(path.join(output_dir, column+"_"+ko+".png"))
                    img_byte_arr = io.BytesIO()
                    img.save(img_byte_arr, format='PNG')
                    img_byte_arr = base64.b64encode(img_byte_arr.getvalue()).decode("utf-8")

                    csv_path = os.path.join(output_dir, prefix + ".csv")
                    
                    output = pd.DataFrame(data.loc[:,[column, ko]]).groupby(column).mean()
                    output.to_csv(csv_path)

                    output["ko"] = '<a href="'+ko_url+ko+'">'+ko+'</a>'
                else:
                    data["category"] = "ko1"             
                    data2["category"] = "ko2"
                    conc = pd.concat([data.loc[:, [column, ko,"category"]],
                        data2.loc[:,[column, ko, "category"]]], axis=0)
                    
                    plt.figure()
                    g = sns.FacetGrid(conc, col="category")
                    g.map(sns.boxplot, column, ko, order=levels)
                    g.savefig(path.join(output_dir, prefix + ".png"))
                    plt.close()
                    
                    img = Image.open(path.join(output_dir, column+"_"+ko+".png"))
                    img_byte_arr = io.BytesIO()
                    img.save(img_byte_arr, format='PNG')
                    img_byte_arr = base64.b64encode(img_byte_arr.getvalue()).decode("utf-8")

                    csv_path = os.path.join(output_dir, prefix + ".csv")
                    
                    corr = ko_table.loc[all_samples, ko].corr(ko_table2.loc[all_samples, ko], method=method)
                    
                    output = pd.DataFrame(conc.loc[:,[column, ko, "category"]]).groupby("category").apply(lambda x: x.groupby(column).mean(ko))
                    output["corr"] = corr
                    output.to_csv(csv_path)

                    output["ko"] = '<a href="'+ko_url+ko+'">'+ko+'</a>'
            
                ## Save the image
                jsonp = prefix + ".jsonp"
                with open(os.path.join(output_dir, jsonp), 'w') as fh:
                    fh.write('load_data("%s",' % prefix)
                    json.dump(img_byte_arr, fh)
                    fh.write(",'")
                    table = q2templates.df_to_html(output, escape=False)
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
