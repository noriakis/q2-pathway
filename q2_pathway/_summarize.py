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
    tables: pd.DataFrame,
    metadata: qiime2.Metadata,
    first: int = 100,
    tss: bool = False,
    method: str = "pearson") -> None:

    tbl_len = len(tables)

    if tss:
        tables = [ko_table.apply(lambda x: x / sum(x), axis=1) for ko_table in tables]

    strat = False

    # all_cols = ko_table.columns.values
    filenames = []
    
    kos = [ko_table.columns.values for ko_table in tables]
    all_kos = list(set.intersection(*map(set, kos)))
    samples = [ko_table.index.values for ko_table in tables]
    all_samples = list(set.intersection(*map(set, samples)))

    ## Assuming unstratified output for ko_table
    ## Filter columns
    metadata = metadata.filter_ids(all_samples)
    metadata = metadata.filter_columns(column_type='categorical')
    metadata = metadata.filter_columns(
        drop_all_unique=True, drop_zero_variance=True, drop_all_missing=True)
    filenames = []

    ## save out metadata for download in viz
    metadata.save(os.path.join(output_dir, 'metadata.tsv'))

    ## Use in script
    metadata_df = metadata.to_dataframe()
    ## First table will be used for subset
    mean_val = tables[0].loc[:, all_kos].apply(lambda x: np.nanmean(x), axis=0).sort_values(ascending=False)
    all_cols = mean_val.head(first).index.values

    ## Loose checking (tax1_K00001)
    if "_" in all_cols[0]:
        strat = True
        ko_table = tables[0]
        
    ## As PICRUSt2 does not produce stratified table by default
    if strat:
        all_taxs = list(set([i.split("_")[0] for i in all_cols]))
        all_kos = list(set([i.split("_")[1] for i in all_cols]))
        
        
        for column in metadata_df.columns:
            metadata_df_filt = metadata_df[metadata_df[column].notna()]
            
            data = pd.concat([ko_table, metadata_df_filt],
                             axis=1, join='inner')
            
            for ko in all_kos:
                prefix = column+"_"+ko

                candidate_columns = [i for i in all_cols if ko in i]
                candidate_columns.append(column)
                output = data.loc[:, candidate_columns]
                output["sample"] = output.index.values
                output = pd.melt(output, id_vars=["sample", column])
                output["variable"] = output["variable"].apply(lambda x: x.split("_")[0])
                
                plt.figure()
                bp = sns.boxplot(data=output, x="variable", y="value")
                fig = bp.get_figure()
                fig.savefig(path.join(output_dir, prefix + ".png"))
                plt.close()
                    
                img = Image.open(path.join(output_dir, column+"_"+ko+".png"))
                img_byte_arr = io.BytesIO()
                img.save(img_byte_arr, format='PNG')
                img_byte_arr = base64.b64encode(img_byte_arr.getvalue()).decode("utf-8")

                csv_path = os.path.join(output_dir, prefix + ".csv")
                    
                output = output.groupby(column).apply(lambda x: x.groupby("variable").mean("value"))
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
                ## Per-KO stratified abundance per species
                
    else: ## If un-stratified
        all_kos = list(set(all_cols))
        for column in metadata_df.columns:
            metadata_df_filt = metadata_df[metadata_df[column].notna()]
            levels = metadata_df_filt[column].unique()
            concs = []
            for e, table in enumerate(tables):
                tmp = pd.concat([table, metadata_df_filt], axis=1, join='inner')
                tmp["category"] = "data"+str(e)
                concs.append(tmp)
                
            corrs = []
            for ko in all_kos:
                prefix = column+"_"+ko            
                conc = pd.concat([data.loc[:, [column, ko,"category"]] for data in concs], axis=0)
                    
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
                                
                output = pd.DataFrame(conc.loc[:,[column, ko, "category"]]).groupby("category").apply(lambda x: x.groupby(column).mean(ko))
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
                
                ## Correlation output per KO
                corrtbl = pd.concat([table.loc[all_samples, ko] for table in tables], axis=1)
                corrtbl.columns = ["data"+str(e) for e, i in enumerate(tables)]
                corr = corrtbl.corr(method=method)

                plt.figure()
                g = sns.heatmap(corr, annot=True)
                fig = g.get_figure()
                fig.savefig(path.join(output_dir, prefix + "_heatmap.png"))
                plt.close()
                
                img = Image.open(path.join(output_dir, column+"_"+ko+"_heatmap.png"))
                img_byte_arr = io.BytesIO()
                img.save(img_byte_arr, format='PNG')
                img_byte_arr = base64.b64encode(img_byte_arr.getvalue()).decode("utf-8")

                csv_path = os.path.join(output_dir, prefix + "_corr.csv")                                
                corr.to_csv(csv_path)
            
            
                
                ## Save the image
                prefix = prefix + "_corr"
                jsonp = prefix + ".jsonp"
                with open(os.path.join(output_dir, jsonp), 'w') as fh:
                    fh.write('load_data("%s",' % prefix)
                    json.dump(img_byte_arr, fh)
                    fh.write(",'")
                    table = q2templates.df_to_html(corr, escape=False)
                    fh.write(table.replace('\n', '').replace("'", "\\'"))
                    fh.write("','")
                    fh.write(quote(prefix+"_corr.csv"))
                    fh.write("');")
                filenames.append(jsonp)
                
                corr = corr.where(np.triu(np.ones(corr.shape)).astype(bool))
                corr = corr.stack().reset_index()
                corr.columns = ['d1','d2','value']
                corrs.append(corr)
        
        all_cor = pd.concat(corrs)
        corsum = all_cor.groupby("d1").apply(lambda x: x.groupby("d2").mean("value"))
        csv_path = os.path.join(output_dir, "whole_corr.csv")                                
        corsum.to_csv(csv_path)
        
        all_cor["label"] = all_cor.d1.map(str) + " - " + all_cor.d2
        plt.figure()
        g = sns.boxplot(all_cor, x="label", y="value")
        fig = g.get_figure()
        fig.savefig(path.join(output_dir, "whole_corr.png"))
        plt.close()
        
        img = Image.open(path.join(output_dir, "whole_corr.png"))
        img_byte_arr = io.BytesIO()
        img.save(img_byte_arr, format='PNG')
        img_byte_arr = base64.b64encode(img_byte_arr.getvalue()).decode("utf-8")

        jsonp = "whole.jsonp"
        with open(os.path.join(output_dir, jsonp), 'w') as fh:
            fh.write('load_data("%s",' % "whole")
            json.dump(img_byte_arr, fh)
            fh.write(",'")
            table = q2templates.df_to_html(corsum, escape=False)
            fh.write(table.replace('\n', '').replace("'", "\\'"))
            fh.write("','")
            fh.write(quote("whole_corr.csv"))
            fh.write("');")
        filenames.append(jsonp)
        
    index = os.path.join(TEMPLATES, 'index.html')
    q2templates.render(index, output_dir, context={
        'columns': [quote(fn) for fn in filenames]
    })
    shutil.copytree(
        os.path.join(TEMPLATES, 'dist'),
        os.path.join(output_dir, 'dist'))
