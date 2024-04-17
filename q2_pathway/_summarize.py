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
import scipy
import matplotlib.pyplot as plt
from itertools import combinations
import pkg_resources
TEMPLATES = pkg_resources.resource_filename('q2_pathway', 'assets')
ko_url = "https://www.genome.jp/dbget-bin/www_bget?"


def summarize(output_dir: str,
    tables: pd.DataFrame,
    metadata: qiime2.Metadata,
    convert_table: qiime2.Metadata = None,
    first: int = 100,
    candidate_pathway: str = None,
    candidate: str = None,
    tss: bool = False,
    split_str: str = None,
    use_p: bool = False,
    method: str = "spearman",
    tables_name: str = None,
    cor_fig_width: int = 12,
    cor_thresh: float = None,
    map_ko: bool = False) -> None:


    if map_ko:
        ## Obtain KO description from KEGG REST API
        
        koname = pd.read_csv("https://rest.kegg.jp/list/ko", sep="\t", header=None, index_col=0)
        kodic = koname[1].to_dict()

    if convert_table is not None:
        ## Converting table should be artifact
        ## with the column name of "converted".
        ## So if shotgun profiled qza is to be read,
        ## the index should be shotgun-ID and "converted" column
        ## corresponds to 16S-ID.
        mapping = convert_table.to_dataframe()
        # mapping = pd.read_csv(convert, sep="\t", header=None, index_col=0)
        change = mapping["converted"].to_dict()

    tbl_len = len(tables)

    if tss:
        tables = [ko_table.apply(lambda x: x / sum(x), axis=1) for ko_table in tables]

    strat = False

    # all_cols = ko_table.columns.values
    filenames = []
    
    kos = [ko_table.columns.values for ko_table in tables]
    all_kos = list(set.intersection(*map(set, kos)))
    samples = [ko_table.index.values for ko_table in tables]

    if split_str is not None:
        samples = [[i.split(quote(split_str))[0] for i in j] for j in samples]

    if convert_table is not None:
        samples = [[change[i] if i in change.keys() else i for i in j] for j in samples]

    all_samples = list(set.intersection(*map(set, samples)))

    if len(all_samples) == 0:
        raise ValueError("No common samples present across datasets")

    if len(all_kos) == 0:
        raise ValueError("No common genes present across datasets")


    ## Assuming unstratified output for ko_table
    ## Filter columns
    metadata = metadata.filter_ids(all_samples)
    metadata = metadata.filter_columns(column_type='categorical')
    metadata = metadata.filter_columns(
        drop_all_unique=False, drop_zero_variance=False, drop_all_missing=True)

    filenames = []

    ## save out metadata for download in viz
    metadata.save(os.path.join(output_dir, 'metadata.tsv'))

    ## Use in script
    metadata_df = metadata.to_dataframe()
    
    if candidate is None:
        ## First table will be used for subset
        ## By default, the genes are sorted based on mean abundance and `p-first` KOs 
        ## are subset for summarization.
        mean_val = tables[0].loc[:, all_kos].apply(lambda x: np.nanmean(x), axis=0).sort_values(ascending=False)
        all_cols = mean_val.head(first).index.values
    else:
        all_cols = [candidate]
    
    ## Override `first` and `candidate` option if pathway is specified        
    if candidate_pathway is not None:
        ## Download pathway info
        kop = pd.read_csv("https://rest.kegg.jp/link/ko/pathway", sep="\t", header=None)
        kop = kop[kop[0].apply(lambda x: "path:ko" in x)]
        kop[0] = kop[0].apply(lambda x: x.split(":")[1])
        kop[1] = kop[1].apply(lambda x: x.split(":")[1])
        all_cols = kop[kop[0].isin([candidate_pathway])][1].tolist()
        all_cols = list(set(all_kos) & set(all_cols))

    ## Loose checking (assuming stratified table is depicted as tax1_K00001)
    if "_" in all_cols[0]:
        strat = True
        ko_table = tables[0]
        
    ## As q2-picrust2 does not produce stratified table by default
    ## What should be the format for ths stratified output?
    ## Currently, the correlation will not be produced for the stratified output.

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
                if len(candidate_columns) == 0:
                    continue
                candidate_columns.append(column)
                output = data.loc[:, candidate_columns]
                output["sample"] = output.index.values
                output = pd.melt(output, id_vars=["sample", column])
                output["variable"] = output["variable"].apply(lambda x: x.split("_")[0])
                
                plt.figure()
                bp = sns.boxplot(data=output, x="variable", y="value")
                fig = bp.get_figure()
                fig.savefig(path.join(output_dir, prefix + ".png"))
                plt.clf()
                plt.close(bp.fig)
                    
                img = Image.open(path.join(output_dir, column+"_"+ko+".png"))
                img_byte_arr = io.BytesIO()
                img.save(img_byte_arr, format='PNG')
                img_byte_arr = base64.b64encode(img_byte_arr.getvalue()).decode("utf-8")

                csv_path = os.path.join(output_dir, prefix + ".csv")
                    
                output = output.groupby(column).apply(lambda x: x.groupby("variable").mean("value"))
                
                if map_ko:
                    if ko in kodic.keys():
                        output["ko_description"] = kodic[ko]
                
                output.to_csv(csv_path)

                if ko.startswith("K"):
                    output["ko"] = '<a href="'+ko_url+ko+'">'+ko+'</a>'
                else:
                    output["ko"] = ko
                    
                
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
                ## Rename sample id based on the parameter
                if split_str is not None:
                    table.index = [i.split(quote(split_str))[0] for i in table.index.values]
                if convert_table is not None:
                    table.index = [change[i] if i in change.keys() else i for i in table.index.values]
                if not use_p:
                    tmp = pd.concat([table, metadata_df_filt], axis=1, join='inner')
                    if tables_name is not None:
                        tmp["category"] = tables_name[e]
                    else:
                        tmp["category"] = "data"+str(e)

                    concs.append(tmp)
                else:
                    ## We only need KO counts
                    concs.append(table)

            ## If Wilcoxon-based correlation, we ignore the 
            ## candidate, candidate_pathway, and first option
            if use_p:
                base = list(combinations(levels, 2))
                for pair in base:
                    level1 = pair[0]
                    level2 = pair[1]
                    prefix = level1 + " - " + level2
                    a = metadata_df_filt[metadata_df_filt[column] == level1].index.tolist()
                    b = metadata_df_filt[metadata_df_filt[column] == level2].index.tolist()
                    vals = []
                    for tbl in concs:
                        val = pd.Series(tbl.columns.map(
                            lambda x: scipy.stats.mannwhitneyu(tbl.loc[a, x], tbl.loc[b, x]).pvalue
                        ))
                        val.index = tbl.columns
                        vals.append(val)
                    kos = [val.index.values for val in vals]
                    all_kos = list(set.intersection(*map(set, kos)))
                    df = pd.DataFrame([val.loc[all_kos] for val in vals]).T
                    if tables_name is not None:
                        df.columns = tables_name
                    else:
                        df.columns = ["data"+str(e) for e, i in enumerate(tables)]
                    corr = df.corr()

                    plt.figure()
                    fig = sns.heatmap(corr, annot=True)
                    figs = fig.get_figure()
                    figs.savefig(path.join(output_dir, prefix + "_heatmap.png"))
                    plt.close()
                    
                    img = Image.open(path.join(output_dir, prefix+"_heatmap.png"))
                    img_byte_arr = io.BytesIO()
                    img.save(img_byte_arr, format='PNG')
                    img_byte_arr = base64.b64encode(img_byte_arr.getvalue()).decode("utf-8")

                    csv_path = os.path.join(output_dir, prefix + "_corr.csv")                                
                    corr.to_csv(csv_path)

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

            else: ## End of summarize for using the p-values
                corrs = []
                for ko in all_kos:
                    prefix = column+"_"+ko
                    if ko not in concs[0].columns.values:
                        continue
                    conc = pd.concat([data.loc[:, [column, ko, "category"]] for data in concs], axis=0)
                        
                    # plt.figure()
                    g = sns.FacetGrid(conc, col="category")
                    g.map(sns.boxplot, column, ko, order=levels)
                    g.savefig(path.join(output_dir, prefix + ".png"))
                    # plt.clf()
                    plt.close(g.fig)
                        
                    img = Image.open(path.join(output_dir, column+"_"+ko+".png"))
                    img_byte_arr = io.BytesIO()
                    img.save(img_byte_arr, format='PNG')
                    img_byte_arr = base64.b64encode(img_byte_arr.getvalue()).decode("utf-8")

                    csv_path = os.path.join(output_dir, prefix + ".csv")
                                    
                    output = pd.DataFrame(conc.loc[:,[column, ko, "category"]]).groupby("category").apply(lambda x: x.groupby(column).mean(ko))
                    if map_ko:
                        if ko in kodic.keys():
                            output["ko_description"] = kodic[ko]
                    
                    output.to_csv(csv_path)

                    if ko.startswith("K"):
                        output["ko"] = '<a href="'+ko_url+ko+'">'+ko+'</a>'
                    else:
                        output["ko"] = ko
                
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
                    
                    csv_path = os.path.join(output_dir, prefix + "_corr.csv")                                

                    ## Correlation output per KO per dataset
                    if tbl_len > 1:
                        ## If multiple tables, append scatterplot
                        corrtbl = pd.concat([table.loc[all_samples, ko] for table in tables], axis=1)
                        if tables_name is not None:
                            corrtbl.columns = tables_name
                        else:
                            corrtbl.columns = ["data"+str(e) for e, i in enumerate(tables)]

                        corr = corrtbl.corr(method=method)
                        base = list(combinations(corrtbl.columns.values, 2))

                        fig, ax =plt.subplots(1,1+len(base), figsize=(cor_fig_width, 4))
                        sns.heatmap(corr, annot=True, ax=ax[0])
                        for e, i in enumerate(base):                
                            sns.scatterplot(corrtbl, x=i[0], y=i[1], ax=ax[e+1])
                        plt.tight_layout()

                        figs = fig.get_figure()
                        figs.savefig(path.join(output_dir, prefix + "_heatmap.png"))
                        # plt.clf()
                        plt.close(fig)
                    else:
                        corrtbl = pd.concat([table.loc[all_samples, ko] for table in tables], axis=1)
                        if tables_name is not None:
                            corrtbl.columns = tables_name
                        else:
                            corrtbl.columns = ["data"+str(e) for e, i in enumerate(tables)]

                        corr = corrtbl.corr(method=method)
                        plt.figure()
                        fig = sns.heatmap(corr, annot=True)
                        figs = fig.get_figure()
                        figs.savefig(path.join(output_dir, prefix + "_heatmap.png"))
                        # plt.clf()
                        plt.close(fig.fig)
                    
                    corrtbl.to_csv(csv_path)

                    img = Image.open(path.join(output_dir, column+"_"+ko+"_heatmap.png"))
                    img_byte_arr = io.BytesIO()
                    img.save(img_byte_arr, format='PNG')
                    img_byte_arr = base64.b64encode(img_byte_arr.getvalue()).decode("utf-8")

                
                    
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
                        fh.write(quote(prefix+".csv"))
                        fh.write("');")
                    filenames.append(jsonp)
                    
                    corr = corr.where(np.triu(np.ones(corr.shape)).astype(bool))
                    corr = corr.stack().reset_index()
                    corr.columns = ['d1','d2','value']
                    corr["ko"] = ko
                    corrs.append(corr)

                ## Correlation summarization by boxplot for 
                ## every dataset pairs
                all_cor = pd.concat(corrs)
                csv_path = os.path.join(output_dir, "whole_corr.csv")                                
                all_cor.to_csv(csv_path)


                ## Correlation statistics
                corsum = all_cor.groupby("d1").apply(lambda x: x.groupby("d2").agg({'value': ['mean', 'median', 'min', 'max']}))

                all_cor["label"] = all_cor.d1.map(str) + " - " + all_cor.d2
                plt.figure(figsize=(12,10))
                g = sns.boxplot(all_cor, x="label", y="value")
                plt.xticks(rotation=45)
                plt.tight_layout()
                fig = g.get_figure()
                fig.savefig(path.join(output_dir, "whole_corr.png"))
                # plt.clf()
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
                ## End of summarize for the raw abundances

                if cor_thresh is not None:
                    ## Subset to between dataset correlation
                    subset_cor = all_cor[all_cor.d1 != all_cor.d2]
                    kostat = subset_cor.groupby("ko").agg({'value': ['mean', 'median', 'min', 'max']})
                    kostat["ko"] = kostat.index.values
                    kostat["mean"] = kostat[("value","mean")]

                    outp = kostat.loc[:,["ko","mean"]]
                    outp = outp[outp["mean"] > cor_thresh]
                    csv_path = os.path.join(output_dir, "cor_thresh_"+str(cor_thresh)+".csv")
                    outp.to_csv(csv_path)

                    plt.figure(figsize=(12,10))
                    g = sns.barplot(kostat, x="ko", y="mean")
                    plt.xticks(rotation=45)
                    plt.tight_layout()
                    fig = g.get_figure()
                    fig.savefig(path.join(output_dir, "corr_thresh.png"))
                    plt.close()
                    
                    img = Image.open(path.join(output_dir, "corr_thresh.png"))
                    img_byte_arr = io.BytesIO()
                    img.save(img_byte_arr, format='PNG')
                    img_byte_arr = base64.b64encode(img_byte_arr.getvalue()).decode("utf-8")

                    pref = "cor_thresh_"+str(cor_thresh)
                    jsonp = pref + ".jsonp"
                    with open(os.path.join(output_dir, jsonp), 'w') as fh:
                        fh.write('load_data("%s",' % pref)
                        json.dump(img_byte_arr, fh)
                        fh.write(",'")
                        table = q2templates.df_to_html(outp, escape=False)
                        fh.write(table.replace('\n', '').replace("'", "\\'"))
                        fh.write("','")
                        fh.write(quote(pref+".csv"))
                        fh.write("');")
                    filenames.append(jsonp)
        
    index = os.path.join(TEMPLATES, 'index.html')
    q2templates.render(index, output_dir, context={
        'columns': [quote(fn) for fn in filenames]
    })
    shutil.copytree(
        os.path.join(TEMPLATES, 'dist'),
        os.path.join(output_dir, 'dist'))
