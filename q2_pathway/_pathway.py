import shutil
import base64
import io
import json
import pykegg
from os import path
import pandas as pd
import qiime2
import q2templates
import numpy as np
import os
from PIL import Image
import pkg_resources
import scipy
from urllib.parse import quote
from itertools import combinations
import subprocess

ko_url = "https://www.genome.jp/dbget-bin/www_bget?"
TEMPLATES = pkg_resources.resource_filename('q2_pathway', 'assets')


def kegg(output_dir: str, ko_table: pd.DataFrame, metadata: qiime2.Metadata, pathway_id: str,
    map_ko: bool = False, low_color: str = "blue", high_color: str = "red", tss: bool = False) -> None:
    
    if tss:
        ko_table = ko_table.apply(lambda x: x/sum(x), axis=1)    	
    
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

    ## Obtain KGML and nodes DataFrame
    graph = pykegg.KGML_graph(pid=pathway_id)
    nodes = graph.get_nodes()
    nodes = nodes[nodes.original_type == "ortholog"]

    ## KO description
    if map_ko:
        koname = pd.read_csv("https://rest.kegg.jp/list/ko", sep="\t", header=None)
        koname.index = koname[0].apply(lambda x: "ko:"+x)
        kodic = koname[1].to_dict()

    for column in metadata_df.columns:
        metadata_df_filt = metadata_df[metadata_df[column].notna()]
        uniq_level = metadata_df_filt[column].unique()
        base = list(combinations(uniq_level, 2))
        prefixes = []
        valdics = dict()
        for comb in base:
            ## Sort the level
            comb = sorted(comb)

            ## User control for ordering these values
            level1 = comb[0]
            level2 = comb[1]
            prefix = column+"_"+level1+"_vs_"+level2
            prefixes.append(prefix)

            ## T-stats
            a = metadata_df_filt[metadata_df_filt[column] == level1].index.tolist()
            b = metadata_df_filt[metadata_df_filt[column] == level2].index.tolist()
            val = pd.Series(ko_table.columns.map(lambda x: scipy.stats.ttest_ind(ko_table.loc[a, x], ko_table.loc[b, x], equal_var=False).statistic))
            val.index = ["ko:"+i for i in ko_table.columns]
            valdic = val.to_dict()
            valdics[prefix] = valdic

            nodes = pykegg.append_colors_continuous_values(nodes, valdic,
                    new_color_column=prefix, two_slope=False, node_name_column="name", delim=" ",
                    colors=["blue","red"], orig_value="tstat_"+prefix)

            output = nodes.loc[:,["name","tstat_"+prefix]]
            if map_ko:
                output["description"] = output["name"].apply(lambda x: " ".join([kodic[q] for q in x.split(" ")]))

            nodes[prefix] = nodes[prefix].fillna("#808080")
            kegg_map_image1 = pykegg.overlay_opencv_image(nodes, pid=pathway_id, fill_color=prefix,
                             transparent_colors=["#BFBFFF", "#FFFFFF", "#BFFFBF"])

            kegg_map_image2 = pykegg.append_legend(kegg_map_image1,
                                     min_value=np.nanmin([valdic[i] for i in valdic.keys()]),
                                     max_value=np.nanmax([valdic[i] for i in valdic.keys()]),
                                     two_slope=False, colors=["blue","red"],
                                     width=4, label="T-statistic")
            kegg_map_image = Image.fromarray(kegg_map_image2)

            img_byte_arr = io.BytesIO()
            kegg_map_image.save(img_byte_arr, format='PNG')
            img_byte_arr = base64.b64encode(img_byte_arr.getvalue()).decode("utf-8")

            csv_path = os.path.join(output_dir, prefix + ".csv")
            output.to_csv(csv_path)
            
            output["name"] = output["name"].apply(lambda x: " ".join(['<a href="'+ko_url+i+'">'+i+'</a>' for i in x.split(" ")]))
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

        ## Concat all the values
        minval = np.nanmin(nodes.loc[:, ["tstat_"+p for p in prefixes]].values)
        maxval = np.nanmax(nodes.loc[:, ["tstat_"+p for p in prefixes]].values)
        for p in prefixes:
            nodes = pykegg.append_colors_continuous_values(nodes, valdics[p],
                    new_color_column=p, two_slope=False, node_name_column="name", delim=" ",
                    colors=["blue","red"], orig_value="tstat_"+p, fix_min=minval, fix_max=maxval)
        
        output = nodes.loc[:,["name"]+["tstat_"+p for p in prefixes]]
        csv_path = os.path.join(output_dir, column + ".csv")
        output.to_csv(csv_path)

        if map_ko:
            output["description"] = output["name"].apply(lambda x: " ".join([kodic[q] for q in x.split(" ")]))

        output["name"] = output["name"].apply(lambda x: " ".join(['<a href="'+ko_url+i+'">'+i+'</a>' for i in x.split(" ")]))
        nodes["concat"] = nodes.loc[:,prefixes].values.tolist()
        nodes["concat"] = nodes["concat"].apply(lambda x: [i if isinstance(i, str) and i.startswith("#") else "#808080" for i in x])

        kegg_map_image1 = pykegg.overlay_opencv_image(nodes, pid=pathway_id, fill_color="concat",
                             transparent_colors=["#BFBFFF", "#FFFFFF", "#BFFFBF"])

        kegg_map_image2 = pykegg.append_legend(kegg_map_image1,
                                     min_value=minval,
                                     max_value=maxval,
                                     two_slope=False, colors=["blue","red"],
                                     width=4, label="T-statistic")
        kegg_map_image = Image.fromarray(kegg_map_image2)

        img_byte_arr = io.BytesIO()
        kegg_map_image.save(img_byte_arr, format='PNG')
        img_byte_arr = base64.b64encode(img_byte_arr.getvalue()).decode("utf-8")
             
        ## Save the whole condition image
        jsonp = column + ".jsonp"
        with open(os.path.join(output_dir, jsonp), 'w') as fh:
            fh.write('load_data("%s",' % column)
            json.dump(img_byte_arr, fh)
            fh.write(",'")
            table = q2templates.df_to_html(output, escape=False)
            fh.write(table.replace('\n', '').replace("'", "\\'"))
            fh.write("','")
            fh.write(quote(column+".csv"))
            fh.write("');")
        filenames.append(jsonp)

    ## Output
    index = os.path.join(TEMPLATES, 'index.html')
    q2templates.render(index, output_dir, context={
        'columns': [quote(fn) for fn in filenames]
    })

    shutil.copytree(
        os.path.join(TEMPLATES, 'dist'),
        os.path.join(output_dir, 'dist'))


def gsea(output_dir: str, ko_table: pd.DataFrame, metadata: qiime2.Metadata, tss: bool = False):
    
    if tss:
        ko_table = ko_table.apply(lambda x: x/sum(x), axis=1)

    ## Filter columns
    metadata = metadata.filter_ids(ko_table.index)
    metadata = metadata.filter_columns(column_type='categorical')
    metadata = metadata.filter_columns(
        drop_all_unique=True, drop_zero_variance=True, drop_all_missing=True)
    filenames = []
    
    ## Download pathway info for GSEA and output
    kop = pd.read_csv("https://rest.kegg.jp/link/ko/pathway", sep="\t", header=None)
    kop = kop[kop[0].apply(lambda x: "path:ko" in x)]
    kop.to_csv(os.path.join(output_dir, "pathway_map.tsv"), sep="\t", index=False, header=None)
    
    ## save out metadata for download in viz
    metadata.save(os.path.join(output_dir, 'metadata.tsv'))

    ## Use in script
    metadata_df = metadata.to_dataframe()

    for column in metadata_df.columns:
        metadata_df_filt = metadata_df[metadata_df[column].notna()]
        uniq_level = metadata_df_filt[column].unique()
        base = list(combinations(uniq_level, 2))
        prefixes = []
        valdics = dict()
        for comb in base:
            ## Sort the level
            comb = sorted(comb)

            ## User control for ordering these values
            level1 = comb[0]
            level2 = comb[1]
            prefix = column+"_"+level1+"_vs_"+level2
            prefixes.append(prefix)

            ## T-stats
            a = metadata_df_filt[metadata_df_filt[column] == level1].index.tolist()
            b = metadata_df_filt[metadata_df_filt[column] == level2].index.tolist()
            val = pd.Series(ko_table.columns.map(lambda x: scipy.stats.ttest_ind(ko_table.loc[a, x], ko_table.loc[b, x], equal_var=False).statistic))
            val.index = ["ko:"+i for i in ko_table.columns]
            val.to_csv(os.path.join(output_dir, "values_"+prefix+".tsv"), sep="\t", header=None)
            
            ## Perform GSEA
            cmd = ["Rscript", path.join(TEMPLATES, "perform_gsea.R"), output_dir, prefix]
            try:
                res = subprocess.run(cmd, check=True)
            except subprocess.CalledProcessError as e:
                raise ValueError("GSEA cannot be performed")
            
            ## Read output image
            gsea_image = Image.open(path.join(output_dir, "gsea_res_"+prefix+".png"))
            img_byte_arr = io.BytesIO()
            gsea_image.save(img_byte_arr, format='PNG')
            img_byte_arr = base64.b64encode(img_byte_arr.getvalue()).decode("utf-8")
            
            gseares = pd.read_csv(path.join(output_dir, "gsea_res_"+prefix+".tsv"), sep="\t", index_col=0, header=0)
            gseares["pathway"] = gseares["pathway"].apply(lambda x: '<a href="https://www.kegg.jp/entry/pathway+'+x.split(":")[1]+'">'+x+'</a>')
            jsonp = prefix + ".jsonp"
            
            ## The same structure as kegg
            with open(os.path.join(output_dir, jsonp), 'w') as fh:
                fh.write('load_data("%s",' % prefix)
                json.dump(img_byte_arr, fh)
                fh.write(",'")
                table = q2templates.df_to_html(gseares, escape=False)
                fh.write(table.replace('\n', '').replace("'", "\\'"))
                fh.write("','")
                fh.write(quote("gsea_res_"+prefix+".tsv"))
                fh.write("');")
            filenames.append(jsonp)
            
    index = os.path.join(TEMPLATES, 'index.html')
    q2templates.render(index, output_dir, context={
        'columns': [quote(fn) for fn in filenames]
    })
    shutil.copytree(
        os.path.join(TEMPLATES, 'dist'),
        os.path.join(output_dir, 'dist'))

                
