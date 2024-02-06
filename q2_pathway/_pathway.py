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

TEMPLATES = pkg_resources.resource_filename('q2_pathway', 'assets')

def kegg(output_dir: str, ko_table: pd.DataFrame, metadata: qiime2.Metadata, pathway_id: str,
    map_ko: bool = False, low_color: str = "blue", high_color: str = "red") -> None:
    
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
    metadata_df.index

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

            output.to_csv()
            ## Save the image
            jsonp = prefix + ".jsonp"
            with open(os.path.join(output_dir, jsonp), 'w') as fh:
                fh.write('load_data("%s",' % prefix)
                json.dump(img_byte_arr, fh)
                fh.write(",'")
                table = q2templates.df_to_html(output)
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
            table = q2templates.df_to_html(output)
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
