import shutil
import json
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
TEMPLATES = pkg_resources.resource_filename("q2_pathway", "assets")


def kegg(
    output_dir: str,
    ko_table: pd.DataFrame,
    metadata: qiime2.Metadata,
    pathway_id: str,
    highlight: qiime2.Metadata = None,
    map_ko: bool = False,
    low_color: str = "blue",
    high_color: str = "red",
    tss: bool = False,
    method: str = "t",
    rank: str = None,
    mc_samples: int = 128,
) -> None:
    
    """
    Specify ranking to be plotted,
    if not specified
    """
    if rank is None:
        if method == "aldex2":
            rank = "effect"
        if method == "deseq2":
            rank = "log2FoldChange"

    if highlight is not None:
        highlight_df = highlight.to_dataframe()
        if "thresholded" not in highlight_df.columns:
            raise ValueError("thresholded column not in highlight metadata")

    try:
        import pykegg
    except ImportError:
        raise ImportError("Please install `pykegg` for this module by running `pip install pykegg`")
        
    if tss:
        ko_table = ko_table.apply(lambda x: x / sum(x), axis=1)
        if method == "aldex2":
            raise ValueError("ALDEx2 is for the raw count data.")
        if method == "deseq2":
            raise ValueError("DESeq2 is for the raw count data.")

    """
    Filter columns
    """
    metadata = metadata.filter_ids(ko_table.index)
    metadata = metadata.filter_columns(column_type="categorical")
    metadata = metadata.filter_columns(
        drop_all_unique=True, drop_zero_variance=True, drop_all_missing=True
    )
    filenames = []

    metadata.save(os.path.join(output_dir, "metadata.tsv"))

    metadata_df = metadata.to_dataframe()

    """
    Obtain KGML and nodes DataFrame
    """
    graph = pykegg.KGML_graph(pid=pathway_id)
    nodes = graph.get_nodes()
    nodes = nodes[nodes.original_type == "ortholog"]

    """
    Obtain KO description
    """
    if map_ko:
        koname = pd.read_csv("https://rest.kegg.jp/list/ko", sep="\t", header=None)
        koname.index = koname[0].apply(lambda x: "ko:" + x)
        kodic = koname[1].to_dict()

    for column in metadata_df.columns:
        metadata_df_filt = metadata_df[metadata_df[column].notna()]
        uniq_level = metadata_df_filt[column].unique()
        base = list(combinations(uniq_level, 2))
        prefixes = []
        valdics = dict()

        """
        If DESeq2, first calculate the DESeq object
        """
        if method == "deseq2":
            metapath = path.join(output_dir, "meta_deseq2.tsv")
            kotablepath = path.join(output_dir, "ko_table.tsv")
            metadata_df_filt.to_csv(metapath, sep="\t")
            ko_table.to_csv(kotablepath, sep="\t")
            # First make DESeq2 object and store it.
            cmd = [
                "Rscript",
                path.join(TEMPLATES, "make_deseq2.R"),
                kotablepath,
                metapath,
                column,
                output_dir,
            ]
            try:
                res = subprocess.run(cmd, check=True)
            except subprocess.CalledProcessError as _:
                raise ValueError(
                    "DESeq2 cannot be performed due to the error in R script."
                )

        for comb in base:

            """
            User control for ordering these levels in the metadata columns
            """
            comb = sorted(comb)
            level1 = comb[0]
            level2 = comb[1]
            prefix = column + "_" + level1 + "_vs_" + level2

            if method == "t":
                ## T-stats
                prefixes.append(prefix)
                a = metadata_df_filt[metadata_df_filt[column] == level1].index.tolist()
                b = metadata_df_filt[metadata_df_filt[column] == level2].index.tolist()
                val = pd.Series(
                    ko_table.columns.map(
                        lambda x: scipy.stats.ttest_ind(
                            ko_table.loc[a, x], ko_table.loc[b, x], equal_var=False
                        ).statistic
                    )
                )
                val.index = ["ko:" + i for i in ko_table.columns]
                valdic = val.to_dict()
                valdics[prefix] = valdic
            elif method == "aldex2":
                ## We should invert the levels when using effect
                sort_col = "wi.ep"
                prefix = column + "_" + level2 + "_vs_" + level1
                prefixes.append(prefix)
                ## ALDEx2
                ## This will take time if you have many KOs across many metadata
                metapath = path.join(output_dir, "meta_aldex2.tsv")
                kotablepath = path.join(output_dir, "ko_table.tsv")
                aldexpath = path.join(output_dir, "aldex2_res_" + prefix + ".tsv")
                aldeximagepath = path.join(output_dir, "aldex2_res_" + prefix + ".png")

                ## Make two conditions
                metadata_df_tmp = metadata_df_filt[
                    (metadata_df_filt[column] == level1)
                    | (metadata_df_filt[column] == level2)
                ]
                metadata_df_tmp.to_csv(metapath, sep="\t")
                ko_table.to_csv(kotablepath, sep="\t")
                cmd = [
                    "Rscript",
                    path.join(TEMPLATES, "perform_aldex2.R"),
                    kotablepath,
                    metapath,
                    column,
                    aldexpath,
                    aldeximagepath,
                    str(mc_samples),
                    output_dir,
                ]
                try:
                    res = subprocess.run(cmd, check=True)
                except subprocess.CalledProcessError as _:
                    raise ValueError(
                        "ALDEx2 cannot be performed due to the error in R script."
                    )
                res = pd.read_csv(aldexpath, sep="\t", index_col=0, header=0)

                ## Better to select by the  users
                val = res[rank]
                val.index = ["ko:" + i for i in val.index.values]
                valdic = val.to_dict()
                valdics[prefix] = valdic

                aldex2_out = prefix + "_aldex2_res"
                ## Save the aldex2 out
                jsonp = aldex2_out + ".jsonp"

                res = res.sort_values(by=sort_col).head(50)

                with open(os.path.join(output_dir, jsonp), "w") as fh:
                    fh.write('load_data("%s",' % aldex2_out)
                    json.dump("aldex2_res_" + prefix + ".png", fh)
                    fh.write(",'")
                    table = q2templates.df_to_html(res, escape=False)
                    fh.write(table.replace("\n", "").replace("'", "\\'"))
                    fh.write("','")
                    fh.write(quote("aldex2_res_" + prefix + ".tsv"))
                    fh.write("');")
                filenames.append(jsonp)
            elif method == "deseq2":

                """
                DESeq2 processing
                Probably use PyDESeq2: consider the version compatibility for conda installation
                Consider prefiltering option, but it is better to handle in `filter-feature`
                """
                sort_col = "padj"
                prefix = column + "_" + level1 + "_vs_" + level2
                prefixes.append(prefix)

                deseq2path = path.join(output_dir, "deseq2_res_" + prefix + ".tsv")
                deseq2imagepath = path.join(
                    output_dir, "deseq2_res_" + prefix + ".png"
                )

                
                """
                Make two conditions
                """
                cmd = [
                    "Rscript",
                    path.join(TEMPLATES, "perform_deseq2.R"),
                    kotablepath,
                    metapath,
                    column,
                    deseq2path,
                    deseq2imagepath,
                    level1,
                    level2,
                    output_dir,
                ]
                try:
                    res = subprocess.run(cmd, check=True)
                except subprocess.CalledProcessError as _:
                    raise ValueError(
                        "DESeq2 cannot be performed due to the error in R script."
                    )
                res = pd.read_csv(deseq2path, sep="\t", index_col=0, header=0)

                val = res[rank]
                val.index = ["ko:" + i for i in val.index.values]
                valdic = val.to_dict()
                valdics[prefix] = valdic

                val.to_csv(
                    os.path.join(output_dir, "values_" + prefix + ".tsv"),
                    sep="\t",
                    header=None,
                )

                deseq2_out = prefix + "_deseq2_res"
                ## Save the aldex2 out
                jsonp = deseq2_out + ".jsonp"

                res = res.sort_values(by=sort_col).head(50)

                with open(os.path.join(output_dir, jsonp), "w") as fh:
                    fh.write('load_data("%s",' % deseq2_out)
                    json.dump("deseq2_res_" + prefix + ".png", fh)
                    fh.write(",'")
                    table = q2templates.df_to_html(res, escape=False)
                    fh.write(table.replace("\n", "").replace("'", "\\'"))
                    fh.write("','")
                    fh.write(quote("deseq2_res_" + prefix + ".tsv"))
                    fh.write("');")
                filenames.append(jsonp)

            else:
                raise ValueError("Method should be set to t, deseq2, or aldex2")

            nodes = pykegg.append_colors_continuous_values(
                nodes,
                valdic,
                new_color_column=prefix,
                two_slope=False,
                node_name_column="name",
                delim=" ",
                colors=[low_color, high_color],
                orig_value="stat_" + prefix,
            )

            output = nodes.loc[:, ["name", "stat_" + prefix]]
            if map_ko:
                output["description"] = output["name"].apply(
                    lambda x: " ".join([kodic[q] for q in x.split(" ")])
                )

            nodes[prefix] = nodes[prefix].fillna("#808080")

            if highlight is not None:
                highlight_nodes = ["ko:"+i for i in highlight_df[highlight_df.thresholded=="True"].index.values]
                highlight_value = []
                for node in nodes["name"]:
                    in_node = [i.replace("...","") for i in node.split(" ")]
                    intersect = set(in_node) & set(highlight_nodes)
                    if len(intersect) > 0:
                        highlight_value.append(True)
                    else:
                        highlight_value.append(False)

                nodes["highlight"] = highlight_value
                nodes = nodes[~nodes["x"].isna()]

                kegg_map_image1 = pykegg.overlay_opencv_image(
                    nodes,
                    pid=pathway_id,
                    fill_color=prefix,
                    transparent_colors=["#BFBFFF", "#FFFFFF", "#BFFFBF"],
                    highlight_nodes="highlight"
                )
            else:
                kegg_map_image1 = pykegg.overlay_opencv_image(
                    nodes,
                    pid=pathway_id,
                    fill_color=prefix,
                    transparent_colors=["#BFBFFF", "#FFFFFF", "#BFFFBF"],
                )


            kegg_map_image2 = pykegg.append_legend(
                kegg_map_image1,
                min_value=np.nanmin([valdic[i] for i in valdic.keys()]),
                max_value=np.nanmax([valdic[i] for i in valdic.keys()]),
                two_slope=False,
                colors=[low_color, high_color],
                width=4,
                label="Statistics",
            )
            kegg_map_image = Image.fromarray(kegg_map_image2)
            kegg_map_image.save(
                os.path.join(output_dir, prefix + "_" + pathway_id + ".png")
            )

            csv_path = os.path.join(output_dir, prefix + ".csv")
            output.to_csv(csv_path)

            output["name"] = output["name"].apply(
                lambda x: " ".join(
                    ['<a href="' + ko_url + i + '">' + i + "</a>" for i in x.split(" ")]
                )
            )

            jsonp = prefix + ".jsonp"
            with open(os.path.join(output_dir, jsonp), "w") as fh:
                fh.write('load_data("%s",' % prefix)
                json.dump(prefix + "_" + pathway_id + ".png", fh)
                fh.write(",'")
                table = q2templates.df_to_html(output, escape=False)
                fh.write(table.replace("\n", "").replace("'", "\\'"))
                fh.write("','")
                fh.write(quote(prefix + ".csv"))
                fh.write("');")
            filenames.append(jsonp)

        """
        Concat all the values
        """
        minval = np.nanmin(nodes.loc[:, ["stat_" + p for p in prefixes]].values)
        maxval = np.nanmax(nodes.loc[:, ["stat_" + p for p in prefixes]].values)
        for p in prefixes:
            nodes = pykegg.append_colors_continuous_values(
                nodes,
                valdics[p],
                new_color_column=p,
                two_slope=False,
                node_name_column="name",
                delim=" ",
                colors=[low_color, high_color],
                orig_value="stat_" + p,
                fix_min=minval,
                fix_max=maxval,
            )

        output = nodes.loc[:, ["name"] + ["stat_" + p for p in prefixes]]
        csv_path = os.path.join(output_dir, column + ".csv")
        output.to_csv(csv_path)

        if map_ko:
            output["description"] = output["name"].apply(
                lambda x: " ".join([kodic[q] for q in x.split(" ")])
            )

        output["name"] = output["name"].apply(
            lambda x: " ".join(
                ['<a href="' + ko_url + i + '">' + i + "</a>" for i in x.split(" ")]
            )
        )
        nodes["concat"] = nodes.loc[:, prefixes].values.tolist()
        nodes["concat"] = nodes["concat"].apply(
            lambda x: [
                i if isinstance(i, str) and i.startswith("#") else "#808080" for i in x
            ]
        )

        """
        Mark high correlation if specified
        """
        if highlight is not None:
            highlight_nodes = ["ko:"+i for i in highlight_df[highlight_df.thresholded=="True"].index.values]
            highlight_value = []
            for node in nodes["name"]:
                in_node = [i.replace("...","") for i in node.split(" ")]
                intersect = set(in_node) & set(highlight_nodes)
                if len(intersect) > 0:
                    highlight_value.append(True)
                else:
                    highlight_value.append(False)
            nodes["highlight"] = highlight_value
            kegg_map_image1 = pykegg.overlay_opencv_image(
                nodes,
                pid=pathway_id,
                fill_color="concat",
                transparent_colors=["#BFBFFF", "#FFFFFF", "#BFFFBF"],
                highlight_nodes="highlight"
            )
        else:
            kegg_map_image1 = pykegg.overlay_opencv_image(
                nodes,
                pid=pathway_id,
                fill_color="concat",
                transparent_colors=["#BFBFFF", "#FFFFFF", "#BFFFBF"],
            )



        kegg_map_image2 = pykegg.append_legend(
            kegg_map_image1,
            min_value=minval,
            max_value=maxval,
            two_slope=False,
            colors=[low_color, high_color],
            width=4,
            label="Statistics",
        )
        kegg_map_image = Image.fromarray(kegg_map_image2)
        kegg_map_image.save(
            os.path.join(output_dir, column + "_" + pathway_id + ".png")
        )

        """
        Save the pathway with whole comparison
        """
        jsonp = column + ".jsonp"
        with open(os.path.join(output_dir, jsonp), "w") as fh:
            fh.write('load_data("%s",' % column)
            json.dump(column + "_" + pathway_id + ".png", fh)
            fh.write(",'")
            table = q2templates.df_to_html(output, escape=False)
            fh.write(table.replace("\n", "").replace("'", "\\'"))
            fh.write("','")
            fh.write(quote(column + ".csv"))
            fh.write("');")
        filenames.append(jsonp)

    index = os.path.join(TEMPLATES, "index.html")
    q2templates.render(
        index, output_dir, context={"columns": [quote(fn) for fn in filenames]}
    )

    shutil.copytree(os.path.join(TEMPLATES, "dist"), os.path.join(output_dir, "dist"))


def trunc(x):
    tmp = x.split("/")
    if len(tmp) >= 5:
        truncstr = "/".join(tmp[0:5]) + str("/truncated")
    else:
        truncstr = "/".join(tmp)
    return truncstr


def gsea(
    output_dir: str,
    tables: pd.DataFrame,
    metadata: qiime2.Metadata,
    tss: bool = False,
    methods: list = ["t"],
    mc_samples: int = 128,
    module: bool = False,
    map_pathway: bool = False,
    tables_name: str = None,
    bg: str = "all",
    rank: str = None,
    same: bool = False,
    min_size: int = 0,
    max_size: int = 500,

):
    """
    Perform gene set enrichment analysis based on pathway to KO relationship
    obtained from KEGG API. T-statistics and ALDEx2 and DESeq2 statistics can
    be used for ranking the genes.

    [TODO] More controlling options for fgsea parameters
    
    Parameters:
    -----------
    rank: str
        Which column to use for ranking in ALDEx2 and DESeq2.
    same: bool
        Use same gene set across the gene family tables.
    """
    if len(methods) > 1:
        if len(methods) != len(tables):
            raise ValueError("Method length and table length do not match")
    else:
        methods = [methods[0] for i in range(len(tables))]

    if same:
        kos = [ko_table.columns.values for ko_table in tables]
        common_kos = list(set.intersection(*map(set, kos)))
        tables = [table.loc[:, common_kos] for table in tables]

    if tss:
        tables = [table.apply(lambda x: x / sum(x), axis=1) for table in tables]
        if "aldex2" in methods:
            raise ValueError("ALDEx2 is for the count data.")
        if "deseq2" in methods:
            raise ValueError("DESeq2 is for the raw count data.")

    filenames = []

    if map_pathway:
        kon = pd.read_csv("https://rest.kegg.jp/list/pathway", sep="\t", header=None)
        kon[0] = kon[0].apply(lambda x: x.replace("map", "path:ko"))
        kon.to_csv(
            os.path.join(output_dir, "pathway_names.tsv"),
            sep="\t",
            index=False,
            header=None,
        )
        konflag = 1
    else:
        konflag = 0

    """
    Download pathway (or module) info for GSEA and output
    """
    if module:
        kop = pd.read_csv("https://rest.kegg.jp/link/ko/module", sep="\t", header=None)
        kop.to_csv(
            os.path.join(output_dir, "pathway_map.tsv"),
            sep="\t",
            index=False,
            header=None,
        )
    else:
        kop = pd.read_csv("https://rest.kegg.jp/link/ko/pathway", sep="\t", header=None)
        kop = kop[kop[0].apply(lambda x: "path:ko" in x)]
        kop.to_csv(
            os.path.join(output_dir, "pathway_map.tsv"),
            sep="\t",
            index=False,
            header=None,
        )

    """
    save out metadata for download in viz
    """
    metadata.save(os.path.join(output_dir, "metadata.tsv"))


    for e, ko_table in enumerate(tables):
        method = methods[e]
        print("Dataset " + str(e) + ": method " + method)
        if rank is None:
            if method == "aldex2":
                rank = "effect"
            if method == "deseq2":
                rank = "log2FoldChange"

        if tables_name is not None:
            dataset_name = tables_name[e]
        else:
            dataset_name = "data" + str(e)

        """
        Filter columns
        """
        tmp_metadata = metadata.filter_ids(ko_table.index)
        tmp_metadata = metadata.filter_columns(column_type="categorical")
        tmp_metadata = metadata.filter_columns(
            drop_all_unique=True, drop_zero_variance=True, drop_all_missing=True
        )

        metadata_df = tmp_metadata.to_dataframe()

        for column in metadata_df.columns:
            metadata_df_filt = metadata_df[metadata_df[column].notna()]
            uniq_level = metadata_df_filt[column].unique()
            base = list(combinations(uniq_level, 2))

            if method == "deseq2":
                metapath = path.join(output_dir, "meta_deseq2.tsv")
                kotablepath = path.join(output_dir, "ko_table.tsv")
                metadata_df_filt.to_csv(metapath, sep="\t")
                ko_table.to_csv(kotablepath, sep="\t")

                """
                First make DESeq2 object and store it.
                """
                cmd = [
                    "Rscript",
                    path.join(TEMPLATES, "make_deseq2.R"),
                    kotablepath,
                    metapath,
                    column,
                    output_dir,
                ]
                try:
                    res = subprocess.run(cmd, check=True)
                except subprocess.CalledProcessError as _:
                    raise ValueError(
                        "DESeq2 cannot be performed due to the error in R script."
                    )

            for comb in base:
                """
                Sort the level
                """
                comb = sorted(comb)

                ## User control for ordering these values
                level1 = comb[0]
                level2 = comb[1]

                if method == "t":
                    """
                    T-statistics from scipy.stats.ttest_ind
                    """
                    prefix = (
                        dataset_name + "_" + column + "_" + level1 + "_vs_" + level2
                    )

                    a = metadata_df_filt[
                        metadata_df_filt[column] == level1
                    ].index.tolist()
                    b = metadata_df_filt[
                        metadata_df_filt[column] == level2
                    ].index.tolist()
                    val = pd.Series(
                        ko_table.columns.map(
                            lambda x: scipy.stats.ttest_ind(
                                ko_table.loc[a, x], ko_table.loc[b, x], equal_var=False
                            ).statistic
                        )
                    )
                    val.index = ["ko:" + i for i in ko_table.columns]
                    val.to_csv(
                        os.path.join(output_dir, "values_" + prefix + ".tsv"),
                        sep="\t",
                        header=None,
                    )
                elif method == "aldex2":
                    """
                    We should invert the levels when using effect
                    """
                    prefix = (
                        dataset_name + "_" + column + "_" + level2 + "_vs_" + level1
                    )

                    metapath = path.join(output_dir, "meta_aldex2.tsv")
                    kotablepath = path.join(output_dir, "ko_table.tsv")
                    aldexpath = path.join(output_dir, "aldex2_res_" + prefix + ".tsv")
                    aldeximagepath = path.join(
                        output_dir, "aldex2_res_" + prefix + ".png"
                    )

                    """
                    Make two conditions
                    """
                    metadata_df_tmp = metadata_df_filt[
                        (metadata_df_filt[column] == level1)
                        | (metadata_df_filt[column] == level2)
                    ]
                    metadata_df_tmp.to_csv(metapath, sep="\t")
                    ko_table.to_csv(kotablepath, sep="\t")
                    cmd = [
                        "Rscript",
                        path.join(TEMPLATES, "perform_aldex2.R"),
                        kotablepath,
                        metapath,
                        column,
                        aldexpath,
                        aldeximagepath,
                        str(mc_samples),
                        output_dir,
                    ]
                    try:
                        res = subprocess.run(cmd, check=True)
                    except subprocess.CalledProcessError as _:
                        raise ValueError("ALDEx2 cannot be performed")
                    res = pd.read_csv(aldexpath, sep="\t", index_col=0, header=0)

                    val = res[rank]
                    val.index = ["ko:" + i for i in val.index.values]
                    val.to_csv(
                        os.path.join(output_dir, "values_" + prefix + ".tsv"),
                        sep="\t",
                        header=None,
                    )

                    aldex2_out = prefix + "_aldex2_res"
                    ## Save the aldex2 out
                    jsonp = aldex2_out + ".jsonp"

                    sort_col = "wi.ep"
                    res = res.sort_values(by=sort_col).head(50)

                    with open(os.path.join(output_dir, jsonp), "w") as fh:
                        fh.write('load_data("%s",' % aldex2_out)
                        json.dump("aldex2_res_" + prefix + ".png", fh)
                        fh.write(",'")
                        table = q2templates.df_to_html(res, escape=False)
                        fh.write(table.replace("\n", "").replace("'", "\\'"))
                        fh.write("','")
                        fh.write(quote("aldex2_res_" + prefix + ".tsv"))
                        fh.write("');")
                    filenames.append(jsonp)

                elif method == "deseq2":
                    """
                    DESeq2 processing

                    Probably use PyDESeq2: consider the version compatibility for conda installation
                    Consider prefiltering option, but it is better to handle in `filter-feature` in QIIME 2
                    for the better reproduciblity.
                    """
                    sort_col = "padj"
                    prefix = (
                        dataset_name + "_" + column + "_" + level1 + "_vs_" + level2
                    )


                    deseq2path = path.join(output_dir, "deseq2_res_" + prefix + ".tsv")
                    deseq2imagepath = path.join(
                        output_dir, "deseq2_res_" + prefix + ".png"
                    )

                    ## Make two conditions
                    cmd = [
                        "Rscript",
                        path.join(TEMPLATES, "perform_deseq2.R"),
                        kotablepath,
                        metapath,
                        column,
                        deseq2path,
                        deseq2imagepath,
                        level1,
                        level2,
                        output_dir,
                    ]
                    try:
                        res = subprocess.run(cmd, check=True)
                    except subprocess.CalledProcessError as _:
                        raise ValueError(
                            "DESeq2 cannot be performed due to the error in R script."
                        )
                    res = pd.read_csv(deseq2path, sep="\t", index_col=0, header=0)

                    val = res[rank]
                    val.index = ["ko:" + i for i in val.index.values]

                    val.to_csv(
                        os.path.join(output_dir, "values_" + prefix + ".tsv"),
                        sep="\t",
                        header=None,
                    )

                    deseq2_out = prefix + "_deseq2_res"
                    jsonp = deseq2_out + ".jsonp"

                    res = res.sort_values(by=sort_col).head(50)

                    with open(os.path.join(output_dir, jsonp), "w") as fh:
                        fh.write('load_data("%s",' % deseq2_out)
                        json.dump("deseq2_res_" + prefix + ".png", fh)
                        fh.write(",'")
                        table = q2templates.df_to_html(res, escape=False)
                        fh.write(table.replace("\n", "").replace("'", "\\'"))
                        fh.write("','")
                        fh.write(quote("deseq2_res_" + prefix + ".tsv"))
                        fh.write("');")
                    filenames.append(jsonp)
                else:
                    raise ValueError("Method should be set to t, deseq2, or aldex2")

                """
                Perform GSEA
                """
                if bg != "all":
                    pd.Series(["ko:" + i for i in ko_table.columns]).to_csv(
                        os.path.join(output_dir, prefix + "_all_KO.txt"),
                        sep="\t",
                        index=False,
                        header=False,
                    )
                cmd = [
                    "Rscript",
                    path.join(TEMPLATES, "perform_gsea.R"),
                    output_dir,
                    prefix,
                    str(konflag),
                    bg,
                    str(min_size),
                    str(max_size)
                ]
                try:
                    res = subprocess.run(cmd, check=True)
                except subprocess.CalledProcessError as _:
                    raise ValueError("GSEA cannot be performed")

                gseares = pd.read_csv(
                    path.join(output_dir, "gsea_res_" + prefix + ".tsv"),
                    sep="\t",
                    index_col=0,
                    header=0,
                )
                if module:
                    gseares["pathway"] = gseares["pathway"].apply(
                        lambda x: '<a href="https://www.kegg.jp/entry/module+'
                        + x.split(":")[1]
                        + '">'
                        + x
                        + "</a>"
                    )
                else:
                    gseares["pathway"] = gseares["pathway"].apply(
                        lambda x: '<a href="https://www.kegg.jp/entry/pathway+'
                        + x.split(":")[1]
                        + '">'
                        + x
                        + "</a>"
                    )
                jsonp = prefix + ".jsonp"
                gseares["leadingEdge"] = gseares["leadingEdge"].apply(
                    lambda x: trunc(x)
                )

                """
                The same structure as kegg
                """
                with open(os.path.join(output_dir, jsonp), "w") as fh:
                    fh.write('load_data("%s",' % prefix)
                    json.dump("gsea_res_" + prefix + ".png", fh)
                    fh.write(",'")
                    table = q2templates.df_to_html(gseares, escape=False)
                    fh.write(table.replace("\n", "").replace("'", "\\'"))
                    fh.write("','")
                    fh.write(quote("gsea_res_" + prefix + ".tsv"))
                    fh.write("');")
                filenames.append(jsonp)

    index = os.path.join(TEMPLATES, "index.html")
    q2templates.render(
        index, output_dir, context={"columns": [quote(fn) for fn in filenames]}
    )
    shutil.copytree(os.path.join(TEMPLATES, "dist"), os.path.join(output_dir, "dist"))
