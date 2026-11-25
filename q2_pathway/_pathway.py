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
import zipfile
import re
import textwrap

from ._consensus import summarize_gsea_directory
from ._stats import (
    harmonize_samples_for_visualization,
    summarize_direction_statistics,
    validate_table_annotations,
    welch_t_statistics,
)

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
        raise ImportError(
            "Please install `pykegg` for this module by running `pip install pykegg`"
        )

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
                """
                Calculation for t-statistics
                """

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
                """
                We should invert the levels when using effect in ALDEx2
                """
                sort_col = "wi.ep"
                prefix = column + "_" + level2 + "_vs_" + level1
                prefixes.append(prefix)

                """
                ALDEx2
                This will take time if you have many KOs across many metadata
                """
                metapath = path.join(output_dir, "meta_aldex2.tsv")
                kotablepath = path.join(output_dir, "ko_table.tsv")
                aldexpath = path.join(output_dir, "aldex2_res_" + prefix + ".tsv")
                aldeximagepath = path.join(output_dir, "aldex2_res_" + prefix + ".png")

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

                val = res[rank]
                val.index = ["ko:" + i for i in val.index.values]
                valdic = val.to_dict()
                valdics[prefix] = valdic

                aldex2_out = prefix + "_aldex2_res"
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
                deseq2imagepath = path.join(output_dir, "deseq2_res_" + prefix + ".png")

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
                highlight_nodes = [
                    "ko:" + i
                    for i in highlight_df[
                        highlight_df.thresholded == "True"
                    ].index.values
                ]
                highlight_value = []
                for node in nodes["name"]:
                    in_node = [i.replace("...", "") for i in node.split(" ")]
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
                    highlight_nodes="highlight",
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
            highlight_nodes = [
                "ko:" + i
                for i in highlight_df[highlight_df.thresholded == "True"].index.values
            ]
            highlight_value = []
            for node in nodes["name"]:
                in_node = [i.replace("...", "") for i in node.split(" ")]
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
                highlight_nodes="highlight",
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


def _safe_output_label(value: str) -> str:
    """Return a stable filename component without changing the display label."""
    label = re.sub(r"[^A-Za-z0-9._-]+", "_", str(value)).strip("._")
    return label or "unnamed"


def _append_method_order(image: Image.Image, table_names: list[str]) -> Image.Image:
    """Add the left-to-right node subdivision order below a KEGG image."""
    label = "Node subdivisions (left to right): " + " | ".join(
        f"{i + 1}. {name}" for i, name in enumerate(table_names)
    )
    lines = textwrap.wrap(label, width=105) or [label]
    from PIL import ImageDraw, ImageFont

    try:
        font = ImageFont.truetype("DejaVuSans.ttf", 14)
    except OSError:
        font = ImageFont.load_default()
    line_height = 20
    footer_height = 12 + line_height * len(lines)
    canvas = Image.new("RGB", (image.width, image.height + footer_height), "white")
    canvas.paste(image.convert("RGB"), (0, 0))
    draw = ImageDraw.Draw(canvas)
    for i, line in enumerate(lines):
        draw.text((10, image.height + 6 + i * line_height), line, fill="black", font=font)
    return canvas


def _r_method_statistics(
    table: pd.DataFrame,
    metadata_df: pd.DataFrame,
    column: str,
    level1: str,
    level2: str,
    method: str,
    rank: str,
    mc_samples: int,
    work_dir: str,
) -> pd.Series:
    """Run the existing R differential-statistic helpers for one KO table."""
    os.makedirs(work_dir, exist_ok=True)
    selected = metadata_df[
        metadata_df[column].astype(str).isin([str(level1), str(level2)])
    ][[column]].copy()
    selected = selected.loc[selected.index.intersection(table.index)]
    selected = selected.loc[table.index.intersection(selected.index)]
    selected_table = table.loc[selected.index]
    if selected[column].astype(str).value_counts().reindex(
        [str(level1), str(level2)], fill_value=0
    ).min() < 2:
        raise ValueError(
            "Differential statistics require at least two samples in each group."
        )

    table_path = os.path.join(work_dir, "ko_table.tsv")
    metadata_path = os.path.join(work_dir, "metadata.tsv")
    result_path = os.path.join(work_dir, "statistics.tsv")
    image_path = os.path.join(work_dir, "diagnostic.png")
    selected_table.to_csv(table_path, sep="\t")
    selected.to_csv(metadata_path, sep="\t")

    if method == "deseq2":
        subprocess.run(
            [
                "Rscript",
                path.join(TEMPLATES, "make_deseq2.R"),
                table_path,
                metadata_path,
                column,
                work_dir,
            ],
            check=True,
        )
        subprocess.run(
            [
                "Rscript",
                path.join(TEMPLATES, "perform_deseq2.R"),
                table_path,
                metadata_path,
                column,
                result_path,
                image_path,
                str(level1),
                str(level2),
                work_dir,
            ],
            check=True,
        )
    elif method == "aldex2":
        subprocess.run(
            [
                "Rscript",
                path.join(TEMPLATES, "perform_aldex2.R"),
                table_path,
                metadata_path,
                column,
                result_path,
                image_path,
                str(mc_samples),
                work_dir,
            ],
            check=True,
        )
    else:
        raise ValueError("method must be t, deseq2, or aldex2.")

    result = pd.read_csv(result_path, sep="\t", index_col=0)
    if rank not in result.columns:
        raise ValueError(f"Statistic column {rank!r} was not produced by {method}.")
    return pd.to_numeric(result[rank], errors="coerce")


def kegg_multi(
    output_dir: str,
    tables: list[pd.DataFrame],
    metadata: qiime2.Metadata,
    pathway_id: str,
    tables_name: list[str] = None,
    metadata_column: str = None,
    level1: str = None,
    level2: str = None,
    map_ko: bool = False,
    low_color: str = "blue",
    mid_color: str = "white",
    high_color: str = "red",
    missing_color: str = "#808080",
    tss: bool = False,
    method: str = "t",
    rank: str = None,
    mc_samples: int = 128,
    direction_threshold: float = 1.0,
    min_methods: int = 2,
) -> None:
    """Map method-specific KO statistics side by side on each KEGG node.

    Every node subdivision corresponds to one input table in ``tables_name``
    order.  Samples are harmonized globally, but KO catalogues are not
    intersected: a KO absent from a method is rendered with ``missing_color``.
    """
    if len(tables) < 2:
        raise ValueError("kegg_multi requires at least two KO tables.")
    names = (
        [f"data{i}" for i in range(len(tables))]
        if tables_name is None
        else [str(name).strip() for name in tables_name]
    )
    if len(names) != len(tables):
        raise ValueError("tables_name must have one value per input table.")
    if any(not name for name in names) or len(set(names)) != len(names):
        raise ValueError("tables_name values must be non-empty and unique.")

    method = str(method).lower()
    if method not in {"t", "deseq2", "aldex2"}:
        raise ValueError("method must be t, deseq2, or aldex2.")
    if tss and method in {"deseq2", "aldex2"}:
        raise ValueError("tss cannot be combined with DESeq2 or ALDEx2.")
    if rank is None:
        rank = {"t": "statistic", "deseq2": "log2FoldChange", "aldex2": "effect"}[
            method
        ]

    metadata_df = metadata.to_dataframe()
    tables, metadata_df, harmonization = harmonize_samples_for_visualization(
        tables, metadata_df, names
    )
    if tss:
        tables = [table.div(table.sum(axis=1), axis=0) for table in tables]

    harmonization.to_csv(
        os.path.join(output_dir, "sample_harmonization.tsv"), sep="\t", index=False
    )
    pd.DataFrame(
        {
            "subdivision": range(1, len(names) + 1),
            "table_name": names,
        }
    ).to_csv(os.path.join(output_dir, "method_order.tsv"), sep="\t", index=False)

    if metadata_column is not None:
        if metadata_column not in metadata_df.columns:
            raise ValueError(f"Metadata column {metadata_column!r} was not found.")
        columns = [metadata_column]
    else:
        columns = [
            column
            for column in metadata_df.columns
            if not pd.api.types.is_numeric_dtype(metadata_df[column])
            and 1 < metadata_df[column].nunique(dropna=True) < len(metadata_df)
        ]
    if not columns:
        raise ValueError("No categorical metadata column is available for comparison.")
    if (level1 is None) != (level2 is None):
        raise ValueError("level1 and level2 must be provided together.")
    if level1 is not None and metadata_column is None:
        raise ValueError("metadata_column is required when selecting levels.")

    try:
        import pykegg
    except ImportError as error:
        raise ImportError(
            "Please install pykegg for KEGG pathway rendering: pip install pykegg"
        ) from error

    graph = pykegg.KGML_graph(pid=pathway_id)
    base_nodes = graph.get_nodes()
    base_nodes = base_nodes[base_nodes.original_type == "ortholog"].copy()
    pathway_kos = {
        ko.replace("...", "")
        for node in base_nodes["name"].astype(str)
        for ko in node.split(" ")
    }
    ko_descriptions = None
    if map_ko:
        ko_names = pd.read_csv("https://rest.kegg.jp/list/ko", sep="\t", header=None)
        ko_names.index = ko_names[0].apply(lambda value: "ko:" + value)
        ko_descriptions = ko_names[1].to_dict()

    filenames = []
    for column in columns:
        available_levels = sorted(
            metadata_df[column].dropna().astype(str).unique().tolist()
        )
        if level1 is not None:
            if str(level1) not in available_levels or str(level2) not in available_levels:
                raise ValueError("The selected levels are absent from the metadata column.")
            contrasts = [(str(level1), str(level2))]
        else:
            contrasts = list(combinations(available_levels, 2))

        for first_level, second_level in contrasts:
            contrast_label = f"{column}_{first_level}_vs_{second_level}"
            safe_contrast = _safe_output_label(contrast_label)
            method_stats: dict[str, pd.Series] = {}
            for index, (name, table) in enumerate(zip(names, tables)):
                if method == "t":
                    statistic = welch_t_statistics(
                        table,
                        metadata_df[column],
                        first_level,
                        second_level,
                    )
                else:
                    statistic = _r_method_statistics(
                        table,
                        metadata_df,
                        column,
                        first_level,
                        second_level,
                        method,
                        rank,
                        mc_samples,
                        os.path.join(
                            output_dir,
                            "method_diagnostics",
                            f"{index + 1}_{_safe_output_label(name)}",
                            safe_contrast,
                        ),
                    )
                statistic.index = statistic.index.astype(str).str.removeprefix("ko:")
                method_stats[name] = statistic

            union_kos = sorted(
                set().union(*(set(statistic.index) for statistic in method_stats.values()))
            )
            ko_statistics = pd.DataFrame(index=union_kos)
            for name in names:
                ko_statistics[name] = method_stats[name].reindex(union_kos)
            support = summarize_direction_statistics(
                ko_statistics,
                direction_threshold=direction_threshold,
                min_methods=min_methods,
            )
            support.index.name = "KO"
            support_path = f"ko_support_{safe_contrast}.tsv"
            support.to_csv(os.path.join(output_dir, support_path), sep="\t")

            finite_pathway_values = []
            for statistic in method_stats.values():
                for ko_id, value in statistic.items():
                    prefixed = ko_id if ko_id.startswith("ko:") else "ko:" + ko_id
                    if prefixed in pathway_kos and np.isfinite(value):
                        finite_pathway_values.append(float(value))
            if not finite_pathway_values:
                raise ValueError(
                    f"No finite KO statistics map to {pathway_id} for {contrast_label}."
                )
            color_limit = max(abs(value) for value in finite_pathway_values)
            if color_limit == 0:
                color_limit = 1.0

            nodes = base_nodes.copy()
            color_columns = []
            statistic_columns = []
            for index, name in enumerate(names):
                color_column = f"method_color_{index + 1}"
                statistic_column = f"statistic_{index + 1}"
                values = {
                    (ko if ko.startswith("ko:") else "ko:" + ko): float(value)
                    for ko, value in method_stats[name].items()
                    if np.isfinite(value)
                }
                nodes = pykegg.append_colors_continuous_values(
                    nodes,
                    values,
                    new_color_column=color_column,
                    two_slope=True,
                    center_value=0,
                    node_name_column="name",
                    delim=" ",
                    colors=[low_color, mid_color, high_color],
                    orig_value=statistic_column,
                    fix_min=-color_limit,
                    fix_center=0,
                    fix_max=color_limit,
                )
                nodes[color_column] = nodes[color_column].fillna(missing_color)
                color_columns.append(color_column)
                statistic_columns.append(statistic_column)
            nodes["method_colors"] = nodes[color_columns].values.tolist()

            rendered = pykegg.overlay_opencv_image(
                nodes,
                pid=pathway_id,
                fill_color="method_colors",
                transparent_colors=["#BFBFFF", "#FFFFFF", "#BFFFBF"],
            )
            rendered = pykegg.append_legend(
                rendered,
                min_value=-color_limit,
                max_value=color_limit,
                center_value=0,
                two_slope=True,
                colors=[low_color, mid_color, high_color],
                width=4,
                label=f"{rank}: {first_level} - {second_level}",
            )
            rendered_image = _append_method_order(Image.fromarray(rendered), names)
            image_name = f"multi_{safe_contrast}_{_safe_output_label(pathway_id)}.png"
            rendered_image.save(os.path.join(output_dir, image_name))

            node_output = nodes.loc[:, ["name"] + statistic_columns].copy()
            node_output.columns = ["KEGG_node"] + names
            if ko_descriptions is not None:
                node_output["description"] = node_output["KEGG_node"].apply(
                    lambda value: " | ".join(
                        ko_descriptions.get(ko, "") for ko in value.split(" ")
                    )
                )
            node_path = f"node_statistics_{safe_contrast}.tsv"
            node_output.to_csv(os.path.join(output_dir, node_path), sep="\t", index=False)
            node_output["KEGG_node"] = node_output["KEGG_node"].apply(
                lambda value: " ".join(
                    f'<a href="{ko_url}{ko}">{ko}</a>' for ko in value.split(" ")
                )
            )

            jsonp = f"multi_{safe_contrast}.jsonp"
            with open(os.path.join(output_dir, jsonp), "w") as handle:
                handle.write(f'load_data("{safe_contrast}",')
                json.dump(image_name, handle)
                handle.write(",'" )
                table_html = q2templates.df_to_html(node_output, escape=False)
                handle.write(table_html.replace("\n", "").replace("'", "\\'"))
                handle.write("','")
                handle.write(quote(node_path))
                handle.write("');")
            filenames.append(jsonp)

    index_path = os.path.join(TEMPLATES, "index.html")
    q2templates.render(
        index_path,
        output_dir,
        context={"columns": [quote(filename) for filename in filenames]},
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
    split_str: str = None,
    table_roles: list = None,
    algorithm_families: list = None,
    top_set_alpha: float = 0.05,
    consensus_alpha: float = 0.05,
    direction_threshold: float = 0.8,
    min_families: int = 2,
    p_adjust_method: str = "bh",
):
    """
    Perform gene set enrichment analysis based on pathway to KO relationship
    obtained from KEGG API. T-statistics and ALDEx2 and DESeq2 statistics can
    be used for ranking the genes.

    [TODO] More controlling options for fgsea parameters.
    Currently, min and max gene set size can be controlled.

    Parameters:
    -----------
    rank: str
        Which column to use for ranking in ALDEx2 and DESeq2.
    same: bool
        Use same gene set across the gene family tables.
    methods: list of str
        Methods to rank the gene families per tables.
        If one method is specified, the method will be used for all the tables.
    """
    table_names, table_roles, algorithm_families = validate_table_annotations(
        len(tables),
        tables_name,
        table_roles,
        algorithm_families,
        min_predictions=1,
    )

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

    if split_str is not None:
        tables = [
            table.set_axis(
                [i.split(quote(split_str))[0] for i in table.index.values.tolist()]
            )
            for table in tables
        ]

    filenames = []
    gsea_run_records = []

    if map_pathway:
        if module:
            kon = pd.read_csv("https://rest.kegg.jp/list/module", sep="\t", header=None)
            kon[0] = kon[0].apply(lambda x: str("md:") + str(x))
            kon.to_csv(
                os.path.join(output_dir, "pathway_names.tsv"),
                sep="\t",
                index=False,
                header=None,
            )
            konflag = 1
        else:
            kon = pd.read_csv(
                "https://rest.kegg.jp/list/pathway", sep="\t", header=None
            )
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
        dataset_rank = rank
        if dataset_rank is None:
            if method == "aldex2":
                dataset_rank = "effect"
            if method == "deseq2":
                dataset_rank = "log2FoldChange"

        dataset_name = table_names[e]

        """
        Filter columns
        """
        tmp_metadata = metadata.filter_ids(ko_table.index)
        tmp_metadata = tmp_metadata.filter_columns(column_type="categorical")
        tmp_metadata = tmp_metadata.filter_columns(
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
                contrast_id = f"{column}__{level1}_vs_{level2}"

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

                    val = res[dataset_rank]
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

                    val = res[dataset_rank]
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
                    str(max_size),
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
                gsea_run_records.append(
                    {
                        "dataset": dataset_name,
                        "role": table_roles[e],
                        "algorithm_family": algorithm_families[e],
                        "method": method,
                        "column": column,
                        "level1": level1,
                        "level2": level2,
                        "contrast_id": contrast_id,
                        "prefix": prefix,
                    }
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
    downloads = []

    ## Finally, concatenate gsea results and summarize
    pd.DataFrame(gsea_run_records).to_csv(
        os.path.join(output_dir, "gsea_runs.tsv"), sep="\t", index=False
    )
    try:
        summarize_gsea_directory(
            output_dir,
            output_dir,
            top_set_alpha=top_set_alpha,
            consensus_alpha=consensus_alpha,
            direction_threshold=direction_threshold,
            min_families=min_families,
            p_adjust_method=p_adjust_method,
        )
    except (ValueError, FileNotFoundError) as error:
        raise ValueError(
            "GSEA summary (correlation/consensus) cannot be performed: " + str(error)
        ) from error
    zip_name = "GSEA_final_outputs.zip"
    zip_path = os.path.join(output_dir, zip_name)

    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as z:
        for fn in os.listdir(output_dir):
            if fn != zip_name:
                print(fn)
                fp = os.path.join(output_dir, fn)
                if os.path.exists(fp):
                    z.write(fp, arcname=fn)

    downloads.append({"href": zip_name, "label": zip_name})
    q2templates.render(
        index,
        output_dir,
        context={
            "columns": [quote(fn) for fn in filenames],
            "downloads": downloads,
        },
    )
    # shutil.copytree(os.path.join(TEMPLATES, "dist"), os.path.join(output_dir, "dist"))
