from os import path
import pandas as pd
import qiime2
import shutil
import q2templates
import numpy as np
import os
import json
import seaborn as sns
from urllib.parse import quote
import scipy
import matplotlib.pyplot as plt
from itertools import combinations
from qiime2 import Metadata
import pkg_resources

TEMPLATES = pkg_resources.resource_filename("q2_pathway", "assets")
ko_url = "https://www.genome.jp/dbget-bin/www_bget?"


def output_json(prefix, output_dir, output):
    with open(os.path.join(output_dir, prefix + ".jsonp"), "w") as fh:
        fh.write('load_data("%s",' % prefix)
        json.dump(prefix + ".png", fh)
        fh.write(",'")
        table = q2templates.df_to_html(output, escape=False)
        fh.write(table.replace("\n", "").replace("'", "\\'"))
        fh.write("','")
        fh.write(quote(prefix + ".csv"))
        fh.write("');")


def summarize(
    output_dir: str,
    tables: pd.DataFrame,
    metadata: qiime2.Metadata = None,
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
    map_ko: bool = False,
    skip: bool = False,
) -> None:
    if map_ko:
        koname = pd.read_csv(
            "https://rest.kegg.jp/list/ko", sep="\t", header=None, index_col=0
        )
        kodic = koname[1].to_dict()

    if convert_table is not None:
        """
        Converting table should be Metadata with the column name of "converted".
        If shotgun profiled qza is to be read, the index should be shotgun-ID and `converted` column
        corresponds to 16S-ID.
        """
        mapping = convert_table.to_dataframe()
        change = mapping["converted"].to_dict()

    tbl_len = len(tables)

    if tss:
        tables = [ko_table.apply(lambda x: x / sum(x), axis=1) for ko_table in tables]

    strat = False

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

    if metadata is not None:
        metadata = metadata.filter_ids(all_samples)
        metadata = metadata.filter_columns(column_type="categorical")
        metadata = metadata.filter_columns(
            drop_all_unique=False, drop_zero_variance=False, drop_all_missing=True
        )
        metadata.save(os.path.join(output_dir, "metadata.tsv"))
        metadata_df = metadata.to_dataframe()
    else:
        metadata_df = pd.DataFrame(index=all_samples)
        metadata_df["All"] = "All"
        metadata_df.to_csv(os.path.join(output_dir, "metadata.tsv"), sep="\t")

    filenames = []


    if candidate is None:
        """
        First table will be used for subset
        By default, the genes are sorted based on mean abundance and `p-first` KOs
        are subset for summarization.
        """
        mean_val = (
            tables[0]
            .loc[:, all_kos]
            .apply(lambda x: np.nanmean(x), axis=0)
            .sort_values(ascending=False)
        )
        all_cols = mean_val.head(first).index.values
    else:
        all_cols = [candidate]

    """
    Override `first` and `candidate` option if pathway is specified
    """
    if candidate_pathway is not None:
        kop = pd.read_csv("https://rest.kegg.jp/link/ko/pathway", sep="\t", header=None)
        kop = kop[kop[0].apply(lambda x: "path:ko" in x)]
        kop[0] = kop[0].apply(lambda x: x.split(":")[1])
        kop[1] = kop[1].apply(lambda x: x.split(":")[1])
        all_cols = kop[kop[0].isin([candidate_pathway])][1].tolist()
        all_cols = list(set(all_kos) & set(all_cols))

    """
    Full (stratified output) processing
    

    Loose checking for stratified (assuming stratified table is shown as K00001|taxonomy1)
    As q2-picrust2 does not produce stratified table by default, what should be the format for ths stratified output?
    Currently, the correlation will not be produced for the stratified output but stratified abundances will be plotted.
    Also, only the first table will be summarized if multiple tables are given.
    """
    if "|" in tables[0].columns[0]:
        strat = True
        ko_table = tables[0]
        all_cols = ko_table.columns

    if strat:
        # all_taxs = list(set([i.split("_")[0] for i in all_cols]))
        if candidate is None:
            all_kos = list(set([i.split("|")[0] for i in all_cols]))
        else:
            all_kos = [candidate]

        for column in metadata_df.columns:
            metadata_df_filt = metadata_df[metadata_df[column].notna()]
            data = pd.concat([ko_table, metadata_df_filt], axis=1, join="inner")

            for ko in all_kos:
                prefix = column + "_" + ko

                candidate_columns = [i for i in all_cols if ko in i]
                if len(candidate_columns) == 0:
                    continue
                candidate_columns.append(column)
                output = data.loc[:, candidate_columns]
                output["sample"] = output.index.values
                output = pd.melt(output, id_vars=["sample", column])
                output["variable"] = output["variable"].apply(
                    lambda x: "|".join(x.split("|")[1:])
                )

                plt.figure()
                bp = sns.boxplot(data=output, x="variable", y="value")
                fig = bp.get_figure()
                fig.savefig(path.join(output_dir, prefix + ".png"))
                plt.close()

                csv_path = os.path.join(output_dir, prefix + ".csv")

                output = output.groupby(column).apply(
                    lambda x: x.groupby("variable").mean("value")
                )

                if map_ko:
                    if ko in kodic.keys():
                        output["ko_description"] = kodic[ko]

                output.to_csv(csv_path)

                if ko.startswith("K"):
                    output["ko"] = '<a href="' + ko_url + ko + '">' + ko + "</a>"
                else:
                    output["ko"] = ko

                output_json(prefix, output_dir, output)
                filenames.append(prefix + ".jsonp")
                """
                End of Per-KO stratified abundance
                """

    else: 
        """
        Destratified output processing
        """
        all_kos = all_cols

        for column in metadata_df.columns:
            metadata_df_filt = metadata_df[metadata_df[column].notna()]
            levels = metadata_df_filt[column].unique()
            concs = []
            for e, table in enumerate(tables):
                """
                Rename sample id based on the parameter
                """
                if split_str is not None:
                    table.index = [
                        i.split(quote(split_str))[0] for i in table.index.values
                    ]
                if convert_table is not None:
                    table.index = [
                        change[i] if i in change.keys() else i
                        for i in table.index.values
                    ]
                if not use_p:
                    tmp = pd.concat([table, metadata_df_filt], axis=1, join="inner")
                    if tables_name is not None:
                        tmp["category"] = tables_name[e]
                    else:
                        tmp["category"] = "data" + str(e)

                    concs.append(tmp)
                else:
                    ## We only need KO counts
                    concs.append(table)

            """
            If Wilcoxon-based correlation (Sun et al. 2020.),
            we ignore the `candidate`, `candidate_pathway`, and `first` option.
            """
            def signp(x, y):
                pv = scipy.stats.ranksums(x, y).pvalue
                return np.log10(pv) * np.sign(np.mean(x) - np.mean(y))

            if use_p:
                base = list(combinations(levels, 2))
                for pair in base:
                    level1 = pair[0]
                    level2 = pair[1]
                    prefix = level1 + " - " + level2
                    a = metadata_df_filt[
                        metadata_df_filt[column] == level1
                    ].index.tolist()
                    b = metadata_df_filt[
                        metadata_df_filt[column] == level2
                    ].index.tolist()
                    vals = []
                    for tbl in concs:
                        val = pd.Series(
                            tbl.columns.map(
                                lambda x: signp(tbl.loc[a, x], tbl.loc[b, x])
                            )
                        )
                        val.index = tbl.columns
                        vals.append(val)
                    kos = [val.index.values for val in vals]
                    all_kos = list(set.intersection(*map(set, kos)))
                    df = pd.DataFrame([val.loc[all_kos] for val in vals]).T
                    if tables_name is not None:
                        df.columns = tables_name
                    else:
                        df.columns = ["data" + str(e) for e, i in enumerate(tables)]
                    corr = df.corr(method=method)

                    prefix = prefix + "_corr"

                    plt.figure()
                    fig = sns.heatmap(corr, annot=True)
                    figs = fig.get_figure()
                    figs.savefig(path.join(output_dir, prefix + ".png"))
                    plt.close()

                    csv_path = os.path.join(output_dir, prefix + ".csv")
                    corr.to_csv(csv_path)

                    output_json(prefix, output_dir, corr)
                    filenames.append(prefix + ".jsonp")

            else:
                corrs = []
                for ko in all_kos:
                    prefix = column + "_" + ko
                    if ko not in concs[0].columns.values:
                        continue
                    conc = pd.concat(
                        [data.loc[:, [column, ko, "category"]] for data in concs],
                        axis=0,
                    )
                    if not skip:
                        g = sns.FacetGrid(conc, col="category")
                        g.map(sns.boxplot, column, ko, order=levels)
                        g.savefig(path.join(output_dir, prefix + ".png"))
                        plt.close(g.fig)

                        csv_path = os.path.join(output_dir, prefix + ".csv")

                        output = (
                            pd.DataFrame(conc.loc[:, [column, ko, "category"]])
                            .groupby("category")
                            .apply(lambda x: x.groupby(column).mean(ko))
                        )
                        if map_ko:
                            if ko in kodic.keys():
                                output["ko_description"] = kodic[ko]

                        output.to_csv(csv_path)

                        if ko.startswith("K"):
                            output["ko"] = (
                                '<a href="' + ko_url + ko + '">' + ko + "</a>"
                            )
                        else:
                            output["ko"] = ko

                        
                        """
                        Save the image
                        """
                        output_json(prefix, output_dir, output)
                        filenames.append(prefix + ".jsonp")

                    """
                    Correlation output per KO per dataset
                    """
                    prefix = prefix + "_corr"
                    if tbl_len > 1:
                        
                        """
                        If multiple tables, append scatterplot
                        This will enlarge the width of the figures
                        """
                        corrtbl = pd.concat(
                            [table.loc[all_samples, ko] for table in tables], axis=1
                        )
                        if tables_name is not None:
                            corrtbl.columns = tables_name
                        else:
                            corrtbl.columns = [
                                "data" + str(e) for e, i in enumerate(tables)
                            ]

                        corr = corrtbl.corr(method=method)

                        if not skip:
                            base = list(combinations(corrtbl.columns.values, 2))

                            fig, ax = plt.subplots(
                                1, 1 + len(base), figsize=(cor_fig_width, 4)
                            )
                            sns.heatmap(corr, annot=True, ax=ax[0])
                            for e, i in enumerate(base):
                                sns.scatterplot(corrtbl, x=i[0], y=i[1], ax=ax[e + 1])
                            plt.tight_layout()

                            figs = fig.get_figure()
                            figs.savefig(path.join(output_dir, prefix + ".png"))
                            plt.close(fig)
                    else:
                        ## Only heatmap
                        corrtbl = pd.concat(
                            [table.loc[all_samples, ko] for table in tables], axis=1
                        )
                        if tables_name is not None:
                            corrtbl.columns = tables_name
                        else:
                            corrtbl.columns = [
                                "data" + str(e) for e, i in enumerate(tables)
                            ]

                        corr = corrtbl.corr(method=method)
                        if not skip:
                            plt.figure()
                            fig = sns.heatmap(corr, annot=True)
                            figs = fig.get_figure()
                            figs.savefig(path.join(output_dir, prefix + ".png"))
                            plt.close(fig.fig)

                    if not skip:
                        csv_path = os.path.join(output_dir, prefix + ".csv")
                        corrtbl.to_csv(csv_path)

                        ## Save the image
                        output_json(prefix, output_dir, corr)
                        filenames.append(prefix + ".jsonp")

                    ## Keep the correlation values
                    corr = corr.where(np.triu(np.ones(corr.shape)).astype(bool))
                    corr = corr.stack().reset_index()
                    corr.columns = ["d1", "d2", "value"]
                    corr["ko"] = ko
                    corrs.append(corr)

                ## Correlation summarization by boxplot for every dataset pairs
                all_cor = pd.concat(corrs)
                csv_path = os.path.join(output_dir, "whole_corr.csv")
                all_cor.to_csv(csv_path)

                ## Correlation statistics
                corsum = all_cor.groupby("d1").apply(
                    lambda x: x.groupby("d2").agg(
                        {"value": ["mean", "median", "min", "max"]}
                    )
                )

                all_cor["label"] = all_cor.d1.map(str) + " - " + all_cor.d2
                plt.figure(figsize=(12, 10))
                g = sns.boxplot(all_cor, x="label", y="value")
                plt.xticks(rotation=45)
                plt.tight_layout()
                fig = g.get_figure()
                fig.savefig(path.join(output_dir, "whole_corr.png"))
                plt.close()

                prefix = "whole_corr"
                output_json(prefix, output_dir, corsum)
                filenames.append(prefix + ".jsonp")

                if cor_thresh is not None:
                    """
                    Subset to between dataset correlation
                    """
                    subset_cor = all_cor[all_cor.d1 != all_cor.d2]
                    kostat = subset_cor.groupby("ko").agg(
                        {"value": ["mean", "median", "min", "max"]}
                    )
                    kostat["ko"] = kostat.index.values
                    kostat["mean"] = kostat[("value", "mean")]

                    outp = kostat.loc[:, ["ko", "mean"]]
                    csv_path = os.path.join(
                        output_dir, "cor_thresh_" + str(cor_thresh) + ".csv"
                    )
                    outp.index.name = "feature-id"
                    outp = outp.reset_index()
                    outp.columns = ["feature-id","ko","mean"]
                    outp["thresholded"] = (outp["mean"] > cor_thresh).apply(lambda x: str(x))
                    outp.index = outp["feature-id"]
                    outp = outp.drop("feature-id", axis=1)

                    """
                    Output metadata that can be used in the `kegg` module
                    """
                    Metadata(outp).save(csv_path)

                    outp = outp[outp["mean"] > cor_thresh]

                    prefix = "cor_thresh_" + str(cor_thresh)

                    plt.figure(figsize=(12, 10))
                    g = sns.barplot(kostat, x="ko", y="mean")
                    plt.xticks(rotation=45)
                    plt.tight_layout()
                    fig = g.get_figure()
                    fig.savefig(path.join(output_dir, prefix + ".png"))
                    plt.close()

                    output_json(prefix, output_dir, outp)
                    filenames.append(prefix + ".jsonp")

    index = os.path.join(TEMPLATES, "index.html")
    q2templates.render(
        index, output_dir, context={"columns": [quote(fn) for fn in filenames]}
    )
    shutil.copytree(os.path.join(TEMPLATES, "dist"), os.path.join(output_dir, "dist"))


def contribute(
    output_dir: str,
    table: pd.DataFrame,
    candidate: str,
    metadata: qiime2.Metadata = None,
    fig_height: int = 16,
) -> None:
    filenames = []

    test = table.columns[0]
    if not test.startswith("K"):
        raise ValueError("No KO in the column")
    if "|" not in test:
        raise ValueError("Seems like not stratified output of `infer`")

    all_samples = [i for i in table.index.values]

    if metadata is not None:
        metadata = metadata.filter_ids(all_samples)
        metadata = metadata.filter_columns(column_type="categorical")
        metadata = metadata.filter_columns(
            drop_all_unique=False, drop_zero_variance=False, drop_all_missing=True
        )

        ## save out metadata for download in viz
        metadata.save(os.path.join(output_dir, "metadata.tsv"))

        ## Use in script
        metadata_df = metadata.to_dataframe()
    else:
        metadata_df = pd.DataFrame(index=all_samples)
        metadata_df["All"] = "All"
        metadata_df.to_csv(os.path.join(output_dir, "metadata.tsv"), sep="\t")

    for column in metadata_df.columns:
        metadata_df_filt = metadata_df[metadata_df[column].notna()]
        data = pd.concat([table, metadata_df_filt], axis=1, join="inner")

        prefix = column + "_" + candidate

        candidate_columns = [i for i in table.columns if candidate in i]
        if len(candidate_columns) == 0:
            raise ValueError("Candidate ID is not present in the specified table.")

        candidate_columns.append(column)
        output = data.loc[:, candidate_columns]
        output["sample"] = output.index.values
        output = pd.melt(output, id_vars=["sample", column])
        output["variable"] = output["variable"].apply(
            lambda x: "|".join(x.split("|")[1:])
        )

        plt.figure(figsize=(8, fig_height))
        bp = sns.boxplot(data=output, x="value", y="variable", hue=column)
        plt.tight_layout()
        fig = bp.get_figure()
        fig.savefig(path.join(output_dir, prefix + ".png"))
        plt.close()

        csv_path = os.path.join(output_dir, prefix + ".csv")

        outputstat = output.groupby(column).apply(
            lambda x: x.groupby("variable").agg(
                {"value": ["mean", "median", "min", "max"]}
            )
        )

        outputstat.to_csv(csv_path)

        output_json(prefix, output_dir, outputstat)
        filenames.append(prefix + ".jsonp")

        # Per-taxon sum across samples for the candidate KO
        if column == "All":
            prefix = column + "_" + candidate + "_sum"
            candidate_columns = [i for i in table.columns if candidate in i]
            sumtable = pd.DataFrame(
                table.loc[:, candidate_columns].apply(lambda x: sum(x))
            )
            sumtable.columns = ["Sum"]
            sumtable["Taxonomy"] = sumtable.index.values
            sumtable = sumtable.sort_values(by="Sum", ascending=False)

            plt.figure(figsize=(8, fig_height))
            bp = sns.barplot(sumtable, x="Sum", y="Taxonomy")
            plt.tight_layout()
            fig = bp.get_figure()
            fig.savefig(path.join(output_dir, prefix + ".png"))
            plt.close()

            csv_path = os.path.join(output_dir, prefix + ".csv")
            sumtable.to_csv(csv_path)

            output_json(prefix, output_dir, sumtable)
            filenames.append(prefix + ".jsonp")

    index = os.path.join(TEMPLATES, "index.html")
    q2templates.render(
        index, output_dir, context={"columns": [quote(fn) for fn in filenames]}
    )
    shutil.copytree(os.path.join(TEMPLATES, "dist"), os.path.join(output_dir, "dist"))



def scipy_cor(x, y, method):
    """
    Avoid internal alignment in Series.corr
    """
    if method=="pearson":
        return(scipy.stats.pearsonr(x,y)[0])
    elif method=="spearman":
        return(scipy.stats.spearmanr(x,y)[0])
    else:
        return(scipy.stats.kendalltau(x,y)[0])

def permute(
    output_dir: str,
    tables: pd.DataFrame,
    nperm: int = 100,
    convert_table: qiime2.Metadata = None,
    candidate_pathway: str = None,
    candidate: str = None,
    tss: bool = False,
    split_str: str = None,
    method: str = "spearman",
    tables_name: str = None,
) -> None:
    """
    Assess correlation based on permutation
    Accept only two tables, one is ground truth
    and the other is inference
    """
    if len(tables)!=2:
        raise ValueError("Two tables should be provided to this module")
    
    if convert_table is not None:
        """
        Converting table should be Metadata with the column name of "converted".
        If shotgun profiled qza is to be read, the index should be shotgun-ID and `converted` column
        corresponds to 16S-ID.
        """
        mapping = convert_table.to_dataframe()
        change = mapping["converted"].to_dict()

    if tss:
        tables = [ko_table.apply(lambda x: x / sum(x), axis=1) for ko_table in tables]

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

    filenames = []


    if candidate is None:
        all_cols = all_kos
    else:
        all_cols = [candidate]

    """
    Override `first` and `candidate` option if pathway is specified
    """
    if candidate_pathway is not None:
        kop = pd.read_csv("https://rest.kegg.jp/link/ko/pathway", sep="\t", header=None)
        kop = kop[kop[0].apply(lambda x: "path:ko" in x)]
        kop[0] = kop[0].apply(lambda x: x.split(":")[1])
        kop[1] = kop[1].apply(lambda x: x.split(":")[1])
        all_cols = kop[kop[0].isin([candidate_pathway])][1].tolist()
        all_cols = list(set(all_kos) & set(all_cols))
        
    concs = []
    for e, table in enumerate(tables):
        """
        Rename sample id based on the parameter
        """
        if split_str is not None:
            table.index = [
                i.split(quote(split_str))[0] for i in table.index.values
            ]
        if convert_table is not None:
            table.index = [
                change[i] if i in change.keys() else i
                for i in table.index.values
            ]

        if tables_name is not None:
            table["category"] = tables_name[e]
        else:
            table["category"] = "data" + str(e)
    
        concs.append(table)
    
    permuted = []
    unpermuted = []
    for ko in all_cols:
        x = concs[0].loc[all_samples, ko].tolist()
        y = concs[1].loc[all_samples, ko]
        unpermuted.append([scipy_cor(x, y.tolist(), method=method), ko])
        for i in range(nperm):
            yp = y.sample(frac=1, random_state=i).tolist()
            permuted.append([scipy_cor(x, yp, method=method), ko])
            
    perm_df = pd.DataFrame(permuted)
    unperm_df = pd.DataFrame(unpermuted)
    perm_df["group"] = "Permuted"
    unperm_df["group"] = "Unpermuted"
    output = pd.concat([perm_df, unperm_df])
    output.columns = ["value","ko","group"]
    
    """
    Plot and save
    """
    prefix = "permute"
    print(output)
    
    plt.figure()
    bp = sns.boxplot(data=output, x="group", y="value")
    fig = bp.get_figure()
    fig.savefig(path.join(output_dir, prefix + ".png"))
    plt.close()

    csv_path = os.path.join(output_dir, prefix + ".csv")
    output.to_csv(csv_path)
    output = output.groupby("group").apply(
        lambda x: x.groupby("ko").mean("value")
    )



    output_json(prefix, output_dir, output)
    filenames.append(prefix + ".jsonp")

    index = os.path.join(TEMPLATES, "index.html")
    q2templates.render(
        index, output_dir, context={"columns": [quote(fn) for fn in filenames]}
    )
    shutil.copytree(os.path.join(TEMPLATES, "dist"), os.path.join(output_dir, "dist"))
