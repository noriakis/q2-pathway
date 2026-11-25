import os
import subprocess
from importlib import resources
from typing import List

import pandas as pd
import qiime2

from ._stats import (
    closure_multiplicative_clr,
    family_balanced_diagnostics,
    harmonize_tables,
    pairwise_spearman_diagnostics,
    validate_table_annotations,
)


def summarize(
    output_dir: str,
    tables: pd.DataFrame,
    metadata: qiime2.Metadata = None,
    convert_table: qiime2.Metadata = None,
    tss: bool = False,
    split_str: str = None,
    tables_name: List[str] = None,
    table_roles: List[str] = None,
    algorithm_families: List[str] = None,
    zero_replacement: float = 1e-6,
) -> None:

    os.makedirs(output_dir, exist_ok=True)

    if isinstance(tables, pd.DataFrame):
        tables = [tables]
    elif not isinstance(tables, (list, tuple)):
        raise TypeError("`tables` must be a DataFrame or a list of DataFrames.")
    tables = list(tables)
    names, roles, families = validate_table_annotations(
        len(tables), tables_name, table_roles, algorithm_families
    )

    mapping = None
    if convert_table is not None:
        mapping_df = convert_table.to_dataframe()
        if "converted" not in mapping_df.columns:
            raise ValueError("convert_table must have a 'converted' column.")
        if mapping_df.index.has_duplicates:
            raise ValueError("convert_table contains duplicate source sample IDs.")
        mapping = mapping_df["converted"].astype(str).to_dict()

    processed_tables = []
    for table in tables:
        frame = table.copy()
        frame.index = frame.index.astype(str)
        frame.columns = frame.columns.astype(str)
        if split_str is not None:
            frame.index = [sample.split(split_str)[0] for sample in frame.index]
        if mapping is not None:
            frame.index = [mapping.get(sample, sample) for sample in frame.index]
        if frame.index.has_duplicates:
            duplicates = frame.index[frame.index.duplicated()].unique().tolist()[:5]
            raise ValueError(
                "Sample conversion produced duplicate sample IDs: "
                + ", ".join(duplicates)
            )
        processed_tables.append(frame)

    harmonized, common_samples, common_kos, coverage = harmonize_tables(
        processed_tables
    )
    coverage.insert(1, "dataset", names)
    coverage.insert(2, "role", roles)
    coverage.insert(3, "algorithm_family", families)
    coverage.to_csv(
        os.path.join(output_dir, "catalog_coverage.tsv"), sep="\t", index=False
    )
    pd.DataFrame(
        {
            "dataset": names,
            "role": roles,
            "algorithm_family": families,
            "tss_requested": bool(tss),
            "zero_replacement": zero_replacement,
        }
    ).to_csv(os.path.join(output_dir, "table_manifest.tsv"), sep="\t", index=False)

    clr_tables = [
        closure_multiplicative_clr(table, zero_replacement) for table in harmonized
    ]
    per_sample, per_feature = pairwise_spearman_diagnostics(clr_tables, names)
    annotations = pd.DataFrame(
        {
            "dataset": names,
            "role": roles,
            "algorithm_family": families,
        }
    )
    per_sample = per_sample.merge(
        annotations.add_suffix("1"),
        left_on="method1",
        right_on="dataset1",
        how="left",
    ).merge(
        annotations.add_suffix("2"),
        left_on="method2",
        right_on="dataset2",
        how="left",
    )
    per_feature = per_feature.merge(
        annotations.add_suffix("1"),
        left_on="method1",
        right_on="dataset1",
        how="left",
    ).merge(
        annotations.add_suffix("2"),
        left_on="method2",
        right_on="dataset2",
        how="left",
    )
    per_sample.to_csv(
        os.path.join(output_dir, "abundance_agreement_per_sample.tsv"),
        sep="\t",
        index=False,
    )
    per_feature.to_csv(
        os.path.join(output_dir, "abundance_agreement_per_feature.tsv"),
        sep="\t",
        index=False,
    )
    (
        family_per_sample,
        family_per_feature,
        family_summary_sample,
        family_summary_feature,
    ) = family_balanced_diagnostics(clr_tables, names, roles, families)
    family_per_sample.to_csv(
        os.path.join(output_dir, "family_agreement_per_sample.tsv"),
        sep="\t",
        index=False,
    )
    family_per_feature.to_csv(
        os.path.join(output_dir, "family_agreement_per_feature.tsv"),
        sep="\t",
        index=False,
    )
    family_summary_sample.to_csv(
        os.path.join(output_dir, "family_balanced_summary_per_sample.tsv"),
        sep="\t",
        index=False,
    )
    family_summary_feature.to_csv(
        os.path.join(output_dir, "family_balanced_summary_per_feature.tsv"),
        sep="\t",
        index=False,
    )

    metadata_written = False
    if metadata is not None:
        md = metadata.filter_ids(common_samples)
        md.to_dataframe().to_csv(
            os.path.join(output_dir, "metadata_for_R.tsv"), sep="\t"
        )
        metadata_written = True

    script_path = (
        resources.files("q2_pathway").joinpath("assets").joinpath("summarize_ko.R")
    )
    subprocess.run(["Rscript", str(script_path), output_dir, "KO-summary"], check=True)

    links = [
        ("KO_agreement_per_sample.pdf", "KO agreement per sample (PDF)"),
        ("KO_agreement_per_feature.pdf", "KO agreement per feature (PDF)"),
        ("abundance_agreement_per_sample.tsv", "Per-sample correlations (TSV)"),
        ("abundance_agreement_per_feature.tsv", "Per-feature correlations (TSV)"),
        ("catalog_coverage.tsv", "Catalog and intersection sizes (TSV)"),
        ("table_manifest.tsv", "Table roles and algorithm families (TSV)"),
        (
            "family_agreement_per_sample.tsv",
            "Family-balanced pairwise sample correlations (TSV)",
        ),
        (
            "family_agreement_per_feature.tsv",
            "Family-balanced pairwise KO correlations (TSV)",
        ),
        (
            "family_balanced_summary_per_sample.tsv",
            "Family-balanced sample summary (TSV)",
        ),
        (
            "family_balanced_summary_per_feature.tsv",
            "Family-balanced KO summary (TSV)",
        ),
    ]
    if metadata_written:
        links.append(("metadata_for_R.tsv", "Metadata used for diagnostics (TSV)"))
    link_html = "\n".join(
        f'      <li><a href="{path}">{label}</a></li>' for path, label in links
    )
    index_html = f"""<!DOCTYPE html>
<html>
  <head><meta charset="utf-8"><title>q2-pathway KO summary</title></head>
  <body>
    <h1>q2-pathway KO-level diagnostics</h1>
    <p>Tables were restricted to shared samples and KOs, closed to relative
    abundance, zero-replaced multiplicatively, and CLR-transformed. Missing
    catalog entries were not imputed as biological zeros.</p>
    <p>{len(common_samples)} common samples and {len(common_kos)} common KOs.</p>
    <ul>
{link_html}
    </ul>
  </body>
</html>
"""
    with open(os.path.join(output_dir, "index.html"), "w", encoding="utf-8") as fh:
        fh.write(index_html)
