from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import pandas as pd
from scipy.stats import ttest_ind


def validate_table_annotations(
    n_tables: int,
    table_names: Sequence[str] | None,
    table_roles: Sequence[str] | None = None,
    algorithm_families: Sequence[str] | None = None,
    *,
    min_predictions: int = 2,
) -> tuple[list[str], list[str], list[str]]:

    names = (
        [f"data{i}" for i in range(n_tables)]
        if table_names is None
        else [str(x).strip() for x in table_names]
    )
    if len(names) != n_tables:
        raise ValueError("tables_name must have one value per input table.")
    if any(not name for name in names):
        raise ValueError("Table names must be non-empty.")
    if len(set(names)) != len(names):
        raise ValueError("Table names must be unique.")

    roles = (
        ["prediction"] * n_tables
        if table_roles is None
        else [str(x).strip().lower() for x in table_roles]
    )
    if len(roles) != n_tables:
        raise ValueError("table_roles must have one value per input table.")
    unknown_roles = sorted(set(roles) - {"prediction", "reference"})
    if unknown_roles:
        raise ValueError(
            "table_roles values must be 'prediction' or 'reference'; got "
            + ", ".join(unknown_roles)
        )
    if roles.count("prediction") < min_predictions:
        raise ValueError(
            f"At least {min_predictions} prediction table(s) are required."
        )

    families = (
        names.copy()
        if algorithm_families is None
        else [str(x).strip() for x in algorithm_families]
    )
    if len(families) != n_tables:
        raise ValueError("algorithm_families must have one value per input table.")
    if any(not family for family in families):
        raise ValueError("Algorithm-family names must be non-empty.")

    return names, roles, families


def harmonize_tables(
    tables: Sequence[pd.DataFrame],
) -> tuple[list[pd.DataFrame], list[str], list[str], pd.DataFrame]:

    if len(tables) < 2:
        raise ValueError("At least two KO tables are required.")

    normalized: list[pd.DataFrame] = []
    coverage_records: list[dict[str, int]] = []
    for i, table in enumerate(tables):
        if not isinstance(table, pd.DataFrame):
            raise TypeError("Every KO table must be a pandas DataFrame.")
        frame = table.copy()
        frame.index = frame.index.astype(str)
        frame.columns = frame.columns.astype(str)
        if frame.index.has_duplicates:
            dup = frame.index[frame.index.duplicated()].unique().tolist()[:5]
            raise ValueError(f"Duplicate sample IDs in table {i}: {dup}")
        if frame.columns.has_duplicates:
            dup = frame.columns[frame.columns.duplicated()].unique().tolist()[:5]
            raise ValueError(f"Duplicate KO IDs in table {i}: {dup}")
        values = frame.to_numpy(dtype=float)
        if not np.isfinite(values).all():
            raise ValueError(f"Non-finite abundance in table {i}.")
        if (values < 0).any():
            raise ValueError(f"Negative abundance in table {i}.")
        coverage_records.append(
            {
                "table_index": i,
                "samples_before": frame.shape[0],
                "kos_before": frame.shape[1],
                "samples_nonzero": int((frame.sum(axis=1) > 0).sum()),
            }
        )
        normalized.append(frame.loc[frame.sum(axis=1) > 0].copy())

    common_samples = sorted(
        set.intersection(*(set(frame.index) for frame in normalized))
    )
    common_kos = sorted(set.intersection(*(set(frame.columns) for frame in normalized)))
    if not common_samples:
        raise ValueError("No common samples are present across datasets.")
    if not common_kos:
        raise ValueError("No common KOs are present across datasets.")

    harmonized = [frame.loc[common_samples, common_kos].copy() for frame in normalized]
    for record in coverage_records:
        record["samples_common"] = len(common_samples)
        record["kos_common"] = len(common_kos)

    return harmonized, common_samples, common_kos, pd.DataFrame(coverage_records)


def harmonize_samples_for_visualization(
    tables: Sequence[pd.DataFrame],
    metadata: pd.DataFrame,
    table_names: Sequence[str] | None = None,
) -> tuple[list[pd.DataFrame], pd.DataFrame, pd.DataFrame]:

    if not tables:
        raise ValueError("At least one KO table is required.")
    names = (
        [f"data{i}" for i in range(len(tables))]
        if table_names is None
        else [str(name).strip() for name in table_names]
    )
    if len(names) != len(tables):
        raise ValueError("table_names must have one value per input table.")
    if any(not name for name in names) or len(set(names)) != len(names):
        raise ValueError("Table names must be non-empty and unique.")

    meta = metadata.copy()
    meta.index = meta.index.astype(str)
    if meta.index.has_duplicates:
        raise ValueError("Metadata sample IDs must be unique.")

    normalized: list[pd.DataFrame] = []
    records: list[dict[str, int | str]] = []
    for name, table in zip(names, tables):
        if not isinstance(table, pd.DataFrame):
            raise TypeError("Every KO table must be a pandas DataFrame.")
        frame = table.copy()
        frame.index = frame.index.astype(str)
        frame.columns = frame.columns.astype(str)
        if frame.index.has_duplicates:
            raise ValueError(f"Duplicate sample IDs in {name}.")
        if frame.columns.has_duplicates:
            raise ValueError(f"Duplicate KO IDs in {name}.")
        values = frame.to_numpy(dtype=float)
        if not np.isfinite(values).all():
            raise ValueError(f"Non-finite abundance in {name}.")
        if (values < 0).any():
            raise ValueError(f"Negative abundance in {name}.")
        normalized.append(frame)
        records.append(
            {
                "table": name,
                "samples_before": frame.shape[0],
                "kos_before": frame.shape[1],
                "zero_total_samples_before": int((frame.sum(axis=1) <= 0).sum()),
            }
        )

    common_samples = sorted(
        set(meta.index).intersection(
            *[set(frame.index) for frame in normalized]
        )
    )
    if not common_samples:
        raise ValueError("No samples are shared by all KO tables and metadata.")

    retained_samples = [
        sample
        for sample in common_samples
        if all(frame.loc[sample].sum() > 0 for frame in normalized)
    ]
    if not retained_samples:
        raise ValueError("No non-zero samples remain across all KO tables.")

    harmonized = [frame.loc[retained_samples].copy() for frame in normalized]
    harmonized_meta = meta.loc[retained_samples].copy()
    for record, frame in zip(records, harmonized):
        record["samples_shared_before_zero_filter"] = len(common_samples)
        record["samples_retained"] = len(retained_samples)
        record["kos_retained"] = frame.shape[1]

    return harmonized, harmonized_meta, pd.DataFrame(records)


def welch_t_statistics(
    table: pd.DataFrame,
    groups: pd.Series,
    level1: str,
    level2: str,
) -> pd.Series:

    aligned_groups = groups.reindex(table.index)
    first = aligned_groups[aligned_groups.astype(str) == str(level1)].index
    second = aligned_groups[aligned_groups.astype(str) == str(level2)].index
    if len(first) < 2 or len(second) < 2:
        raise ValueError(
            "Welch statistics require at least two samples in each selected group."
        )
    values = [
        ttest_ind(
            table.loc[first, feature].to_numpy(dtype=float),
            table.loc[second, feature].to_numpy(dtype=float),
            equal_var=False,
            nan_policy="omit",
        ).statistic
        for feature in table.columns
    ]
    return pd.Series(values, index=table.columns.astype(str), dtype=float)


def summarize_direction_statistics(
    statistics: pd.DataFrame,
    direction_threshold: float = 1.0,
    min_methods: int = 2,
) -> pd.DataFrame:

    if not 0.5 <= direction_threshold <= 1.0:
        raise ValueError("direction_threshold must be between 0.5 and 1.0.")
    if min_methods < 1:
        raise ValueError("min_methods must be at least 1.")
    values = statistics.apply(pd.to_numeric, errors="coerce")
    positive = (values > 0).sum(axis=1)
    negative = (values < 0).sum(axis=1)
    directional = positive + negative
    available = values.notna().sum(axis=1)
    concordance = pd.Series(np.nan, index=values.index, dtype=float)
    has_direction = directional > 0
    concordance.loc[has_direction] = (
        pd.concat([positive, negative], axis=1).max(axis=1).loc[has_direction]
        / directional.loc[has_direction]
    )
    dominant = pd.Series("none", index=values.index, dtype=object)
    dominant.loc[positive > negative] = "positive"
    dominant.loc[negative > positive] = "negative"
    dominant.loc[(positive == negative) & has_direction] = "tie"

    output = values.copy()
    output["n_methods_available"] = available
    output["n_directional_methods"] = directional
    output["n_positive"] = positive
    output["n_negative"] = negative
    output["dominant_direction"] = dominant
    output["direction_concordance"] = concordance
    output["direction_consistent"] = (
        (directional >= min_methods) & (concordance >= direction_threshold)
    )
    return output


def closure_multiplicative_clr(
    table: pd.DataFrame,
    zero_replacement: float = 1e-6,
) -> pd.DataFrame:
    """
    Calculate compositions, replace zeros, and apply CLR.

    """
    if not 0 < zero_replacement < 1:
        raise ValueError("zero_replacement must be between 0 and 1.")

    values = table.to_numpy(dtype=float, copy=True)
    if not np.isfinite(values).all():
        raise ValueError("KO table contains non-finite values.")
    if (values < 0).any():
        raise ValueError("KO table contains negative values.")

    row_sums = values.sum(axis=1)
    if (row_sums <= 0).any():
        bad = table.index[row_sums <= 0].astype(str).tolist()[:5]
        raise ValueError(f"Samples with zero total abundance: {bad}")
    compositions = values / row_sums[:, None]

    for i in range(compositions.shape[0]):
        row = compositions[i]
        zero_mask = row == 0
        n_zeros = int(zero_mask.sum())
        if n_zeros == 0:
            continue
        positive = row[~zero_mask]
        if positive.size == 0:
            raise ValueError(f"Sample {table.index[i]} has no positive KOs.")
        delta = min(
            zero_replacement,
            float(positive.min()) / 2.0,
            0.5 / n_zeros,
        )
        zero_mass = delta * n_zeros
        row[~zero_mask] *= 1.0 - zero_mass
        row[zero_mask] = delta
        compositions[i] = row / row.sum()

    logged = np.log(compositions)
    clr = logged - logged.mean(axis=1, keepdims=True)
    return pd.DataFrame(clr, index=table.index, columns=table.columns)


def pairwise_spearman_diagnostics(
    clr_tables: Sequence[pd.DataFrame],
    names: Sequence[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:

    if len(clr_tables) != len(names):
        raise ValueError("names must have one value per CLR table.")
    per_sample: list[dict[str, object]] = []
    per_feature: list[dict[str, object]] = []
    for i in range(len(clr_tables)):
        for j in range(i + 1, len(clr_tables)):
            left = clr_tables[i]
            right = clr_tables[j]
            for sample_id in left.index:
                rho = left.loc[sample_id].corr(right.loc[sample_id], method="spearman")
                per_sample.append(
                    {
                        "sample": sample_id,
                        "rho": rho,
                        "method1": names[i],
                        "method2": names[j],
                    }
                )
            for ko_id in left.columns:
                rho = left[ko_id].corr(right[ko_id], method="spearman")
                per_feature.append(
                    {
                        "KO": ko_id,
                        "rho": rho,
                        "method1": names[i],
                        "method2": names[j],
                    }
                )
    return (
        pd.DataFrame(per_sample, columns=["sample", "rho", "method1", "method2"]),
        pd.DataFrame(per_feature, columns=["KO", "rho", "method1", "method2"]),
    )


def family_balanced_diagnostics(
    clr_tables: Sequence[pd.DataFrame],
    names: Sequence[str],
    roles: Sequence[str],
    families: Sequence[str],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:

    lengths = {len(clr_tables), len(names), len(roles), len(families)}
    if len(lengths) != 1:
        raise ValueError("Tables, names, roles, and families must have equal length.")

    grouped: dict[tuple[str, str], list[pd.DataFrame]] = {}
    for table, role, family in zip(clr_tables, roles, families):
        grouped.setdefault((role, family), []).append(table)

    family_tables = []
    family_names = []
    annotations = []
    for (role, family), members in grouped.items():
        stacked = np.stack([member.to_numpy(dtype=float) for member in members])
        collapsed = pd.DataFrame(
            np.median(stacked, axis=0),
            index=members[0].index,
            columns=members[0].columns,
        )
        label = f"{role}:{family}"
        family_tables.append(collapsed)
        family_names.append(label)
        annotations.append(
            {
                "method": label,
                "role": role,
                "algorithm_family": family,
                "n_outputs_collapsed": len(members),
            }
        )

    pair_sample, pair_feature = pairwise_spearman_diagnostics(
        family_tables, family_names
    )
    annotation = pd.DataFrame(annotations).set_index("method")
    for frame in (pair_sample, pair_feature):
        for suffix in ("1", "2"):
            frame[f"role{suffix}"] = frame[f"method{suffix}"].map(annotation["role"])
            frame[f"algorithm_family{suffix}"] = frame[f"method{suffix}"].map(
                annotation["algorithm_family"]
            )
            frame[f"n_outputs_collapsed{suffix}"] = frame[f"method{suffix}"].map(
                annotation["n_outputs_collapsed"]
            )

    def summarize(frame: pd.DataFrame, item: str) -> pd.DataFrame:
        prediction_pairs = frame[
            (frame["role1"] == "prediction") & (frame["role2"] == "prediction")
        ]
        heldout_pairs = frame[
            (frame["role1"] != frame["role2"])
            & frame[["role1", "role2"]].isin(["prediction", "reference"]).all(axis=1)
        ]
        source_items = (
            family_tables[0].index if item == "sample" else family_tables[0].columns
        )
        items = pd.Index(source_items, name=item)
        output = pd.DataFrame(index=items)
        output["prediction_family_concordance"] = prediction_pairs.groupby(item)[
            "rho"
        ].median()
        output["heldout_reference_agreement"] = heldout_pairs.groupby(item)[
            "rho"
        ].median()
        output["n_prediction_families"] = len(
            {family for role, family in grouped if role == "prediction"}
        )
        output["n_reference_families"] = len(
            {family for role, family in grouped if role == "reference"}
        )
        return output.reset_index()

    return (
        pair_sample,
        pair_feature,
        summarize(pair_sample, "sample"),
        summarize(pair_feature, "KO"),
    )
