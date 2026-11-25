from __future__ import annotations

from collections import defaultdict
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


@dataclass(frozen=True)
class DatasetAnnotation:
    role: str = "prediction"
    family: str | None = None


def adjust_pvalues(pvalues: pd.Series, method: str = "bh") -> pd.Series:
    method_key = method.strip().lower().replace("_", "-")
    result = pd.Series(np.nan, index=pvalues.index, dtype=float)
    valid = pvalues.dropna().astype(float).clip(0, 1)
    if valid.empty:
        return result
    if method_key in {"bonferroni", "bonf"}:
        result.loc[valid.index] = np.minimum(valid * len(valid), 1.0)
        return result
    if method_key not in {"bh", "fdr-bh", "benjamini-hochberg"}:
        raise ValueError("p_adjust_method must be 'bh' or 'bonferroni'.")

    ordered = valid.sort_values()
    ranks = np.arange(1, len(ordered) + 1, dtype=float)
    raw = ordered.to_numpy() * len(ordered) / ranks
    corrected = np.minimum.accumulate(raw[::-1])[::-1]
    result.loc[ordered.index] = np.minimum(corrected, 1.0)
    return result


def cauchy_combine(pvalues: list[float] | np.ndarray) -> float:
    p = np.asarray(pvalues, dtype=float)
    p = p[np.isfinite(p)]
    if p.size == 0:
        return np.nan
    p = np.clip(p, 1e-15, 1 - 1e-15)
    statistic = np.mean(np.tan((0.5 - p) * np.pi))
    combined = 0.5 - np.arctan(statistic) / np.pi
    return float(np.clip(combined, 0.0, 1.0))


def fisher_combine(pvalues: list[float] | np.ndarray) -> float:
    p = np.asarray(pvalues, dtype=float)
    p = p[np.isfinite(p)]
    if p.size == 0:
        return np.nan
    p = np.clip(p, 1e-15, 1.0)
    statistic = -2.0 * np.log(p).sum()
    return float(stats.chi2.sf(statistic, 2 * p.size))


def _standardize_gsea(frame: pd.DataFrame) -> pd.DataFrame:
    renamed = {str(col).lower(): col for col in frame.columns}
    pathway_col = renamed.get("pathway", renamed.get("term"))
    p_col = renamed.get("pval", renamed.get("pvalue"))
    q_col = renamed.get("padj", renamed.get("fdr", p_col))
    nes_col = renamed.get("nes", renamed.get("enrichment_score"))
    if pathway_col is None or p_col is None or nes_col is None:
        raise ValueError("GSEA table must contain pathway, p-value, and NES columns.")
    result = pd.DataFrame(
        {
            "pathway": frame[pathway_col].astype(str).str.strip(),
            "pval": pd.to_numeric(frame[p_col], errors="coerce"),
            "padj": pd.to_numeric(frame[q_col], errors="coerce"),
            "nes": pd.to_numeric(frame[nes_col], errors="coerce"),
        }
    )
    result = result.drop_duplicates("pathway", keep="first").set_index("pathway")
    return result


def _majority_sign(values: list[float]) -> int:
    signs = np.sign(np.asarray(values, dtype=float))
    signs = signs[np.isfinite(signs) & (signs != 0)]
    if signs.size == 0:
        return 0
    positives = int((signs > 0).sum())
    negatives = int((signs < 0).sum())
    if positives == negatives:
        return 0
    return 1 if positives > negatives else -1


def build_consensus(
    results: Mapping[str, pd.DataFrame],
    annotations: Mapping[str, DatasetAnnotation],
    *,
    top_set_alpha: float = 0.05,
    consensus_alpha: float = 0.05,
    direction_threshold: float = 0.8,
    min_families: int = 2,
    p_adjust_method: str = "bh",
) -> tuple[pd.DataFrame, pd.DataFrame]:

    if not 0 < top_set_alpha < 1 or not 0 < consensus_alpha < 1:
        raise ValueError("Alpha values must be between 0 and 1.")
    if not 0.5 <= direction_threshold <= 1:
        raise ValueError("direction_threshold must be between 0.5 and 1.")
    if min_families < 1:
        raise ValueError("min_families must be positive.")
    if set(results) != set(annotations):
        raise ValueError("Every result must have exactly one annotation.")

    standardized = {name: _standardize_gsea(df) for name, df in results.items()}
    prediction_names = [
        name for name, ann in annotations.items() if ann.role.lower() == "prediction"
    ]
    reference_names = [
        name for name, ann in annotations.items() if ann.role.lower() == "reference"
    ]
    if len(prediction_names) < 2:
        raise ValueError("At least two prediction outputs are required.")

    family_members: dict[str, list[str]] = defaultdict(list)
    for name in prediction_names:
        family = annotations[name].family or name
        family_members[family].append(name)
    if len(family_members) < min_families:
        raise ValueError("Fewer prediction families than min_families.")

    pathways = sorted(
        set().union(*(set(standardized[name].index) for name in prediction_names))
    )
    family_detail: list[dict[str, object]] = []

    for pathway in pathways:
        for family, members in family_members.items():
            member_p: list[float] = []
            member_nes: list[float] = []
            for name in members:
                table = standardized[name]
                if pathway not in table.index:
                    continue
                row = table.loc[pathway]
                if np.isfinite(row["pval"]):
                    member_p.append(float(row["pval"]))
                if np.isfinite(row["nes"]):
                    member_nes.append(float(row["nes"]))

            family_p = cauchy_combine(member_p)
            family_sign = _majority_sign(member_nes)
            family_nes = float(np.median(member_nes)) if member_nes else np.nan
            family_detail.append(
                {
                    "pathway": pathway,
                    "family": family,
                    "members": "|".join(members),
                    "family_p": family_p,
                    "family_nes": family_nes,
                    "family_sign": family_sign,
                    "n_family_members_present": len(member_p),
                }
            )

    family_detail_frame = pd.DataFrame(family_detail)
    family_detail_frame["family_q"] = family_detail_frame.groupby("family", sort=False)[
        "family_p"
    ].transform(lambda values: adjust_pvalues(values, p_adjust_method))
    family_detail_frame["family_significant"] = (
        family_detail_frame["family_q"] < top_set_alpha
    ) & (family_detail_frame["family_sign"] != 0)

    records: list[dict[str, object]] = []
    for pathway in pathways:
        detail = family_detail_frame[family_detail_frame["pathway"] == pathway]
        finite_p = detail["family_p"].notna()
        family_pvalues = detail.loc[finite_p, "family_p"].astype(float).tolist()
        family_signs = (
            detail.loc[detail["family_sign"] != 0, "family_sign"].astype(int).tolist()
        )
        family_nes_values = (
            detail.loc[detail["family_nes"].notna(), "family_nes"]
            .astype(float)
            .tolist()
        )
        significant_family_signs = (
            detail.loc[detail["family_significant"], "family_sign"].astype(int).tolist()
        )
        output_pvalues = [
            float(standardized[name].at[pathway, "pval"])
            for name in prediction_names
            if pathway in standardized[name].index
            and np.isfinite(standardized[name].at[pathway, "pval"])
        ]

        majority_sign = _majority_sign(family_signs)
        direction_concordance = (
            float(np.mean(np.asarray(family_signs) == majority_sign))
            if family_signs and majority_sign != 0
            else np.nan
        )
        significant_majority_sign = _majority_sign(significant_family_signs)
        n_sig_same_direction = (
            sum(sign == significant_majority_sign for sign in significant_family_signs)
            if significant_majority_sign != 0
            else 0
        )
        n_families = len(family_pvalues)
        n_direction_families = len(family_signs)
        strict_intersection = bool(
            n_families >= min_families
            and len(significant_family_signs) == n_families
            and len(set(significant_family_signs)) == 1
        )
        majority_vote = bool(
            n_families >= min_families
            and n_sig_same_direction >= (n_families // 2 + 1)
        )
        records.append(
            {
                "pathway": pathway,
                "n_prediction_outputs": len(output_pvalues),
                "n_families": n_families,
                "n_direction_families": n_direction_families,
                "family_cauchy_p": cauchy_combine(family_pvalues),
                "naive_fisher_p": fisher_combine(output_pvalues),
                "consensus_nes": (
                    float(np.median(family_nes_values)) if family_nes_values else np.nan
                ),
                "majority_sign": majority_sign,
                "direction_concordance": direction_concordance,
                "strict_intersection": strict_intersection,
                "majority_vote": majority_vote,
                "prediction_outputs": "|".join(prediction_names),
                "reference_outputs_excluded": "|".join(reference_names),
            }
        )

    consensus = pd.DataFrame(records).set_index("pathway")
    consensus["family_cauchy_q"] = adjust_pvalues(
        consensus["family_cauchy_p"], p_adjust_method
    )
    consensus["naive_fisher_q"] = adjust_pvalues(
        consensus["naive_fisher_p"], p_adjust_method
    )
    consensus["reportable"] = (
        (consensus["family_cauchy_q"] < consensus_alpha)
        & (consensus["n_families"] >= min_families)
        & (consensus["n_direction_families"] >= min_families)
        & (consensus["direction_concordance"] >= direction_threshold)
    )
    return consensus.reset_index(), family_detail_frame


def evaluate_against_references(
    results: Mapping[str, pd.DataFrame],
    annotations: Mapping[str, DatasetAnnotation],
    consensus: pd.DataFrame,
    *,
    alpha: float = 0.05,
) -> pd.DataFrame:

    standardized = {name: _standardize_gsea(df) for name, df in results.items()}
    reference_names = [
        name for name, ann in annotations.items() if ann.role.lower() == "reference"
    ]
    prediction_names = [
        name for name, ann in annotations.items() if ann.role.lower() == "prediction"
    ]
    if not reference_names:
        return pd.DataFrame()

    cons = consensus.set_index("pathway")
    candidates: dict[str, pd.DataFrame] = {
        name: standardized[name].rename(
            columns={"padj": "q", "nes": "effect", "pval": "p"}
        )[["p", "q", "effect"]]
        for name in prediction_names
    }
    family_cauchy = cons.rename(
        columns={
            "family_cauchy_p": "p",
            "family_cauchy_q": "q",
            "consensus_nes": "effect",
        }
    )[["p", "q", "effect"]]
    candidates["family_cauchy_q_only"] = family_cauchy.copy()
    family_cauchy.loc[~cons["reportable"].astype(bool), "q"] = 1.0
    candidates["family_cauchy"] = family_cauchy
    candidates["naive_fisher"] = cons.rename(
        columns={
            "naive_fisher_p": "p",
            "naive_fisher_q": "q",
            "consensus_nes": "effect",
        }
    )[["p", "q", "effect"]]
    any_method_significant = pd.concat(
        [
            (standardized[name]["padj"].reindex(cons.index) < alpha).rename(name)
            for name in prediction_names
        ],
        axis=1,
    ).any(axis=1)
    any_method_union = cons[["consensus_nes"]].rename(
        columns={"consensus_nes": "effect"}
    )
    any_method_union["p"] = np.where(any_method_significant, 0.0, 1.0)
    any_method_union["q"] = any_method_union["p"]
    candidates["any_method_union"] = any_method_union[["p", "q", "effect"]]
    for strategy, flag in (
        ("strict_intersection", "strict_intersection"),
        ("majority_vote", "majority_vote"),
    ):
        tmp = cons[["consensus_nes", flag]].rename(columns={"consensus_nes": "effect"})
        tmp["p"] = np.where(tmp[flag], 0.0, 1.0)
        tmp["q"] = tmp["p"]
        candidates[strategy] = tmp[["p", "q", "effect"]]

    selection_rules = {
        **{name: "method adjusted p < alpha" for name in prediction_names},
        "family_cauchy": (
            "family Cauchy adjusted p < alpha, minimum family count, and "
            "direction-concordance threshold"
        ),
        "family_cauchy_q_only": "family Cauchy adjusted p < alpha only",
        "naive_fisher": "naive Fisher adjusted p < alpha (sensitivity only)",
        "any_method_union": (
            "adjusted p < alpha in at least one prediction output; direction from "
            "the family-balanced consensus effect"
        ),
        "strict_intersection": "all algorithm families significant in one direction",
        "majority_vote": "majority of algorithm families significant in one direction",
    }

    records: list[dict[str, object]] = []
    for reference_name in reference_names:
        reference = standardized[reference_name].rename(
            columns={"pval": "p_ref", "padj": "q_ref", "nes": "effect_ref"}
        )[["p_ref", "q_ref", "effect_ref"]]
        common_pathways = set(reference.index[reference["effect_ref"].notna()])
        for prediction_name in prediction_names:
            prediction = standardized[prediction_name]
            common_pathways &= set(prediction.index[prediction["nes"].notna()])
        common_index = pd.Index(sorted(common_pathways), name="pathway")
        if common_index.empty:
            continue
        scoped_reference = reference.reindex(common_index)
        for strategy, candidate in candidates.items():
            joined = (
                candidate.reindex(common_index)
                .join(scoped_reference, how="inner")
                .dropna(subset=["effect", "effect_ref"])
            )
            if joined.empty:
                continue
            signed_candidate = np.sign(joined["effect"]) * -np.log10(
                joined["p"].clip(1e-300, 1)
            )
            signed_reference = np.sign(joined["effect_ref"]) * -np.log10(
                joined["p_ref"].clip(1e-300, 1)
            )
            predicted = set(joined.index[joined["q"] < alpha])
            observed = set(joined.index[joined["q_ref"] < alpha])
            intersection = predicted & observed
            union = predicted | observed
            direction_set = [
                pathway
                for pathway in intersection
                if np.sign(joined.at[pathway, "effect"])
                == np.sign(joined.at[pathway, "effect_ref"])
            ]
            records.append(
                {
                    "reference": reference_name,
                    "strategy": strategy,
                    "selection_rule": selection_rules[strategy],
                    "pathway_universe": "global_intersection",
                    "n_common_pathways": len(joined),
                    "spearman_nes": joined["effect"].corr(
                        joined["effect_ref"], method="spearman"
                    ),
                    "spearman_signed_logp": signed_candidate.corr(
                        signed_reference, method="spearman"
                    ),
                    "n_predicted_significant": len(predicted),
                    "n_reference_significant": len(observed),
                    "precision": (
                        len(intersection) / len(predicted) if predicted else np.nan
                    ),
                    "recall": len(intersection) / len(observed) if observed else np.nan,
                    "directional_precision": (
                        len(direction_set) / len(predicted) if predicted else np.nan
                    ),
                    "directional_recall": (
                        len(direction_set) / len(observed) if observed else np.nan
                    ),
                    "jaccard": len(intersection) / len(union) if union else np.nan,
                    "direction_agreement_significant": (
                        len(direction_set) / len(intersection)
                        if intersection
                        else np.nan
                    ),
                }
            )
    return pd.DataFrame(records)


def pairwise_gsea_diagnostics(
    results: Mapping[str, pd.DataFrame],
    *,
    top_set_alpha: float = 0.05,
) -> pd.DataFrame:

    standardized = {name: _standardize_gsea(df) for name, df in results.items()}
    names = list(standardized)
    records: list[dict[str, object]] = []
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            left_name, right_name = names[i], names[j]
            joined = standardized[left_name].join(
                standardized[right_name],
                how="inner",
                lsuffix="_left",
                rsuffix="_right",
            )
            signed_left = np.sign(joined["nes_left"]) * -np.log10(
                joined["pval_left"].clip(1e-300, 1)
            )
            signed_right = np.sign(joined["nes_right"]) * -np.log10(
                joined["pval_right"].clip(1e-300, 1)
            )
            left_top = set(joined.index[joined["padj_left"] < top_set_alpha])
            right_top = set(joined.index[joined["padj_right"] < top_set_alpha])
            overlap = left_top & right_top
            union = left_top | right_top
            direction = [
                pathway
                for pathway in overlap
                if np.sign(joined.at[pathway, "nes_left"])
                == np.sign(joined.at[pathway, "nes_right"])
            ]
            records.append(
                {
                    "method1": left_name,
                    "method2": right_name,
                    "n_common_pathways": len(joined),
                    "spearman_nes": joined["nes_left"].corr(
                        joined["nes_right"], method="spearman"
                    ),
                    "spearman_signed_logp": signed_left.corr(
                        signed_right, method="spearman"
                    ),
                    "jaccard": len(overlap) / len(union) if union else np.nan,
                    "direction_concordance": (
                        len(direction) / len(overlap) if overlap else np.nan
                    ),
                }
            )
    return pd.DataFrame(records)


def summarize_gsea_directory(
    gsea_dir: str | Path,
    out_dir: str | Path,
    *,
    top_set_alpha: float = 0.05,
    consensus_alpha: float = 0.05,
    direction_threshold: float = 0.8,
    min_families: int = 2,
    p_adjust_method: str = "bh",
) -> list[str]:

    gsea_path = Path(gsea_dir)
    output_path = Path(out_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    manifest_path = gsea_path / "gsea_runs.tsv"
    if not manifest_path.exists():
        raise FileNotFoundError(f"gsea_runs.tsv not found: {manifest_path}")
    runs = pd.read_csv(manifest_path, sep="\t")
    required = {
        "dataset",
        "contrast_id",
        "prefix",
        "role",
        "algorithm_family",
    }
    missing = required - set(runs.columns)
    if missing:
        raise ValueError(
            "gsea_runs.tsv is missing columns: " + ", ".join(sorted(missing))
        )

    generated: list[str] = []
    skipped: list[dict[str, object]] = []
    for contrast_id, subset in runs.groupby("contrast_id", sort=True):
        results: dict[str, pd.DataFrame] = {}
        annotations: dict[str, DatasetAnnotation] = {}
        for row in subset.itertuples(index=False):
            result_path = gsea_path / f"gsea_res_{row.prefix}.tsv"
            if not result_path.exists():
                raise FileNotFoundError(f"GSEA result not found: {result_path}")
            if row.dataset in results:
                raise ValueError(
                    f"Duplicate dataset '{row.dataset}' for contrast {contrast_id}."
                )
            results[row.dataset] = pd.read_csv(result_path, sep="\t")
            annotations[row.dataset] = DatasetAnnotation(
                role=str(row.role), family=str(row.algorithm_family)
            )

        prediction_annotations = [
            annotation
            for annotation in annotations.values()
            if annotation.role.lower() == "prediction"
        ]
        prediction_families = {
            annotation.family for annotation in prediction_annotations
        }
        if len(prediction_annotations) < 2 or len(prediction_families) < min_families:
            skipped.append(
                {
                    "contrast_id": contrast_id,
                    "n_prediction_outputs": len(prediction_annotations),
                    "n_prediction_families": len(prediction_families),
                    "required_families": min_families,
                    "reason": "insufficient independent prediction outputs/families",
                }
            )
            continue

        consensus, family_detail = build_consensus(
            results,
            annotations,
            top_set_alpha=top_set_alpha,
            consensus_alpha=consensus_alpha,
            direction_threshold=direction_threshold,
            min_families=min_families,
            p_adjust_method=p_adjust_method,
        )
        evaluation = evaluate_against_references(
            results, annotations, consensus, alpha=consensus_alpha
        )
        pairwise = pairwise_gsea_diagnostics(results, top_set_alpha=top_set_alpha)
        safe_id = str(contrast_id).replace("/", "_")
        reportable = consensus[consensus["reportable"]].sort_values(
            ["family_cauchy_q", "direction_concordance"],
            ascending=[True, False],
        )
        reportable = reportable.copy()
        reportable.insert(0, "selection_rank", range(1, len(reportable) + 1))
        outputs = {
            f"GSEA_{safe_id}_consensus.tsv": consensus,
            f"GSEA_{safe_id}_reportable.tsv": reportable,
            f"GSEA_{safe_id}_family_detail.tsv": family_detail,
            f"GSEA_{safe_id}_heldout_evaluation.tsv": evaluation,
            f"GSEA_{safe_id}_pairwise_diagnostics.tsv": pairwise,
        }
        for filename, frame in outputs.items():
            frame.to_csv(output_path / filename, sep="\t", index=False)
            generated.append(filename)

    settings = pd.DataFrame(
        [
            {
                "top_set_alpha": top_set_alpha,
                "consensus_alpha": consensus_alpha,
                "direction_threshold": direction_threshold,
                "min_families": min_families,
                "p_adjust_method": p_adjust_method,
                "primary_combination": "family-weighted Cauchy",
                "sensitivity_combination": "naive Fisher",
            }
        ]
    )
    settings.to_csv(output_path / "consensus_settings.tsv", sep="\t", index=False)
    generated.append("consensus_settings.tsv")
    if skipped:
        pd.DataFrame(skipped).to_csv(
            output_path / "consensus_skipped.tsv", sep="\t", index=False
        )
        generated.append("consensus_skipped.tsv")
    return generated
