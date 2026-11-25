from qiime2.plugin import Plugin, Str, Int, Bool, Float, Citations, Metadata, List
from q2_types.feature_table import FeatureTable, Frequency
from q2_types.feature_data import FeatureData, Sequence
import q2_pathway

from q2_pathway import (
    PiphillinDatabase,
    PiphillinDatabaseFormat,
    CnTsvGzFormat,
    GFCnTsvGzFormat,
    FastaGzFormat,
    DatabaseManifestFormat,
)

import importlib

citations = Citations.load("citations.bib", package="q2_pathway")

plugin = Plugin(
    name="pathway",
    version=q2_pathway.__version__,
    website="https://github.com/noriakis/q2-pathway",
    package="q2_pathway",
    description=(
        "QIIME 2 plugin for comparing and analyzing biological pathway information based on gene family abundances."
    ),
    short_description="Compare and analyze gene family and biological pathway information",
    citations=[
        citations["Sato2023Bioinformatics"],
        citations["Kanehisa2000NAR"],
        citations["Fernandes2013PLOSOne"],
        citations["Love2014GenomeBiol"],
        citations["Narayan2020BMCGenomics"],
        citations["Wemheuer2020EnvMicro"],
        citations["Sun2020Microbiome"],
    ],
)

# Register semantic types
plugin.register_semantic_types(PiphillinDatabase)

# Register formats
plugin.register_formats(
    PiphillinDatabaseFormat,
    CnTsvGzFormat,
    GFCnTsvGzFormat,
    FastaGzFormat,
    DatabaseManifestFormat,
)

# Define and register new ArtifactClass
plugin.register_artifact_class(
    PiphillinDatabase,
    PiphillinDatabaseFormat,
    description="A Piphillin database format.",
)

plugin.visualizers.register_function(
    function=q2_pathway.kegg,
    inputs={"ko_table": FeatureTable[Frequency]},
    input_descriptions={"ko_table": "table containing KO abundance per sample"},
    parameters={
        "metadata": Metadata,
        "pathway_id": Str,
        "map_ko": Bool,
        "low_color": Str,
        "high_color": Str,
        "highlight": Metadata,
        "tss": Bool,
        "method": Str,
        "mc_samples": Int,
        "rank": Str,
    },
    parameter_descriptions={
        "pathway_id": "pathway to visualize, should start with ko.",
        "map_ko": "insert the KO description on the output table",
        "low_color": "color for low value",
        "high_color": "color for high value",
        "highlight": "highlight the nodes based on metadata. Should have 'thresholded' column",
        "tss": "total-sum scaling per sample before all the analysis",
        "method": "which value to show in the image",
        "mc_samples": "parameter for ALDEx2::aldex",
        "rank": "Which column to use for ranking in ALDEx2 and DESeq2.",
    },
    name="Plot KEGG PATHWAY",
    description=(
        "Plot the statistics of KO abundances between group on KEGG PATHWAY image."
    ),
)

plugin.visualizers.register_function(
    function=q2_pathway.kegg_multi,
    inputs={"tables": List[FeatureTable[Frequency]]},
    input_descriptions={
        "tables": "method-specific KO abundance tables to display side by side"
    },
    parameters={
        "metadata": Metadata,
        "pathway_id": Str,
        "tables_name": List[Str],
        "metadata_column": Str,
        "level1": Str,
        "level2": Str,
        "map_ko": Bool,
        "low_color": Str,
        "mid_color": Str,
        "high_color": Str,
        "missing_color": Str,
        "tss": Bool,
        "method": Str,
        "rank": Str,
        "mc_samples": Int,
        "direction_threshold": Float,
        "min_methods": Int,
    },
    parameter_descriptions={
        "pathway_id": "KEGG pathway identifier, for example ko00280.",
        "tables_name": "Labels and left-to-right node subdivision order, one per KO table.",
        "metadata_column": "Categorical metadata column to compare. If omitted, all eligible categorical columns are used.",
        "level1": "First group; the displayed direction is level1 minus level2.",
        "level2": "Second group; must be supplied together with level1.",
        "map_ko": "Include KEGG KO descriptions in the node-level output table.",
        "low_color": "Color for negative statistics.",
        "mid_color": "Color for zero statistics.",
        "high_color": "Color for positive statistics.",
        "missing_color": "Color for a KO not available from a prediction output.",
        "tss": "Apply total-sum scaling before Welch statistics.",
        "method": "KO-level statistic: t, deseq2, or aldex2.",
        "rank": "Statistic column for DESeq2 or ALDEx2; defaults to log2FoldChange or effect.",
        "mc_samples": "Monte Carlo samples used by ALDEx2.",
        "direction_threshold": "Required fraction of non-zero method statistics sharing one direction.",
        "min_methods": "Minimum number of directional method statistics needed to label a KO direction-consistent.",
    },
    name="Plot multiple prediction outputs on one KEGG pathway",
    description=(
        "Display one colored subdivision per prediction output in every KEGG node, "
        "using a shared color scale and exporting the underlying KO statistics."
    ),
)

plugin.visualizers.register_function(
    function=q2_pathway.gsea,
    inputs={"tables": List[FeatureTable[Frequency]]},
    input_descriptions={"tables": "list of tables containing KO abundance per sample"},
    parameters={
        "metadata": Metadata,
        "tss": Bool,
        "methods": List[Str],
        "mc_samples": Int,
        "map_pathway": Bool,
        "module": Bool,
        "tables_name": List[Str],
        "bg": Str,
        "rank": Str,
        "same": Bool,
        "min_size": Int,
        "max_size": Int,
        "split_str": Str,
        "table_roles": List[Str],
        "algorithm_families": List[Str],
        "top_set_alpha": Float,
        "consensus_alpha": Float,
        "direction_threshold": Float,
        "min_families": Int,
        "p_adjust_method": Str,
    },
    parameter_descriptions={
        "tss": "total-sum scaling per sample before all the analysis",
        "methods": "which value to show in the image",
        "mc_samples": "parameter for ALDEx2::aldex",
        "map_pathway": "map pathway name",
        "module": "If specified, perform GSEA based on module - KO relationship. default to False",
        "tables_name": "table name for the output, must be the same length as the specified table list",
        "bg": "Background KOs, default to `all`. If other option is specified, subset for the KOs within the corresponding table. ORA option.",
        "rank": "Which column to use for ranking in ALDEx2 and DESeq2.",
        "same": "Use same gene set across the gene family tables",
        "min_size": "parameter in fgsea: Minimal size of a gene set to test",
        "max_size": "parameter in fgsea: Maximal size of a gene set to test",
        "split_str": "Split the sample name and take the first argument, effective for all tables",
        "table_roles": "One role per table: prediction or reference. Reference tables are held out from consensus scoring.",
        "algorithm_families": "One algorithm family per table. Database variants from the same algorithm should share a family name.",
        "top_set_alpha": "Adjusted-p threshold used for overlap, intersection, and majority-vote diagnostics.",
        "consensus_alpha": "Adjusted-p threshold used to report consensus pathways.",
        "direction_threshold": "Minimum fraction of algorithm families agreeing on NES direction.",
        "min_families": "Minimum number of prediction algorithm families required for consensus.",
        "p_adjust_method": "Multiple-testing correction for combined pathway p-values: bh or bonferroni.",
    },
    name="Perform GSEA by the R package fgsea",
    description=("Perform GSEA by the R package fgsea"),
)

plugin.visualizers.register_function(
    function=q2_pathway.summarize,
    inputs={"tables": List[FeatureTable[Frequency]]},
    input_descriptions={"tables": "list of tables containing KO abundance per sample"},
    parameters={
        "metadata": Metadata,
        "tss": Bool,
        "split_str": Str,
        "convert_table": Metadata,
        "tables_name": List[Str],
        "table_roles": List[Str],
        "algorithm_families": List[Str],
        "zero_replacement": Float,
    },
    parameter_descriptions={
        "tss": (
            "Compatibility flag. KO diagnostics always close each sample before "
            "zero replacement and CLR, so this flag does not change the result."
        ),
        "split_str": 'Split the string of column names and takes the first argument, e.g. if XXXX_1234_5678 and "_" is specified, the column will be converted to XXX',
        "convert_table": 'Converting table artifact (Column name "converted" will be used)',
        "tables_name": "table name for the output, must be the same length as the specified table list",
        "table_roles": "One role per table: prediction or reference.",
        "algorithm_families": "One algorithm family per table; database variants should share a family.",
        "zero_replacement": "Zero-replacement value applied after closure and before CLR.",
    },
    name="Summarize the output of functional prediction.",
    description=("Summarize the output of functional prediction."),
)

plugin.visualizers.register_function(
    function=q2_pathway.summarize_detail,
    inputs={"tables": List[FeatureTable[Frequency]]},
    input_descriptions={"tables": "list of tables containing KO abundance per sample"},
    parameters={
        "metadata": Metadata,
        "first": Int,
        "tss": Bool,
        "method": Str,
        "candidate": Str,
        "candidate_pathway": Str,
        "split_str": Str,
        "convert_table": Metadata,
        "map_ko": Bool,
        "cor_fig_width": Int,
        "cor_thresh": Float,
        "use_p": Bool,
        "tables_name": List[Str],
        "skip": Bool,
    },
    parameter_descriptions={
        "tss": "Performs total-sum scaling per sample before all the analysis",
        "first": "If `candidate` or `candidate_pathway` is not specified, The `first` genes sorted by the average abundance of the first table will be summarized.",
        "method": "Correlation method, default to `spearman`",
        "candidate": "Candidate KO to summarize",
        "candidate_pathway": "Candidate pathway ID in KEGG PATHWAY to summarize",
        "split_str": 'Split the string of column names and takes the first argument, e.g. if XXXX_1234_5678 and "_" is specified, the column will be converted to XXX',
        "convert_table": 'Converting table artifact (Column name "converted" will be used)',
        "map_ko": "Insert the KO description in the output",
        "cor_fig_width": "Correlation figure width",
        "cor_thresh": "If specified, make additional visualization using only the KOs above the specified threshold, default to None",
        "use_p": "Use p-values from Wilcoxon rank sum tests between conditions in the metadata for evaulation (Sun et al. 2020)",
        "tables_name": "table name for the output, must be the same length as the specified table list",
        "skip": "Skip the producing of the correlation and heatmap image (useful when whole gene families are to be evaluated and the resulting QZV does not load)",
    },
    name="Summarize the output of functional prediction.",
    description=("Summarize the output of functional prediction."),
)

plugin.methods.register_function(
    function=q2_pathway.infer,
    inputs={
        "sequences": FeatureData[Sequence],
        "seq_table": FeatureTable[Frequency],
        "database": PiphillinDatabase,
    },
    input_descriptions={
        "sequences": "Representative sequences to be profiled",
        "seq_table": "Sequence count table",
        "database": "Piphillin database artifact",
    },
    outputs=[("table", FeatureTable[Frequency])],
    parameters={
        "threads": Int,
        "full": Bool,
        "full_id": Str,
        "full_ko_id": Str,
        "pct_id": Float,
    },
    parameter_descriptions={
        "threads": "The number of threads",
        "full": "Output the full stratified table",
        "full_id": "When the full stratified table is calculated, retain only this reference-genome ID.",
        "full_ko_id": "When the full stratified table is calculated, retain only this KO ID.",
        "pct_id": "Percent of identity, default to 0.99",
    },
    name="Run Piphillin algorithm",
    description=("Run Piphillin algorithm"),
)

plugin.methods.register_function(
    function=q2_pathway.infer_t4f2,
    inputs={
        "sequences": FeatureData[Sequence],
        "seq_table": FeatureTable[Frequency],
    },
    input_descriptions={
        "sequences": "Representative sequences to be profiled",
        "seq_table": "Sequence count table",
    },
    outputs=[("table", FeatureTable[Frequency])],
    parameters={
        "threads": Int,
        "pct_id": Float,
        "database_mode": Str,
        "resource_dir": Str,
        "auto_install": Bool,
    },
    parameter_descriptions={
        "threads": "The number of threads",
        "pct_id": "Percent of identity, default to 0.99",
        "database_mode": "Ref99NR or Ref100NR",
        "resource_dir": "User-writable directory containing Tax4Fun2_ReferenceData_v2.",
        "auto_install": "Download verified reference data when absent.",
    },
    name="Run Tax4Fun2 algorithm",
    description=("Run Tax4Fun2 algorithm"),
)

plugin.methods.register_function(
    function=q2_pathway.aggregate,
    inputs={"table": FeatureTable[Frequency]},
    input_descriptions={"table": "table containing KO abundance per sample"},
    outputs=[("out_table", FeatureTable[Frequency])],
    parameters={"method": Str, "module": Bool},
    parameter_descriptions={
        "method": "how to aggregate the abundance",
        "module": "If specified, perform GSEA based on module - KO relationship. default to False",
    },
    name="Aggregate family abundance to high order abundance",
    description=("Aggregate family abundance to high order abundance"),
)

plugin.visualizers.register_function(
    function=q2_pathway.contribute,
    inputs={"table": FeatureTable[Frequency]},
    input_descriptions={"table": "table containing KO abundance per sample"},
    parameters={
        "metadata": Metadata,
        "candidate": Str,
        "fig_height": Int,
    },
    parameter_descriptions={
        "candidate": "KO to be evaluated",
        "fig_height": "Figure height",
    },
    name="Plot per-taxon abundance of specified gene family",
    description=("Plot per-taxon abundance of specified gene family."),
)


plugin.visualizers.register_function(
    function=q2_pathway.permute,
    inputs={"tables": List[FeatureTable[Frequency]]},
    input_descriptions={
        "tables": "list of two tables containing KO abundance per sample"
    },
    parameters={
        "tss": Bool,
        "method": Str,
        "candidate": Str,
        "candidate_pathway": Str,
        "split_str": Str,
        "convert_table": Metadata,
        "tables_name": List[Str],
        "nperm": Int,
    },
    parameter_descriptions={
        "tss": "Performs total-sum scaling per sample before all the analysis",
        "method": "Correlation method, default to `spearman`",
        "candidate": "Candidate KO to summarize",
        "candidate_pathway": "Candidate pathway ID in KEGG PATHWAY to summarize",
        "split_str": 'Split the string of column names and takes the first argument, e.g. if XXXX_1234_5678 and "_" is specified, the column will be converted to XXX',
        "convert_table": 'Converting table artifact (Column name "convert" will be used)',
        "tables_name": "table name for the output, must be the same length as the specified table list",
        "nperm": "Number of permutations",
    },
    name="Assess correlation of the output of functional prediction based on permutation.",
    description=("Summarize the output of functional prediction."),
)

importlib.import_module("q2_pathway._transformers")
