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
)

import importlib

citations = Citations.load("citations.bib", package="q2_pathway")

plugin = Plugin(
    name="pathway",
    version="2024.5",
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
    ],
)

# Register semantic types
plugin.register_semantic_types(PiphillinDatabase)

# Register formats
plugin.register_formats(PiphillinDatabaseFormat, CnTsvGzFormat, GFCnTsvGzFormat, FastaGzFormat)

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
        "rank": Str
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
        "rank": "Which column to use for ranking in ALDEx2 and DESeq2."
    },
    name="Plot KEGG PATHWAY",
    description=(
        "Plot the statistics of KO abundances between group on KEGG PATHWAY image."
    ),
)

plugin.visualizers.register_function(
    function=q2_pathway.gsea,
    inputs={"tables": List[FeatureTable[Frequency]]},
    input_descriptions={"tables": "list of tables containing KO abundance per sample"},
    parameters={
        "metadata": Metadata,
        "tss": Bool,
        "method": Str,
        "mc_samples": Int,
        "map_pathway": Bool,
        "module": Bool,
        "tables_name": List[Str],
        "bg": Str,
        "rank": Str,
        "same": Bool,
        "min_size": Int,
        "max_size": Int
    },
    parameter_descriptions={
        "tss": "total-sum scaling per sample before all the analysis",
        "method": "which value to show in the image",
        "mc_samples": "parameter for ALDEx2::aldex",
        "map_pathway": "map pathway name",
        "module": "If specified, perform GSEA based on module - KO relationship. default to False",
        "tables_name": "table name for the output, must be the same length as the specified table list",
        "bg": "Background KOs, default to `all`. If other option is specified, subset for the KOs within the corresponding table. ORA option.",
        "rank": "Which column to use for ranking in ALDEx2 and DESeq2.",
        "same": "Use same gene set across the gene family tables",
        "min_size": "parameter in fgsea: Minimal size of a gene set to test",
        "max_size": "parameter in fgsea: Maximal size of a gene set to test"
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
        "convert_table": 'Converting table artifact (Column name "convert" will be used)',
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
        "pct_id": Float,
    },
    parameter_descriptions={
        "threads": "The number of threads",
        "full": "Output the full stratified table",
        "full_id": "When the full stratified table is calculated, only this ID will be remained.",
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
    parameters={"threads": Int, "pct_id": Float, "database_mode": Str},
    parameter_descriptions={
        "threads": "The number of threads",
        "pct_id": "Percent of identity, default to 0.99",
        "database_mode": "Ref99NR or Ref100NR",
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

importlib.import_module("q2_pathway._transformers")
