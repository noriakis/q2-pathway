from qiime2.plugin import Plugin, Str, Int, Bool, Float, Citations, Metadata, List
from q2_types.feature_table import FeatureTable, Frequency
from q2_types.feature_data import FeatureData, Sequence
import q2_pathway

from q2_pathway import (
    T4F2Database,
    T4F2DatabaseFormat)
# import importlib


citations = Citations.load("citations.bib", package="q2_pathway")

plugin = Plugin(
    name="pathway",
    version="2024.2",
    website="https://github.com/noriakis/q2-pathway",
    package="q2_pathway",
    description=(
        "QIIME2 plugin for analyzing biological pathway information based on gene family abundances."
    ),
    short_description="Analyze biological pathway information",
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
plugin.register_semantic_types(T4F2Database)

# Register formats
plugin.register_formats(T4F2DatabaseFormat)

# Define and register new ArtifactClass
plugin.register_artifact_class(T4F2Database,
                               T4F2DatabaseFormat,
                               description="A Tax4Fun2 database directory.")

## kegg
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
        "tss": Bool,
        "method": Str,
        "mc_samples": Int,
    },
    parameter_descriptions={
        "pathway_id": "pathway to visualize, should start with ko.",
        "map_ko": "insert the KO description on the output table",
        "low_color": "color for low value",
        "high_color": "color for high value",
        "tss": "total-sum scaling per sample before all the analysis",
        "method": "which value to show in the image",
        "mc_samples": "parameter for ALDEx2::aldex",
    },
    name="Plot KEGG PATHWAY",
    description=(
        "Plot the statistics of KO abundances between group on KEGG PATHWAY image."
    ),
)


## gsea
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
    },
    parameter_descriptions={
        "tss": "total-sum scaling per sample before all the analysis",
        "method": "which value to show in the image",
        "mc_samples": "parameter for ALDEx2::aldex",
        "map_pathway": "map pathway name",
        "module": "If specified, perform GSEA based on module - KO relationship. default to False",
        "tables_name": "table name for the output, must be the same length as the specified table list",
        "bg": "Background KOs, default to `all`. If other option is specified, subset for the KOs within the corresponding table. ORA option.",
    },
    name="Perform GSEA by the R package fgsea (experimental)",
    description=("Perform GSEA by the R package fgsea (experimental)"),
)

## summarize
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
    },
    name="Summarize the output of functional prediction.",
    description=("Summarize the output of functional prediction."),
)


## infer
plugin.methods.register_function(
    function=q2_pathway.infer,
    inputs={"sequences": FeatureData[Sequence], "seq_table": FeatureTable[Frequency]},
    input_descriptions={
        "sequences": "Representative sequences to be profiled",
        "seq_table": "Sequence count table",
    },
    outputs=[("table", FeatureTable[Frequency])],
    parameters={
        "threads": Int,
        # 'reference_sequences': Str,
        # 'cn_table': Str,
        # 'cn_16s_table': Str,
        "full": Bool,
        "pct_id": Float
    },
    parameter_descriptions={
        "threads": "The number of threads",
        # 'reference_sequences': '16S reference sequences, default to the preset database.',
        # 'cn_table': 'gene copy number table, default to the preset database.',
        # 'cn_16s_table': '16S gene copy number table, default to the preset database.',
        "full": "Output the full stratified table",
        "pct_id": "Percent of identity, default to 0.99"
    },
    name="Run Piphillin algorithm",
    description=("Run Piphillin algorithm"),
)

## infer_t4f2
plugin.methods.register_function(
    function=q2_pathway.infer_t4f2,
    inputs={"sequences": FeatureData[Sequence], "seq_table": FeatureTable[Frequency], "database": T4F2Database},
    input_descriptions={
        "sequences": "Representative sequences to be profiled",
        "seq_table": "Sequence count table",
        "database": "Tax4Fun2 default database artifact"
    },
    outputs=[("table", FeatureTable[Frequency])],
    parameters={
        "threads": Int,
        "pct_id": Float
    },
    parameter_descriptions={
        "threads": "The number of threads",
        "pct_id": "Percent of identity, default to 0.99",
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

# importlib.import_module('q2_pathway._transformers')
