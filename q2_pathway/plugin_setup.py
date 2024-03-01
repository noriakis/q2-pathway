

from qiime2.plugin import (Plugin, Str, Choices, Int, Bool, Range, Float, Citations, Metadata, List)
from q2_types.feature_table import FeatureTable, Frequency
from q2_types.feature_data import FeatureData, Sequence
import q2_pathway

citations = Citations.load('citations.bib', package='q2_pathway')

plugin = Plugin(
    name='pathway',
    version="2024.2",
    website='https://github.com/noriakis/q2-pathway',
    package='q2_pathway',
    description=('QIIME2 plugin for visualizing and analyzing pathway information based on gene abundances.'),
    short_description='Visualize KEGG PATHWAY.',
    citations=[citations['Sato2023Bioinformatics'], citations['Kanehisa2000NAR'], citations["Narayan2020BMCGenomics"],
               citations['Fernandes2013PLOSOne']]
)

## kegg
plugin.visualizers.register_function(
    function=q2_pathway.kegg,
    inputs={'ko_table': FeatureTable[Frequency]},
    input_descriptions={'ko_table': 'table containing KO abundance per sample'},
    parameters={
        'metadata': Metadata,
        'pathway_id': Str,
        'map_ko': Bool,
        'low_color': Str,
        'high_color': Str,
        'tss': Bool,
        'method': Str,
        'mc_samples': Int
    },
    parameter_descriptions={
        'pathway_id': 'pathway to visualize, should start with ko.',
        'map_ko': 'insert the KO description on the output table',
        'low_color': 'color for low value',
        'high_color': 'color for high value',
        'tss': 'total-sum scaling per sample before all the analysis',
        'method': 'which value to show in the image',
        'mc_samples': 'parameter for ALDEx2::aldex'    
    },
    name="Plot KEGG PATHWAY",
    description=("Plot the statistics of KO abundances between group on KEGG PATHWAY image.")
)


## gsea
plugin.visualizers.register_function(
    function=q2_pathway.gsea,
    inputs={'ko_table': FeatureTable[Frequency]},
    input_descriptions={'ko_table': 'table containing KO abundance per sample'},
    parameters={
        'metadata': Metadata,
        'tss': Bool,
        'method': Str,
        'mc_samples': Int,
        'module': Bool
    },
    parameter_descriptions={
        'tss': 'total-sum scaling per sample before all the analysis',
        'method': 'which value to show in the image',
        'mc_samples': 'parameter for ALDEx2::aldex',
        'module': 'If specified, perform GSEA based on module - KO relationship. default to False'   
    },
    name="Perform GSEA by the R package fgsea (experimental)",
    description=("Perform GSEA by the R package fgsea (experimental)")
)

## summarize
plugin.visualizers.register_function(
    function=q2_pathway.summarize,
    inputs={
        'tables': List[FeatureTable[Frequency]]
    },
    input_descriptions={'tables': 'list of tables containing KO abundance per sample'},
    parameters={
        'metadata': Metadata,
        'first': Int,
        'tss': Bool,
        'method': Str,
        'candidate': Str,
        'candidate_pathway': Str,
        'split_str': Str,
        'convert': Str,
        'map_ko': Bool,
        'cor_fig_width': Int,
        'tables_name': List[Str]
    },
    parameter_descriptions={
        'tss': 'total-sum scaling per sample before all the analysis',
        'first': 'If candidate or candidate_pathway is not specified, The `first` genes sorted by average abundance will be summarized.',
        'method': 'correlation method, default to spearman',
        'candidate': 'candidate KO',
        'candidate_pathway': 'candidate pathway ID in KEGG PATHWAY',
        'split_str': 'split the string of column names and takes the first argument, e.g. if XXX_1234 and "_" is specified, the column will be converted to XXX',
        'convert': 'Converting table (first column the original sample ID and second column the converted ID)',
        'map_ko': 'Insert the KO description in the output',
        'cor_fig_width': 'Correlation figure width',
        'tables_name': 'table name for the output, must be the same length as the specified table list'
    },
    name="Summarize the output of functional prediction.",
    description=("Summarize the output of functional prediction.")
)


## infer
plugin.methods.register_function(
    function=q2_pathway.infer,
    inputs={
        'sequences': FeatureData[Sequence],
        'seq_table': FeatureTable[Frequency]
    },
    input_descriptions={
        'sequences': 'representative sequences to be profiled',
        'seq_table': 'sequence count table'
    },
    outputs=[('table', FeatureTable[Frequency])],
    parameters={
        'threads': Int,
        'reference_sequences': Str,
        'cn_table': Str,
        'cn_16s_table': Str,
        'full': Bool,
        'pct_id': Float
    },
    parameter_descriptions={
        'threads': 'The number of threads',
        'reference_sequences': '16S reference sequences, default to the preset database.',
        'cn_table': 'gene copy number table, default to the preset database.',
        'cn_16s_table': '16S gene copy number table, default to the preset database.',
        'full': 'Output the full stratified table',
        'pct_id': 'Percent of identity, default to 0.99'
    },
    name="Run Piphillin algorithm",
    description=("Run Piphillin algorithm")
)
