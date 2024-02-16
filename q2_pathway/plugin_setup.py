

from qiime2.plugin import (Plugin, Str, Choices, Int, Bool, Range, Float, Citations, Metadata)
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
    citations=[citations['Sato2023Bioinformatics'], citations['Kanehisa2000NAR'], citations["Narayan2020BMCGenomics"]]
)

## kegg
plugin.visualizers.register_function(
    function=q2_pathway.kegg,
    inputs={'ko_table': FeatureTable[Frequency]},
    parameters={'metadata': Metadata, 'pathway_id': Str, 'map_ko': Bool, 'low_color': Str, 'high_color': Str},
    name="Plot KEGG PATHWAY",
    description=("Plot the statistics of KO abundances between group on KEGG PATHWAY image.")
)


## gsea
plugin.visualizers.register_function(
    function=q2_pathway.gsea,
    inputs={'ko_table': FeatureTable[Frequency]},
    parameters={'metadata': Metadata},
    name="Perform GSEA by the R package fgsea (experimental)",
    description=("Perform GSEA by the R package fgsea (experimental)")
)

## summarize
plugin.visualizers.register_function(
    function=q2_pathway.summarize,
    inputs={'ko_table': FeatureTable[Frequency]},
    parameters={'metadata': Metadata},
    name="Summarize the stratified output of functional prediction.",
    description=("Summarize the stratified output of functional prediction.")
)


## infer
plugin.methods.register_function(
    function=q2_pathway.infer,
    inputs={'sequences': FeatureData[Sequence],
        'seq_table': FeatureTable[Frequency]
    },
    outputs=[('table', FeatureTable[Frequency])],
    parameters={'threads': Int, 'reference_sequences': Str,
    'cn_table': Str, 'cn_16s_table': Str, 'full': Bool, 'pct_id': Float},
    name="Run Piphillin algorithm",
    description=("Run Piphillin algorithm")
)
