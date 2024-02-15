
from setuptools import setup, find_packages

setup(
    name="q2-pathway",
    version="2024.2",
    packages=find_packages(),
    author="Noriaki Sato",
    author_email="nori@hgc.jp",
    description="QIIME2 plugin for visualizing and analyzing pathway information based on gene abundances.",
    license='BSD-3-Clause',
    url="https://github.com/noriakis/q2-pathway",
    entry_points={
        'qiime2.plugins': ['q2-pathway=q2_pathway.plugin_setup:plugin']
    },
    package_data={
        'q2_pathway': ['citations.bib', 'assets/index.html', 'assets/dist/*', 'assets/perform_gsea.R']
    },
    zip_safe=False,
)
