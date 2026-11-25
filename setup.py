from setuptools import setup, find_packages

setup(
    name="q2-pathway",
    version="2026.10.0rc3",
    packages=find_packages(),
    author="Noriaki Sato",
    author_email="nori@hgc.jp",
    description="QIIME 2 plugin for analyzing biological pathway information based on gene family abundances.",
    license="BSD-3-Clause",
    url="https://github.com/noriakis/q2-pathway",
    entry_points={
        "qiime2.plugins": ["q2-pathway=q2_pathway.plugin_setup:plugin"],
        "console_scripts": [
            "q2-pathway-build-db=q2_pathway.database_builder:main",
            "q2-pathway-install-tax4fun2-data=q2_pathway._infer:install_t4f2_cli",
        ],
    },
    package_data={
        "q2_pathway": [
            "citations.bib",
            "assets/index.html",
            "assets/dist/*",
            "assets/*.R",
            "tests/data/*",
            "tests/data/test-db/*",
        ]
    },
    zip_safe=False,
)
