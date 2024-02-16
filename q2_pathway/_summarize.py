from os import path
import pandas as pd
import qiime2
import q2templates
import numpy as np
import os
import subprocess
from tempfile import TemporaryDirectory
import pkg_resources
TEMPLATES = pkg_resources.resource_filename('q2_pathway', 'assets')

def summarize(output_dir: str, ko_table: pd.DataFrame, ko_table2: pd.DataFrame,
	metadata: qiime2.Metadata) -> None:
    return(1)