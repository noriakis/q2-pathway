from os import path
import pandas as pd
import qiime2
import shutil
import q2templates
import numpy as np
import os
import json
import base64
import subprocess
from tempfile import TemporaryDirectory
import seaborn as sns
from PIL import Image
from urllib.parse import quote
import io
import scipy
import matplotlib.pyplot as plt
from itertools import combinations
import pkg_resources
TEMPLATES = pkg_resources.resource_filename('q2_pathway', 'assets')

def aggregate(table: pd.DataFrame,
    method: str = "sum",
    module: bool = False
    ) -> pd.DataFrame:

    ## Download pathway (or module) info
    if module:
        kop = pd.read_csv("https://rest.kegg.jp/link/ko/module", sep="\t", header=None)
    else:
        kop = pd.read_csv("https://rest.kegg.jp/link/ko/pathway", sep="\t", header=None)
        kop = kop[kop[0].apply(lambda x: "path:ko" in x)]
    
    normtbl = kop.groupby(1).apply(lambda x: len(x[0].unique().tolist()))

    ttable = table.T
    all_path = kop[0].unique()
    kos = ["ko:"+i for i in ttable.index.to_list()]
    outs = []
    for i in all_path:
        inkos = kop[kop[0] == i][1].tolist()
        inkos = list(set(kos) & set(inkos))
        if len(inkos) >= 1:
            out = pd.DataFrame(ttable.loc[[i.split(":")[1] for i in inkos], :].apply(lambda x: sum(x), axis=0))
            out.columns = [i]
            outs.append(out)
        else:
            pass
    conc = pd.concat(outs, axis=1)
    return(conc)