import pandas as pd
import pkg_resources

TEMPLATES = pkg_resources.resource_filename("q2_pathway", "assets")


def aggregate(
    table: pd.DataFrame, method: str = "sum", module: bool = False
) -> pd.DataFrame:
    ## Download pathway (or module) info
    if module:
        kop = pd.read_csv("https://rest.kegg.jp/link/ko/module", sep="\t", header=None)
    else:
        kop = pd.read_csv("https://rest.kegg.jp/link/ko/pathway", sep="\t", header=None)
        kop = kop[kop[0].apply(lambda x: "path:ko" in x)]

    ttable = table.T
    all_path = kop[0].unique()
    kos = ["ko:" + i for i in ttable.index.to_list()]
    outs = []
    for i in all_path:
        inkos = kop[kop[0] == i][1].tolist()
        inkos = list(set(kos) & set(inkos))
        if len(inkos) >= 1:
            if method == "sum":
                out = pd.DataFrame(
                    ttable.loc[[i.split(":")[1] for i in inkos], :].apply(
                        lambda x: sum(x), axis=0
                    )
                )
            elif method == "median":
                out = pd.DataFrame(
                    ttable.loc[[i.split(":")[1] for i in inkos], :].apply(
                        lambda x: x.median(), axis=0
                    )
                )
            elif method == "mean":
                out = pd.DataFrame(
                    ttable.loc[[i.split(":")[1] for i in inkos], :].apply(
                        lambda x: x.mean(), axis=0
                    )
                )
            else:
                raise ValueError("Please specify mean, median or sum in method.")

            out.columns = [i]
            outs.append(out)
        else:
            pass
    conc = pd.concat(outs, axis=1)
    return conc
