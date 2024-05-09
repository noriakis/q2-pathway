from os import path
import pandas as pd
import subprocess
from tempfile import TemporaryDirectory
import pkg_resources
import requests

TEMPLATES = pkg_resources.resource_filename("q2_pathway", "assets")

def install(
    output_dir: str,
    type: str = "tax4fun2",
    ) -> None:
    """Currently, `infer` module has a parameter `reference_database` for specifying
    Tax4Fun2 database path, but the path specification in the parameter should be avoided.
    This function downloads the archive and database from Zenodo, and install the package
    and places the database on plugin directory, which can subsequently be used within 
    the infer function without specifying the path. This is currently implemented as 
    visualizer. Probably output the files as QZV and uses the file for the subsequent analysis?
    """
    print("1. Downloading Tax4Fun2 archive...")
    t4f2_url = "https://zenodo.org/records/10035668/files/Tax4Fun2_1.1.5.tar.gz"
    fn = path.join(TEMPLATES, "Tax4Fun2_1.1.5.tar.gz")
    res = requests.get(t4f2_url, stream=True)
    if res.status_code == 200:
        with open(fn, 'wb') as file:
            for chunk in res:
                file.write(chunk)
                
    ## Install Tax4Fun2
    print("2. Installing Tax4Fun2...")
    cmd = ["R", "CMD", "INSTALL", fn]
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as _:
        raise ValueError("Error installing tax4fun2.")

    print("3. Downloading Tax4Fun2 DB...")        
    t4f2_db_url = "https://zenodo.org/records/10035668/files/Tax4Fun2_ReferenceData_v2.tar.gz"
    fn = path.join(TEMPLATES, "Tax4Fun2_ReferenceData_v2.tar.gz")
    res = requests.get(t4f2_db_url, stream=True)
    if res.status_code == 200:
        with open(fn, 'wb') as file:
            for chunk in res:
                file.write(chunk)


def infer(
    sequences: pd.Series,
    seq_table: pd.DataFrame,
    threads: int = 1,
    full: bool = False,
    pct_id: float = 0.99,
    algorithm: str = "piphillin",
    reference_database: str = None,
) -> pd.DataFrame:
    ## [Idea] Option to use q2-gcn-norm for normalizing by rrnDB
    reference_sequences = path.join(TEMPLATES, "16S_seqs.fasta.gz")
    cn_table = path.join(TEMPLATES, "ko_copynum.tsv.gz")
    cn_16s_table = path.join(TEMPLATES, "16S_cn.tsv.gz")

    with TemporaryDirectory() as temp_dir:
        repseq = path.join(temp_dir, "rep_seqs.fna")

        with open(repseq, "w") as fna:
            for seqname, sequence in sequences.items():
                print(">" + str(seqname) + "\n" + str(sequence), file=fna)
        seq_table.T.to_csv(path.join(temp_dir, "seqtab.txt"), sep="\t")
        if algorithm == "piphillin":
            cmd = [
                "vsearch",
                "--usearch_global",
                repseq,
                "--db",
                reference_sequences,
                "--id",
                str(pct_id),
                "--top_hits_only",
                "--maxaccepts",
                "0",
                "--maxrejects",
                "0",
                "--uc_allhits",
                "--blast6out",
                path.join(temp_dir, "blast_out.txt"),
                "--threads",
                str(threads),
            ]
            try:
                subprocess.run(cmd, check=True)
            except subprocess.CalledProcessError as _:
                raise ValueError("Error running vsearch.")

            cmd = [
                "Rscript",
                path.join(TEMPLATES, "perform_piphillin.R"),
                temp_dir,
                cn_table,
                cn_16s_table,
                str(full),
            ]
            try:
                subprocess.run(cmd, check=True)
            except subprocess.CalledProcessError as _:
                raise ValueError("Error running piphillin algorithm.")
            ko = pd.read_csv(path.join(temp_dir, "ko_table.tsv"), sep="\t", header=0)

            if full:
                ko.index = ko.iloc[:, 0].map(str) + "_" + ko.iloc[:, 1]
                ko = ko.drop(ko.columns[0], axis=1)
                ko = ko.drop(ko.columns[0], axis=1)
                ## This will output tax1_K00001 type index

            return ko.T
        elif algorithm == "tax4fun2":
            ## Although the input filepath is discouraged, the Tax4Fun2 database is structured
            ## inside the archive and cannot be properly converted to QZA.
            if reference_database is None:
                raise ValueError(
                    "Please provide Tax4Fun2 default database path to `reference_database`."
                )
            cmd = [
                "Rscript",
                path.join(TEMPLATES, "perform_tax4fun2.R"),
                temp_dir,
                "rep_seqs.fna",  ## rep-seqs
                "seqtab.txt",
                reference_database,
                str(pct_id),
                str(threads),
            ]
            try:
                subprocess.run(cmd, check=True)
            except subprocess.CalledProcessError as _:
                raise ValueError("Error running tax4fun2 algorithm.")
            ko = pd.read_csv(
                path.join(temp_dir, "functional_prediction.txt"),
                sep="\t",
                header=0,
                index_col=0,
            )
            ko = ko.drop("description", axis=1)
            return ko.T

        else:
            raise ValueError("Please specify appropriate algorithm name")
