from os import path
import os
import pandas as pd
import subprocess
from tempfile import TemporaryDirectory
import pkg_resources

TEMPLATES = pkg_resources.resource_filename("q2_pathway", "assets")

def infer(
    sequences: pd.Series,
    seq_table: pd.DataFrame,
    threads: int = 1,
    full: bool = False,
    pct_id: float = 0.99
) -> pd.DataFrame:
    """
    [TODO] Option to use q2-gcn-norm for normalizing by rrnDB
    Should we separate the reference files to artifact and distribute, rather than
    putting them in conda release?
    """
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
                    "Please provide Tax4Fun2 default database artifact to `reference_database`."
                )
                
            db = qiime2.Artifact.load(reference_database)
            d = db.view(T4F2DirectoryFormat)
            db_path = d.path
            
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


def infer_t4f2(
    sequences: pd.Series,
    seq_table: pd.DataFrame,
    database: str,
    threads: int = 1,
    database_mode: str = "Ref99NR",
    pct_id: float = 0.99,
) -> pd.DataFrame:
    """
    We need a prior installation of Tax4Fun2 in the R in QIIME 2 environment.
    Download https://zenodo.org/records/10035668/files/Tax4Fun2_1.1.5.tar.gz, and run
    `R CMD INSTALL Tax4Fun2_1.1.5.tar.gz`.
    """

    if threads < 1:
        raise ValueError("Thread number should be positive.")

    with TemporaryDirectory() as temp_dir:
        repseq = path.join(temp_dir, "rep_seqs.fna")

        with open(repseq, "w") as fna:
            for seqname, sequence in sequences.items():
                print(">" + str(seqname) + "\n" + str(sequence), file=fna)
        seq_table.T.to_csv(path.join(temp_dir, "seqtab.txt"), sep="\t")

        ## Needs to extract database every time (use cache)
        excmd = [
            "tar", "-zxf",
            database + "/Tax4Fun2_ReferenceData_v2.tar.gz",
            "-C", temp_dir
            ]
        try:
            subprocess.run(excmd, check=True)
        except subprocess.CalledProcessError as _:
            raise ValueError("Error extracting database.")
   
        cmd = [
            "Rscript",
            path.join(TEMPLATES, "perform_tax4fun2.R"),
            temp_dir,
            "rep_seqs.fna",  ## rep-seqs
            "seqtab.txt",
            path.join(temp_dir, "Tax4Fun2_ReferenceData_v2"),
            str(pct_id),
            str(threads),
            str(database_mode)
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
