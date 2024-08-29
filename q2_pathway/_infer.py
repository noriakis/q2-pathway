from os import path
import pandas as pd
import subprocess
from tempfile import TemporaryDirectory
import pkg_resources
import requests
import hashlib

TEMPLATES = pkg_resources.resource_filename("q2_pathway", "assets")


def infer(
    sequences: pd.Series,
    seq_table: pd.DataFrame,
    database: str,
    threads: int = 1,
    full: bool = False,
    full_id: str = None,
    pct_id: float = 0.99,
) -> pd.DataFrame:
    """
    Inference using Piphillin algorithm

    """
    reference_sequences = path.join(database, "16S_seqs.fasta.gz")
    cn_table = path.join(database, "ko_copynum.tsv.gz")
    cn_16s_table = path.join(database, "16S_cn.tsv.gz")

    with TemporaryDirectory() as temp_dir:
        repseq = path.join(temp_dir, "rep_seqs.fna")

        with open(repseq, "w") as fna:
            for seqname, sequence in sequences.items():
                print(">" + str(seqname) + "\n" + str(sequence), file=fna)
        seq_table.T.to_csv(path.join(temp_dir, "seqtab.txt"), sep="\t")
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
            str(full_id),
        ]
        try:
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as _:
            raise ValueError("Error running piphillin algorithm.")

        if full:
            ko = pd.read_csv(path.join(temp_dir, "ko_table.txt"), sep="\t", header=0)
            ko["ind"] = ko.KO + "|" + ko.ID
            # ko = ko.iloc[:, 3].map(str) + "_" + ko.iloc[:, 1]
            ko = ko.pivot_table(index="ind", columns="sample", values="value").fillna(0)
            return ko.T
            ## This will output tax1_K00001 type index
        else:
            ko = pd.read_csv(path.join(temp_dir, "ko_table.tsv"), sep="\t", header=0)
            return ko.T


def install_t4f2():
    correct = "5c2edff2f597e5f061f79a5782a27567"

    fn = path.join(TEMPLATES, "Tax4Fun2_ReferenceData_v2.tar.gz")
    if path.isfile(fn):
        print("Tax4Fun2 already installed; proceeding.")
        print("If error was occurred in downloading, please remove " + fn)
        with open(fn, "rb") as file:
            fileData = file.read()
            db_hash = hashlib.md5(fileData).hexdigest()
        if db_hash != correct:
            print("Hash incorrect, proceeding anyway.")
        else:
            print("Correct hash for database.")
        return(0)

    fns = []
    urls = []
    print("1. Downloading Tax4Fun2 archive...")
    t4f2_url = "https://zenodo.org/records/10035668/files/Tax4Fun2_1.1.5.tar.gz"
    fn = path.join(TEMPLATES, "Tax4Fun2_1.1.5.tar.gz")
    if path.isfile(fn):
        raise ValueError("The file is already downloaded, please remove "+ fn +" if try to proceed.")
    fns.append(fn)
    urls.append(t4f2_url)

    res = requests.get(t4f2_url, stream=True)
    if res.status_code == 200:
        with open(fn, "wb") as file:
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
    t4f2_db_url = (
        "https://zenodo.org/records/10035668/files/Tax4Fun2_ReferenceData_v2.tar.gz"
    )
    fn = path.join(TEMPLATES, "Tax4Fun2_ReferenceData_v2.tar.gz")
    if path.isfile(fn):
        raise ValueError("The file is already downloaded, please remove "+ fn +" if try to proceed.")
    fns.append(fn)
    urls.append(t4f2_db_url)

    res = requests.get(t4f2_db_url, stream=True)
    if res.status_code == 200:
        with open(fn, "wb") as file:
            for chunk in res:
                file.write(chunk)

    hashes = []
    for fn in fns:
        with open(fn, "rb") as file:
            fileData = file.read()
            hashes.append(hashlib.md5(fileData).hexdigest())

    """
    Check database hash
    """
    if hashes[1] != correct:
        print("Hash check for reference database failed. " + fns[1])

    print("4. Extracting archive...")
    cmd = ["tar", "-zxf", fn[1], "-C", TEMPLATES]
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as _:
        raise ValueError("Error extracting tax4fun2 database archive.")

    return(0)



def infer_t4f2(
    sequences: pd.Series,
    seq_table: pd.DataFrame,
    threads: int = 1,
    database_mode: str = "Ref99NR",
    pct_id: float = 0.99,
) -> pd.DataFrame:
    """
    We need a prior installation of Tax4Fun2 in the R in QIIME 2 environment.
    Download https://zenodo.org/records/10035668/files/Tax4Fun2_1.1.5.tar.gz, and run
    `R CMD INSTALL Tax4Fun2_1.1.5.tar.gz`.
    The artifact of the database is available under GNU General Public License v3.0 or later.
    The original files are downloaded from: https://zenodo.org/records/10035668.
    """

    """
    Install the software and database in assets at the first try
    """
    install_t4f2()

    if threads < 1:
        raise ValueError("Thread number should be positive.")

    with TemporaryDirectory() as temp_dir:
        repseq = path.join(temp_dir, "rep_seqs.fna")

        with open(repseq, "w") as fna:
            for seqname, sequence in sequences.items():
                print(">" + str(seqname) + "\n" + str(sequence), file=fna)
        seq_table.T.to_csv(path.join(temp_dir, "seqtab.txt"), sep="\t")

        cmd = [
            "Rscript",
            path.join(TEMPLATES, "perform_tax4fun2.R"),
            temp_dir,
            "rep_seqs.fna",
            "seqtab.txt",
            path.join(TEMPLATES, "Tax4Fun2_ReferenceData_v2"),
            str(pct_id),
            str(threads),
            str(database_mode),
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
