from os import path
import os
import pandas as pd
import subprocess
from tempfile import TemporaryDirectory
import requests
import hashlib
from pathlib import Path
import tarfile
import argparse

TEMPLATES = str(Path(__file__).resolve().parent / "assets")


def _tokenize_sample_index(
    seq_table: pd.DataFrame,
) -> tuple[pd.DataFrame, dict[str, str]]:
    """
    for the safe read in R environment
    """
    frame = seq_table.copy()
    sample_ids = frame.index.astype(str)
    if sample_ids.has_duplicates:
        duplicates = sample_ids[sample_ids.duplicated()].unique().tolist()[:5]
        raise ValueError(f"Duplicate sample IDs: {duplicates}")
    tokens = [f"q2sample{i:08d}" for i in range(len(sample_ids))]
    frame.index = tokens
    return frame, dict(zip(tokens, sample_ids))


def _restore_sample_index(
    table: pd.DataFrame, token_to_sample: dict[str, str]
) -> pd.DataFrame:
    result = table.copy()
    result.index = result.index.astype(str)
    unknown = sorted(set(result.index) - set(token_to_sample))
    if unknown:
        raise ValueError(
            "Inference output contains unrecognized temporary sample IDs: "
            + ", ".join(unknown[:5])
        )
    result.index = [token_to_sample[token] for token in result.index]
    return result


def infer(
    sequences: pd.Series,
    seq_table: pd.DataFrame,
    database: str,
    threads: int = 1,
    full: bool = False,
    full_id: str = None,
    full_ko_id: str = None,
    pct_id: float = 0.99,
    max_accepts: int = 0,
    max_rejects: int = 0,
) -> pd.DataFrame:
    """Infer KO profiles with Piphillin.

    ``full_id`` filters stratified output by reference-genome ID, whereas
    ``full_ko_id`` filters by KO ID and avoids expanding unrelated KOs.
    """
    reference_sequences = path.join(database, "16S_seqs.fasta.gz")
    cn_table = path.join(database, "ko_copynum.tsv.gz")
    cn_16s_table = path.join(database, "16S_cn.tsv.gz")
    tokenized_table, token_to_sample = _tokenize_sample_index(seq_table)

    with TemporaryDirectory() as temp_dir:
        repseq = path.join(temp_dir, "rep_seqs.fna")

        with open(repseq, "w") as fna:
            for seqname, sequence in sequences.items():
                print(">" + str(seqname) + "\n" + str(sequence), file=fna)
        tokenized_table.T.to_csv(path.join(temp_dir, "seqtab.txt"), sep="\t")
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
            str(max_accepts),
            "--maxrejects",
            str(max_rejects),
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
            str(full_ko_id),
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
            return _restore_sample_index(ko.T, token_to_sample)
            ## This will output tax1_K00001 type index
        else:
            ko = pd.read_csv(path.join(temp_dir, "ko_table.tsv"), sep="\t", header=0)
            return _restore_sample_index(ko.T, token_to_sample)


TAX4FUN2_DB_URL = (
    "https://zenodo.org/records/10035668/files/Tax4Fun2_ReferenceData_v2.tar.gz"
)
TAX4FUN2_DB_MD5 = "5c2edff2f597e5f061f79a5782a27567"


def _default_t4f2_resource_dir() -> Path:
    configured = os.environ.get("Q2_PATHWAY_TAX4FUN2_DIR")
    if configured:
        return Path(configured).expanduser().resolve()
    return Path.home() / ".cache" / "q2-pathway" / "tax4fun2"


def _md5(path_to_file: Path) -> str:
    digest = hashlib.md5()
    with path_to_file.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def install_t4f2(resource_dir: str = None) -> Path:
    """Download verified Tax4Fun2 reference data to a user-writable cache.

    The installed Python/R package directory is never modified.  Tax4Fun2 itself
    must already be installed in the active R environment.
    """
    target = (
        Path(resource_dir).expanduser().resolve()
        if resource_dir is not None
        else _default_t4f2_resource_dir()
    )
    target.mkdir(parents=True, exist_ok=True)
    database_dir = target / "Tax4Fun2_ReferenceData_v2"
    if database_dir.is_dir():
        return database_dir

    archive = target / "Tax4Fun2_ReferenceData_v2.tar.gz"
    if not archive.is_file():
        with requests.get(TAX4FUN2_DB_URL, stream=True, timeout=(10, 120)) as response:
            response.raise_for_status()
            with archive.open("wb") as handle:
                for chunk in response.iter_content(chunk_size=1024 * 1024):
                    if chunk:
                        handle.write(chunk)
    digest = _md5(archive)
    if digest != TAX4FUN2_DB_MD5:
        archive.unlink(missing_ok=True)
        raise ValueError("Tax4Fun2 reference database checksum verification failed.")

    with tarfile.open(archive, "r:gz") as tar:
        root = target.resolve()
        for member in tar.getmembers():
            if member.issym() or member.islnk():
                raise ValueError("Tax4Fun2 database archive contains a link.")
            destination = (target / member.name).resolve()
            if root not in destination.parents and destination != root:
                raise ValueError("Unsafe path in Tax4Fun2 database archive.")
        tar.extractall(target)
    if not database_dir.is_dir():
        raise ValueError("Tax4Fun2 archive did not contain the expected directory.")
    return database_dir


def install_t4f2_cli() -> None:
    parser = argparse.ArgumentParser(
        description="Install verified Tax4Fun2 reference data to a user cache."
    )
    parser.add_argument("--resource-dir")
    args = parser.parse_args()
    print(install_t4f2(args.resource_dir))


def infer_t4f2(
    sequences: pd.Series,
    seq_table: pd.DataFrame,
    threads: int = 1,
    database_mode: str = "Ref99NR",
    pct_id: float = 0.99,
    resource_dir: str = None,
    auto_install: bool = False,
) -> pd.DataFrame:
    """
    We need a prior installation of Tax4Fun2 in the R in QIIME 2 environment.
    Download https://zenodo.org/records/10035668/files/Tax4Fun2_1.1.5.tar.gz, and run
    `R CMD INSTALL Tax4Fun2_1.1.5.tar.gz`.
    The artifact of the database is available under GNU General Public License v3.0 or later.
    The original files are downloaded from: https://zenodo.org/records/10035668.
    """

    resource_root = (
        Path(resource_dir).expanduser().resolve()
        if resource_dir is not None
        else _default_t4f2_resource_dir()
    )
    database_dir = resource_root / "Tax4Fun2_ReferenceData_v2"
    if not database_dir.is_dir():
        if auto_install:
            database_dir = install_t4f2(str(resource_root))
        else:
            raise ValueError(
                "Tax4Fun2 reference data were not found. Run "
                "q2-pathway-install-tax4fun2-data or pass resource_dir."
            )

    if threads < 1:
        raise ValueError("Thread number should be positive.")

    tokenized_table, token_to_sample = _tokenize_sample_index(seq_table)

    with TemporaryDirectory() as temp_dir:
        repseq = path.join(temp_dir, "rep_seqs.fna")

        with open(repseq, "w") as fna:
            for seqname, sequence in sequences.items():
                print(">" + str(seqname) + "\n" + str(sequence), file=fna)
        tokenized_table.T.to_csv(path.join(temp_dir, "seqtab.txt"), sep="\t")

        cmd = [
            "Rscript",
            path.join(TEMPLATES, "perform_tax4fun2.R"),
            temp_dir,
            "rep_seqs.fna",
            "seqtab.txt",
            str(database_dir),
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
        return _restore_sample_index(ko.T, token_to_sample)
