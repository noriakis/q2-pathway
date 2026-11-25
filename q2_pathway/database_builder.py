from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import shutil


REQUIRED_OUTPUTS = {
    "sequences": "16S_seqs.fasta.gz",
    "copy_numbers": "16S_cn.tsv.gz",
    "ko_copy_numbers": "ko_copynum.tsv.gz",
}


def _is_gzip(path: Path) -> bool:
    with path.open("rb") as handle:
        return handle.read(2) == b"\x1f\x8b"


def _open_text(path: Path):
    if _is_gzip(path):
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("rt", encoding="utf-8", newline="")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _copy_as_gzip(source: Path, target: Path) -> None:
    if _is_gzip(source):
        shutil.copyfile(source, target)
        return
    with source.open("rb") as input_handle, target.open("wb") as raw_output:
        with gzip.GzipFile(
            filename="", mode="wb", compresslevel=6, fileobj=raw_output, mtime=0
        ) as output_handle:
            shutil.copyfileobj(input_handle, output_handle)


def _read_fasta_ids(path: Path) -> list[str]:
    identifiers = []
    with _open_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith(">"):
                identifier = line[1:].strip().split()[0]
                if not identifier:
                    raise ValueError(f"Empty FASTA identifier at line {line_number}.")
                identifiers.append(identifier)
    if not identifiers:
        raise ValueError("No FASTA records were found.")
    if len(identifiers) != len(set(identifiers)):
        raise ValueError("16S FASTA identifiers must be unique.")
    return identifiers


def _read_copy_number_ids(path: Path) -> list[str]:
    identifiers = []
    with _open_text(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        for line_number, row in enumerate(reader, start=1):
            if len(row) != 2:
                raise ValueError(
                    f"16S copy-number row {line_number} must have two columns."
                )
            try:
                value = float(row[1])
            except ValueError as error:
                raise ValueError(
                    f"Invalid 16S copy number at row {line_number}: {row[1]}"
                ) from error
            if value <= 0:
                raise ValueError(
                    f"16S copy number must be positive at row {line_number}."
                )
            identifiers.append(row[0])
    if len(identifiers) != len(set(identifiers)):
        raise ValueError("16S copy-number identifiers must be unique.")
    return identifiers


def _read_ko_table(path: Path) -> tuple[list[str], list[str]]:
    reference_ids = []
    with _open_text(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration as error:
            raise ValueError("KO copy-number table is empty.") from error
        if len(header) < 2:
            raise ValueError("KO copy-number table must contain KOs as columns.")
        kos = [value.strip() for value in header[1:]]
        if len(kos) != len(set(kos)):
            raise ValueError("KO identifiers in the header must be unique.")
        for line_number, row in enumerate(reader, start=2):
            if len(row) != len(header):
                raise ValueError(
                    f"KO copy-number row {line_number} has {len(row)} columns; "
                    f"expected {len(header)}."
                )
            try:
                values = [float(value) for value in row[1:]]
            except ValueError as error:
                raise ValueError(
                    f"Non-numeric KO copy number at row {line_number}."
                ) from error
            if any(value < 0 for value in values):
                raise ValueError(f"Negative KO copy number at row {line_number}.")
            reference_ids.append(row[0])
    if len(reference_ids) != len(set(reference_ids)):
        raise ValueError("KO table reference identifiers must be unique.")
    return reference_ids, kos


def build_database(
    sequences: str | Path,
    copy_numbers: str | Path,
    ko_copy_numbers: str | Path,
    output_dir: str | Path,
    *,
    database_name: str,
    source_version: str,
    source_url: str | None = None,
) -> Path:

    sources = {
        "sequences": Path(sequences).resolve(),
        "copy_numbers": Path(copy_numbers).resolve(),
        "ko_copy_numbers": Path(ko_copy_numbers).resolve(),
    }
    for label, source in sources.items():
        if not source.is_file():
            raise FileNotFoundError(f"{label} input not found: {source}")
    target = Path(output_dir).resolve()
    if target.exists() and any(target.iterdir()):
        raise FileExistsError(f"Output directory is not empty: {target}")
    target.mkdir(parents=True, exist_ok=True)

    fasta_ids = _read_fasta_ids(sources["sequences"])
    cn_ids = _read_copy_number_ids(sources["copy_numbers"])
    ko_reference_ids, kos = _read_ko_table(sources["ko_copy_numbers"])
    fasta_cn_overlap = len(set(fasta_ids) & set(cn_ids))
    cn_ko_overlap = len(set(cn_ids) & set(ko_reference_ids))
    if fasta_cn_overlap == 0:
        raise ValueError("No reference IDs overlap between FASTA and 16S copy numbers.")
    if cn_ko_overlap == 0:
        raise ValueError(
            "No reference IDs overlap between 16S and KO copy-number tables."
        )

    output_hashes = {}
    source_hashes = {}
    for label, filename in REQUIRED_OUTPUTS.items():
        source_hashes[label] = _sha256(sources[label])
        destination = target / filename
        _copy_as_gzip(sources[label], destination)
        output_hashes[label] = _sha256(destination)

    manifest = {
        "schema_version": 1,
        "database_name": database_name,
        "source_version": source_version,
        "source_url": source_url,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "files": {
            label: {
                "source": str(sources[label]),
                "source_sha256": source_hashes[label],
                "artifact_filename": REQUIRED_OUTPUTS[label],
                "artifact_sha256": output_hashes[label],
            }
            for label in REQUIRED_OUTPUTS
        },
        "statistics": {
            "fasta_sequences": len(fasta_ids),
            "copy_number_rows": len(cn_ids),
            "ko_reference_rows": len(ko_reference_ids),
            "distinct_kos": len(kos),
            "fasta_copy_number_id_overlap": fasta_cn_overlap,
            "copy_number_ko_reference_id_overlap": cn_ko_overlap,
        },
    }
    with (target / "database_manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return target


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build a validated q2-pathway Piphillin database directory."
    )
    parser.add_argument("--sequences", required=True)
    parser.add_argument("--copy-numbers", required=True)
    parser.add_argument("--ko-copy-numbers", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--database-name", required=True)
    parser.add_argument("--source-version", required=True)
    parser.add_argument("--source-url")
    args = parser.parse_args()
    output = build_database(
        args.sequences,
        args.copy_numbers,
        args.ko_copy_numbers,
        args.output_dir,
        database_name=args.database_name,
        source_version=args.source_version,
        source_url=args.source_url,
    )
    print(output)


if __name__ == "__main__":
    main()
