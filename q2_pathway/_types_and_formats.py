from qiime2.plugin import SemanticType, model, ValidationError
import gzip
import json


class FastaGzFormat(model.BinaryFileFormat):
    """
    A gzipped fasta file.

    """

    def _validate_(self, level):
        with self.open() as fh:
            if fh.peek(2)[:2] != b"\x1f\x8b":
                raise ValidationError("File is uncompressed")
        with gzip.open(str(self), "rt", encoding="utf-8") as f:
            first = f.readline()
            if not first.startswith(">"):
                raise ValidationError(
                    "First FASTA line must be a description starting with '>'."
                )


class CnTsvGzFormat(model.BinaryFileFormat):
    """
    A gzipped TSV file containing two-column 16S copy numbers.

    """

    def _validate_(self, level):
        with self.open() as fh:
            if fh.peek(2)[:2] != b"\x1f\x8b":
                raise ValidationError("File is uncompressed")

        with gzip.open(str(self), "rt", encoding="utf-8") as f:
            first = f.readline().rstrip("\n").split("\t")
            if len(first) != 2:
                raise ValidationError(
                    "The 16S copy-number table must have two columns."
                )
            try:
                value = float(first[1])
            except ValueError as error:
                raise ValidationError(
                    "The second 16S copy-number column must be numeric."
                ) from error
            if value <= 0:
                raise ValidationError("16S copy numbers must be positive.")


class GFCnTsvGzFormat(model.BinaryFileFormat):
    """
    A gzipped TSV file containing KO (or other gene famimlies) copy numbers.

    """

    def _validate_(self, level):
        with self.open() as fh:
            if fh.peek(2)[:2] != b"\x1f\x8b":
                raise ValidationError("File is uncompressed")
        with gzip.open(str(self), "rt", encoding="utf-8") as f:
            header = f.readline().rstrip("\n").split("\t")
            first = f.readline().rstrip("\n").split("\t")
            if len(header) < 2 or len(first) != len(header):
                raise ValidationError(
                    "The KO copy-number table must have a reference-ID column "
                    "followed by numeric KO columns."
                )
            try:
                values = [float(value) for value in first[1:]]
            except ValueError as error:
                raise ValidationError("KO copy numbers must be numeric.") from error
            if any(value < 0 for value in values):
                raise ValidationError("KO copy numbers must be non-negative.")


class DatabaseManifestFormat(model.TextFileFormat):
    """JSON provenance and coverage summary created by the database builder."""

    def _validate_(self, level):
        try:
            with self.open() as handle:
                manifest = json.load(handle)
        except (OSError, json.JSONDecodeError) as error:
            raise ValidationError("Database manifest is not valid JSON.") from error
        required = {"schema_version", "database_name", "source_version", "files"}
        missing = required - set(manifest)
        if missing:
            raise ValidationError(
                "Database manifest is missing: " + ", ".join(sorted(missing))
            )


PiphillinDatabase = SemanticType("PiphillinDatabase")


class PiphillinDatabaseFormat(model.DirectoryFormat):
    seqs = model.File(r"16S_seqs.fasta.gz", format=FastaGzFormat)
    cn = model.File(r"16S_cn.tsv.gz", format=CnTsvGzFormat)
    kocn = model.File(r"ko_copynum.tsv.gz", format=GFCnTsvGzFormat)
    manifest = model.File(
        r"database_manifest.json", format=DatabaseManifestFormat, optional=True
    )

    def _validate_(self, level):
        ## If all the files are set, no further validation is needed.
        pass
