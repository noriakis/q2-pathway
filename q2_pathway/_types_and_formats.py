from qiime2.plugin import SemanticType, model, ValidationError, BinaryFileFormat
import hashlib

T4F2Database = SemanticType("T4F2Database")


class T4F2DatabaseFileFormat(BinaryFileFormat):
    def _validate_(self, level):
        correct = "5c2edff2f597e5f061f79a5782a27567"
        filehash = hashlib.md5(open(str(self), "rb").read()).hexdigest()
        if correct == filehash:
            pass
        else:
            raise ValidationError("Seems like it is not the correct artifact")


T4F2DatabaseFormat = model.SingleFileDirectoryFormat(
    "T4F2DatabaseFormat", "Tax4Fun2_ReferenceData_v2.tar.gz", T4F2DatabaseFileFormat
)


class FastaGzFormat(model.BinaryFileFormat):
    """
    A gzipped fasta file.

    """

    def _validate_(self, level):
        with self.open() as fh:
            if fh.peek(2)[:2] != b"\x1f\x8b":
                raise ValidationError("File is uncompressed")


class TsvGzFormat(model.BinaryFileFormat):
    """
    A gzipped TSV file.

    """

    def _validate_(self, level):
        with self.open() as fh:
            if fh.peek(2)[:2] != b"\x1f\x8b":
                raise ValidationError("File is uncompressed")


PiphillinDatabase = SemanticType("PiphillinDatabase")


class PiphillinDatabaseFormat(model.DirectoryFormat):
    seqs = model.File(r"16S_seqs.fasta.gz", format=FastaGzFormat)
    cn = model.File(r"16S_cn.tsv.gz", format=TsvGzFormat)
    kocn = model.File(r"ko_copynum.tsv.gz", format=TsvGzFormat)

    def _validate_(self, level):
        pass
