from qiime2.plugin import SemanticType, model, ValidationError, BinaryFileFormat
import gzip
import hashlib

class FastaGzFormat(model.BinaryFileFormat):
    """
    A gzipped fasta file.

    """

    def _validate_(self, level):
        with self.open() as fh:
            if fh.peek(2)[:2] != b"\x1f\x8b":
                raise ValidationError("File is uncompressed")
        with gzip.open(str(self), "rt", "utf-8") as f:
            lines = f.readlines()            
            if lines[0][0] != '>':
                    raise ValidationError("First line of file is not a valid "
                                          "description. Descriptions must "
                                          "start with '>'")


class CnTsvGzFormat(model.BinaryFileFormat):
    """
    A gzipped TSV file containing two-column 16S copy numbers.

    """
    def _validate_(self, level):
        with self.open() as fh:
            if fh.peek(2)[:2] != b"\x1f\x8b":
                raise ValidationError("File is uncompressed")

        with gzip.open(str(self), "rt", "utf-8") as f:
            lines = f.readlines()
            first = lines[0].split("\t")
            if len(first)!=2:
                raise ValidationError("The 16S copy number table shuold be two-column layout.")
            if not first[1].replace("\n","").isnumeric():
                raise ValidationError("The second column of 16S copy number table shuold be numeric.")

class GFCnTsvGzFormat(model.BinaryFileFormat):
    """
    A gzipped TSV file containing KO (or other gene famimlies) copy numbers.

    """

    def _validate_(self, level):
        with self.open() as fh:
            if fh.peek(2)[:2] != b"\x1f\x8b":
                raise ValidationError("File is uncompressed")
        with gzip.open(str(self), 'rt', 'utf-8') as f:
            lines = f.readlines()
            first = lines[1].split("\t")
            if type(first[0]) is not str:
                raise ValidationError("The first column of gene family copy number table shuold be string.")
            for col in first[1:]:
                if not col.replace(".","").replace("\n","").isnumeric():
                    raise ValidationError("The columns of gene family copy number table shuold be numeric.")


PiphillinDatabase = SemanticType("PiphillinDatabase")


class PiphillinDatabaseFormat(model.DirectoryFormat):
    seqs = model.File(r"16S_seqs.fasta.gz", format=FastaGzFormat)
    cn = model.File(r"16S_cn.tsv.gz", format=CnTsvGzFormat)
    kocn = model.File(r"ko_copynum.tsv.gz", format=GFCnTsvGzFormat)

    def _validate_(self, level):
        ## If all the files are set, no further validation is needed.
        pass
