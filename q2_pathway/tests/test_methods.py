from qiime2.plugin.testing import TestPluginBase
from qiime2.plugin.util import transform
from q2_pathway._infer import infer
from q2_pathway._summarize import summarize
import pandas as pd
from q2_types.feature_data import DNAFASTAFormat
from qiime2 import Metadata
from tempfile import TemporaryDirectory

class InferTests(TestPluginBase):
    """
    Test the functionality of infer using small database, sequence and count table
    """
    package = 'q2_pathway.tests'
    def test_infer(self):
        sequence1 = transform(
            self.get_data_path('test-seqs.fasta'),
            from_type=DNAFASTAFormat,
            to_type=pd.Series)
        table1 = pd.read_csv(self.get_data_path("test-feature-table.tsv"),
            sep="\t", index_col=0, header=0) 
        db1 = self.get_data_path("test-db")
        test = infer(sequence1, table1, db1)



class SummarizeTests(TestPluginBase):
    """
    Test the functionality of summarize
    """
    package = 'q2_pathway.tests'
    def test_summarize(self):
        sequence1 = transform(
            self.get_data_path('test-seqs.fasta'),
            from_type=DNAFASTAFormat,
            to_type=pd.Series)
        table1 = pd.read_csv(self.get_data_path("test-feature-table.tsv"),
            sep="\t", index_col=0, header=0) 
        db1 = self.get_data_path("test-db")
        test = infer(sequence1, table1, db1)
        test2 = infer(sequence1, table1, db1)
        meta = Metadata.load(self.get_data_path('test-metadata.tsv'))
        with TemporaryDirectory() as temp_dir:
            summarize(temp_dir, [test, test2], meta)
