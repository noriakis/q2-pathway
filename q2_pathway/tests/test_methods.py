from qiime2.plugin.testing import TestPluginBase
from qiime2.plugin.util import transform
from q2_pathway._infer import infer
import pandas as pd
from q2_types.feature_data import DNAFASTAFormat

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

