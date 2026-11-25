from qiime2.plugin.testing import TestPluginBase
from qiime2.plugin.util import transform
from q2_pathway._infer import infer
from q2_pathway._summarize_ko import summarize
import pandas as pd
from q2_types.feature_data import DNAFASTAFormat
from qiime2 import Metadata
from tempfile import TemporaryDirectory
import os


class InferTests(TestPluginBase):
    """
    Test the functionality of infer using small database, sequence and count table
    """

    package = "q2_pathway.tests"

    def test_infer(self):
        sequence1 = transform(
            self.get_data_path("test-seqs.fasta"),
            from_type=DNAFASTAFormat,
            to_type=pd.Series,
        )
        table1 = pd.read_csv(
            self.get_data_path("test-feature-table.tsv"),
            sep="\t",
            index_col=0,
            header=0,
        )
        db1 = self.get_data_path("test-db")
        observed = infer(sequence1, table1, db1)
        self.assertEqual(observed.shape, (34, 10426))
        self.assertEqual(set(observed.index), set(table1.index.astype(str)))
        self.assertAlmostEqual(observed.at["L3S360", "K00943"], 11.0)
        self.assertAlmostEqual(observed.at["L6S20", "K03111"], 122.6666667)

    def test_infer_single_sample(self):
        sequence1 = transform(
            self.get_data_path("test-seqs.fasta"),
            from_type=DNAFASTAFormat,
            to_type=pd.Series,
        )
        table1 = pd.read_csv(
            self.get_data_path("test-feature-table.tsv"),
            sep="\t",
            index_col=0,
            header=0,
        )
        single_sample = table1.loc[["L3S360"]]
        db1 = self.get_data_path("test-db")

        observed = infer(sequence1, single_sample, db1)

        self.assertEqual(observed.shape, (1, 10426))
        self.assertEqual(observed.index.tolist(), single_sample.index.tolist())
        self.assertAlmostEqual(observed.at["L3S360", "K00943"], 11.0)

    def test_infer_full(self):
        sequence1 = transform(
            self.get_data_path("test-seqs.fasta"),
            from_type=DNAFASTAFormat,
            to_type=pd.Series,
        )
        table1 = pd.read_csv(
            self.get_data_path("test-feature-table.tsv"),
            sep="\t",
            index_col=0,
            header=0,
        )
        db1 = self.get_data_path("test-db")
        observed = infer(sequence1, table1, db1, full=True)
        self.assertEqual(observed.shape, (34, 31278))
        self.assertEqual(set(observed.index), set(table1.index.astype(str)))
        l3s360 = observed.loc["L3S360"]
        collapsed = l3s360.groupby(
            l3s360.index.str.split("|", regex=False).str[0]
        ).sum()
        self.assertAlmostEqual(collapsed.at["K00943"], 11.0)
        self.assertAlmostEqual(collapsed.at["K03111"], 55.0)

        filtered = infer(
            sequence1,
            table1,
            db1,
            full=True,
            full_id="RS_GCF_004923255.1",
        )
        self.assertEqual(filtered.shape, (34, 10426))
        self.assertTrue(
            all(
                column.endswith("|RS_GCF_004923255.1")
                for column in filtered.columns
            )
        )

        filtered_ko = infer(
            sequence1,
            table1,
            db1,
            full=True,
            full_ko_id="K00943",
        )
        self.assertGreater(filtered_ko.shape[1], 0)
        self.assertEqual(set(filtered_ko.index), set(table1.index.astype(str)))
        self.assertTrue(
            all(column.startswith("K00943|") for column in filtered_ko.columns)
        )
        expected_ko = observed.loc[
            :, [column for column in observed.columns if column.startswith("K00943|")]
        ]
        pd.testing.assert_frame_equal(
            filtered_ko.sort_index(axis=1),
            expected_ko.sort_index(axis=1),
            check_dtype=False,
        )


class SummarizeTests(TestPluginBase):
    """
    Test the functionality of summarize
    """

    package = "q2_pathway.tests"

    def test_summarize(self):
        sequence1 = transform(
            self.get_data_path("test-seqs.fasta"),
            from_type=DNAFASTAFormat,
            to_type=pd.Series,
        )
        table1 = pd.read_csv(
            self.get_data_path("test-feature-table.tsv"),
            sep="\t",
            index_col=0,
            header=0,
        )
        db1 = self.get_data_path("test-db")
        test = infer(sequence1, table1, db1)
        test2 = infer(sequence1, table1, db1)
        meta = Metadata.load(self.get_data_path("test-metadata.tsv"))
        with TemporaryDirectory() as temp_dir:
            summarize(
                temp_dir,
                [test, test2],
                meta,
                tables_name=["Piphillin-A", "Piphillin-B"],
                table_roles=["prediction", "prediction"],
                algorithm_families=["Piphillin", "Piphillin"],
            )
            self.assertTrue(
                os.path.exists(
                    os.path.join(temp_dir, "abundance_agreement_per_feature.tsv")
                )
            )
            self.assertTrue(
                os.path.exists(
                    os.path.join(temp_dir, "family_balanced_summary_per_feature.tsv")
                )
            )
