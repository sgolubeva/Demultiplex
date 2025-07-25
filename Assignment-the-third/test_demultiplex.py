import unittest
import os
import gzip
from demultiplex import rev_compliment, process_record, write_file

class TestDemultiplex(unittest.TestCase):

    def test_rev_compliment(self):
        self.assertEqual(rev_compliment("AGCT"), "AGCT")
        self.assertEqual(rev_compliment("GATTACA"), "TGTAATC")
        self.assertEqual(rev_compliment("NNNN"), "NNNN")

    def test_process_record(self):
        with open("test.fq", "w") as f:
            f.write("@header\n")
            f.write("AGCT\n")
            f.write("+\n")
            f.write("!!!!\n")
        with open("test.fq", "r") as f:
            self.assertEqual(process_record(f), ["@header", "AGCT", "+", "!!!!"])
        os.remove("test.fq")

    def test_write_file(self):
        r1 = ["@header1", "AGCT", "+", "!!!!"]
        r4 = ["@header4", "TCGA", "+", "!!!!"]
        with open("test_R1.fq", "w") as f1, open("test_R4.fq", "w") as f4:
            write_file((f1, f4), "AGCT", "TCGA", r1, r4)

        with open("test_R1.fq", "r") as f1, open("test_R4.fq", "r") as f4:
            self.assertEqual(f1.read(), "@header1 AGCT-TCGA\nAGCT\n+\n!!!!\n")
            self.assertEqual(f4.read(), "@header4 AGCT-TCGA\nTCGA\n+\n!!!!\n")

        os.remove("test_R1.fq")
        os.remove("test_R4.fq")

    def test_main_logic(self):
        # Create dummy input files
        with gzip.open("test_R1.fq.gz", "wt") as f1, \
             gzip.open("test_R2.fq.gz", "wt") as f2, \
             gzip.open("test_R3.fq.gz", "wt") as f3, \
             gzip.open("test_R4.fq.gz", "wt") as f4:

            # Matched
            f1.write("@h1\nAGCT\n+\n!!!!\n")
            f2.write("@h2\nGCTT\n+\n!!!!\n")
            f3.write("@h3\nAAGC\n+\n!!!!\n")
            f4.write("@h4\nTCGA\n+\n!!!!\n")

            # Index-hopped
            f1.write("@h1\nAAAA\n+\n!!!!\n")
            f2.write("@h2\nGGGG\n+\n!!!!\n")
            f3.write("@h3\nCCCC\n+\n!!!!\n")
            f4.write("@h4\nTTTT\n+\n!!!!\n")

            # Unknown
            f1.write("@h1\nNNNN\n+\n!!!!\n")
            f2.write("@h2\nNNNN\n+\n!!!!\n")
            f3.write("@h3\nNNNN\n+\n!!!!\n")
            f4.write("@h4\nNNNN\n+\n!!!!\n")

        with open("indexes.txt", "w") as f:
            f.write("GCTT\t\n")
            f.write("GGGG\t\n")

        # Run the script
        os.system("python Assignment-the-third/demultiplex.py -r1 test_R1.fq.gz -r2 test_R2.fq.gz -r3 test_R3.fq.gz -r4 test_R4.fq.gz")

        # Check output files
        self.assertTrue(os.path.exists("R1_GCTT_matching.fq"))
        self.assertTrue(os.path.exists("R1_index_hopped.fq"))
        self.assertTrue(os.path.exists("R1_unknown.fq"))

        # Clean up files
        for f in os.listdir("."):
            if f.endswith(".fq") or f.endswith(".fq.gz") or f.endswith(".tsv") or f.endswith(".txt"):
                os.remove(f)


if __name__ == '__main__':
    unittest.main()
