#!/usr/bin/env python
"""
version 0.2
minor revisions:
 -still needs lots of work, but this accounts for version 0.4
"""
#TODO implement disabling of sortTestMethodsUsing
import unittest
import os
from Bio import Entrez
import sys


sys.dont_write_bytecode = True


sys.path.append(os.path.join(
    os.path.dirname(os.path.dirname(__file__)), "pyutilsnrw"))

from pyutilsnrw.get_genomes import parse_accession_list, fetch_and_write_seqs


class GetGenomesTestCase(unittest.TestCase):
    """Tests for `get_genomes.py`."""
    def setUp(self):
        self.test_dir = os.path.join(os.path.dirname(__file__),
                                     "output_getgenomes_tests")
        self.acc_file = os.path.join(os.path.dirname(__file__),
                                     "references",
                                     "accessions.txt")
        self.grouped_file = os.path.join(os.path.dirname(__file__),
                                         "references",
                                         "grouped_output.fasta")

    def test_len_accessionlist_is_2(self):
        """Make sure this parses out 2 accessions"""
        self.assertTrue(2 == len(parse_accession_list(self.acc_file)))

    def test_entrez_canreach_ncbi(self):
        """test that ncbi can be reached ok
        """
        def make_urlcall(acc="LC075482.1"):
            try:
                return(Entrez.efetch(db="nucleotide",
                                     id=acc, rettype="fasta", retmode="text"))
            except Exception as e:
                print(e)
                return(1)

        test = make_urlcall()
        self.assertTrue(test != 1)

    def test_fetchandwriteseqs_makes_output(self):
        """Does the output have the right number of lines?
        """
        with open(self.grouped_file, "r") as doc:
            test_doc = doc.readlines()
        self.assertTrue(len(test_doc) == 56)  # determined emperically

    # def test_zcleanup_if_done(self):
    #     print(outputpath)
    #     os.remove(outputpath)  # bye bye temp file from tests. Begins with z cause of the sorting :(. see sortTestMethodsUsing if you are feeling ambitions
    #     self.assertTrue(not os.path.isfile(outputpath))


if __name__ == '__main__':
    unittest.main()
