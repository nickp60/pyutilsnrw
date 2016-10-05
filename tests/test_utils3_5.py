# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 08:57:31 2016
@author: nicholas
The Goal of this is to have a unified place to put the useful
python 3.5 functions or templates

how I got the fastq file
# seqtk sample -s 27 ~/GitHub/FA/pseudochromosome/data/20150803_Abram1/ \
    reads/3123-1_1_trimmed.fastq .0005

bam file was from a riboseed mapping; md5: 939fbf2c282091aec0dfa278b05e94ec

mapped bam was made from bam file with the following command
 samtools view -Bh -F 4 /home/nicholas/GitHub/FB/Ecoli_comparative_genomics/
    scripts/riboSeed_pipeline/batch_coli_unpaired/map/
    mapping_20160906_region_7_riboSnag/
    test_smalt4_20160906_region_7_riboSnagS.bam >
     ~/GitHub/pyutilsnrw/tests/test_mapped.sam
md5: 27944249bf064ba54576be83053e82b0

"""
__version__ = "0.0.3"
import time
import sys
import shutil
import logging
import subprocess
import os
import unittest
import hashlib
import glob
import argparse
sys.dont_write_bytecode = True

from pyutilsnrw.utils3_5 import make_output_prefix, check_installed_tools,\
    copy_file, get_ave_read_len_from_fastq, get_number_mapped,\
    extract_mapped_and_mappedmates, keep_only_first_contig, md5,\
    combine_contigs, clean_temp_dir, get_genbank_record, get_fasta_lengths,\
    file_len


# def get_args():
#     parser = argparse.ArgumentParser(
#         description="test suite for pyutilsnrw repo")
#     parser.add_argument("-k", "--keep_temps", dest='keep_temps',
#                         action="store_true",
#                         help="set if you want to inspect the output files",
#                         default=False)
#     args = parser.parse_args()
#     return(args)


@unittest.skipIf((sys.version_info[0] != 3) or (sys.version_info[1] < 5),
                 "Subprocess.call among otherthings wont run if you try this" +
                 " with less than python 3.5")
class utils3_5TestCase(unittest.TestCase):
    def setUp(self):
        # self.genbank_filename = "references/n_equitans.gbk"
        self.genbank_filename = os.path.join(os.path.dirname(__file__),
                                             str("references" + os.path.sep +
                                                 'n_equitans.gbk'))
        self.multigenbank_filename = os.path.join(os.path.dirname(__file__),
                                             str("references" + os.path.sep +
                                                 'uams1_rs.gb'))
        # curdir = os.getcwd()
        self.testdirname = os.path.join(os.path.dirname(__file__),
                                        "output_utils3_5_tests")
        self.test_fastq_file = os.path.join(os.path.dirname(__file__),
                                            str("references" + os.path.sep +
                                                'reads_reference.fastq'))
        self.test_bam_file = os.path.join(os.path.dirname(__file__),
                                          str("references" + os.path.sep +
                                              "mapping_reference.bam"))
        self.test_sam_mapped_file = os.path.join(os.path.dirname(__file__),
                                                 str("references" +
                                                     os.path.sep +
                                                     "mapping_reference_mapped.sam"))
        self.test_multifasta = os.path.join(os.path.dirname(__file__),
                                            str("references" + os.path.sep +
                                                "test_multiseqs_reference.fasta"))
        self.test_singlefasta = os.path.join(os.path.dirname(__file__),
                                             str("references" + os.path.sep +
                                                 "test_only_first_reference.fasta"))
        self.test_combined = os.path.join(os.path.dirname(__file__),
                                          str("references" + os.path.sep +
                                              "combined_contigs_reference.fa"))
        self.test_md5s_prefix = os.path.join(os.path.dirname(__file__),
                                             str("references" + os.path.sep +
                                                 "md5"))
        self.samtools_exe = "samtools"

    def test_file_len(self):
        self.assertEqual(file_len(self.test_combined), 32)

    def test_make_testing_dir(self):
        if not os.path.exists(self.testdirname):
            os.makedirs(self.testdirname)
        self.assertTrue(os.path.exists(self.testdirname))

    def test_clean_temp_dir(self):
        """ I tried to do something like
        @unittest.skipUnless(clean_temp, "temporary files were retained")
        but couldnt get the variabel to be passed through.
        """
        if not os.path.exists(os.path.join(self.testdirname, "test_subdir")):
            os.makedirs(os.path.join(self.testdirname, "test_subdir"))
        clean_temp_dir(self.testdirname)

    def test_make_output_prefix(self):
        test_prefix = make_output_prefix(self.testdirname, "utils_3.5")
        self.assertEqual(test_prefix,
                         "".join([self.testdirname, os.path.sep, "utils_3.5"]))

    def test_check_installed_tools(self):
        """is pwd on all mac/linux systems?
        #TODO replace with better passing test
        """
        check_installed_tools(["pwd"])
        # test fails properly
        with self.assertRaises(SystemExit):
            check_installed_tools(["thisisnotapathtoanactualexecutable"])

    def test_md5_strings(self):
        """ minimal md5 examples
        """
        self.assertEqual(md5("thisstringisidenticalto", string=True),
                         md5("thisstringisidenticalto", string=True))
        self.assertNotEqual(md5("thisstringisntidenticalto", string=True),
                            md5("thisstringisnotidenticalto", string=True))

    def test_references_md5(self):
        """ is this paranoia, as well as bad testing?
        """
        test_pairs = [["3ba332f8a3b5d935ea6c4e410ccdf44b",
                       "references/combined_contigs_reference.fa"],
                      ["939fbf2c282091aec0dfa278b05e94ec",
                       "references/mapping_reference.bam"],
                      ["27944249bf064ba54576be83053e82b0",
                       "references/mapping_reference_mapped.sam"],
                      ["ac80c75f468011ba11e72ddee8560b33",
                       "references/md5_a.txt"],
                      ["ac80c75f468011ba11e72ddee8560b33",
                       "references/md5_b.txt"],
                      ["92fc8592819173343a75a40874d86144",
                       "references/md5_fail.txt"],
                      ["d6b0e5b28d0b4de431f10a03042ff37b",
                       "references/reads_reference.fastq"],
                      ["40ac496ec5b221636db81ce09e04c1d9",
                       "references/test_multiseqs_reference.fasta"],
                      ["920b5c9dc69fb2a9fed50b18f3e92895",
                       "references/test_only_first_reference.fasta"]]
        for i in test_pairs:
            self.assertEqual(i[0],
                             md5(os.path.join(os.path.dirname(__file__),
                                              i[1])))

    def test_md5_files(self):
        md5_a = md5(str(self.test_md5s_prefix + "_a.txt"))
        md5_b = md5(str(self.test_md5s_prefix + "_b.txt"))
        md5_fail = md5(str(self.test_md5s_prefix + "_fail.txt"))
        self.assertEqual(md5_a, md5_b)
        self.assertNotEqual(md5_a, md5_fail)

    def test_copy_file(self):
        if not os.path.exists(self.test_fastq_file):
            raise("test file is gone!  where is test_reads.fastq ?")
        new_path = copy_file(current_file=self.test_fastq_file,
                             dest_dir=self.testdirname,
                             name="newname.fastq", overwrite=False)
        # test path to copied file is constructed properly
        self.assertEqual(new_path, os.path.join(self.testdirname, "newname.fastq"))
        # test identity of files
        self.assertEqual(md5(new_path),
                         md5(os.path.join(self.testdirname, "newname.fastq")))
        # test overwrite exit
        with self.assertRaises(SystemExit):
            new_path = copy_file(current_file=self.test_fastq_file,
                                 dest_dir=self.testdirname,
                                 name="newname.fastq", overwrite=False)
        os.remove(new_path)

    def test_get_ave_read_len_from_fastq(self):
        """this probably could/should be refined to have a better test
        """
        mean_read_len = get_ave_read_len_from_fastq(self.test_fastq_file, N=5)
        self.assertEqual(217.8, mean_read_len)

    def test_get_number_mapped(self):
        """
        """
        result = get_number_mapped(self.test_bam_file, self.samtools_exe)
        reference = "151 + 0 mapped (0.56% : N/A)"
        self.assertEqual(result, reference)

    def test_extract_mapped_and_mappedmates(self):
        """ dont trust this if  make_output_prefix test fails
        some help from PSS on SO:
        http://stackoverflow.com/questions/16874598/
            how-do-i-calculate-the-md5-checksum-of-a-file-in-python
        """
        # copy files
        test_bam_dup = copy_file(current_file=self.test_bam_file,
                                 dest_dir=self.testdirname,
                                 name="", overwrite=False)
        test_mapped_dup = copy_file(current_file=self.test_sam_mapped_file,
                                    dest_dir=self.testdirname,
                                    name="", overwrite=False)
        ref_dir = os.path.join(self.testdirname)
        prefix = make_output_prefix(output_dir=ref_dir,
                                    name="mapping_reference")
        extract_mapped_and_mappedmates(map_results_prefix=prefix,
                                       fetch_mates=False,
                                       keep_unmapped=False,
                                       samtools_exe=self.samtools_exe)
        # reference mapping md5
        mapped_md5 = "27944249bf064ba54576be83053e82b0"
        md5_returned = md5(str(prefix + "_mapped.sam"))
        # Finally compare original MD5 with freshly calculated
        self.assertEqual(mapped_md5, md5_returned)
        # delete files created
        files_created = ["_mapped.bam",
                         "_mapped.bam.bai",
                         "_mapped.sam",
                         ".bam"]
        if mapped_md5 == md5_returned:
            for i in files_created:
                os.remove(str(prefix + i))

    def test_keep_only_first_contig(self):
        """copy_file
        """
        # copy to test dir
        copy_file(current_file=self.test_multifasta,
                  dest_dir=self.testdirname,
                  name='duplicated_multifasta.fasta', overwrite=False,
                  logger=None)
        path_to_dup = os.path.join(self.testdirname, "duplicated_multifasta.fasta")
        keep_only_first_contig(path_to_dup, newname="contig1")
        self.assertEqual(md5(path_to_dup), md5(self.test_singlefasta))
        os.remove(path_to_dup)

    def test_combine_contigs(self):
        duplicated_multifasta = copy_file(current_file=self.test_multifasta,
                                          dest_dir=self.testdirname,
                                          name='multifasta_test_combine.fasta',
                                          overwrite=False,
                                          logger=None)
        for_first_contig = copy_file(current_file=self.test_multifasta,
                                     dest_dir=self.testdirname,
                                     name='single_fasta_test_combine.fasta',
                                     overwrite=False,
                                     logger=None)
        keep_only_first_contig(for_first_contig, newname="contig1")
        combined_contigs = combine_contigs(self.testdirname,
                                           pattern="*test_combine",
                                           contigs_name="combined_contigs.fa",
                                           ext=".fasta",
                                           verbose=False)
        self.assertEqual(md5(self.test_combined), md5(combined_contigs))
        for i in [duplicated_multifasta, for_first_contig, combined_contigs]:
            os.remove(i)

    def test_get_genbank_record(self):
        """Reads records from a GenBank file.
        """
        records = get_genbank_record(self.genbank_filename)
        assert(type(records) == list)
        multirecords = get_genbank_record(self.multigenbank_filename,
                                          first_only=False)
        assert(type(multirecords) == list)

    def test_get_fasta_lengths(self):
        self.assertEqual(get_fasta_lengths(self.test_singlefasta), [169])
        self.assertEqual(get_fasta_lengths(self.test_multifasta),
                         [169, 161, 159, 159, 151, 133, 128])

    def tearDown(self):
        pass


if __name__ == '__main__':
    # Commented out because processing command-line arguments in this way
    # breaks expected behaviour of unittests in Python, see:
    # https://docs.python.org/3/library/unittest.html
    # e.g. `tests.py -v` should enable verbose unit test output, but catching
    # the cmd-line like this prevents the framework seeing the arguments.
    #args = get_args()
    # curdir = os.getcwd()
    # testdirname = os.path.join(os.path.dirname(__file__),
    #                            "output_utils3_5_tests")
    # test_fastq_file = os.path.join(os.path.dirname(__file__),
    #                                str("references" + os.path.sep +
    #                                    'reads_reference.fastq'))
    # test_bam_file = os.path.join(os.path.dirname(__file__),
    #                              str("references" + os.path.sep +
    #                                  "mapping_reference.bam"))
    # test_sam_mapped_file = os.path.join(os.path.dirname(__file__),
    #                                     str("references" + os.path.sep +
    #                                         "mapping_reference_mapped.sam"))
    # test_multifasta = os.path.join(os.path.dirname(__file__),
    #                                str("references" + os.path.sep +
    #                                    "test_multiseqs_reference.fasta"))
    # test_singlefasta = os.path.join(os.path.dirname(__file__),
    #                                 str("references" + os.path.sep +
    #                                     "test_only_first_reference.fasta"))
    # test_combined = os.path.join(os.path.dirname(__file__),
    #                              str("references" + os.path.sep +
    #                                  "combined_contigs_reference.fa"))
    # test_md5s_prefix = os.path.join(os.path.dirname(__file__),
    #                                 str("references" + os.path.sep +
    #                                     "md5"))
    # samtools_exe = "samtools"
    unittest.main()
