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
import sys
import os
import unittest
import logging

from pyutilsnrw.utils3_5 import make_output_prefix, check_installed_tools,\
    copy_file, get_ave_read_len_from_fastq, get_number_mapped,\
    extract_mapped_and_mappedmates, keep_only_first_contig, md5,\
    combine_contigs, clean_temp_dir, get_genbank_record, get_fasta_lengths,\
    file_len, multisplit, check_version_from_init, check_version_from_cmd

sys.dont_write_bytecode = True

logger = logging


@unittest.skipIf((sys.version_info[0] != 3) or (sys.version_info[1] < 5),
                 "Subprocess.call among otherthings wont run if you try this" +
                 " with less than python 3.5")
class utils3_5TestCase(unittest.TestCase):
    """ Test the utils3_5 collection of functions
    """
    def setUp(self):
        self.genbank_filename = os.path.join(os.path.dirname(__file__),
                                             str("references" + os.path.sep +
                                                 'n_equitans.gbk'))
        self.multigenbank_filename = os.path.join(os.path.dirname(__file__),
                                                  str("references" +
                                                      os.path.sep +
                                                      'uams1_rs.gb'))
        self.testdirname = os.path.join(os.path.dirname(__file__),
                                        "output_utils3_5_tests")
        self.test_fastq_file = os.path.join(os.path.dirname(__file__),
                                            str("references" + os.path.sep +
                                                'reads_reference.fastq'))
        self.test_empty_file = os.path.join(os.path.dirname(__file__),
                                            str("references" + os.path.sep +
                                                'empty.fasta'))
        self.test_bam_file = os.path.join(os.path.dirname(__file__),
                                          str("references" + os.path.sep +
                                              "mapping_reference.bam"))
        self.test_sam_mapped_file = os.path.join(
            os.path.dirname(__file__),
            str("references" +
                os.path.sep +
                "mapping_reference_mapped.sam"))
        self.test_multifasta = os.path.join(
            os.path.dirname(__file__),
            str("references" + os.path.sep +
                "test_multiseqs_reference.fasta"))
        self.test_singlefasta = os.path.join(
            os.path.dirname(__file__),
            str("references" + os.path.sep +
                "test_only_first_reference.fasta"))
        self.test_combined = os.path.join(
            os.path.dirname(__file__),
            str("references" + os.path.sep +
                "combined_contigs_reference.fa"))
        self.test_md5s_prefix = os.path.join(
            os.path.dirname(__file__),
            str("references" + os.path.sep + "md5"))
        self.samtools_exe = "samtools"

    def test_check_init(self):
        """this checks the version number in an init file
        """
        initf = os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            "pyutilsnrw", "__init__.py")
        print(initf)
        self.assertTrue(isinstance(check_version_from_init(
            init_file=initf, min_version="0.0.0"), str))
        with self.assertRaises(FileNotFoundError):
            check_version_from_init(
                init_file="notheinitfile", min_version="0.0.0")
        with self.assertRaises(ValueError):
            check_version_from_init(init_file=initf, min_version="10.0.0")

    def test_cmd_version(self):
        """This isnt a great test, as it has to be changed when I
        """
        # samtools_verison = check_version_from_cmd(
        #     cmd='samtools',
        #     line=3,
        #     pattern=r"\s*Version: (?P<version>[^(]+)", where='stderr',
        #                        min_version="0.0.0")
        pip_version = check_version_from_cmd(
            exe='pip',
            cmd=' --version',
            line=1,
            pattern=r"pip (?P<version>[^from]+)",
            # logger=logger,
            where='stdout',
            min_version="0.0.0")
        self.assertTrue(pip_version > '7.0.0')

    def test_file_len(self):
        """ test against file of known length
        """
        self.assertEqual(file_len(self.test_combined), 32)

    def test_make_testing_dir(self):
        """ make a place for the temporary files
        """
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
        """check against known install (python) and a non_executable
        """
        check_installed_tools("python")
        # test fails properly
        nonex = "thisisnotapathtoanactualexecutable"
        with self.assertRaises(SystemExit):
            check_installed_tools(nonex)
        self.assertFalse(check_installed_tools(nonex,
                                               hard=False))

    def test_md5_strings(self):
        """ minimal md5 examples with strings
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
        """ test file contests identiy
        """
        md5_a = md5(str(self.test_md5s_prefix + "_a.txt"))
        md5_b = md5(str(self.test_md5s_prefix + "_b.txt"))
        md5_fail = md5(str(self.test_md5s_prefix + "_fail.txt"))
        self.assertEqual(md5_a, md5_b)
        self.assertNotEqual(md5_a, md5_fail)

    def test_copy_file(self):
        """ make sure copied files are the same
        """
        if not os.path.exists(self.test_fastq_file):
            raise FileNotFoundError("test file is gone!  " +
                                    "where is test_reads.fastq ?")
        new_path = copy_file(current_file=self.test_fastq_file,
                             dest_dir=self.testdirname,
                             name="newname.fastq", overwrite=False)
        # test path to copied file is constructed properly
        self.assertEqual(new_path, os.path.join(self.testdirname,
                                                "newname.fastq"))
        # test identity of files
        self.assertEqual(md5(new_path),
                         md5(os.path.join(self.testdirname, "newname.fastq")))
        # test overwrite exit
        with self.assertRaises(SystemExit):
            new_path = copy_file(current_file=self.test_fastq_file,
                                 dest_dir=self.testdirname,
                                 name="newname.fastq", overwrite=False)
        os.remove(new_path)

    def test_average_read_len(self):
        """ tests get_ave_read_len_from_fastq
        this probably could/should be refined to have a better test
        """
        mean_read_len = get_ave_read_len_from_fastq(self.test_fastq_file, N=5)
        self.assertEqual(217.8, mean_read_len)

    def test_get_number_mapped(self):
        """ checks flagstat
        """
        result = get_number_mapped(self.test_bam_file, self.samtools_exe)
        reference = "151 + 0 mapped (0.56% : N/A)"
        self.assertEqual(result, reference)

    def test_extraction(self):
        """ tests extract_mapped_and_mappedmates
        dont trust this if  make_output_prefix test fails
        some help from PSS on SO:
        http://stackoverflow.com/questions/16874598/
            how-do-i-calculate-the-md5-checksum-of-a-file-in-python
        """
        # copy files
        test_bam_dup = copy_file(current_file=self.test_bam_file,
                                 dest_dir=self.testdirname,
                                 name="", overwrite=False)
        self.assertTrue(os.path.exists(test_bam_dup))
        test_mapped_dup = copy_file(current_file=self.test_sam_mapped_file,
                                    dest_dir=self.testdirname,
                                    name="", overwrite=False)
        self.assertTrue(os.path.exists(test_mapped_dup))
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
        path_to_dup = os.path.join(self.testdirname,
                                   "duplicated_multifasta.fasta")
        keep_only_first_contig(path_to_dup, newname="contig1")
        self.assertEqual(md5(path_to_dup), md5(self.test_singlefasta))
        os.remove(path_to_dup)

    def test_combine_contigs(self):
        """ compine two files, compare lengths, check construction
        """
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
        assert isinstance(records, list)
        multirecords = get_genbank_record(self.multigenbank_filename,
                                          first_only=False)
        assert isinstance(multirecords, list)

    def test_get_fasta_lengths(self):
        """ get the lengths of the multifasta entries
        """
        self.assertEqual(get_fasta_lengths(self.test_singlefasta), [169])
        self.assertEqual(get_fasta_lengths(self.test_multifasta),
                         [169, 161, 159, 159, 151, 133, 128])

    def test_multisplit(self):
        """ split a string that has multiple delimiters
        """
        test_string = "look_this+is+a locus_that_is+multi-delimited"
        list_of_things = multisplit(["-", "_", "+", " "], test_string)
        test_other_string = "look_this+is+a\faillocus_that_is+multi-delimited"
        list_of_other_things = multisplit(["-", "_", "+", " "],
                                          test_other_string)
        self.assertEqual(list_of_things, ["look", "this", "is", "a", "locus",
                                          "that", "is", "multi", "delimited"])
        self.assertNotEqual(list_of_other_things, ["look", "this", "is", "a",
                                                   "locus", "that", "is",
                                                   "multi", "delimited"])

    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()

###
### below are the functions to still write tests for
###
# def run_quast(contigs, output, quast_exe, ref="", threads=1, logger=None):
#     """Reference is optional. This is, honestly, a pretty dumb feature
#     requires sys, subprocess, (system install of quast)
#     """
# def setup_protein_blast(input_file, input_type="fasta", dbtype="prot",
#                         title="blastdb", out="blastdb",
#                         makeblastdb_exe='', logger=None):
#     """
#     This runs make blast db with the given parameters
#     requires logging, os, subprocess, shutil
#     """
#     if makeblastdb_exe == '':
#         makeblastdb_exe = shutil.which("makeblastdb")
#     makedbcmd = str("{0} -in {1} -input_type {2} -dbtype {3} " +
#                     "-title {4} -out {5}").format(makeblastdb_exe,
#                                                   input_file,
#                                                   input_type,
#                                                   dbtype, title, out)
#     if logger:
#         logger.info("Making blast db: {0}".format(makedbcmd))
#     try:
#         subprocess.run(makedbcmd, shell=sys.platform != "win32",
#                        stdout=subprocess.PIPE,
#                        stderr=subprocess.PIPE, check=True)
#         logging.debug("BLAST database '{0}' created here: {1}".format(
#             title, out))
#         return(0)
#     except:
#         if logger:
#             logging.error("Something bad happened when trying to make " +
#                           "a blast database")
#         sys.exit(1)


# def run_blastp(input_file, database_name, outfmt, blastp_exe='', logger=None):
#     """
#     requires logging subprocess, os, shutil
#     """
#     output_file = os.path.join(os.path.split(input_file)[0],
#                                str(os.path.splitext(
#                                    os.path.basename(input_file))[0] +
#                                    "_blast_hits.tab"))
#     if blastp_exe == '':
#         blastp_exe = shutil.which("blastp")
#     blastpcmd = str("{0} -db {1} -query {2} -out {3} -outfmt " +
#                     "{4}").format(blastp_exe, database_name, input_file,
#                                   output_file, outfmt)
#     if logger:
#         logger.info("Running blastp: {0}".format(blastpcmd))
#     try:
#         subprocess.run(blastpcmd, shell=sys.platform != "win32",
#                        stdout=subprocess.PIPE,
#                        stderr=subprocess.PIPE, check=True)
#         if logger:
#             logger.debug("Results from BLASTing {0} are here: {1}".format(
#                 input_file, output_file))
#         return(0)
#     except:
#         if logger:
#             logger.error("Something bad happened when running blast")
#         sys.exit(1)


# #def merge_outfiles(filelist, outfile_name):
# def merge_blast_tab_outfiles(filelist, outfile_name, logger=None):
#     """
#     #TODO needs a test for headers
#     for combining tab formated blast output format 6
#     returns 0 if successful
#     requires logging
#     """
#     # only grab .tab files, ie, the blast output
#     # logger=logging.getLogger()
#     filelist = [i for i in filelist if i.split(".")[-1:] == ['tab']]
#     if len(filelist) == 1:
#         if logger:
#             logger.warning("only one file found! no merging needed")
#         return(0)
#     elif len(filelist) == 0:
#         if logger:
#             logger.error("filelist empt; cannot perform merge!")
#         return(1)
#     else:
#         if logger:
#             logger.info("merging all the blast results to %s" % outfile_name)
#         nfiles = len(filelist)
#         fout = open(outfile_name, "a")
#         # first file:
#         for line in open(filelist[0]):
#             fout.write(line)
#         #  now the rest:
#         for num in range(1, nfiles):
#             f = open(filelist[num])
#             for line in f:
#                 fout.write(line)
#             f.close()  # not really needed
#         fout.close()
#         return(0)


# def cleanup_output_to_csv(infile,
#                           accession_pattern='(?P<accession>[A-Z _\d]*\.\d*)',
#                           logger=None):
#     """
#     given .tab from merge_blast_tab_outfiles, assign pretty column names,
#     """
#     # if logger:
#     #     logger=logging.getLogger(name=None)
#     print("cleaning up the csv output")
#     colnames = ["query_id", "subject_id", "identity_perc", "alignment_length",
#                 "mismatches", "gap_opens", "q_start", "q_end", "s_start",
#                 "s_end", "evalue", "bit_score"]
#     csv_results = pd.read_csv(open(infile), comment="#", sep="\t",
#                               names=colnames)
#     #This default regex will probably break things eventually...
#     # it looks for capital letter and numbers, dot, number, ie SHH11555JJ8.99
#     csv_results["accession"] = csv_results.query_id.str.extract(accession_pattern)
#     # write out results with new headers or with new headers and merged metadat from accessions.tab
#     genes = open(genelist, "r")
#     genedf = pd.read_csv(genes, sep=",")
#     output_path_csv = str(os.path.splitext(infile)[0] + ".csv")
#     results_annotated = pd.merge(csv_results, genedf, how="left",
#                                  on="accession")
#     results_annotated.to_csv(open(output_path_csv, "w"))
#     print("wrote final csv to %s" % output_path_csv)
# #%%



# def check_single_scaffold(input_genome_path):
#     """Test for single scaffold. from genbank
#     """


# def get_genbank_seq(input_genome_path, first_only=False):
#     """Get all sequences from genbank, return a list, unless first only
#     get the sequence from the FIRST record only in a genbank file
#     """
