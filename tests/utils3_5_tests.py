# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 08:57:31 2016
@author: nicholas
The Goal of this is to have a unified place to put the useful
python 3.5 functions or templates

how I got the fastq file
# seqtk sample -s 27 ~/GitHub/FA/pseudochromosome/data/20150803_Abram1/reads/3123-1_1_trimmed.fastq .0005 

bam file was from a riboseed mapping; md5: 939fbf2c282091aec0dfa278b05e94ec

mapped bam was made from bam file with the following command
 samtools view -Bh -F 4 /home/nicholas/GitHub/FB/Ecoli_comparative_genomics/scripts/riboSeed_pipeline/batch_coli_unpaired/map/mapping_20160906_region_7_riboSnag/test_smalt4_20160906_region_7_riboSnagS.bam > ~/GitHub/pyutilsnrw/tests/test_mapped.sam
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

sys.dont_write_bytecode = True

from pyutilsnrw import utils3_5
from pyutilsnrw.utils3_5 import make_output_prefix, check_installed_tools,\
    copy_file, get_ave_read_len_from_fastq, get_number_mapped,\
    extract_mapped_and_mappedmates, keep_only_first_contig, md5

    

@unittest.skipIf((sys.version_info[0] != 3) or (sys.version_info[1] < 5), "Subprocess.call, among otherthings wont run if you try this with less than python 3.5")
class utils3_5TestCase(unittest.TestCase):
    def setUp(self):
        pass

    def test_make_testing_dir(self):
        if not os.path.exists(testdirname):
            os.makedirs(testdirname)
        self.assertTrue(os.path.exists(testdirname))

    def test_make_output_prefix(self):
        test_prefix = make_output_prefix(testdirname, "utils_3.5")
        self.assertEqual(test_prefix,
                         "".join([testdirname,os.path.sep, "utils_3.5"]))
    
    def test_check_installed_tools(self):
        """is pwd on all mac/linux systems?
        #TODO replace with better passing test
        """
        check_installed_tools(["pwd"])
        # test fails properly
        with self.assertRaises(SystemExit):
            check_installed_tools(["thisisnotapathtoanactualexecutable"])

    def test_copy_file(self):
        if not os.path.exists(test_fastq_file):
            raise("test file is gone! There should be a file called test_reads.fastq")
        new_path = copy_file(current_file=test_fastq_file, dest_dir=testdirname,
                  name="newname.fastq", overwrite=False)
        self.assertEqual(new_path, os.path.join(testdirname,"newname.fastq"))
        # test overwrite exit
        with self.assertRaises(SystemExit):
            new_path = copy_file(current_file=test_fastq_file, 
                                 dest_dir=testdirname,
                                 name="newname.fastq", overwrite=False)
        os.remove(new_path)
    
    def test_get_ave_read_len_from_fastq(self):
        """this probably could/should be refined to have a better test
        """
        mean_read_len = get_ave_read_len_from_fastq(test_fastq_file, N=5)
        self.assertEqual(217.8, mean_read_len)
    # def test_this_fails(self):
    #     self.assertEqual("pinecone", 42)

    # def  tearDown(self): # and test clean temp dir
    #     """ I tried to do something like
    #     @unittest.skipUnless(clean_temp, "temporary files were retained")
    #     but couldnt get the variabel to be passed through.  For now, just make
    #     sure this is the last
    #     """
    #     print("running teardown method")
    #     clean_temp_dir(testdirname)

    def test_get_number_mapped(self):
        """This is bad cause I hardcoded the path for samtools
        """
        result = get_number_mapped(test_bam_file, samtools_exe)
        reference = "151 + 0 mapped (0.56% : N/A)"
        self.assertEqual(result, reference)

    def test_extract_mapped_and_mappedmates(self):
        """ dont trust this if  make_output_prefix test fails
        some help from PSS on SO:
        http://stackoverflow.com/questions/16874598/how-do-i-calculate-the-md5-checksum-of-a-file-in-python
        """
        prefix = make_output_prefix(output_dir=os.path.join(os.path.dirname(__file__),
                                                            "references"), 
                                    name="pyutilsnrw_sample")
        extract_mapped_and_mappedmates(map_results_prefix=prefix, 
                                       fetch_mates=False, 
                                       keep_unmapped=False,
                                       samtools_exe=samtools_exe)
        # reference mapping md5
        mapped_md5 = "27944249bf064ba54576be83053e82b0"  
        md5_returned = md5(str(prefix+"_mapped.sam"))

        # Finally compare original MD5 with freshly calculated
        self.assertEqual(mapped_md5, md5_returned)
        # delete files created
        files_created = ["_mapped.bam",
                         "_mapped.bam.bai",
                         "_mapped.sam"]
        if mapped_md5 == md5_returned:
            for i in files_created:
                os.remove(str(prefix+i))

    def test_keep_only_first_contig(self):
        """copy_file
        """
        copy_file(current_file=test_multifasta,
                  dest_dir=os.path.dirname(test_multifasta),
                  name='duplicated_multifasta.fasta', overwrite=False, 
                  logger=None)
        keep_only_first_contig(test_multifasta, newname="contig1")
        self.assertEqual(md5(test_multifasta), md5(test_singlefasta))
        copy_file(current_file=os.path.join(os.path.dirname(test_multifasta),
                               "duplicated_multifasta.fasta"),
                  dest_dir=os.path.dirname(test_multifasta),
                  name=os.path.basename(test_multifasta), overwrite=True, 
                  logger=None)
        os.remove(os.path.join(os.path.dirname(test_multifasta),
                               "duplicated_multifasta.fasta"))

def run_spades(output, ref, ref_as_contig, pe1_1='', pe1_2='', pe1_s='',
               as_paired=True, keep_best=True, prelim=False,
               k="21,33,55,77,99", seqname=''):
    """wrapper for common spades setting for long illumina reads
        ref_as_contig should be either blank, 'trusted', or 'untrusted'
        prelim flag is True, only assembly is run, and without coverage correction
        #TODO
        the seqname variable is used only for renaming the resulting contigs
        during iterative assembly.  It would be nice to inheirit from "ref",
        but that is changed with each iteration. This should probably be addressed
        before next major version change
    """
    if seqname == '':
        seqname = ref
    kmers = k  # .split[","]
    success = False
    #  prepare reference, if being used
    if not ref_as_contig == "":
        alt_contig = str("--%s-contigs %s" % (ref_as_contig, ref))
    else:
        alt_contig = ''
    # prepare read types, etc
    if as_paired and pe1_s != "":  # for libraries with both
        singles = str("--pe1-s %s " % pe1_s)
        pairs = str("--pe1-1 %s --pe1-2 %s " % (pe1_1, pe1_2))
    elif as_paired and pe1_s == "":  # for libraries with just PE
        singles = ""
        pairs = str("--pe1-1 %s --pe1-2 %s " % (pe1_1, pe1_2))
    elif pe1_s == "":  # for libraries treating paired ends as two single-end libs
        singles = ''
        pairs = str("--pe1-s %s --pe2-s %s " % (pe1_1, pe1_2))
    else:  # for 3 single end libraries
        singles = str("--pe1-s %s " % pe1_s)
        pairs = str("--pe2-s %s --pe3-s %s " % (pe1_1, pe1_2))
    reads = str(pairs+singles)
#    spades_cmds=[]
    if prelim:
        prelim_cmd =\
            str(args.spades_exe + " --only-assembler --cov-cutoff off --sc --careful -k {0}" +
                " {1} {2} -o {3}").format(kmers, reads, alt_contig, output)
        logger.info("Running the following command:\n{0}".format(prelim_cmd))
        subprocess.run(prelim_cmd,
                       shell=sys.platform != "win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE, check=True)
        success = output_from_subprocess_exists(os.path.join(output,
                                                             "contigs.fasta"))
        if prelim and keep_best and success:
            logger.info("reserving first contig")
            keep_only_first_contig(str(os.path.join(output, "contigs.fasta")),
                                   newname=
                                       os.path.splitext(os.path.basename(seqname))[0])
    else:
        spades_cmd = str(args.spades_exe + " --careful -k {0} {1} {2} -o " +
                         "{3}").format(kmers, reads, alt_contig, output)
        logger.info("Running the following command:\n{0}".format(spades_cmd))
        subprocess.run(spades_cmd,
                       shell=sys.platform != "win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE)
        # not check=True; dont know spades return codes
        success = output_from_subprocess_exists(os.path.join(output,
                                                             "contigs.fasta"))
    return("{0}contigs.fasta".format(os.path.join(output, "")), success)


def run_quast(contigs, output, ref=""):
    """Reference is optional. This is, honestly, a pretty dumb feature
    """
    if not ref == "":
        ref = str("-R %s" % ref)
    quast_cmd = str(args.quast_exe + " {0} {1} -o " +
                    "{2}").format(contigs, ref, output)
    logger.info("Running quast as follows: {0}".format(quast_cmd))
    subprocess.run(quast_cmd,
                   shell=sys.platform != "win32",
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE)




def combine_contigs(mauve_path, contigs_name="riboSeedContigs_aggregated"):
    """changed over to biopython
    combine all *.fasta in dir, return path to concatenated file
    requires Bio.SeqIO, glob, os
    """
    output = os.path.join(mauve_path, str(contigs_name+".fasta"))
    fastas = glob.glob(str(mauve_path+"*.fasta"))
    with open(output, 'w') as w_file:
        for filen in fastas:
            with open(filen, 'rU') as o_file:
                seq_records = SeqIO.parse(o_file, 'fasta')
                SeqIO.write(seq_records, w_file, 'fasta')

    return(output)


def setup_protein_blast(input_file, input_type="fasta", dbtype="prot",
                        title="blastdb", out="blastdb",
                        makeblastdb_exe=''):
    """
    This runs make blast db with the given parameters
    requires logging, os, subprocess, shutil
    """
    logger = logging.getLogger(__name__)
    #logging.getLogger(name=None)
    logger.debug("TESTING I 2 3!")
    if makeblastdb_exe == '':
        makeblastdb_exe = shutil.which("makeblastdb")
    makedbcmd = str("{0} -in {1} -input_type {2} -dbtype {3} " +
                    "-title {4} -out {5}").format(makeblastdb_exe, 
                        input_file, input_type, dbtype, title, out)
    logger.info("Making blast db: {0}".format(makedbcmd))
    try:
        subprocess.run(makedbcmd, shell=sys.platform != "win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE, check=True)
        logging.debug("BLAST database '{0}' created here: {1}".format(
            title, out))
        return(0)
    except:
        logging.error("Something bad happened when trying to make " +
                     "a blast database")
        sys.exit(1)


def run_blastp(input_file, database_name, outfmt, blastp_exe=''):
    """
    requires logging subprocess, os, shutil
    """
    #logger = logging.getLogger(name=None)
    logger = logging.getLogger(__name__)
    output_file = os.path.join(os.path.split(input_file)[0],
                               str(os.path.splitext(
                                   os.path.basename(input_file))[0] +
                                   "_blast_hits.tab"))
    if blastp_exe == '':
        blastp_exe = shutil.which("blastp")
    blastpcmd = str("{0} -db {1} -query {2} -out {3} -outfmt " +
                    "{4}").format(blastp_exe, database_name, input_file,
                                  output_file, outfmt)
    logger.info("Running blastp: {0}".format(blastpcmd))
    try:
        subprocess.run(blastpcmd, shell=sys.platform != "win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE, check=True)
        logger.debug("Results from BLASTing {0} are here: {1}".format(
            input_file, output_file))
        return(0)
    except:
        logger.error("Something bad happened when running blast")
        sys.exit(1)


#def merge_outfiles(filelist, outfile_name):
def merge_blast_tab_outfiles(filelist, outfile_name):
    """
    #TODO needs a test for headers
    for combining tab formated blast output format 6
    returns 0 if successful
    requires logging
    """
    # only grab .tab files, ie, the blast output
    logger=logging.getLogger()
    filelist = [i for i in filelist if i.split(".")[-1:] == ['tab']]
    if len(filelist) == 1:
        logger.warning("only one file found! no merging needed")
        return(0)
    elif len(filelist) == 0:
        logger.error("filelist empt; cannot perform merge!")
        return(1)
    else:
        logger.info("merging all the blast results to %s" % outfile_name)
        nfiles = len(filelist)
        fout = open(outfile_name, "a")
        # first file:
        for line in open(filelist[0]):
            fout.write(line)
        #  now the rest:
        for num in range(1, nfiles):
            f = open(filelist[num])
            for line in f:
                fout.write(line)
            f.close()  # not really needed
        fout.close()
        return(0)


def cleanup_output_to_csv(infile, accession_pattern='(?P<accession>[A-Z _\d]*\.\d*)'):
    """
    given .tab from merge_blast_tab_outfiles, assign pretty column names,
    """
    logger=logging.getLogger(name=None)
    print("cleaning up the csv output")
    colnames = ["query_id", "subject_id", "identity_perc", "alignment_length", "mismatches",
                "gap_opens", "q_start", "q_end", "s_start", "s_end", "evalue", "bit_score"]
    csv_results = pd.read_csv(open(infile), comment="#", sep="\t", names=colnames)
    #This regex will probably break things rather badly before too long...
    # it looks for capital letter and numbers, dot, number, ie SHH11555JJ8.99
    csv_results["accession"] = csv_results.query_id.str.extract(accession_pattern)
    # write out results with new headers or with new headers and merged metadat from accessions.tab
    genes = open(genelist, "r")
    genedf = pd.read_csv(genes, sep=",")
    output_path_csv = str(os.path.splitext(infile)[0]+".csv")
    results_annotated = pd.merge(csv_results, genedf, how="left",  on="accession")
    results_annotated.to_csv(open(output_path_csv, "w"))
    print("wrote final csv to %s" % output_path_csv)
#%%


if __name__ == '__main__':
    # logger = utils3_5.set_up_logging(verbosity=1, name=__name__,
    #                                  outfile=os.path.join(curdir, 
    #                                                       "tests_logger.txt"))
    curdir=os.getcwd()
    samtools_exe="/usr/bin/samtools"
    testdirname=os.path.join(os.path.dirname(__file__), "utils3_5tests")
    test_fastq_file = os.path.join(os.path.dirname(__file__), 
                                   str("references"+os.path.sep+'test_reads.fastq'))
    test_bam_file = os.path.join(os.path.dirname(__file__),
                                   str("references"+os.path.sep+
                                       "pyutilsnrw_sample.bam"))
    test_bam_mapped_file = os.path.join(os.path.dirname(__file__), 
                                   str("references"+os.path.sep+
                                       "pyutilsnrw_sample_mapped.bam"))
    test_multifasta = os.path.join(os.path.dirname(__file__), 
                                   str("references"+os.path.sep+
                                       "test_multiseqs_reference.fasta"))
    test_singlefasta = os.path.join(os.path.dirname(__file__), 
                                   str("references"+os.path.sep+
                                       "test_only_first_reference.fasta"))
    clean_temp = True
    unittest.main()
