# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 08:57:31 2016
@author: nicholas
The Goal of this is to have a unified place to put the useful
python 3.5 functions or templates
"""
__version__ = "0.0.2"
import time
import sys
import shutil
import logging
import subprocess
import os
import unittest
sys.dont_write_bytecode = True


@unittest.skipIf((sys.system_info[0] != 3) or (sys.system_info[1] < 5), "Subprocess.call, among otherthings wont run if you try this with les than python 3.5")
from pyutilsnrw import utils3_5
class utils3_5TestCase(unittest.TestCase):
    def test_set_up_root_logging(self):
        """
        This checks that set_up_root_logging ouputs a file with the
        given verbosity. 
        """
        logfile = "test_log.txt"
        test_logger = utils3_5.set_up_root_logging(3, logfile)
        self.assertEqual('hello', 'world')
        # check that s.split fails when the separator is not a string
        with self.assertRaises(TypeError):
            s.split(2)

def set_up_root_logging(verbosity, outfile):
    """derived from set_up_logging; had problem where functions in modules
    could only log to root. If you cant lick 'em, join 'em.
    requires logging, os, sys, time
    logs debug level to file, and [verbosity] level to stderr
    return a logger object
     """
    import logging
    if (verbosity*10) not in range(10, 60, 10):
        raise ValueError('Invalid log level: %s' % verbosity)
    try:
        logging.basicConfig(level=logging.DEBUG,
                            format="%(asctime)s - %(levelname)s - %(message)s",
                            datefmt='%m-%d %H:%M:%S',
                            filename=outfile,
                            filemode='w')
    except:
        logger.error("Could not open {0} for logging".format(outfile))
        sys.exit(1)
    # define a Handler which writes INFO messages or higher to the sys.stderr
    console = logging.StreamHandler(sys.stderr)
    console.setLevel(level=(verbosity *10))
    # set a format which is simpler for console use
    formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    # tell the handler to use this format
    console.setFormatter(formatter)
    # add the handler to the root logger
    logging.getLogger('').addHandler(console)
    # Now, we can log to the root logger, or any other logger. First the root...
    logging.info("Initializing logger")
    logging.debug("logging at level {0}".format(verbosity))
    logger = logging.getLogger()    
    return(logger)

def last_exception():
    """ From Pyani
    Returns last exception as a string, or use in logging.
    """
    exc_type, exc_value, exc_traceback = sys.exc_info()
    return(''.join(traceback.format_exception(exc_type, exc_value,
                                              exc_traceback)))

# def make_outdir_carful(dirname):
def make_outdir(dirname):
    """ From Pyani
    Make the output directory, if required.
    This is a little involved.  If the output directory already exists,
    we take the safe option by default, and stop with an error.  We can,
    however, choose to force the program to go on, in which case we can
    either clobber the existing directory, or not.  The options turn out
    as the following, if the directory exists:
    DEFAULT: stop and report the collision
    FORCE: continue, and remove the existing output directory
    NOCLOBBER+FORCE: continue, but do not remove the existing output
    """
    if os.path.exists(dirname):
        if not args.force:
            logger.error("Output directory %s would " % dirname +
                         "overwrite existing files (exiting)")
            sys.exit(1)
        else:
            logger.info("Removing directory %s and everything below it" %
                        dirname)
            if args.noclobber:
                logger.warning("NOCLOBBER: not actually deleting directory")
            else:
                shutil.rmtree(args.output)
    logger.info("Creating directory %s" % dirname)
    try:
        os.makedirs(dirname)   # We make the directory recursively
        # Depending on the choice of method, a subdirectory will be made for
        # alignment output files
    except OSError:
        # This gets thrown if the directory exists. If we've forced overwrite/
        # delete and we're not clobbering, we let things slide
        if args.noclobber and args.force:
            logger.info("NOCLOBBER+FORCE: not creating directory")
        else:
            logger.error(last_exception)
            sys.exit(1)


#def make_output_prefix(map_output_dir, exp_name):
def make_output_prefix(output_dir, name):
    """ makes output prefix from output directory and name.
    requires os, logger
    """
    logger.debug(" output_prefix: %s" % os.path.join(output_dir, name))
    return(os.path.join(output_dir, name))


#def copy_ref_to_temp(current_file, dest_dir, overwrite=False):
def copy_file(current_file, dest_dir, overwrite=False):
    """Copy reference fasta file to dest_dir to avoid making messy
    indexing files everywhere (generated during mapping).  Could all be done
    with shutil.copyfile, but I cant figure a way to handle the errors as well..
    This uses a system call to "rm x -f"; is this safe? I have to have the
    -f flag because otherwise it prompts, which disrumpt the flow.
    require shutil, logger, subprocess, sys, os, subprocess
    """
    new_ref = os.path.join(dest_dir, os.path.basename(current_file))
    if os.path.exists(new_ref):
        if overwrite:
            try:
                logger.debug("removing {0} to be overwritten with {1}.".format(
                             new_ref, current_file))
                rm_cmd = "rm -f {0}".format(new_ref)
                logger.debug(rm_cmd)
                subprocess.run(rm_cmd, shell=sys.platform != "win32",
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE, check=True)
            except:
                logger.error("couldn't overwrite file with copy_ref_to_temp")
                sys.exit(1)
        else:
            logger.error("cannot overwrite {0} to {1}".format(new_ref, current_file))
            sys.exit(1)
    else:
        logger.error
    logger.info("copying fasta from %s to %s" % (current_file, new_ref))
    cmd = str("cp %s %s" % (current_file, new_ref))
    subprocess.run(cmd, shell=sys.platform != "win32",
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE, check=True)
    return(new_ref)


def check_installed_tools(list_of_tools):
    """given a list of executables (or things that should be in PATH,
    raise error if executable for a tool is not found
    requires shutil, logger, sys
    """
    for i in list_of_tools:
        if not shutil.which(i):
            logger.error("Must have {0} installed in PATH!".format(i))
            sys.exit(1)
        else:
            logger.debug("{0} executable found".format(i))




def get_ave_read_len_from_fastq(fastq1, N=50):
    """from LP; return average read length in fastq1 file from first N reads
    """
    count, tot = 0, 0
    if os.path.splitext(fastq1)[-1] in ['.gz', '.gzip']:
        open_fun = gzip.open
    else:
        open_fun = open
    data = SeqIO.parse(open_fun(fastq1, "rt"), "fastq")
    for read in data:
        count += 1
        tot += len(read)
        if count >= N:
            break
    logger.info("From the first {0} reads in {1}, mean length is {2}".format(
                N, os.path.basename(fastq1), float(tot/count)))
    return(float(tot/count))


#def get_number_mapped(bam):
def get_number_mapped(bam, samtools_exe):
    """use samtools flagstats to retrieve total mapped reads as a diagnostic
    returns a string to be printed, the 4th line of flagstat
    requires subprocess, sys, logger (SAMtools)
    """
    flagstats = subprocess.run(str("{0} flagstat {1}").format(samtools_exe, bam),
                               shell=sys.platform != "win32",
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               check=True)
    try:
        printout = flagstats.stdout.decode("utf-8").split("\n")[4]
    except IndexError:
        logger.error("Error reading {0} to determine number of mapped reads".format(bam))
        sys.exit(1)
    return(printout)


def extract_mapped_and_mappedmates(map_results_prefix, fetch_mates, keep_unmapped):
    """
    IF fetch_mates is true, mapped reads are extracted, and mates are feteched with the LC_ALL line.
    If not, that part is skipped, and just the mapped reads are extracted.
    Setting keep_unmapped to true will output a bam file with all the remaining reads.  This could
    be used if you are really confident there are no duplicate mappings you are interested in.
     -F 4 option selects mapped reads
    Note that the umapped output includes reads whose pairs were mapped.  This is
    to try to catch the stragglers.
    LC_ALL=C  call from pierre lindenbaum. No idea how it does, but its magic
    """
    extract_cmds = []
    logger.info("Extracting the reads of interest")
    # Either get nates or ignore mates
    if fetch_mates:
        samview = str(args.samtools_exe + "  view -h -F 4 {0}.bam  | cut -f1 > " +
                      "{0}_mappedIDs.txt").format(map_results_prefix)
        lc_cmd = str("LC_ALL=C grep -w -F -f {0}_mappedIDs.txt  < {0}.sam >  " +
                     "{0}_mapped.sam").format(map_results_prefix)
        extract_cmds.extend([samview, lc_cmd])
    else:
        samview = str(args.samtools_exe + " view -hS -F 4 {0}.bam > " +
                      "{0}_mapped.sam").format(map_results_prefix)
        extract_cmds.extend([samview])
    samsort = str(args.samtools_exe + " view -bhS {0}_mapped.sam " +
                  "| samtools sort - > {0}_mapped.bam").format(map_results_prefix)
    samindex = str(args.samtools_exe + " index {0}_mapped.bam").format(map_results_prefix)
    extract_cmds.extend([samsort, samindex])
    if keep_unmapped:
        samviewU = str(args.samtools_exe + "  view -f 4 {0}.bam  | cut -f1 > " +
                       "{0}_unmappedIDs.txt").format(map_results_prefix)
        lc_cmdU = str("LC_ALL=C grep -w -F -f {0}_unmappedIDs.txt  < " +
                      "{0}.sam > {0}_unmapped.sam").format(map_results_prefix)
        samindexU = str(args.samtools_exe + " view -bhS {0}_unmapped.sam | samtools " +
                        "sort - -o {0}_unmapped.bam && samtools index " +
                        "{0}_unmapped.bam").format(map_results_prefix)
        extract_cmds.extend([samviewU, lc_cmdU, samindexU])
    logger.debug("running the following commands to extract reads:")
    for i in extract_cmds:
        logger.debug(i)
        subprocess.run(i, shell=sys.platform != "win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE, check=True)


#def clean_temp_dir(map_output_dir):
def clean_temp_dir(temp_dir):
    """ from http://stackoverflow.com/questions/185936/delete-folder-contents-in-python
        this should fail on read-only files
        requires shutil, os
    """
    for the_file in os.listdir(temp_dir):
        file_path = os.path.join(temp_dir, the_file)
        try:
            if os.path.isfile(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)
        except Exception as e:
            print(e)


def output_from_subprocess_exists(output):
    """just what it says: given an output path or dir, return True if it exists
    requires os
    """
    if os.path.exists(output):
        return(True)
    else:
        return(False)


def keep_only_first_contig(ref, newname="contig1"):
    #TODO make a biopython version
    """ 
    given a multi fasta from SPAdes, extract first entry,
    rename "NODE_1" with newname, overwrite file
    requires os, re
    """
    temp = os.path.join(os.path.dirname(ref), "temp.fasta")
    with open(temp, "w") as outfile:
        with open(ref, "r") as file_handle:
            lines = file_handle.readlines()
            if lines[0][0] != ">":
                raise ValueError("Something wrong with spades output contig! Not a valid fasta!")
            new_header = str(re.sub("NODE_1", newname, str(lines[0])))
            outfile.write(new_header)
            for line in range(1, len(lines)):
                if lines[line][0] == ">":
                    break  # stop after first entry
                outfile.write(str(lines[line]))
    os.remove(ref)
    os.rename(temp, ref)


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
    unittest.main()
