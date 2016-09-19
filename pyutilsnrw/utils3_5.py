# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 08:57:31 2016
@author: nicholas
The Goal of this is to have a unified place to put the useful
python 3.5 functions or templates

Minor version changes:
- added explicit if logger:
- added hashlib, re, various bits and bobs
- found bug in combine_contigs that ignored dir if didnt end in
  path sep
"""
__version__ = "0.0.3"
import time
import sys
import shutil
import logging
import subprocess
import os
from Bio import SeqIO
import hashlib
import re
import glob
import gzip
import errno
import argparse

logger = logging.getLogger('root')


def get_args_template():
    """
    Template for argument parsing; dont import this, by the way
    requires argparse
    """
    parser = argparse.ArgumentParser(
        description="Given regions from riboSnag, assembles the mapped reads")
    parser.add_argument("seed_dir", action="store",
                        help="path to roboSnag results directory")
    parser.add_argument("-F", "--fastq1", dest='fastq1', action="store",
                        help="forward fastq reads, can be compressed",
                        type=str, default="")
    parser.add_argument("-R", "--fastq2", dest='fastq2', action="store",
                        help="reverse fastq reads, can be compressed",
                        type=str, default="")
    parser.add_argument("-S", "--fastq_single", dest='fastqS', action="store",
                        help="single fastq reads", type=str, default="")
    parser.add_argument("-n", "--experiment_name", dest='exp_name',
                        action="store",
                        help="prefix for results files; default: %(default)s",
                        default="riboSeed", type=str)
    parser.add_argument("-m", "--method_for_map", dest='method',
                        action="store",
                        help="availible mappers: smalt; default: %(default)s",
                        default='smalt', type=str)
    parser.add_argument("-c", "--cores", dest='cores', action="store",
                        default=1, type=int,
                        help="cores for multiprocessing workers" +
                        "; default: %(default)s")
    parser.add_argument("-r", "--reference_genome", dest='reference_genome',
                        action="store", default='', type=str,
                        help="fasta reference genome, used for estimating " +
                        "insert sizes, QUAST, and SPAdes")
    parser.add_argument("-o", "--output", dest='output', action="store",
                        help="output directory; " +
                        "default: %(default)s", default=os.getcwd(), type=str)
    parser.add_argument("--paired_inference", dest='paired_inference',
                        action="store_true", default=False,
                        help="if --paired_inference, mapped read's " +
                        "pairs are included; default: %(default)s")
    parser.add_argument("--subtract", dest='subtract', action="store_true",
                        default=False, help="if --subtract, reads aligned " +
                        "to each reference will not be aligned to future " +
                        "iterations.  Probably you shouldnt do this" +
                        "unless you really happen to want to")
    parser.add_argument("--keep_unmapped", dest='keep_unmapped',
                        action="store_true", default=False,
                        help="if --keep_unmapped fastqs are generated " +
                        "containing the unmapped reads; default: %(default)s")
    parser.add_argument("--ref_as_contig", dest='ref_as_contig',
                        action="store",
                        default="", type=str,
                        help="if 'trusted', SPAdes will  use the seed " +
                        "sequences as a --trusted-contig; if 'untrusted', " +
                        "SPAdes will treat as --untrusted-contig. if '', " +
                        "seeds will not be used during assembly. " +
                        "See SPAdes docs; default: %(default)s")
    parser.add_argument("--temps", dest='temps', action="store_true",
                        default=False,
                        help="if --temps, intermediate files will be " +
                        "kept; default: %(default)s")
    parser.add_argument("-i", "--iterations", dest='iterations',
                        action="store",
                        default=2, type=int,
                        help="if  > 1, repeated iterationss will occur after\
                        assembly of seed regions ; default: %(default)s")
    parser.add_argument("-v", "--verbosity", dest='verbosity', action="store",
                        default=2, type=int,
                        help="1 = debug(), 2 = info(), 3 = warning(), " +
                        "4 = error() and 5 = critical(); default: %(default)s")
    parser.add_argument("--DEBUG", dest='DEBUG', action="store_true",
                        default=False,
                        help="if --DEBUG, test data will be " +
                        "used; default: %(default)s")
    parser.add_argument("--force", dest='force', action="store_true",
                        default=False,
                        help="if --force, existing results dirs will be " +
                        "used; default: %(default)s")
    ##TODO  Make these check a config file
    parser.add_argument("--spades_exe", dest="spades_exe",
                        action="store", default="spades.py",
                        help="Path to spades executable; default: %(default)s")
    parser.add_argument("--samtools_exe", dest="samtools_exe",
                        action="store", default="samtools",
                        help="Path to bwa executable; default: %(default)s")
    parser.add_argument("--smalt_exe", dest="smalt_exe",
                        action="store", default="smalt",
                        help="Path to smalt executable; default: %(default)s")
    parser.add_argument("--quast_exe", dest="quast_exe",
                        action="store", default="quast.py",
                        help="Path to quast executable; default: %(default)s")
    args = parser.parse_args()
    return(args)

#def set_up_logging_given_numeric_verbosity(verbosity, outfile):
def set_up_logging(verbosity, outfile, name):
    """
    Set up logging a la pyani, with
    a little help from:
    https://aykutakin.wordpress.com/2013/08/06/logging-to-console-and-file-in-python/
    requires logging, os, sys, time
    logs debug level to file, and [verbosity] level to stderr
    return a logger object
    """
    import logging
    # levels = ["debug", "info", "warning", "error", "critical"]
    logger = logging.getLogger(str(name + ": %s" %
                               time.asctime()))
    if (verbosity*10) not in range(10, 60, 10):
        raise ValueError('Invalid log level: %s' % verbosity)
    # setting root level; otherwisse, only logs warining and up
    # setting root level; otherwisse, only logs warining and up
    # logging.basicConfig(level=logging.DEBUG)
    logger.setLevel(logging.DEBUG)
    # create console handler and set level to given verbosity
    console_err = logging.StreamHandler(sys.stderr)
    console_err.setLevel(level = (verbosity * 10))
    console_err_formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    console_err.setFormatter(console_err_formatter)
    logger.addHandler(console_err)
    # create debug file handler and set level to debug
    try:
        logfile_handler = logging.FileHandler(outfile, "w")
        logfile_handler.setLevel(logging.DEBUG)
        logfile_handler_formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
        logfile_handler.setFormatter(logfile_handler_formatter)
        logger.addHandler(logfile_handler)
    except:
        logger.error("Could not open {0} for logging".format(outfile))
        sys.exit(1)
    logger.info("Initializing logger")
    logger.debug("logging at level {0}".format(verbosity))
    return(logger)


def set_up_root_logging(verbosity, outfile):
    """derived from set_up_logging; had problem where functions in modules
    could only log to root. If you cant lick 'em, join 'em.
    requires logging, os, sys, time
    logs debug level to file, and [verbosity] level to stderr
    return a logger object
     """
    import logging
    if (verbosity * 10) not in range(10, 60, 10):
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
    console.setLevel(level = (verbosity *10))
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


def make_outdir(path, logger=None):
    """makes a directory if it doesnt exist
    requires os, errno
    """
    try:
        os.makedirs(path)
        if logger:
            logger.debug("creating new directory: {0}".format(path))
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            if logger:
                logger.error("Cannot make directory {0}!".format(path))
            sys.exit(1)


#def make_output_prefix(map_output_dir, exp_name):
def make_output_prefix(output_dir, name, logger=None):
    """ makes output prefix from output directory and name.
    requires os, logger
    """
    if logger:
        logger.debug(" output_prefix: %s" % os.path.join(output_dir, name))
    return(os.path.join(output_dir, name))


#def copy_ref_to_temp(current_file, dest_dir, overwrite=False):
def copy_file(current_file, dest_dir, name='', overwrite=False, logger=None):
    """Copy reference fasta file to dest_dir to avoid making messy
    indexing files everywhere (generated during mapping).
    There is an option to rename resulting file. Probably could all be done
    with shutil.copyfile, but I cant figure a way to handle the errors as well..
    This uses a system call to "rm x -f"; is this safe? I have to have the
    -f flag because otherwise it prompts, which disrumpt the flow.
    require shutil, logger, subprocess, sys, os, subprocess
    returns new path
    """
    if type(name) is str and name != "":
        new_file_name = str(name)
    else:
        new_file_name = os.path.basename(current_file)
    new_ref = os.path.join(dest_dir, new_file_name)
    if os.path.exists(new_ref):
        if overwrite:
            try:
                rm_cmd = "rm -f {0}".format(new_ref)
                if logger:
                    logger.debug(str("removing {0} to be overwritten " +
                                     "with {1}.").format(
                        new_ref, current_file))
                    logger.debug(rm_cmd)
                subprocess.run(rm_cmd, shell=sys.platform != "win32",
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE, check=True)
            except:
                if logger:
                    logger.error(str("couldn't overwrite file with" +
                                     " copy_ref_to_temp"))
                sys.exit(1)
        else:
            if logger:
                logger.error(str("cannot overwrite {0} " +
                             "to {1}").format(new_ref, current_file))
            sys.exit(1)
    else:
        if logger:
            logger.info(str("copying fasta from {0} to " +
                            "{1}").format(current_file, new_ref))
    cmd = str("cp %s %s" % (current_file, new_ref))
    subprocess.run(cmd, shell=sys.platform != "win32",
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE, check=True)
    return(new_ref)


def check_installed_tools(list_of_tools, logger=None):
    """given a list of executables (or things that should be in PATH,
    raise error if executable for a tool is not found
    requires shutil, logger, sys
    """
    for i in list_of_tools:
        if not shutil.which(i):
            if logger:
                logger.error("Must have {0} installed in PATH!".format(i))
            sys.exit(1)
        else:
            if logger:
                logger.debug("{0} executable found".format(i))


def get_ave_read_len_from_fastq(fastq1, N=50, logger=None):
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
    if logger:
        logger.info(str("From the first {0} reads in {1}, " +
                        "mean length is {2}").format(N,
                                                     os.path.basename(fastq1),
                                                     float(tot / count)))
    return(float(tot / count))


#def get_number_mapped(bam):
def get_number_mapped(bam, samtools_exe, logger=None):
    """use samtools flagstats to retrieve total mapped reads as a diagnostic
    returns a string to be printed, the 4th line of flagstat
    requires subprocess, sys, logger (SAMtools)
    """
    flagstatcmd = str("{0} flagstat {1}").format(samtools_exe, bam)
    flagstats = subprocess.run(flagstatcmd,
                               shell=sys.platform != "win32",
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               check=True)
    try:
        printout = flagstats.stdout.decode("utf-8").split("\n")[4]
        #TODO test for none mapped
    except IndexError:
        if logger:
            logger.error("Error reading {0}".format(bam))
        sys.exit(1)
    return(printout)


def extract_mapped_and_mappedmates(map_results_prefix, fetch_mates,
                                   keep_unmapped, samtools_exe, logger=None):
    """
    Take a prefix for a dir containing your mapped bam file.
    IF fetch_mates is true, mapped reads are extracted,
    and mates are feteched with the LC_ALL line.If not, that part is
    skipped, and just the mapped reads are extracted.
    Setting keep_unmapped to true will output a bam file with
    all the remaining reads.  This could be used if you are really confident
    there are no duplicate mappings you are interested in.
     -F 4 option selects mapped reads
    Note that the umapped output includes reads whose pairs were mapped.
    This is to try to catch the stragglers.
    LC_ALL=C  call from pierre lindenbaum. No idea how it does, but its magic
    """
    extract_cmds = []
    # Either get nates or ignore mates
    if fetch_mates:
        samview = str(samtools_exe + "  view -h -F 4 {0}.bam  | cut -f1 > " +
                      "{0}_mappedIDs.txt").format(map_results_prefix)
        lc_cmd = str("LC_ALL=C grep -w -F -f {0}_mappedIDs.txt  < {0}.sam > " +
                     "{0}_mapped.sam").format(map_results_prefix)
        extract_cmds.extend([samview, lc_cmd])
    else:
        samview = str(samtools_exe + " view -hS -F 4 {0}.bam > " +
                      "{0}_mapped.sam").format(map_results_prefix)
        extract_cmds.extend([samview])
    samsort = str(samtools_exe + " view -bhS {0}_mapped.sam " +
                  "| samtools sort - > " +
                  "{0}_mapped.bam").format(map_results_prefix)
    samindex = " {0} index {1}_mapped.bam".format(samtools_exe,
                                                  map_results_prefix)
    extract_cmds.extend([samsort, samindex])
    if keep_unmapped:
        samviewU = str(samtools_exe + "  view -f 4 {0}.bam  | cut -f1 > " +
                       "{0}_unmappedIDs.txt").format(map_results_prefix)
        lc_cmdU = str("LC_ALL=C grep -w -F -f {0}_unmappedIDs.txt  < " +
                      "{0}.sam > {0}_unmapped.sam").format(map_results_prefix)
        samindexU = str("{1} view -bhS {0}_unmapped.sam | samtools " +
                        "sort - -o {0}_unmapped.bam && samtools index " +
                        "{0}_unmapped.bam").format(map_results_prefix,
                                                   samtools_exe)
        extract_cmds.extend([samviewU, lc_cmdU, samindexU])
    if logger:
        logger.debug("running the following commands to extract reads:")
    for i in extract_cmds:
        if logger:
            logger.debug(i)
        subprocess.run(i, shell=sys.platform != "win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE, check=True)


#def clean_temp_dir(map_output_dir):
def clean_temp_dir(temp_dir, logger=None):
    """ from http://stackoverflow.com/questions/
            185936/delete-folder-contents-in-python
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
            if logger:
                logger.error(e)
            sys.exit(1)


def output_from_subprocess_exists(output):
    """just what it says: given an output path or dir,
    return True if it exists
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
    fasta_sequences = SeqIO.parse(open(ref),'fasta')
    with open(output_file) as out_file:
        for fasta in fasta_sequences:
            name, sequence = fasta.id, fasta.seq.tostring()
            new_sequence = some_function(sequence)
            write_fasta(out_file)

    """
    temp = os.path.join(os.path.dirname(ref), "temp.fasta")
    with open(temp, "w") as outfile:
        with open(ref, "r") as file_handle:
            lines = file_handle.readlines()
            if lines[0][0] != ">":
                raise ValueError(str("Error  with spades output contig!" +
                                     " Not a valid fasta!"))
            new_header = str(re.sub("NODE_\d*", newname, str(lines[0])))
            outfile.write(new_header)
            for line in range(1, len(lines)):
                if lines[line][0] == ">":
                    break  # stop after first entry
                outfile.write(str(lines[line]))
    os.remove(ref)
    os.rename(temp, ref)


# def run_spades(output, ref, ref_as_contig, pe1_1='', pe1_2='', pe1_s='',
#                as_paired=True, keep_best=True, prelim=False,
#                k="21,33,55,77,99", seqname=''):
#     """wrapper for common spades setting for long illumina reads
#         ref_as_contig should be either blank, 'trusted', or 'untrusted'
#         prelim flag is True, only assembly is run, and without coverage correction
#         #TODO
#         the seqname variable is used only for renaming the resulting contigs
#         during iterative assembly.  It would be nice to inheirit from "ref",
#         but that is changed with each iteration. This should probably be addressed
#         before next major version change
#     """
#     if seqname == '':
#         seqname = ref
#     kmers = k  # .split[","]
#     success = False
#     #  prepare reference, if being used
#     if not ref_as_contig == "":
#         alt_contig = str("--%s-contigs %s" % (ref_as_contig, ref))
#     else:
#         alt_contig = ''
#     # prepare read types, etc
#     if as_paired and pe1_s != "":  # for libraries with both
#         singles = str("--pe1-s %s " % pe1_s)
#         pairs = str("--pe1-1 %s --pe1-2 %s " % (pe1_1, pe1_2))
#     elif as_paired and pe1_s == "":  # for libraries with just PE
#         singles = ""
#         pairs = str("--pe1-1 %s --pe1-2 %s " % (pe1_1, pe1_2))
#     elif pe1_s == "":  # for libraries treating paired ends as two single-end libs
#         singles = ''
#         pairs = str("--pe1-s %s --pe2-s %s " % (pe1_1, pe1_2))
#     else:  # for 3 single end libraries
#         singles = str("--pe1-s %s " % pe1_s)
#         pairs = str("--pe2-s %s --pe3-s %s " % (pe1_1, pe1_2))
#     reads = str(pairs+singles)
# #    spades_cmds=[]
#     if prelim:
#         prelim_cmd =\
#             str(args.spades_exe + " --only-assembler --cov-cutoff off --sc --careful -k {0}" +
#                 " {1} {2} -o {3}").format(kmers, reads, alt_contig, output)
#         logger.info("Running the following command:\n{0}".format(prelim_cmd))
#         subprocess.run(prelim_cmd,
#                        shell=sys.platform != "win32",
#                        stdout=subprocess.PIPE,
#                        stderr=subprocess.PIPE, check=True)
#         success = output_from_subprocess_exists(os.path.join(output,
#                                                              "contigs.fasta"))
#         if prelim and keep_best and success:
#             logger.info("reserving first contig")
#             keep_only_first_contig(str(os.path.join(output, "contigs.fasta")),
#                                    newname=
#                                        os.path.splitext(os.path.basename(seqname))[0])
#     else:
#         spades_cmd = str(args.spades_exe + " --careful -k {0} {1} {2} -o " +
#                          "{3}").format(kmers, reads, alt_contig, output)
#         logger.info("Running the following command:\n{0}".format(spades_cmd))
#         subprocess.run(spades_cmd,
#                        shell=sys.platform != "win32",
#                        stdout=subprocess.PIPE,
#                        stderr=subprocess.PIPE)
#         # not check=True; dont know spades return codes
#         success = output_from_subprocess_exists(os.path.join(output,
#                                                              "contigs.fasta"))
#     return("{0}contigs.fasta".format(os.path.join(output, "")), success)


def run_quast(contigs, output, quast_exe, ref="", threads=1, logger=None):
    """Reference is optional. This is, honestly, a pretty dumb feature
    requires sys, subprocess, (system install of quast)
    """
    if not ref == "":
        ref = str("-R %s" % ref)
    quast_cmd = str("{3}  {0} {1} -t {3} -o " +
                    "{2}").format(contigs, ref, output, quast_exe, threads)
    if logger:
        logger.info("Running quast as follows: {0}".format(quast_cmd))
    subprocess.run(quast_cmd,
                   shell=sys.platform != "win32",
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE)


def combine_contigs(contigs_dir, pattern = "*",
                    contigs_name="riboSeedContigs_aggregated",
                    ext=".fasta", verbose=False, logger=None):
    """changed over to biopython
    combine all *ext files in dir, return path to concatenated file
    requires Bio.SeqIO, glob, os
    """
    if contigs_dir[-1] != os.path.sep:
        contigs_dir = str(contigs_dir + os.path.sep)
    output = os.path.join(contigs_dir, str(contigs_name + ext))
    fastas = glob.glob(str(contigs_dir + pattern + ext))
    fastas.sort()
    if verbose:
        print(str("combining the following files matching pattern " +
                  "{0}:{1}".format(pattern, " ".join(fastas))))
    if logger:
        logger.info(str("combining the following files matching pattern " +
                        "{0}:{1}".format(pattern, " ".join(fastas))))
    if len(fastas) == 0:
        if logger:
            logger.error("No matching files to combine found" +
                         " in {0}!".format(contigs_dir))
        sys.exit(1)
    with open(output, 'w') as w_file:
        for filen in fastas:
            with open(filen, 'r') as o_file:
                seq_records = SeqIO.parse(o_file, 'fasta')
                SeqIO.write(seq_records, w_file, 'fasta')

    return(output)


def setup_protein_blast(input_file, input_type="fasta", dbtype="prot",
                        title="blastdb", out="blastdb",
                        makeblastdb_exe='', logger=None):
    """
    This runs make blast db with the given parameters
    requires logging, os, subprocess, shutil
    """
    if makeblastdb_exe == '':
        makeblastdb_exe = shutil.which("makeblastdb")
    makedbcmd = str("{0} -in {1} -input_type {2} -dbtype {3} " +
                    "-title {4} -out {5}").format(makeblastdb_exe,
                                                  input_file,
                                                  input_type,
                                                  dbtype, title, out)
    if logger:
        logger.info("Making blast db: {0}".format(makedbcmd))
    try:
        subprocess.run(makedbcmd, shell=sys.platform != "win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE, check=True)
        logging.debug("BLAST database '{0}' created here: {1}".format(
            title, out))
        return(0)
    except:
        if logger:
            logging.error("Something bad happened when trying to make " +
                          "a blast database")
        sys.exit(1)


def run_blastp(input_file, database_name, outfmt, blastp_exe='', logger=None):
    """
    requires logging subprocess, os, shutil
    """
    #logger = logging.getLogger(name=None)
    # logger = logging.getLogger(__name__)
    output_file = os.path.join(os.path.split(input_file)[0],
                               str(os.path.splitext(
                                   os.path.basename(input_file))[0] +
                                   "_blast_hits.tab"))
    if blastp_exe == '':
        blastp_exe = shutil.which("blastp")
    blastpcmd = str("{0} -db {1} -query {2} -out {3} -outfmt " +
                    "{4}").format(blastp_exe, database_name, input_file,
                                  output_file, outfmt)
    if logger:
        logger.info("Running blastp: {0}".format(blastpcmd))
    try:
        subprocess.run(blastpcmd, shell=sys.platform != "win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE, check=True)
        if logger:
            logger.debug("Results from BLASTing {0} are here: {1}".format(
                input_file, output_file))
        return(0)
    except:
        if logger:
            logger.error("Something bad happened when running blast")
        sys.exit(1)


#def merge_outfiles(filelist, outfile_name):
def merge_blast_tab_outfiles(filelist, outfile_name, logger=None):
    """
    #TODO needs a test for headers
    for combining tab formated blast output format 6
    returns 0 if successful
    requires logging
    """
    # only grab .tab files, ie, the blast output
    # logger=logging.getLogger()
    filelist = [i for i in filelist if i.split(".")[-1:] == ['tab']]
    if len(filelist) == 1:
        if logger:
            logger.warning("only one file found! no merging needed")
        return(0)
    elif len(filelist) == 0:
        if logger:
            logger.error("filelist empt; cannot perform merge!")
        return(1)
    else:
        if logger:
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


def cleanup_output_to_csv(infile,
                          accession_pattern='(?P<accession>[A-Z _\d]*\.\d*)',
                          logger=None):
    """
    given .tab from merge_blast_tab_outfiles, assign pretty column names,
    """
    # if logger:
    #     logger=logging.getLogger(name=None)
    print("cleaning up the csv output")
    colnames = ["query_id", "subject_id", "identity_perc", "alignment_length",
                "mismatches", "gap_opens", "q_start", "q_end", "s_start",
                "s_end", "evalue", "bit_score"]
    csv_results = pd.read_csv(open(infile), comment="#", sep="\t",
                              names=colnames)
    #This default regex will probably break things eventually...
    # it looks for capital letter and numbers, dot, number, ie SHH11555JJ8.99
    csv_results["accession"] = csv_results.query_id.str.extract(accession_pattern)
    # write out results with new headers or with new headers and merged metadat from accessions.tab
    genes = open(genelist, "r")
    genedf = pd.read_csv(genes, sep=",")
    output_path_csv = str(os.path.splitext(infile)[0] + ".csv")
    results_annotated = pd.merge(csv_results, genedf, how="left",
                                 on="accession")
    results_annotated.to_csv(open(output_path_csv, "w"))
    print("wrote final csv to %s" % output_path_csv)
#%%


def md5(fname, string=False):
    """straight from quantumsoup
    http://stackoverflow.com/questions/3431825/
        generating-an-md5-checksum-of-a-file
    updated 20160916 for strings; if true, fname can be a string
    """
    if string:
        return(hashlib.md5(fname.encode('utf-8')).hexdigest())
    else:
        hash_md5 = hashlib.md5()
        with open(fname, "rb") as f:
            for chunk in iter(lambda: f.read(4096), b""):
                hash_md5.update(chunk)
            return(hash_md5.hexdigest())


def check_single_scaffold(input_genome_path):
    """Test for single scaffold.
    """
    print("testing for multiple scaffolds...")
    counter = -1  # if all goes well, returns 0, else returns 1 or more or -1
    for line in open(input_genome_path, "r"):
        if re.search("ORIGIN", line) is not None:
            counter = counter + 1
    return(counter)


def get_genbank_record(input_genome_path, check=False):
    """reads the FIRST record only from a genbank file;
        will probably only work for first scaffold
    """
    print("Reading genbank file...")
    with open(input_genome_path) as input_genome_handle:
        genome_seq_record = next(SeqIO.parse(input_genome_handle, "genbank"))
    # to avoid issues from working multiple times with open file handles, this
    # just reads it in fresh
    if check:
        with open(input_genome_path) as input_genome_handle:
            genome_seq = next(SeqIO.parse(input_genome_handle, "genbank")).seq
            if genome_seq[0: 100] == str("N" * 100):
                print("Careful: the first 100 nucleotides are N's; " +
                      "did you download the a truncated .gb file?")
                sys.exit(1)
    return(genome_seq_record)


def get_genbank_seq(input_genome_path):
    """get the sequence from the FIRST record only in a genbank file
    """
    print("fetching nucleotide sequence from genbank file...")
    with open(input_genome_path) as input_genome_handle:
        genome_seq_record = next(SeqIO.parse(input_genome_handle, "genbank"))
    return(genome_seq_record.seq)
