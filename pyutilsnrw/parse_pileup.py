import os
import sys
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

"""  There is no easy way to make a consunsus from a samtools
pileup with low coverage.  Here is my minimal solution
"""

def check_samtools_pileup(pileup):
    """checks that the file is a valid samtools pileup. outputs a list
    """
    res = []
    try:
        with open(pileup, "r") as file:
            for f in file:
                f = re.split(r'\t+', f.strip())
                if len(f) != 6:
                    print("6 columns not found")
                if not isinstance(int(f[1]), int):
                    print("second column isnt int")
                res.append(f)
    except:
        raise ValueError("Error with reading pileup file")
    return(res)
def reconstruct_seq(refpath, pileup, verbose=True, veryverb=False,
                      logger=None):
    """ This is a bit of a mess, to say the least.  Given a list from
    check samtools pileup, and a reference fasta, this reconstructs ambiguous
    regions
    """
    if verbose and logger:
        log_status = logger.info
    elif verbose or veryverb:
        log_status = print
    else:
        pass
    if verbose:
        log_status(str("reconstucting consensus sequence " +
                       "from {0} and pileup").format(refpath))
    seqfile = SeqIO.parse(open(refpath, "r"), "fasta")
    for i in seqfile:
        ref = str(i.seq)
    ref = "${0}".format(ref) # make 0 based
    new = ""
    skip = 0
    indels = 0
    N_deletions, N_insertions = 0, 0
    insert = re.compile('[,\\.]{0,1}\\+[0-9]+[ACGTNacgtn]+')
    delete = re.compile('[\\.,]{0,1}-[0-9]+[ACGTNacgtn]+')
    j = 0  # counter for pileup
    for i in range(0, len(ref)):  # counter for ref index
        if veryverb and verbose:
            try:
                log_status("ref. index: %i\n pile index: %i" %(i,j))
                log_status(ref[i])
                log_status(pileup[j])
            except:
                pass
        # This is how we handle deletions; decrement skip, and next iteration
        if skip > 0:
            skip = skip - 1
            if verbose:
                log_status("skipping {0}".format(i))
        # if j is greater then length of pileup, go with ref.
        # This should avoid out of range issues
        elif j > len(pileup)-1:
            new = "".join([new, ref[i]])
        #  if index isnt in second col of pileup, skip, filling with ref
        # note that because the reference is now zero base, no correction needed
        elif i != int((pileup[j][1])):
            if verbose:
                log_status("no entry in pileup for %i" % i)
            new = "".join([new, ref[i]])
            # this should keep pileup counter the same when
            j = j - 1
        # if (N in ref, or pilup differs from ref), and
        # (pilup has single value, or all values same) and
        # pileup char isnt $*, go with pileup
        # *(start char is ^W, which is two chars, breaks 3rd line of conditions
        # NOTE: lowercase letters converted to upper, because orientation
        #       is already handled by samtools.
        elif (ref[i] == "N" or pileup[j][4][0] != ref[i]) and \
             (len(pileup[j][4]) == 1 or \
              all(x == pileup[j][4][0] for x in list(pileup[j][4]))) and \
             pileup[j][4][0] not in [",", ".", "^", "$"]:
            new = "".join([new, pileup[j][4][0].upper()]) # append  upper
        # This is tp handle insetions;  could use a lamda?
        elif re.match(insert, pileup[j][4]) is not None and \
            all([hits == re.findall(insert, pileup[j][4])[0] for hits in \
                 re.findall(insert, pileup[j][4])]):
            if verbose:
                log_status("found insert!")
            insert_seq = re.search('[ACGTNacgtn]+',  pileup[j][4]).group(0)
            insert_N = int(re.search('[0-9]+',  pileup[j][4]).group(0))
            if not len(insert_seq) == insert_N:
                raise ValueError("error parsing insert")
            new="".join([new, insert_seq])
            indels = indels + insert_N
            N_insertions = N_insertions + insert_N
        # deletions
        elif re.match(delete, pileup[j][4]) is not None and \
            all([hits == re.findall(insert, pileup[j][4])[0] for hits in \
                 re.findall(insert, pileup[j][4])]):
            if verbose:
                log_status("found deletion! {0}".format(pileup[j][4]))
            delete_N = int(re.search('[0-9]+',  pileup[j][4]).group(0))
            skip = delete_N
            indels = indels + delete_N
            N_deletions = N_deletions + delete_N
        # Most cases fall in this category
        elif pileup[j][4][0] in [",", ".", "^", "$"] or \
            not all(x == pileup[j][4][0] for x in list(pileup[j][4])):
            if verbose:
                log_status("using ref")
            new="".join([new, ref[i]])
        else:
            if verbose:
                log_status("Case Not covered!")
            sys.exit(1)
        j = j + 1  # increment the pileup counter
    if verbose:
        log_status(str("total indels: {0}\n\tdeletions {1}\n\tinsetions: "
                        +"{2}").format(indels, N_deletions, N_insertions))
    else:
        print(str("total indels: {0}\n\tdeletions {1}\n\tinsetions: "
              +"{2}").format(indels, N_deletions, N_insertions))
    return(new[1:])  # [1:] gets rid of starting dollar character


if __name__ == "__main__":
	if sys.argv[1] in ["-h", "h", "--help"]:
		print("USAGE: parse_pileup.py pathToSamtoolsPileup.txt pathToRefSequence.fasta")
		sys.exit(0)
    path = sys.argv[1]
    reference = sys.argv[2]
    print(sys.argv[1:])
    pi = check_samtools_pileup(pileup=path)
    with open (reference, "r") as reffile:
        seq = SeqIO.parse(reffile, "fasta")
        for i in seq:
            refseq = i.seq
    print(reference)
    newseq = reconstruct_seq(refpath=reference, pileup=pi, verbose=True, veryverb=True)
    with open("output.fasta", 'w') as out:
        SeqIO.write(SeqRecord(Seq(newseq,IUPAC.IUPACAmbiguousDNA()),
                              id = "contigs_consensus",
                              description=""), out, 'fasta')
