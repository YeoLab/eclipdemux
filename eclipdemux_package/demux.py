#!/usr/bin/env python

"""
Converts randomer + barcoded fastq files
into something that can be barcode collapsed and mapped
"""

# pylint: disable=fixme

# transitionning to python2/python3 support
# uncomment from this compatibility import list, as py3/py2 support progresses
from __future__ import print_function
from __future__ import division
# from __future__  import absolute_import
# from __future__  import unicode_literals
# from future import standard_library
# from future.builtins import builtins
# from future.builtins import utils

from collections import Counter
# from collections import defaultdict
from collections import OrderedDict
from collections import namedtuple
#from itertools import izip  Python 3 has the zip function built-in and it need not be imported via itertools.
import gzip
import os
from optparse import OptionParser

from future.utils import raise_with_traceback
from future.utils import iteritems


###############################################################################
# This module is organized in 8 sections:
#########################################
#
# 0) exceptions :
#
#     custom exceptions for invalid inputs to some of this module's function
#
# 1) barcode dictionary IOs :
#
#     file reading and data structure utilities to handle
#     a set of barcodes ids and sequences
#
# 2) fastq reads IOs :
#
#      file reading/writing to handle paired end fastq files,
#      data structure for an isolated pair of fastq reads
#      (both with title, sequence and quality strings)
#
# 3) matching and hamming distance functions
#
#      utilities to calculate matching distances between barcodes and reads,
#      to look for and decide best matching barcode or None,
#      sanity checks on fastq reads, handling of nomatch and tooshort
#
# 4) fastqread barcode assigning
#
#      fastqread barcode assigning
#
# 5) fastqread reformatting
#
#      this is the core motivation for this module
#      this includes utilities to :
#      rename paired end reads, adding randomer sequence at the title start
#      trim barcode from front of 1st read in pair
#      trim randomer from front of 2nd read in pair
#
# 6) umi repeats counting
#
#      umi repeats counting
#
# 7) output file handling, core demux functionx
#
#      core demux function
#      output files naming, opened files dict handling
#      handing over integer and string values for next cwl tool
#
# 8) this module's main function
#
#     main function
#     calling from command line
#


###############################################################################
# exceptions
############

class IllegalBarcodeIdValueError(ValueError):
    """raised when a barcode id is found that is a reserved word"""


class IllegalBarcodeSeqValueError(ValueError):
    """raised when a barcode seq is found that is not a string of A C G T N"""


class HammingUnequalLengthsValueError(ValueError):
    """raised when hamming dist called for strings of unequal lengths"""


class SequenceTooShortToTrimValueError(ValueError):
    """raised if sequence checked for trimming is shorter then trim length"""


class PairedEndReadsNameMisMatchValueError(ValueError):
    """raised if two paired end reads do no have identical names"""


class BarcodesFileInvalidFormatValueError(ValueError):
    """raised if barcodes_dict_create is called on non-fasta non-tsv file"""


###############################################################################
# barcode dictionary
####################


def _assert_legal_barcodeid(barcodeid):
    """
    Raise ValueError if barcodeid is a reserved word
    Args:
        barcodeid: a string, barcode id to be checked
    Returns: None
    """
    reserved_barcodeids = ["NIL", "SHO", "unassigned", "tooshort"]
    if any(barcodeid == word for word in reserved_barcodeids):
        raise_with_traceback(IllegalBarcodeIdValueError(barcodeid))


def _assert_legal_barcodeseq(barcodeseq):
    """
    Raise ValueError if barcodeseq is not composed of legal letters
    Args:
        barcodeseq: a string, barcode seq to be checked
    Returns: None
    """
    if any(ltr not in "ACGTN" for ltr in barcodeseq) or len(barcodeseq) == 0:
        raise_with_traceback(IllegalBarcodeSeqValueError(barcodeseq))


def _barcodedict_from_tsvio(tsvio):
    """
    makes a barcode dict from reading a tsv file, (used by legacy eclip)
    Args:
        tsvio: opened tsv file with barcodes,
        columns are: sequence, barcode-id
    Returns: a dict,with key / value being barcode-id / sequence
    """
    odict = OrderedDict()
    for line in tsvio:
        barcodeseq, barcodeid = line.strip().split("\t")
        _assert_legal_barcodeid(barcodeid)
        _assert_legal_barcodeseq(barcodeseq)
        odict[barcodeseq] = barcodeid
    return odict


def _barcodedict_from_fastaio(fastaio):
    """
    makes a barcode dict from reading a fasta file (used by new eclip)
    Args:
        fastaio: opened fasta file with barcodes
    Returns: dict with key/value being barcode sequence/id
    """
    odict = OrderedDict()
    lines = fastaio.readlines()
    for even, odd in zip(lines[0::2], lines[1::2]):
        barcodeid = even.strip()[1:]
        barcodeseq = odd.strip()
        _assert_legal_barcodeid(barcodeid)
        _assert_legal_barcodeseq(barcodeseq)
        odict[barcodeseq] = barcodeid
    return odict


def _barcodedict_add_edgecases(barcodedict):
    """
    add to barcodedicat key/values for unassigned and tooshort edgecases
    Args:
        barcodedict: dict,with key/value being barcode sequence/id
    Returns: None
    """
    # case of fastqread's unassigned to any barcode
    barcodedict[""] = "NIL"
    # case of fastqread's too short to trim (read1's bc, read2's randomer)
    barcodedict[None] = "SHO"


def _barcodedict_reorder(barcodedict):
    """
    reorder barcodes dictionanry,
    so as to start matching attempts with most expected barcodes first
    Args:
        barcodedict:  a dict with keys/values as barcode seq. / barcode id
    Returns: barcodedict, sorted by decr nbr of _ in barcode id
    """
    # barcodes should be sorted by their size
    # ie if they are annotated as present in my barcode file or not
    # if the barcodes are annotated put them first
    # ie they are expected to be demuxed
    # TODO legacy barcodedict reorder sorts by nb of _ in barcode id
    # TODO use expectedbarcodeida, expectedbarcodeidb to sort barcodedict
    barcodedict = OrderedDict(sorted(iteritems(barcodedict),
                                     key=lambda item: len(item[1].split("_")),
                                     reverse=True))
    return barcodedict


def barcodedict_create(barcode_filename):
    """
    creates a barcode dict from either a fasta or a tsv file
    Args:
        barcode_filename: fasta or tsv filename to read
    Returns: barcodedict, a dict with keys/values as barcode seq. / barcode id
    """
    file_ext = barcode_filename.split('.')[-1]
    if file_ext == "fasta":
        with open(barcode_filename, 'r') as fastaio:
            barcodedict = _barcodedict_from_fastaio(fastaio)
    elif file_ext == "tsv":
        with open(barcode_filename, 'r') as tsvio:
            barcodedict = _barcodedict_from_tsvio(tsvio)
    else:
        barcodedict = None  # because "local var might be refd b4 assignmt"
        raise_with_traceback(BarcodesFileInvalidFormatValueError(barcode_filename))
    _barcodedict_reorder(barcodedict)
    #_barcodedict_add_edgecases(barcodedict)
    return barcodedict


###############################################################################
# fastq reads IOs
#################


def _fastqread_from_fastqio(fastqio):
    """
    Args:
        fastqio: a stream (for ex opened file) with the content of a fastq file
    Returns: a dictionary object with attributes title, sequence, quality
    """
    title = fastqio.readline().strip()[1:]  # discard the initial '@' letter
    sequence = fastqio.readline().strip()
    _ = fastqio.readline().strip()  # title may be repeated, but ignore it
    quality = fastqio.readline().strip()

    # if one or more read line was  "", EOF was reached, return None
    if "" in [title, sequence, quality]:
        fastqread = None
    else:
        fastqread = {'title': title, 'sequence': sequence, 'quality': quality}
    return fastqread


def _fastqread_to_fastqio(fastqread, fastqio):
    """
    Args:
        fastqread: fastq dict (with title, sequence and quality properties)
        fastqio: a stream (ex opened file) with the content of a fastq file
    Returns: None
    """
    fastqio.write((
        '@' + fastqread['title'] + '\n' +
        fastqread['sequence'] + '\n' +
        '+' + '\n' +
        fastqread['quality'] + '\n'
    ).encode())                         #File was opened in binary format so cannot directly write strings to it
    return fastqio


###############################################################################
# barcode matching
##################


def _assert_in_alphabet(nucl):
    """
    Raise ValueError if nucl not one of 'A', 'C', 'G', 'T', 'N',
    lowercase not accepted
    Args:
        nucl: one of 'A', 'C', 'G', 'T', 'N'
    Returns: None
    """
    alphabet = ['A', 'C', 'G', 'T', 'N']
    if nucl not in alphabet:
        raise_with_traceback(IllegalBarcodeSeqValueError(nucl))


def _matching_strict(nucl1, nucl2):
    """
    checks if 2 nucleotides strictly match, N matches nothing, even not N
    Args:
        nucl1: one of 'A', 'C', 'G', 'T', 'N'
        nucl2: one of 'A', 'C', ''G', 'T', 'N'
    Returns: boolean, true if matching strictly
    """
    _assert_in_alphabet(nucl1)
    _assert_in_alphabet(nucl2)
    return nucl1 == nucl2 and nucl1 != "N"


def _matching_loose(nucl1, nucl2):
    """
    checks if 2 nucleotides loosely match, allowing N (mask) to match anything
    Args:
        nucl1: one of 'a', 'c', ''g', 't', 'N'
        nucl2: one of 'a', 'c', ''g', 't', 'N'
    Returns: boolean, true if matching loosely
    """
    _assert_in_alphabet(nucl1)
    _assert_in_alphabet(nucl2)
    return nucl1 == "N" or nucl2 == "N" or nucl1 == nucl2


def _hamming(matching_function, word1, word2):
    """
    nb of miss-matches between 2 words, comparing letters via matching_function
    word1 and word2 must have same length
    Args:
        matching_function: a comparison function for a pair of letters
        word1: a sequence string, using alphabet: A, C, G, T and N
        word2: another sequence string using alphabet: A, C, G, T and N
    Returns: integer, nbr of miss-matches
    """
    # TODO clean this
    # only throw exception if first word shorter then second
    # because we should always have bcfound as long as the longest barcode
    if len(word1) < len(word2):
        raise_with_traceback(HammingUnequalLengthsValueError(word1, word2))
    return sum(not matching_function(letter1, letter2)
               for letter1, letter2 in zip(word1, word2))  #Changed it to zip - Python 3 specifications


def hamming_strict(word1, word2):
    """
    hamming distance between two same-length words, with N matching nothing
    This is only for the sake of testing and contrasting with loose_hamming
    Args:
        word1: a sequence string, using alphabet: A, C, G, T and N
        word2: another sequence string using alphabet: A, C, G, T and N
    Returns: an integer, the hamming distance between the 2 sequences,
    total number of positions mismatched, N being a match to any letter
    """
    return _hamming(_matching_strict, word1, word2)


def hamming_loose(word1, word2):
    """
    hamming distance between two same-length words, with N matching anything
    Ns should appear similar to other barcodes as this makes the results more
    stringent, not less
    Args:
        word1: a sequence string, using alphabet: A, C, G, T and N
        word2: another sequence string using alphabet: A, C, G, T and N
    Returns: an integer, the hamming distance between the 2 sequences,
    total number of positions mismatched, N being a match to any letter
    """
    return _hamming(_matching_loose, word1, word2)


def barcodedict_distances(bcdict):
    """
    list strict hamming distances for all pairs of barcodes in a barcodes dict
    This is not used by the main demux function
    It helps checking how different our standard Yeo Lab barcodes are
    Args:
        bcdict: a dict,with key/value being barcode sequence/id
    Returns: a sorted list of increasing distances between barcodes
    """
    distances_triplets = []
    for bcseq1 in bcdict:
        for bcseq2 in bcdict:
            if bcseq1 < bcseq2:
                distances_triplets.append(
                    [hamming_loose(bcseq1, bcseq2),
                     bcseq1, bcdict[bcseq1], bcseq2, bcdict[bcseq2]]
                )
    return sorted(distances_triplets)


###############################################################################
# barcode assigning
###################


def assignbarcode(readsequence, barcodedict, allowed_distance=0):
    """
    Checks if fastqread has one of the barcodes in barcodedict and returns it
    Args:
        fastqread: str, read to check if barcode exists in
        barcodedict: dict of barcode sequence / barcode id
        allowed_distance: integer, max hamming distance between given barcode
        and barcode in read, defaults to 0
    Returns: best barcode match in read and the real barcode sequence,
    none if none found
    """

    ###########################################################################
    # TODO "for each barcode length see if known barcodes appear"   really?
    ###########################################################################
    # see:  (a single end version)
    # https://github.com/YeoLab/gscripts/blob/master/gscripts/clipseq/demultiplex_barcoded_fastq.py
    #
    # but direct anceestor of this script is actually:
    # https://github.com/YeoLab/gscripts/blob/master/gscripts/clipseq/demux_paired_end.py
    ###########################################################################
    # assumes larger barcodes are less likely, and searches for them first
    #
    # barcode_lengths = {len(barcode) for barcode in barcodes.keys()}
    # for barcode_length in sorted(barcode_lengths, reverse=True):
    # cur_barcode =  bla bla bla
    #     if cur_barcode in barcodes:
    #         barcode = cur_barcode
    #         break
    ###########################################################################

    # TODO cleanup
    # print('readsequence: ', readsequence)
    # print('barcodedict: ', barcodedict)
    # print('allowed_distance: ', allowed_distance)

    longest_barcode_length = max([len(bcseq)
                                  for bcseq in barcodedict if bcseq])
    bcfound = readsequence[:longest_barcode_length]

    # Gets min hamming distance between barcode in read and barcode in list,
    # assigned barcode non deterministic if tie on min distance
    distances_to_bcseqs = [(hamming_loose(bcfound, bcseq), bcseq)
                           for bcseq in barcodedict]
    closest_dist, closest_bcseq = min(distances_to_bcseqs, key=lambda x: x[0])
    if closest_dist <= allowed_distance:
        assigned_barcode = closest_bcseq
    else:
        assigned_barcode = ''
    return assigned_barcode  # , barcode_found


###############################################################################
# fastq reads reformatting
##########################

def _assert_same_name(fastqread1, fastqread2):
    """
    Raise PairedEndReadsNameMisMatchValueError if the 2 paired end fastq reads
    have mis-matching names (first part of their title string)
    Args:
        fastqread1: fastq dict (with title, sequence and quality properties)
        fastqread2: fastq dict for other end pairing with fastqread1
    Returns: None
    """
    # print('fastqread1', fastqread1)
    # print('fastqread2', fastqread2)
    fastqread1_name = fastqread1["title"].split()[0]
    fastqread2_name = fastqread2["title"].split()[0]
    if fastqread1_name != fastqread2_name:
        raise_with_traceback(PairedEndReadsNameMisMatchValueError(
            fastqread1_name, fastqread2_name))


def _assert_can_trim(fastqread, trimlength):
    if len(fastqread["sequence"]) <= trimlength:
        raise_with_traceback(
            SequenceTooShortToTrimValueError(fastqread, trimlength))


def _fastqread_trim(fastqread, trimlength):
    """
    trim trimlength letters form the front of both sequence and quality
    string attributes in fastqread dictionary
    Args:
        fastqread: fastq dict (with title, sequence and quality properties)
        trimlength: integer, nbr bases to cut from both sequence and quality

    Returns: fastqread dict with shortened sequence and quality strings
    """
    # _assert_can_trim(fastqread, trimlength)  # Not raising error any more
    fastqread['sequence'] = fastqread['sequence'][trimlength:]
    fastqread['quality'] = fastqread['quality'][trimlength:]
    return fastqread


def _fastqread_titleprefix(fastqread, randomersequence):
    """
    modify title of fastqread to include randomersequence at the beginning
    Args:
        fastqread: a fastq dict (with title, sequence and quality properties)
        randomersequence: a sequence string
    Returns: fastqread dict with modified title
    """
    # first letter in title is normally a '@' but we'll just force one there
    fastqread["title"] = randomersequence + ":" + fastqread["title"]
    return fastqread


def _pairedreads_reformat(fastqread1, fastqread2,
                          barcodelength, randomerlength):
    """
    reformat paired fastqread's r=to include randomesequnece in title
    and cut barcode from front of read1
    and cutting randomer from front of read2
    Args:
        fastqread1: fastq dict (with title, sequence and quality properties)
        fastqread2: fastq dict for other end pairing with fastqread1
        barcodelength: length of matched barcode to trim at front of fastqread1
        randomerlength: length of randomer to trim at front of fastqread2
    Returns: 4-tuple of moddigied paired dastqread's re-titled and trimmed
        fastqread1:
        fastqread2:
        barcodesequencefound: actual barcode sequence found at front of read1
        randomersequence: randomer sequence identified
    """
    bcseqfound = fastqread1["sequence"][0:barcodelength]
    randomerseq = fastqread2["sequence"][0:randomerlength]

    fastqread1 = _fastqread_titleprefix(fastqread1, randomerseq)
    fastqread2 = _fastqread_titleprefix(fastqread2, randomerseq)

    fastqread1 = _fastqread_trim(fastqread1, barcodelength)
    fastqread2 = _fastqread_trim(fastqread2, randomerlength)

    # either read is too short for trimming, and by convention :
    # bcseqfound and bcseqassigned (outside this call) should be set to None
    if len(fastqread1["sequence"]) == 0 or len(fastqread2["sequence"]) == 0:
        bcseqfound = None

    # np barcode sequence assigned:
    # if barcodelength is 0, bcseqfound will be "" and that is fine

    return fastqread1, fastqread2, bcseqfound, randomerseq


###############################################################################
# umi repeats counting
######################


def repeatscounter_create():
    """
    initiate a counter of reads repeats (same barcode/realbarcode/umi)
    umi is aka randomer sequence
    Args:
        #barcodedict: dict of barcode sequence / id
    Returns:
        Object: a 3-tier dict to count repeats
    """
    # each key is a barcode sequence
    # the associated value 3-levels-of-keys dict: [bcassigned][bcfound][umi]
    # and final leaf values are counts of reads
    # umi is aka randomer sequence

    # repeatscounter = {}
    # for bcseq in barcodedict:
    #     repeatscounter[bcseq] = defaultdict(defaultdict(defaultdict(0)))
    # repeatscounter[""] = defaultdict(Counter)         # no barcode assigned
    # repeatscounter[None] = defaultdict(Counter)       # to short to trim
    # return repeatscounter
    return Counter()


def repeatscounter_increment(repeatscounter, barcode, bcfound, umi):
    """
    increment repeatscounter for one additional read seen with the
    combination barcode / bcfound / umi
    Args:
        repeatscounter: data structure used to count repeats
        barcode: barcode
        bcfound: actual barcode found at front of read 1
        umi: randomer sequence
    Returns: None
    """
    # print("incrementing", barcode, bcfound, umi)
    repeatscounter[(barcode, bcfound, umi)] += 1
    return None


def repeatscounter_write(repeatscounter, metrics_file_name):
    """
    write to metricsfilename data from 3-tier dict bc_bcfound_rand_cnts_dict
    Args:
        repeatscounter: 3-tier dict
        metrics_file_name: file to write to
    Returns: None
    """
    with open(metrics_file_name, 'w') as opened_metrics_file:
        opened_metrics_file.write(
            "%s\t%s\t%s\t%s\n" %
            ("assigned", "found", "umi", "count"))
        # TODO do not sort
        for (bcassigned, bcfound, umi), count \
                in sorted(iteritems(repeatscounter)):
            opened_metrics_file.write(
                "%s\t%s\t%s\t%s\n" % (bcassigned, bcfound, umi, count))

        # for bcassigned, bcfounds \
        #         in repeatscounter.items():
        #     for bcfound, randomers in bcfounds.items():
        #         for randomer, count in randomers.items():
        #             opened_metrics_file.write(
        #                 "%s\t%s\t%s\t%s\n" %
        #                 (bcassigned, bcfound, randomer, count))


###############################################################################
# outputs file handling
#######################


def _fqgz_open(fastq_or_fastqgz_filename):
    """
    open a fastq or fastq.gz file with appropriate opener
    Args:
        fastq_or_fastqgz_filename: name of file to open
    Returns: opened file handle
    """
    if os.path.splitext(fastq_or_fastqgz_filename)[1] == ".gz":
        opened_fastq_io = gzip.open(fastq_or_fastqgz_filename)
    else:
        opened_fastq_io = open(fastq_or_fastqgz_filename)
    return opened_fastq_io


def cwlhandover_write(dataset, newname,
                      expectedbarcodeida, expectedbarcodeidb, direct='./'):
    """
    hand over these strings to next cwl tool (parsebarcodes)
    done here instead of in cwl InitialWorkDirRequirement
    because toil does not support InitialWorkDirRequirement yet
    Args:
        dataset:
        newname:
        expectedbarcodeida:
        expectedbarcodeidb:
    Returns: None
    """
    with open(direct + newname, mode='w') as newname_fh:
        newname_fh.write(newname)
    with open(direct + dataset, mode='w') as newname_fh:
        newname_fh.write(dataset)
    with open(direct + expectedbarcodeida, mode='w') as newname_fh:
        newname_fh.write(expectedbarcodeida)
    with open(direct + expectedbarcodeidb, mode='w') as newname_fh:
        newname_fh.write(expectedbarcodeidb)


def _demuxfilename(dataset, newname, r1_or_r2, barcodeid):
    """
    construct filename based on  barcodeid and read parity
    Args:
        dataset:
        newname:
        r1_or_r2:
        barcodeid:
    Returns: file name
    """
    splitfilename = (dataset.split(".") +  # could be a singleton list
                     newname.split(".") +  # could be a singleton list
                     [barcodeid, r1_or_r2, "fq", "gz"])
    jointfilename = ".".join(splitfilename)
    return jointfilename


def _demuxdicts_demuxedfiles_opener(demuxeddict1, demuxeddict2, dataset, newname):
    """
    based on 2 demuxeddicts, make a function that
    takes two arguments of the seq and id for a barcode
    and appends to the 2 demuxeddicts 2 open files dedicated
    to collect reads assigned to that barcode
    Args:
        demuxeddict1:
        demuxeddict2:
        dataset:
        newname:
    Returns: an appending function,
    """
    def demuxedfiles_open(bcseq, bcid):
        """
        Appends 2 open files for saving reads assigned with barcode bcseq/bcid
        This function is dynamically generated and has "clojured" values
        (demuxeddict1, demuxeddict2, dataset, newname)
        Args:
            bcseq: barcode sequence
            bcid:  barcode id
        Returns: None
        """
        filename1 = _demuxfilename(dataset, newname, "r1", bcid)
        demuxeddict1[bcseq] = gzip.open(filename1, 'w')
        filename2 = _demuxfilename(dataset, newname, "r2", bcid)
        demuxeddict2[bcseq] = gzip.open(filename2, 'w')

    return demuxedfiles_open


def _demuxeddicts_create(barcodedict, dataset, newname):
    """
    creates dicts of barcode seq / demuxed filehandles to assign reads to
    Args:
        barcodedict: dict of barcode sequence / id
        dataset:
        newname:
    Returns: 1 dict for read1's and 1 dict for read2's
        keys are each barcode sequence,
        each value is an open file to write read1's assigned to a barcode seq.
    """
    demuxeddict1, demuxeddict2 = dict(), dict()

    demuxedfiles_opener = _demuxdicts_demuxedfiles_opener(
        demuxeddict1, demuxeddict2, dataset, newname)
    for bcseq, bcid in barcodedict.items():
        demuxedfiles_opener(bcseq, bcid)
    # for reads unassigned to any barcode
    demuxedfiles_opener("", "NIL")
    # for reads too short for trimming
    demuxedfiles_opener(None, "SHO")

    return demuxeddict1, demuxeddict2


def demuxeddicts_close(demuxeddict1, demuxeddict2):
    """
    close all opened demuxed files
    Args:
        demuxeddict1: dict of barcode seq / open demuxed read1's file
        demuxeddict2: dict of barcode seq / open demuxed read1's file
    Returns: None
    """
    for fwd_file in demuxeddict1.values():
        fwd_file.close()
    for rev_file in demuxeddict2.values():
        rev_file.close()


###############################################################################
# datastructures setup and wrapup
#################################

DemuxDataStructure = namedtuple(
    'DemuxDataStructure',
    'barcodedict, repeatscounter, demuxeddict1, demuxeddict2')

def datastructure_setup(dataset, newname, bcida, bcidb, barcodesfile):
    """
    this module's datastructures setup function
    Args:
        dataset: dataset name to prefix all demuxed file names with
        newname: root name to use for demuxed files
        bcida: expected barcodeid a
        bcidb: expected barcodeid b
        barcodesfile: reference barcodes file either a fasta or tsv
    Returns: a four tuple of
        barcodedict: dict of barcode sequence / id
        repeatscounter:
        demuxeddict1:
        demuxeddict2:
    """
    cwlhandover_write(dataset, newname, bcida, bcidb)
    barcodedict = barcodedict_create(barcodesfile)
    repeatscounter = repeatscounter_create()
    demuxeddict1, demuxeddict2 = _demuxeddicts_create(barcodedict,
                                                      dataset,
                                                      newname)
    return DemuxDataStructure(
        barcodedict=barcodedict,
        repeatscounter=repeatscounter,
        demuxeddict1=demuxeddict1,
        demuxeddict2=demuxeddict2
        )


def datastructure_wrapup(datastr, metricsfile):
    """
    this module's datastructures wrapup function
    Args:
        datastr: demuxdatastructure
                 ('barcodedict, repeatscounter, demuxeddict1, demuxeddict2')
        metricsfile: filename to save metrics to
    Returns: None
    """
    demuxeddicts_close(datastr.demuxeddict1, datastr.demuxeddict2)
    repeatscounter_write(datastr.repeatscounter, metricsfile)


###############################################################################
# core demux function
#####################


def demux(datastr, fastq1, fastq2, randomer_length, allowed_distance):
    """
    demultiplex paired end fastq files named fastq1 and fastq2
    Args:
        datastr: demuxdatastructure:
                 ('barcodedict, repeatscounter, demuxeddict1, demuxeddict2')
        fastq1: file name of first in paired ends fastq
        fastq2: file name of first in paired ends fastq
        randomer_length: Length of randomer sequence at front of second read
        allowed_distance: Max Hamming distance between read barcode
            and given barcodes to assign a read to a given barcode
    Returns: None
    """
    # TODO cleanup
    # print('now demuxing')
    # print('datastr', datastr)
    # print('fastq1', fastq1)
    # print('fastq2', fastq2)
    # print('randomer_length', randomer_length)
    # print('allowed_distance', allowed_distance)

    # TODO do not use globals DEMUXEDDICT1 /2 in demux
    # TODO demux should return the reformatted fastqdicts ?
    with _fqgz_open(fastq1) as fastqfile1, _fqgz_open(fastq2) as fastqfile2:
        while True:

            fastqread1 = _fastqread_from_fastqio(fastqfile1)
            fastqread2 = _fastqread_from_fastqio(fastqfile2)
            if fastqread1 is None or fastqread2 is None:
                break

            _assert_same_name(fastqread1, fastqread2)

            bcassigned = assignbarcode(fastqread1["sequence"],
                                       datastr.barcodedict,
                                       allowed_distance)

            bc_length = len(bcassigned) if bcassigned else 0

            fastqread1, fastqread2, bcfound, umi = _pairedreads_reformat(
                fastqread1, fastqread2, bc_length, randomer_length)

            # unassigned to any barcode
            # TODO legacy "unassigned" id, but bcassigned/bcfound/umi =?
            if bcfound == "":
                bcassigned, bcfound, umi = "", "", umi

            # too short for trimming
            # TODO legacy "too-short" id, but bcassigned/bcfound/umi =?
            if bcfound is None:
                bcassigned, bcfound, umi = None, None, umi

            repeatscounter_increment(datastr.repeatscounter,
                                     bcassigned, bcfound, umi)

            #  TODO write "tooshort" reads to demuxed "SHO" file ?
            #  write unassigned reads to demuxed "NIL"/"unassigned" file
            _fastqread_to_fastqio(fastqread1,
                                  datastr.demuxeddict1[bcassigned])
            _fastqread_to_fastqio(fastqread2,
                                  datastr.demuxeddict2[bcassigned])



###############################################################################
# this module's main function
#############################


def main():
    """
    demux module's main function
    Returns: None
    """

    usage = "takes raw fastq files and " + \
    "demultiplex inline randomer and adapter sequences"
    parser = OptionParser(usage)

    parser.add_option("--dataset", dest="dataset")
    parser.add_option("--newname", dest="newname")
    # parser.add_option("--out_file_1", dest="out_file_1")
    # parser.add_option("--out_file_2", dest="out_file_2")
    parser.add_option("-m", "--metrics", dest="metricsfile")

    parser.add_option(
        "--fastq_1", dest="fastq1", help="fastq file to barcode demux")
    parser.add_option(
        "--fastq_2", dest="fastq2", help="fastq file to barcode demux")

    parser.add_option("--expectedbarcodeida", dest="expectedbarcodeida")
    parser.add_option("--expectedbarcodeidb", dest="expectedbarcodeidb")
    parser.add_option("-b", "--barcodesfile", dest="barcodesfile",
                      help="reference barcodes fasta or tsv file "
                      "(tsv columns: barcode_sequence/barcode_id)")

    parser.add_option("--length", dest="randomer_length",
                      help="Length of randomer sequence at front of 2nd read",
                      type=int,
                      default=10)  # TODO randomer_length legacy default was 3
    parser.add_option("--max_hamming_distance", dest="allowed_distance",
                      help="Max Hamming distance between read barcode and "
                      "given barcodes to assign a barcode to a read",
                      type=int,
                      default=1)
    parser.add_option("--legacy", dest="legacy",
                      action="store_true", default=False,
                      help="")

    (opts, _) = parser.parse_args()   # args unused

    # BARCODEDICT, REPEATSCOUNTER, DEMUXEDDICT1, DEMUXEDDICT2 = \
    datastr = \
        datastructure_setup(
            opts.dataset, opts.newname,
            opts.expectedbarcodeida, opts.expectedbarcodeidb,
            opts.barcodesfile)

    demux(datastr,
          opts.fastq1, opts.fastq2,
          opts.randomer_length, opts.allowed_distance)

    datastructure_wrapup(datastr,
                         opts.metricsfile)


if __name__ == "__main__":
    main()

###############################################################################
# extras
########

# debugging  barcodedicts

# ------------------------

# import json
# def dump_dict(dict, filename):
#     with open(filename, 'w') as openfile:
#         json.dump(dict, openfile)

# namedtuple stuff
# -----------------

# from collections import namedtuple

# FastqRead = namedtuple('FastqRead', 'title sequence plus quality')

# def next_fastqread_from_fastqio(fastqio):
#     return FastqRead(
#         title=fastqio.readline(),
#         sequence= fastqio.readline(),
#         plus=fastqio.readline(),
#         quality=fastqio.readline()
#         )

# spare docstring stuff
# ----------------------
# metricsfile: file name to write metrics to
# fastq1: fastq file to barcode demux
# fastq2: fastq file to barcode demux
#
#
# randomer_length: Length of randomer sequence at front of second read
# allowed_distance: Max Hamming distance between read barcode
#     and given barcodes to assign a read to a given barcode
