#!/bin/envpython

"""
unit tests for demux.py
"""


# pylint: disable=missing-docstring
# pylint: disable=redefined-outer-name

# what is this error for ?
# pylint: disable=protected-access

from __future__ import print_function

from io import BytesIO
from collections import OrderedDict
from collections import Counter
from distutils import dir_util
import os

import pytest

import eclipdemux_package.demux as dmx


@pytest.fixture()
def barcodes_fasta_file():
    file_content = '\n'.join([
        '>bcAG',
        'AAAGGG',
        '>bcCT',
        'CCCTTT'
    ])
    # fasta_file =
    return BytesIO(file_content)

@pytest.fixture()
def barcodes_tsv_file():
    file_content = '\n'.join([
        'AAGG\tbcAG',
        'CCTT\tbcCT'
    ])
    return BytesIO(file_content)

@pytest.fixture()
def barcodedict_example():
    #print("in barcodedict_example")
    return OrderedDict({"AAAGGG": "bcAG", "CCCTTT":"bcCT"})


@pytest.fixture()
def fastqio_toreadfrom():
    file_content = '\n'.join([
        '@D00611:409:H3HWKBCXY:1:1107:1276:2158 1:N:0',
        'AACATGTGTGAAGGGGACAG',
        '+',
        '<..<.<<<<GAAGIIGAIII',
        '@D00611:409:H3HWKBCXY:1:1107:1305:2159 1:N:0',
        'TCCTCACCTGGCACTGTCCC',
        '+',
        '<<G.<.AGA<GA<AAAA<.I'
    ])
    return BytesIO(file_content)


@pytest.fixture()
def fastqio_toreadfrom_interrupted():
    file_content = '\n'.join([
        '@D00611:409:H3HWKBCXY:1:1107:1276:2158 1:N:0',
        'AACATGTGTGAAGGGGACAG',
        '+',
        '<..<.<<<<GAAGIIGAIII',
        '@D00611:409:H3HWKBCXY:1:1107:1305:2159 1:N:0',
        'TCCTCACCTGGCACTGTCCC',
    ])
    return BytesIO(file_content)


@pytest.fixture()
def fastqio_towriteto():
    return BytesIO('')


@pytest.fixture()
def fastqread_example():
    return  {'title': 'D00611', 'sequence': 'ACGTACGT', 'quality': '<..<.<<<'}


@pytest.fixture()
def repeatscounter():
    counter = Counter()
    counter[('AAA', 'AAA', 'umi_1')] = 1
    counter[('AAA', 'AAA', 'umi_2')] = 2
    counter[('AAA', 'AAA', 'umi_3')] = 3
    counter[('AAA', 'AAC', 'umi_4')] = 1
    counter[('AAA', 'AAC', 'umi_5')] = 2
    counter[('AAA', 'AAC', 'umi_6')] = 3
    counter[('GGG', 'GGG', 'umi_7')] = 1
    counter[('GGG', 'GGG', 'umi_8')] = 2
    counter[('GGG', 'GGG', 'umi_9')] = 3
    counter[('GGG', 'GGT', 'umi_10')] = 1
    counter[('GGG', 'GGT', 'umi_11')] = 2
    counter[('GGG', 'GGT', 'umi_12')] = 3
    return counter

@pytest.fixture()
def datadir(tmpdir, request):
    '''
    Fixture responsible for searching a folder with the same name of test
    module and, if available, moving all contents to a temporary directory so
    tests can use them freely.
    '''
    filename = request.module.__file__
    test_dir, _ = os.path.splitext(filename)

    if os.path.isdir(test_dir):
        dir_util.copy_tree(test_dir, bytes(tmpdir))

    return tmpdir

###############################################################################
# fastq IOs
###########


def test__fastqread_from_fastqio_1(fastqio_toreadfrom):
    obtained = dmx._fastqread_from_fastqio(fastqio_toreadfrom)
    expected = {'title': 'D00611:409:H3HWKBCXY:1:1107:1276:2158 1:N:0',
                'sequence': 'AACATGTGTGAAGGGGACAG',
                'quality': '<..<.<<<<GAAGIIGAIII'}
    assert obtained == expected


def test__fastqread_from_fastqio_2(fastqio_toreadfrom):
    obtained1 = dmx._fastqread_from_fastqio(fastqio_toreadfrom)
    obtained2 = dmx._fastqread_from_fastqio(fastqio_toreadfrom)
    expected1 = {'title': 'D00611:409:H3HWKBCXY:1:1107:1276:2158 1:N:0',
                 'sequence': 'AACATGTGTGAAGGGGACAG',
                 'quality': '<..<.<<<<GAAGIIGAIII'}
    expected2 = {'title': 'D00611:409:H3HWKBCXY:1:1107:1305:2159 1:N:0',
                 'sequence': 'TCCTCACCTGGCACTGTCCC',
                 'quality': '<<G.<.AGA<GA<AAAA<.I'}
    assert obtained1 == expected1
    assert obtained2 == expected2


def test__fastqread_from_fastqio_3(fastqio_toreadfrom_interrupted):
    obtained1 = dmx._fastqread_from_fastqio(fastqio_toreadfrom_interrupted)
    obtained2 = dmx._fastqread_from_fastqio(fastqio_toreadfrom_interrupted)
    expected1 = {'title': 'D00611:409:H3HWKBCXY:1:1107:1276:2158 1:N:0',
                 'sequence': 'AACATGTGTGAAGGGGACAG',
                 'quality': '<..<.<<<<GAAGIIGAIII'}
    expected2 = None
    assert obtained1 == expected1
    assert obtained2 == expected2


def test__fastqread_to_fastqio(tmpdir):
    fastqfile = tmpdir.mkdir("test_fastqread_to_fastqio").join("towriteto.fastq")
    fastqdict = {'title': 'D00611:409:H3HWKBCXY:1:1107:1276:2158 1:N:0',
                 'sequence': 'AACATGTGTGAAGGGGACAG',
                 'quality': '<..<.<<<<GAAGIIGAIII'}
    dmx._fastqread_to_fastqio(fastqdict, fastqfile)
    obtained = fastqfile.read()
    expected = ('@D00611:409:H3HWKBCXY:1:1107:1276:2158 1:N:0\n'
                + 'AACATGTGTGAAGGGGACAG\n'
                + '+\n'
                + '<..<.<<<<GAAGIIGAIII\n')
    assert obtained == expected


###############################################################################
# barcode dictionary IOs
########################


@pytest.mark.xfail(raises=dmx.IllegalBarcodeIdValueError)
def test__assert_legal_barcodeid_1():
    dmx._assert_legal_barcodeid("NIL")


@pytest.mark.xfail(raises=dmx.IllegalBarcodeSeqValueError)
def test__assert_legal_barcodeseq_1():
    dmx._assert_legal_barcodeseq("ACGTNZ")


def test__barcodedict_from_tsvio(barcodes_tsv_file):
    dict_expected = OrderedDict()
    dict_expected["AAGG"] = "bcAG"
    dict_expected["CCTT"] = "bcCT"
    dict_produced = dmx._barcodedict_from_tsvio(barcodes_tsv_file)
    for item_expected, item_obtained in zip(dict_expected.items(),
                                            dict_produced.items()):
        assert item_obtained == item_expected


def test__barcodedict_from_fastaio(barcodes_fasta_file):
    dict_expected = OrderedDict()
    dict_expected["AAAGGG"] = "bcAG"
    dict_expected["CCCTTT"] = "bcCT"
    dict_produced = dmx._barcodedict_from_fastaio(barcodes_fasta_file)
    for item_expected, item_obtained in zip(dict_expected.items(),
                                            dict_produced.items()):
        assert item_obtained == item_expected


def test__barcodedict_add_edgecases():
    dic = OrderedDict()
    dic["AAAGGG"] = "bcAG"
    dic["CCCTTT"] = "bcCT"
    dmx._barcodedict_add_edgecases(dic)
    assert dic[""] == "NIL"
    assert dic[None] == "SHO"

def test__barcodedict_reorder():
    dic = OrderedDict()
    dic["AAAGGG"] = "bcAG_1_2"
    dic["CCCTTT"] = "bcCT_1"
    dic["GGGAAA"] = "bcGA_1_2_3"
    dic["TTTCCC"] = "bcTC"
    dic = dmx._barcodedict_reorder(dic)
    assert dic == OrderedDict([('GGGAAA', 'bcGA_1_2_3'),
                               ('AAAGGG', 'bcAG_1_2'),
                               ('CCCTTT', 'bcCT_1'),
                               ('TTTCCC', 'bcTC')])


def test_barcodedict_create_fasta(tmpdir):
    expected = OrderedDict()
    expected["AAACCC"] = "bcAC"
    expected["GGGTTT"] = "bcGT"
    #expected[''] = 'NIL'
    #expected[None] = 'SHO'
    fastafile = tmpdir.mkdir("test__barcodedict_create").join("barcodes.fasta")
    filecontent = '>bcAC\n' 'AAACCC\n' '>bcGT\n' 'GGGTTT\n'
    fastafile.write(filecontent)
    obtained = dmx.barcodedict_create(str(fastafile))
    assert obtained == expected

def test_barcodedict_create_tsv(tmpdir):
    expected = OrderedDict()
    expected["AAACCC"] = "bcAC"
    expected["GGGTTT"] = "bcGT"
    #expected[''] = 'NIL'
    #expected[None] = 'SHO'
    tsvfile = tmpdir.mkdir("test__barcodedict_create").join("barcodes.tsv")
    filecontent = 'AAACCC\t' 'bcAC\n' 'GGGTTT\t' 'bcGT\n'
    tsvfile.write(filecontent)
    obtained = dmx.barcodedict_create(str(tsvfile))
    assert obtained == expected


@pytest.mark.xfail(raises=dmx.BarcodesFileInvalidFormatValueError)
def test_barcodedict_create_fastq():
    dmx.barcodedict_create('basenamedoesnotmatter.fastq')


@pytest.mark.xfail(raises=dmx.BarcodesFileInvalidFormatValueError)
def test_barcodedict_create_csv():
    dmx.barcodedict_create('basenamedoesnotmatter.csv')


###############################################################################
# matching and hamming distance functions
#########################################


@pytest.mark.xfail(raises=dmx.IllegalBarcodeSeqValueError)
def test__assert_in_alphabet():
    dmx._assert_in_alphabet('Z')


def test__matching_strict():
    assert dmx._matching_strict("A", "A") == 1
    assert dmx._matching_strict("A", "N") == 0
    assert dmx._matching_strict("N", "A") == 0
    assert dmx._matching_strict("N", "N") == 0
    assert dmx._matching_strict("A", "C") == 0


def test__matching_loose():
    assert dmx._matching_loose("A", "A") == 1
    assert dmx._matching_loose("A", "N") == 1
    assert dmx._matching_loose("N", "A") == 1
    assert dmx._matching_loose("N", "N") == 1
    assert dmx._matching_loose("A", "C") == 0#first longer then second is ok


def test__hamming():
    assert dmx._hamming(str.__eq__, "ACGT", "ACGT") == 0
    assert dmx._hamming(str.__eq__, "ACGT", "ACTT") == 1
    assert dmx._hamming(str.__eq__, "ACGT", "TTTT") == 3
    assert dmx._hamming(str.__eq__, "ACGT", "ACGN") == 1
    assert dmx._hamming(str.__eq__, "ACGT", "NNNN") == 4
    assert dmx._hamming(str.__eq__, "ACGT", "TTNN") == 4
    assert dmx._hamming(str.__eq__, "ACNN", "ACNN") == 0


@pytest.mark.xfail(raises=dmx.HammingUnequalLengthsValueError)
def test__hamming_2():
    # first shorter then second is not ok
    dmx._hamming(str.__eq__, "123", "1234")


def test__hamming_3():
    dmx._hamming(str.__eq__, "1234", "1234")


def test__hamming_4():
    #first longer then second is ok
    dmx._hamming(str.__eq__, "1234", "123")


def test_hamming_strict():
    assert dmx.hamming_strict("ACGT", "ACGT") == 0
    assert dmx.hamming_strict("ACGT", "ACTT") == 1
    assert dmx.hamming_strict("ACGT", "TTTT") == 3
    assert dmx.hamming_strict("ACGT", "ACGN") == 1
    assert dmx.hamming_strict("ACGT", "NNNN") == 4
    assert dmx.hamming_strict("ACGT", "TTNN") == 4
    assert dmx.hamming_strict("ACNN", "ACNN") == 2


def test_hamming_loose():
    assert dmx.hamming_loose("ACGT", "ACGT") == 0
    assert dmx.hamming_loose("ACGT", "ACTT") == 1
    assert dmx.hamming_loose("ACGT", "TTTT") == 3
    assert dmx.hamming_loose("ACGT", "ACGN") == 0
    assert dmx.hamming_loose("ACGT", "NNNN") == 0
    assert dmx.hamming_loose("ACGT", "TTNN") == 2
    assert dmx.hamming_loose("ACNN", "ACNN") == 0


def test_barcodedict_distances(barcodedict_example):
    assert dmx.barcodedict_distances(barcodedict_example) \
           == [[6, 'AAAGGG', 'bcAG', 'CCCTTT', 'bcCT']]

def test_assignbarcode():
    dic = OrderedDict()
    dic["AAAGGG"] = "bcAG"
    dic["CCCTTT"] = "bcCT"
    assert dmx.assignbarcode("AAAGGGTTTTTTT", dic, 0) == 'AAAGGG'
    assert dmx.assignbarcode("CCCCCCTTTTTTT", dic, 0) == ''


###############################################################################
# fastq reads updates
#####################


@pytest.mark.xfail(raises=dmx.PairedEndReadsNameMisMatchValueError)
def test__assert_same_name():
    fastqread1 = {'title': 'D00611', 'sequence': 'ACGT', 'quality': '<..<'}
    fastqread2 = {'title': 'D00612', 'sequence': 'TCGA', 'quality': '.<<.'}
    dmx._assert_same_name(fastqread1, fastqread2)


def test__assert_can_trim_pass():
    fastqread1 = {'title': 'D00611', 'sequence': 'ACGT', 'quality': '<..<'}
    dmx._assert_can_trim(fastqread1, 3)


@pytest.mark.xfail(raises=dmx.SequenceTooShortToTrimValueError)
def test__assert_can_trim_fail():
    fastqread1 = {'title': 'D00611', 'sequence': 'ACGT', 'quality': '<..<'}
    dmx._assert_can_trim(fastqread1, 4)


def test__fastqread_trim_1():
    fastqread_input = {'title': 'D00611',
                       'sequence': 'ACGTACGT',
                       'quality': '<..<.<<.'}
    obtained = dmx._fastqread_trim(fastqread_input, 5)
    expected = {'title': 'D00611',
                'sequence': 'CGT',
                'quality': '<<.'}
    assert obtained == expected


def test__fastqread_trim_2():
    fastqread_input = {'title': 'D00611',
                       'sequence': 'ACGT',
                       'quality': '<..<'}
    obtained = dmx._fastqread_trim(fastqread_input, 8)
    expected = {'title': 'D00611',
                'sequence': '',
                'quality': ''}
    assert obtained == expected


def test__fastqread_titleprefix():
    fastqread1 = {'title': 'D00611',
                  'sequence': 'ACGT',
                  'quality': '<..<'}
    randomersequence1 = "acgtacgt"
    obtained = dmx._fastqread_titleprefix(fastqread1, randomersequence1)
    expected = {'title': 'acgtacgt:D00611',
                'sequence': 'ACGT',
                'quality': '<..<'}
    assert obtained == expected


def test__pairedreads_reformat_1():
    fastqread1 = {'title': 'D00611', 'sequence': 'ACGT', 'quality': '<..<'}
    fastqread2 = {'title': 'D00611', 'sequence': 'TCGA', 'quality': '<I><'}
    barcodelength = 1
    randomerlength = 2
    obtained = dmx._pairedreads_reformat(fastqread1, fastqread2,
                                         barcodelength, randomerlength)
    expected = (
        {'title': 'TC:D00611', 'sequence': 'CGT', 'quality': '..<'},
        {'title': 'TC:D00611', 'sequence': 'GA', 'quality': '><'},
        'A',
        'TC'
    )
    assert obtained == expected


def test__pairedreads_reformat_2():
    fastqread1 = {'title': 'D00611', 'sequence': 'ACGT', 'quality': '<..<'}
    fastqread2 = {'title': 'D00611', 'sequence': 'TCGA', 'quality': '<I><'}
    barcodelength = 0
    randomerlength = 2
    obtained = dmx._pairedreads_reformat(fastqread1, fastqread2,
                                         barcodelength, randomerlength)
    expected = (
        {'title': 'TC:D00611', 'sequence': 'ACGT', 'quality': '<..<'},
        {'title': 'TC:D00611', 'sequence': 'GA', 'quality': '><'},
        '',
        'TC'
    )
    assert obtained == expected

def test__pairedreads_reformat_3():
    fastqread1 = {'title': 'D00611', 'sequence': 'ACGT', 'quality': '<..<'}
    fastqread2 = {'title': 'D00611', 'sequence': 'TCGA', 'quality': '<I><'}
    barcodelength = 5
    randomerlength = 2
    obtained = dmx._pairedreads_reformat(fastqread1, fastqread2,
                                         barcodelength, randomerlength)
    expected = (
        {'title': 'TC:D00611', 'sequence': '', 'quality': ''},
        {'title': 'TC:D00611', 'sequence': 'GA', 'quality': '><'},
        None,
        'TC'
    )
    assert obtained == expected


###############################################################################
# fastq reads updates
#####################


# def test_cwlhandover_write(tmpdir):
#     testdir = tmpdir.mkdir("test_cwlhandover_write")
#     dmx.cwlhandover_write('dataset1', 'newname1', 'bcida', 'bcidb', str(testdir))
#     # assert testdir.join('dataset1').read() == 'dataset1'
#     # assert testdir.join('newname1').read() == 'newname1'
#     # assert testdir.join('bcida').read() == 'bcida'
#     # assert testdir.join('bcidb').read() == 'bcidb'


def test_cwlhandover_write(tmpdir):
    cwd = os.getcwd()
    os.chdir(str(tmpdir))
    dmx.cwlhandover_write('dataset1', 'newname1', 'bcida', 'bcidb')
    os.chdir(cwd)
    obtained = [tmpdir.join('dataset1').read(),
                tmpdir.join('newname1').read(),
                tmpdir.join('bcida').read(),
                tmpdir.join('bcidb').read()]
    assert obtained == ['dataset1', 'newname1', 'bcida', 'bcidb']


def test__demuxfilename():
    assert dmx._demuxfilename("s.RBFOX2", "IP.1", "r1", "X1A") \
        == "s.RBFOX2.IP.1.X1A.r1.fq.gz"
    assert dmx._demuxfilename("s.RBFOX2", "IP.1", "r2", "X1A") \
        == "s.RBFOX2.IP.1.X1A.r2.fq.gz"


###############################################################################
# umi repeats counting
######################


def test_repeatscounter_create():
    counter = dmx.repeatscounter_create()
    assert isinstance(counter, Counter)

def test_repeatscounter_write(repeatscounter, tmpdir):
    repeatscountfile = tmpdir.mkdir("test_repeatscounter_write")\
        .join("repeatscountfile.txt")
    dmx.repeatscounter_write(repeatscounter, str(repeatscountfile))
    obtained = repeatscountfile.read()

    expected = ('assigned\tfound\tumi\tcount\n' +
                'AAA\tAAA\tumi_1\t1\n' +
                'AAA\tAAA\tumi_2\t2\n' +
                'AAA\tAAA\tumi_3\t3\n' +
                'AAA\tAAC\tumi_4\t1\n' +
                'AAA\tAAC\tumi_5\t2\n' +
                'AAA\tAAC\tumi_6\t3\n' +
                'GGG\tGGG\tumi_7\t1\n' +
                'GGG\tGGG\tumi_8\t2\n' +
                'GGG\tGGG\tumi_9\t3\n' +
                'GGG\tGGT\tumi_10\t1\n' +
                'GGG\tGGT\tumi_11\t2\n' +
                'GGG\tGGT\tumi_12\t3\n')
    assert obtained == expected


def test_repeatscounter_increment(repeatscounter):
    dmx.repeatscounter_increment(repeatscounter, 'AAA', 'AAC', 'umi_4')
    dmx.repeatscounter_increment(repeatscounter, 'GGG', 'GGG', 'umi_9')
    dmx.repeatscounter_increment(repeatscounter, 'GGG', 'GGT', 'umi_13')
    expected = {
        ('AAA', 'AAA', 'umi_1'): 1,
        ('AAA', 'AAA', 'umi_2'): 2,
        ('AAA', 'AAA', 'umi_3'): 3,
        ('AAA', 'AAC', 'umi_4'): 2,
        ('AAA', 'AAC', 'umi_5'): 2,
        ('AAA', 'AAC', 'umi_6'): 3,
        ('GGG', 'GGG', 'umi_7'): 1,
        ('GGG', 'GGG', 'umi_8'): 2,
        ('GGG', 'GGG', 'umi_9'): 4,
        ('GGG', 'GGT', 'umi_10'): 1,
        ('GGG', 'GGT', 'umi_11'): 2,
        ('GGG', 'GGT', 'umi_12'): 3,
        ('GGG', 'GGT', 'umi_13'): 1
    }
    assert repeatscounter == expected


###############################################################################
# output file handling, core demux function
###########################################


def test_fqgz_open():
    pass


###############################################################################
# core demux function
#####################


def test_demux():
    pass


###############################################################################
# datastructures setup and wrapup
#################################

def test_datastructure_setup(tmpdir):
    tsvfile = tmpdir.mkdir("test__barcodedict_create").join("barcodes.tsv")
    filecontent = 'AAACCC\t' 'bcAC\n' 'GGGTTT\t' 'bcGT\n'
    tsvfile.write(filecontent)

    cwd = os.getcwd()
    os.chdir(str(tmpdir))
    _ = dmx.datastructure_setup(
        'NeedADataStructHere', 'dataset1' 'newname1', 'bcida', 'bcidb',
        str(tsvfile))
    os.chdir(cwd)

    # TODO


def test_datastructure_wrapup():
    # datastr =
    # metricsfile =
    # dmx.datastructure_wrapup(datastr, metricsfile)
    # TODO
    pass


###############################################################################
# module's main function
#############################

def test_main():
    pass


# TUTORIAL
##########
# def setup_module(module):
#     print ("setup_module      module:%s" % module.__name__)
#
# def teardown_module(module):
#     print ("teardown_module   module:%s" % module.__name__)
#
# def setup_function(function):
#     print ("setup_function    function:%s" % function.__name__)
#
# def teardown_function(function):
#     print ("teardown_function function:%s" % function.__name__)


# TUTORIAL
##########
# @pytest.fixture()
# def starting_next_test():
#     print('\nStarting next test')


# from pytest import fixture
# import os


# from __future__ import unicode_literals
# import sys
# sys.path.insert(0, '../eclipdemux/eclipdemux')


# from __future__ import unicode_literals
# from distutils import dir_util


# TUTORIAL
##########
# fixture directory
####################
# work out path for directory where test files ares tored
# FIXTURE_DIR = os.path.join(
#     os.path.dirname(os.path.realpath(__file__)),
#     'testfiles',
#     )


# @fixture
# def datadir(tmpdir, request):
#     '''
#     Fixture responsible for searching a folder with the same name of test
#     module and, if available, moving all contents to a temporary directory so
#     tests can use them freely.
#     '''
#     filename = request.module.__file__
#     test_dir, _ = os.path.splitext(filename)
#
#     if os.path.isdir(test_dir):
#         dir_util.copy_tree(test_dir, bytes(tmpdir))
#
#     return tmpdir

# @pytest.mark.datafiles(
#     os.path.join(FIXTURE_DIR, 'barcodes.fasta')
# )


# TUTORIAL
##########
if __name__ == '__main__':
    pytest.main()
