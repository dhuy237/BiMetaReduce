from mrjob.job import MRJob
from mrjob.protocol import RawValueProtocol
from mrjob.protocol import TextProtocol
from mrjob.protocol import JSONProtocol

import itertools as it
from Bio.Seq import Seq
from multiprocessing import Pool, Array, Value
from gensim import corpora
import numpy as np

import sys

sys.path.append("../")  # Add "../" to utils folder path
from utils import globals


def gen_kmers(klist):
    bases = ["A", "C", "G", "T"]
    kmers_list = []

    for k in klist:
        kmers_list += ["".join(p) for p in it.product(bases, repeat=k)]

    # reduce a half of k-mers due to symmetry
    kmers_dict = dict()
    for myk in kmers_list:
        k_reverse_complement = Seq(myk).reverse_complement()
        if not myk in kmers_dict and not str(k_reverse_complement) in kmers_dict:
            kmers_dict[myk] = 0

    return list(kmers_dict.keys())


def create_document(read, klist):
    """
    Create a set of document from reads, consist of all k-mer in each read
    For example:
    k = [3, 4, 5]
    documents =
    [
        'AAA AAT ... AAAT AAAC ... AAAAT AAAAC' - read 1
        'AAA AAT ... AAAT AAAC ... AAAAT AAAAC' - read 2
        ...
        'AAA AAT ... AAAT AAAC ... AAAAT AAAAC' - read n
    ]
    :param reads:
    :param klist: list of int
    :return: list of str
    """
    # create a set of document
    k_mers_read = []
    for k in klist:
        k_mers_read += [read[j : j + k] for j in range(0, len(read) - k + 1)]
    return k_mers_read


def create_dictionary(klist):
    # create k-mer dictionary
    k_mers_set = [gen_kmers(klist)]
    dictionary = corpora.Dictionary(k_mers_set)
    dictionary.save(globals.DATA_PATH + "dictionary.pkl")


class CreateDocument(MRJob):

    INPUT_PROTOCOL = JSONProtocol

    def mapper(self, _, line):
        """
        Convert the input from string to list
        Example of the input: "['AAC...TTG', '0']"
        read_label[0]: read
        read_label[1]: label
        """
        # read_label = line.strip("']['").split("', '")
        read_label = line[1]

        documents = create_document(str(read_label[0]), klist=globals.LENGTHS_OF_K_MERS)

        yield None, (line[0], read_label[0], read_label[1], documents)

    def reducer(self, key, values):
        # Store documents into a list
        # documents = []
        # for item in values:
        #     documents += [item[2]]
        # yield key, str(documents)

        for value in values:
            yield key, (value[0], value[1], value[2], value[3])


create_dictionary(globals.LENGTHS_OF_K_MERS)
CreateDocument.run()