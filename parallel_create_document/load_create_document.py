from mrjob.job import MRJob
from mrjob.job import MRStep
from mrjob.protocol import RawValueProtocol
from mrjob.protocol import TextProtocol
from mrjob.protocol import RawProtocol

import itertools as it
from Bio.Seq import Seq
from multiprocessing import Pool, Array, Value
from gensim import corpora
import numpy as np

LENGTHS_OF_K_MERS = [4]
N_WORKERS = 30

def gen_kmers( klist=LENGTHS_OF_K_MERS ):
    bases = ['A', 'C', 'G', 'T']
    kmers_list = []
    for k in klist:
        kmers_list += [''.join(p) for p in it.product(bases, repeat=k)]

    # reduce a half of k-mers due to symmetry
    kmers_dict = dict()
    for myk in kmers_list:
        k_reverse_complement = Seq(myk).reverse_complement()
        if not myk in kmers_dict and not str(k_reverse_complement) in kmers_dict:
            kmers_dict[myk] = 0

    return list(kmers_dict.keys())

def create_document( reads, klist = LENGTHS_OF_K_MERS ):
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
    documents = []
    for read in reads:
        k_mers_read = []
        for k in klist:
            k_mers_read += [read[j:j + k] for j in range(0, len(read) - k + 1)]
        documents.append(k_mers_read)
    return documents

class CreateDocument(MRJob):

    def mapper_create(self, _, line):
        # create k-mer dictionary
        k_mers_set = [gen_kmers( line[0] )] #[genkmers(val) for val in klist]
        # dictionary = corpora.Dictionary(k_mers_set)

        yield None, (line[0], line[1])

    def reducer_create(self, key, values):
        pass

if __name__ == '__main__':
    CreateDocument.run()