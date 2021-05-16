from mrjob.job import MRJob
from mrjob.protocol import RawValueProtocol
from mrjob.protocol import TextProtocol
from mrjob.protocol import JSONProtocol

import itertools as it
from Bio.Seq import Seq
from multiprocessing import Pool, Array, Value
from gensim import corpora
import numpy as np
from datetime import datetime

# import sys
# sys.path.append("../")  # Add "../" to utils folder path
# from bimeta.utils import globals

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


class CreateDocument(MRJob):

    INPUT_PROTOCOL = JSONProtocol

    def configure_args(self):
        # Passthrough option
        # That file will be downloaded to each task’s local directory 
        # and the value of the option will magically be changed to its path. 
        super(CreateDocument, self).configure_args()
        self.add_passthru_arg("-k", "--k_mers", help="Lengths of k-mers", default='4', type=str, nargs='?')

    def mapper(self, _, line):
        """
        Convert the input from string to list
        Example of the input: "['AAC...TTG', '0']"
        read_label[0]: read
        read_label[1]: label
        """
        # read_label = line.strip("']['").split("', '")
        read_label = line[1]

        documents = create_document(str(read_label[0]), klist=map(int, self.options.k_mers))

        yield None, (line[0], read_label[0], read_label[1], documents)

    def reducer(self, key, values):
        # Store documents into a list
        # documents = []
        # for item in values:
        #     documents += [item[2]]
        # yield key, str(documents)

        for value in values:
            yield key, (value[0], value[1], value[2], value[3])


CreateDocument.run()