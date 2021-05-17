import itertools as it
from Bio.Seq import Seq
from multiprocessing import Pool, Array, Value
from gensim import corpora
import numpy as np
import argparse
import json
import os
import re

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help = "Input file")
parser.add_argument("-o", "--output", help = "Output file")
parser.add_argument("-k", "--k_mers", help="Lengths of k-mers", default='4', type=str, nargs='?')
args = parser.parse_args()

# LENGTHS_OF_K_MERS = [4]
N_WORKERS = 30


def create_document(reads, klist):
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
    
    
def parallel_create_document(reads, klist, n_workers=2):
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

    documents = []
    reads_str_chunk = [list(item) for item in np.array_split(reads, n_workers)]
    chunks = [(reads_str_chunk[i], klist) for i in range(n_workers)]
    pool = Pool(processes=n_workers)

    result = pool.starmap(create_document, chunks)
    for item in result:
        documents += item

    return documents


def read_file(filename):
    """
    For reading output_1_1 file

    """
    reads = []
    labels = []

    with open(filename) as f:
        content = f.readlines()

    for line in content:
        _, read, label = re.sub('[null\t\n\[\]\"]', '', line).replace(' ', '').split(',')
        reads.append(read)
        labels.append(label)
        
    return reads, labels


def convert2json(reads, labels, documents):
    """
    For saving to output_1_2 file

    """
    result = []
    for i, item in enumerate(reads):
        result.append([i, item, str(labels[i]), documents[i]])
    return result


def save_file(result, output_path):
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, 'w+') as f:
        for item in result:
            f.write("null\t%s\n" % json.dumps(item))


reads, labels = read_file(args.input)

documents = create_document(reads, klist=list(map(int, args.k_mers)))

result = convert2json(reads, labels, documents)
save_file(result, args.output)