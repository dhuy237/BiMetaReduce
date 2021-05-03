import itertools as it
from Bio.Seq import Seq
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


def create_dictionary(klist):
    # create k-mer dictionary
    k_mers_set = [gen_kmers(klist)]
    dictionary = corpora.Dictionary(k_mers_set)
    dictionary.save(globals.DATA_PATH + "dictionary.pkl")


create_dictionary(globals.LENGTHS_OF_K_MERS)