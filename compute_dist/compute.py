from mrjob.job import MRJob
from mrjob.protocol import RawValueProtocol
from mrjob.protocol import TextProtocol
from mrjob.protocol import JSONProtocol

import gensim
import numpy as np
import sys

sys.path.append("../")  # Add "../" to utils folder path
from utils import globals


def compute_dist(dist, groups, seeds, only_seed=True):
    res = []
    if only_seed:
        print(seeds)
        for seednodes in seeds:
            tmp = dist[seednodes, :]
            print(tmp)
            if globals.GROUP_AGGREGATION == "MEAN":
                res += [np.mean(tmp, axis=0)]
            elif globals.GROUP_AGGREGATION == "MEDIAN":
                res += [np.median(tmp, axis=0)]
    else:
        for groupnodes in groups:
            tmp = dist[groupnodes, :]
            if globals.GROUP_AGGREGATION == "MEAN":
                res += [np.mean(tmp, axis=0)]

            elif globals.GROUP_AGGREGATION == "MEDIAN":
                res += [np.median(tmp, axis=0)]
    print("------------------")
    print(np.array(res))
    print("------------------")
    return np.array(res)


# Run command for compute_dist()
# corpus_m = gensim.matutils.corpus2dense(corpus, len(dictionary.keys())).T
# kmer_group_dist = compute_dist(corpus_m, GL, SL, only_seed=False)


class ComputeDist(MRJob):

    INPUT_PROTOCOL = JSONProtocol

    def mapper(self, _, line):
        yield None, line

    def reducer(self, key, values):
        for value in values:
            yield key, value


ComputeDist.run()