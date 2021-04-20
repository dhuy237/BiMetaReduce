from mrjob.job import MRJob
from mrjob.protocol import RawValueProtocol
from mrjob.protocol import TextProtocol
from mrjob.protocol import JSONProtocol
from mrjob.step import MRStep

import itertools as it
import networkx as nx
import nxmetis

import sys

sys.path.append("../")  # Add "../" to utils folder path
from utils import globals


def build_hash_table(reads):
    # Create hash table with q-mers are keys
    # Building hash table
    lmers_dict = dict()
    for idx, r in enumerate(reads):
        for j in range(0, len(r) - globals.LENGTH_OF_Q_MERS + 1):
            lmer = r[j : j + globals.LENGTH_OF_Q_MERS]
            if lmer in lmers_dict:
                lmers_dict[lmer] += [idx]
            else:
                lmers_dict[lmer] = [idx]

    return lmers_dict


def build_edge(lmers_dict):
    # Building edges
    E = dict()
    for lmer in lmers_dict:
        for e in it.combinations(lmers_dict[lmer], 2):
            if e[0] != e[1]:
                e_curr = (e[0], e[1])
            if e_curr in E:
                E[e_curr] += 1  # Number of connected lines between read a and b
            else:
                E[e_curr] = 1

    # Contain pairs of reads that is connected (edge weight >= NUM_SHARED_READS is connected)
    E_Filtered = {kv[0]: kv[1] for kv in E.items() if kv[1] >= globals.NUM_SHARED_READS}

    return E_Filtered


def build_graph(E_Filtered, labels):
    # Building graph
    G = nx.Graph()
    # Adding nodes
    color_map = {
        0: "red",
        1: "green",
        2: "blue",
        3: "yellow",
        4: "darkcyan",
        5: "violet",
        6: "black",
        7: "grey",
        8: "sienna",
        9: "wheat",
        10: "olive",
        11: "lightgreen",
        12: "cyan",
        13: "slategray",
        14: "navy",
        15: "hotpink",
    }
    for i in range(0, len(labels)):
        G.add_node(i, label=labels[i], color=color_map[labels[i]])

    # Adding edges
    for kv in E_Filtered.items():
        G.add_edge(kv[0][0], kv[0][1], weight=kv[1])

    return G


def metis_partition_groups_seeds(G):
    CC = [cc for cc in nx.connected_components(G)]
    GL = []
    for subV in CC:
        if len(subV) > globals.MAXIMUM_COMPONENT_SIZE:
            # use metis to split the graph
            subG = nx.subgraph(G, subV)
            nparts = int(len(subV) / globals.MAXIMUM_COMPONENT_SIZE + 1)
            (_, parts) = nxmetis.partition(subG, nparts, edge_weight="weight")

            # only add connected components
            for p in parts:
                pG = nx.subgraph(G, p)
                GL += [list(cc) for cc in nx.connected_components(pG)]

            # add to group list
            # Maybe ignore this line
            GL += parts
        else:
            GL += [list(subV)]
    SL = []
    for p in GL:
        pG = nx.subgraph(G, p)
        SL += [nx.maximal_independent_set(pG)]

    return GL, SL


class BuildOverlapGraph(MRJob):

    INPUT_PROTOCOL = JSONProtocol

    def mapper_create_hash_table(self, _, line):
        # Create hash table
        # Something wrong with this line when trying to run the code again
        # r = line[1].strip("']['").split("', '")[0]

        r = line[1][0]
        for j in range(0, len(r) - globals.LENGTH_OF_Q_MERS + 1):
            lmer = r[j : j + globals.LENGTH_OF_Q_MERS]
            # line[0]: index of the line
            yield lmer, line[0]

    def reducer_create_hash_table(self, key, values):
        result = []
        for value in values:
            result.append(value)
        yield key, result

    def mapper_create_edge(self, key, value):
        for e in it.combinations(value, 2):
            if e[0] != e[1]:
                e_curr = (e[0], e[1])
                yield e_curr, 1

    def reducer_create_edge(self, key, values):
        yield key, sum(values)

    def steps(self):
        return [
            MRStep(
                mapper=self.mapper_create_hash_table,
                reducer=self.reducer_create_hash_table,
            ),
            MRStep(mapper=self.mapper_create_edge, reducer=self.reducer_create_edge),
        ]


BuildOverlapGraph.run()