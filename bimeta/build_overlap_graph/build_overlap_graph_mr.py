from mrjob.job import MRJob
from mrjob.protocol import RawValueProtocol
from mrjob.protocol import TextProtocol
from mrjob.protocol import JSONProtocol
from mrjob.step import MRStep

import itertools as it

# import sys
# sys.path.append("../")  # Add "../" to utils folder path
# from bimeta.utils import globals

class BuildOverlapGraph(MRJob):

    INPUT_PROTOCOL = JSONProtocol

    def configure_args(self):
        # Passthrough option
        # That file will be downloaded to each task’s local directory 
        # and the value of the option will magically be changed to its path. 
        super(BuildOverlapGraph, self).configure_args()
        self.add_passthru_arg("-qm", "--q_mers", help="Lengths of q-mers", default=30, type=int)
 
    def init_configure(self): 
        self.q_mers = self.options.q_mers

    def mapper_create_hash_table(self, _, line):
        # Create hash table
        # Something wrong with this line when trying to run the code again
        # r = line[1].strip("']['").split("', '")[0]

        r = line[1][0]
        q_mers = self.q_mers

        for j in range(0, len(r) - q_mers + 1):
            lmer = r[j: j + q_mers]
            # line[0]: index of the line
            yield lmer, line[0]

    def reducer_create_hash_table(self, key, values):
        result = []
        for value in values:
            result.append(value)
        yield key, result

    def mapper_create_edge(self, key, values):
        for e in it.combinations(values, 2):
            if e[0] != e[1]:
                e_curr = (e[0], e[1])
                yield e_curr, 1

    def reducer_create_edge(self, key, values):
        yield key, sum(values)

    def steps(self):
        return [
            MRStep(
                mapper_init=self.init_configure,
                mapper=self.mapper_create_hash_table,
                reducer=self.reducer_create_hash_table
            ),
            MRStep(mapper=self.mapper_create_edge,
                   reducer=self.reducer_create_edge)
        ]

BuildOverlapGraph.run()