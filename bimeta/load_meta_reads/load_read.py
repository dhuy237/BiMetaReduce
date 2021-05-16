#!./venv/bin/python
from mrjob.job import MRJob
from mrjob.protocol import RawValueProtocol
from mrjob.protocol import TextProtocol
from mrjob.protocol import RawProtocol

from Bio import SeqIO
from Bio.Seq import Seq
import re
from datetime import datetime
import json
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument("-t", "--time", help = "Output overview file")
# Use parse_known_args() to ignore all the arguments for mrjob class
args, unknown = parser.parse_known_args()

def format_read(read):
    z = re.split("[|={,]+", read.description)
    return read.seq, z[3]


class LoadMetaRead(MRJob):

    # INPUT_PROTOCOL = RawValueProtocol # RawValueProtocol: Default in python3
    # INTERNAL_PROTOCOL = RawValueProtocol
    # OUTPUT_PROTOCOL = RawProtocol
    def configure_args(self):
        super(LoadMetaRead, self).configure_args()
        # Define the argument again to ignore the argparse when running the MR job
        self.add_passthru_arg("-t", "--time", help = "Output overview file")

    def mapper_raw(self, file_path, file_uri):
        from Bio import SeqIO
        from Bio.Seq import Seq

        seqs = list(SeqIO.parse(file_path, "fasta"))

        is_paired_end = False
        if len(seqs) > 2 and seqs[0].id[-1:] != seqs[1].id[-1:]:
            is_paired_end = True

        label_list = dict()
        label_index = 0

        for i in range(0, len(seqs), 2 if is_paired_end else 1):
            read, label = format_read(seqs[i])
            if is_paired_end:
                read2, _ = format_read(seqs[i + 1])
                read += read2

            if label not in label_list:
                label_list[label] = label_index
                label_index += 1

            yield None, (str(read), str(label_list[label]))

    def reducer(self, key, values):
        for i, value in enumerate(values):
            # yield i to know their line position in the dataset
            yield key, (i, value)

    # combiner = reducer

start_time = datetime.now()
LoadMetaRead.run()
execute_time = (datetime.now() - start_time).total_seconds()
# print("Step 1.1:", execute_time)
# sys.stderr.write(str(execute_time))
# data = {"step1":str(execute_time)}
# json.dump(data, sys.stderr)
# data = {}
# data["1.1"] = execute_time
# with open(args.time+'/overview.json', 'w') as outfile:
#     json.dump(data, outfile)