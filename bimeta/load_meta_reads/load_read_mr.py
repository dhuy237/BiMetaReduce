#!./venv/bin/python
from mrjob.job import MRJob
from mrjob.protocol import RawValueProtocol
from mrjob.protocol import TextProtocol
from mrjob.protocol import RawProtocol

from Bio import SeqIO
from Bio.Seq import Seq
import re

def format_read(read):
    z = re.split("[|={,]+", read.description)
    return read.seq, z[3]


class LoadMetaRead(MRJob):

    # INPUT_PROTOCOL = RawValueProtocol # RawValueProtocol: Default in python3
    # INTERNAL_PROTOCOL = RawValueProtocol
    # OUTPUT_PROTOCOL = RawProtocol

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

LoadMetaRead.run()
