from Bio import SeqIO
from Bio.Seq import Seq
import re
import argparse
import json
import os

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help = "Input file")
parser.add_argument("-o", "--output", help = "Output file")
args = parser.parse_args()


def format_read(read):
    z = re.split("[|={,]+", read.description)
    return read.seq, z[3]


def load_meta_reads(filename, type='fasta'):
    try:
        seqs = list(SeqIO.parse(filename, type))
        reads = []
        labels = []

        # Detect for paired-end or single-end reads
        # If the id of two first reads are different (e.g.: .1 and .2), they are paired-end reads
        is_paired_end = False
        if len(seqs) > 2 and seqs[0].id[-1:] != seqs[1].id[-1:]:
            is_paired_end = True

        label_list = dict()
        label_index = 0

        for i in range(0, len(seqs), 2 if is_paired_end else 1):
            read, label = format_read(seqs[i])
            if is_paired_end:
                read2, label2 = format_read(seqs[i + 1])
                read += read2
            reads += [str(read)]

            # Create labels
            if label not in label_list:
                label_list[label] = label_index
                label_index += 1
            labels.append(label_list[label])

        del seqs

        return reads, labels
    except:
        print('Error when loading file {} '.format(filename))
        return []


def convert2json(reads, labels):
    result = []
    for i, item in enumerate(reads):
        result.append([i, [item, str(labels[i])]])
    return result

    
def save_file(result, output_path):
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, 'w+') as f:
        for item in result:
            f.write("null\t%s\n" % json.dumps(item))


reads, labels = load_meta_reads(args.input)
result = convert2json(reads, labels)
save_file(result, args.output)

