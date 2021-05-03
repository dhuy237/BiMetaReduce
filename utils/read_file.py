from Bio import SeqIO
from Bio.Seq import Seq
import re
import json

"""
For getting the overview of the .fna file.
Example of output file reads_summary.json
[
    {
        "species": 1,
        "number": 500,
        "code": "15668172",
        "name": "\"Methanocaldococcus jannaschii DSM 2661 chromosome\" (392b1054a4bf536ea1cc349545ace50120973c3a)"
    },
    ...
]
"""

DATA_PATH = "/home/dhuy237/thesis/code/bimetaReduce/data/R4_medium/R4_medium.fna"
SAVE_PATH = "/home/dhuy237/thesis/code/bimetaReduce/data/R4_medium/"

def format_read(read):
    z = re.split("[|={,]+", read.description)
    return read.seq, z[3], z[-1]

def load_meta_reads(filename, type='fasta'):
    seqs = list(SeqIO.parse(filename, type))
    labels = []

    # Detect for paired-end or single-end reads
    # If the id of two first reads are different (e.g.: .1 and .2), they are paired-end reads
    is_paired_end = False
    if len(seqs) > 2 and seqs[0].id[-1:] != seqs[1].id[-1:]:
        is_paired_end = True

    label_list = dict()
    result = []
    label_index = 0

    for i in range(0, len(seqs), 2 if is_paired_end else 1):
        _, label, name = format_read(seqs[i])

        # Create labels
        if label not in label_list:
            inf = {}
            inf["species"] = label_index
            inf["number"] = 0
            inf["code"] = label
            inf["name"] = name
            result.append(inf)

            label_list[label] = label_index
            label_index += 1
            
        labels.append(label_list[label])
    
    # Count number of occurrences of each species
    # Add the species number with 1 because it start from 0
    for i, item in enumerate(result):
        result[i]["number"] = labels.count(int(item["species"]))
        result[i]["species"] += 1
        
    del seqs

    return result

def convert2json(data, save_path):
    with open(save_path+'reads_summary.json', 'w+', encoding='utf-8') as f:
        json.dump(data, f, ensure_ascii=False, indent=4)

data = load_meta_reads(DATA_PATH)
convert2json(data, SAVE_PATH)