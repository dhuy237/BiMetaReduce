from Bio import SeqIO
from Bio.Seq import Seq
import re
import json
import seaborn as sns

"""
For getting the overview of the .fna file.
Example of output file reads_summary.json
[
    {
        "species": 1,
        "number": 500,
        "code": "15668172",
        "name": "\"Methanocaldococcus jannaschii DSM 2661 chromosome\" (392b1054a4bf536ea1cc349545ace50120973c3a)",
        "color: "#F1F1F1"
    },
    ...
]
"""

DATA_PATH = "/home/dhuy237/thesis/code/bimetaReduce/bimeta/data/R4_medium/R4_medium.fna"
SAVE_PATH = "/home/dhuy237/thesis/code/bimetaReduce/bimeta/data/R4_medium/"
JSON_PATH = "/home/dhuy237/thesis/code/bimetaReduce/bimeta/data/test"
MAXIMUM_SPECIES = 20

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

    # Generate color for each species
    palette = sns.color_palette(None, MAXIMUM_SPECIES).as_hex()
    
    for i in range(0, len(seqs), 2 if is_paired_end else 1):
        _, label, name = format_read(seqs[i])

        # Create labels
        if label not in label_list:
            inf = {}
            inf["species"] = label_index
            inf["number"] = 0
            inf["code"] = label
            inf["name"] = name
            inf["color"] = palette[label_index]
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

def readjson(data_path):
    with open(data_path+'/overview.json', 'r') as outfile:
        data = json.load(outfile)   
        return data
       

data = load_meta_reads(DATA_PATH)
convert2json(data, SAVE_PATH)

data_json = readjson(JSON_PATH)
print(data_json)