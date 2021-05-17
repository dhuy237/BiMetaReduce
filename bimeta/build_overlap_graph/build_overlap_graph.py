import re
import argparse
import os 
import itertools as it

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help = "Input file")
parser.add_argument("-o", "--output", help = "Output file")
parser.add_argument("-qm", "--q_mers", help="Lengths of q-mers", default=30, type=int)
parser.add_argument("-r", "--num_reads", help = "Number of shared reads", default=45, type=int)
args = parser.parse_args()


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


def build_hash_table(reads, q_mers):
    lmers_dict = {}

    for idx, r in enumerate(reads):
        for j in range(0, len(r) - q_mers + 1):
            lmer = r[j:j+q_mers]
            if lmer in lmers_dict:
                lmers_dict[lmer] += [idx]
            else:
                lmers_dict[lmer] = [idx]

    return lmers_dict


def build_edges(lmers_dict, num_reads):
    E = {}

    for lmer in lmers_dict:

        for e in it.combinations(lmers_dict[lmer], 2):
            if e[0] != e[1]:
                e_curr = (e[0], e[1])
            if e_curr in E:
                E[e_curr] += 1 # Number of connected lines between read a and b
            else:
                E[e_curr] = 1
    
    E_Filtered = {kv[0]: kv[1] for kv in E.items() if kv[1] >= num_reads}
    
    return E_Filtered


def save_file(E_Filtered, output_path):
    """
    For saving to output_2_1 file

    """
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, 'w+') as f:
        for k, v in E_Filtered.items():
            f.write("%s\t%s\n" % (list(k), v))


reads, labels = read_file(args.input)

lmers_dict = build_hash_table(reads, args.q_mers)

E_Filtered = build_edges(lmers_dict, args.num_reads)

save_file(E_Filtered, args.output)