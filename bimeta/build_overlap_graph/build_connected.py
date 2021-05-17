import re
import argparse
import os 
import networkx as nx
import json

parser = argparse.ArgumentParser()
parser.add_argument("-v", "--vertices", help = "Input vertices file")
parser.add_argument("-e", "--edges", help = "Input edges file")
parser.add_argument("-o", "--output", help = "Output file")
parser.add_argument("-g", "--output_graph", help = "Output graph file")
args = parser.parse_args()


def build_overlap_graph(E_Filtered, labels):
    G = nx.Graph()

    # Add nodes
    for i in range(0, len(labels)):
        G.add_node(i, label=labels[i])

    # Add edges
    for kv in E_Filtered.items():
        G.add_edge(kv[0][0], kv[0][1], weight=kv[1])
    
    return G

def get_connected_components(G):
    CC = [cc for cc in nx.connected_components(G)]
    GL = []
    for subV in CC:
        GL += [list(subV)]

    return GL


def read_vertices(filename):
    """
    For reading output_1_1 file.
    Get only label attribute.

    """
    labels = []

    with open(filename) as f:
        content = f.readlines()

    for line in content:
        label = re.sub('[null\t\n\[\]\"]', '', line).replace(' ', '').split(',')[2]
        labels.append(label)
        
    return labels


def read_edges(filename):
    E_Filtered = {}

    with open(filename) as f:
        content_edges = f.readlines()

    for line in content_edges:
        # Get number character only -> remove leading and trailing whitespaces -> splits a string into a list (separator: whitespace)
        clean_line = re.sub("[^0-9]", " ", line).strip().split()
        E_Filtered[(int(clean_line[0]), int(clean_line[1]))] = int(clean_line[2])

    return E_Filtered


def save_file(GL, output_path):
    """
    Save output as (key, value) format. Only works for saving the file to local path.
    key: None
    value: each item in the GL list

    Note: Remember to export in json format for MapReduce job be able to read the file
    """
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, 'w+') as f:
        for item in GL:
            f.write("%s\n" % json.dumps(item))


labels = read_vertices(args.vertices)

E_Filtered = read_edges(args.edges)

G = build_overlap_graph(E_Filtered, labels)

GL = get_connected_components(G)

save_file(GL, args.output)