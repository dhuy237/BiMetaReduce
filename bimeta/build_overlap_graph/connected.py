import pyspark
from pyspark.sql import SparkSession
from graphframes import GraphFrame

import json
import networkx as nx
import re
from matplotlib import pyplot as plt
from igraph import Graph, plot
import pandas as pd
from pandas.core.common import flatten
import sys
import argparse
import seaborn as sns
import json

# sys.path.append("../")  # Add "../" to utils folder path
# from utils import globals

# Use only for include utils as .zip file: --py-files utils.zip
# from utils import globals

# FILENAME_VERTICES = globals.DATA_PATH + "output_1_1_2.txt"
# FILENAME_EDGES = globals.DATA_PATH + "output_2_1_2.txt"
# CHECKPOINT_DIR = "/home/dhuy237/graphframes_cps"
# OUTPUT_PATH = globals.DATA_PATH+'temp.txt'

MAXIMUM_SPECIES = 20

parser = argparse.ArgumentParser()
parser.add_argument("-v", "--vertices", help = "Input vertices file")
parser.add_argument("-e", "--edges", help = "Input edges file")
parser.add_argument("-c", "--checkpoint", help = "Checkpoint directory")
parser.add_argument("-o", "--output", help = "Output file")
parser.add_argument("-g", "--output_graph", help = "Output graph file")
parser.add_argument("-r", "--num_reads", help = "Number of shared reads", default=45, type=int)
args, unknown = parser.parse_known_args()

def build_vertices(filename_vertices):
    V = []

    with open(filename_vertices) as f:
        content_vertices = f.readlines()

    for line in content_vertices:
        clean_line = re.sub('[null\t\n\[\]\"]', "",
                            line).replace(" ", "").split(",")
        V.append(clean_line)

    df_vertices = pd.DataFrame(V, columns=["id", "read", "label"])

    return df_vertices


def build_edges(filename_edges, num_reads):
    E = []

    with open(filename_edges) as f:
        content_edges = f.readlines()

    for line in content_edges:
        # Get number character only -> remove leading and trailing whitespaces -> splits a string into a list (separator: whitespace)
        clean_line = re.sub("[^0-9]", " ", line).strip().split()
        E.append([clean_line[0], clean_line[1], clean_line[2]])

    E_Filtered = [kv for kv in E if int(kv[2]) >= num_reads]

    df_edges = pd.DataFrame(E_Filtered, columns=["src", "dst", "weight"])

    return df_edges


def get_connected_components(vertices_path, edges_path, checkpoint_dir, num_reads):
    # Read vertices and edges files
    df_vertices = build_vertices(vertices_path)
    df_edges = build_edges(edges_path, num_reads)

    # Build Graph
    spark = SparkSession.builder.appName("build_graph").getOrCreate()
    vertices = spark.createDataFrame(df_vertices)

    edges = spark.createDataFrame(df_edges)
    g = GraphFrame(vertices, edges)

    # Display Graph
    g.vertices.show()
    g.edges.show()

    # Connected Components
    # Get SparkContext using spark.sparkContext
    spark.sparkContext.setCheckpointDir(dirName=checkpoint_dir)
    result = g.connectedComponents()

    dictionary = {}

    sorted_result = result.select("id", "component").orderBy('component', ascending=False)

    for row in sorted_result.collect():
        if row[1] in dictionary:
            dictionary[row[1]].append(row[0])
        else:
            dictionary[row[1]] = [row[0]]

    GL = []

    for _, value in dictionary.items():
        GL.append(value)

    return GL, spark, g


def save_file_local(GL, path):
    """
    Save output as (key, value) format. Only works for saving the file to local path.
    key: null
    value: each item in the GL list

    Note: Remember to export in json format for MapReduce job be able to read the file
    """
    with open(path, 'w+') as f:
        for item in GL:
            f.write("null\t%s\n" % json.dumps(item))


def save_file_hdfs(GL, session, path):
    """
    Use this to save file to HDFS.
    The saved file will be named "part-00000"
    """
    # First need to convert the list to parallel RDD
    rdd_list = session.sparkContext.parallelize(GL)

    # Use the map function to write one element per line and write all elements to a single file (coalesce)
    rdd_list.coalesce(1).map(lambda row: str(row)).saveAsTextFile(path)


def visualize_graph(color_dict, graphframes, output_path):
    # Create iGraph from graphframes edges
    ig = Graph.TupleList(graphframes.edges.collect(), directed=False)

    # Get label from graphframe vertices
    label = graphframes.vertices.select("label").rdd.flatMap(lambda x: x).collect()

    # Create "color" attribute for vertex label
    ig.vs["color"] = [color_dict[label[int(name)]] for name in ig.vs["name"]]

    # Save plot as .png
    out = plot(ig, vertex_size=10, bbox=(0, 0, 500, 500))
    out.save(output_path+"/graph.png")

palette = sns.color_palette(None, MAXIMUM_SPECIES).as_hex()
color_dict = {str(i):palette[i] for i in range(MAXIMUM_SPECIES)}

GL, spark, g = get_connected_components(args.vertices, args.edges, args.checkpoint, args.num_reads)

# save_file_local(GL, args.output)
save_file_hdfs(GL, spark, args.output)
visualize_graph(color_dict, g, args.output_graph)
