import pyspark
from pyspark.sql import SparkSession
from graphframes import GraphFrame

import argparse
import seaborn as sns
import json
from igraph import Graph, plot
import re
import pandas as pd

# For generating the color palette
MAXIMUM_SPECIES = 20

parser = argparse.ArgumentParser()
parser.add_argument("-v", "--vertices", help = "Input vertices file")
parser.add_argument("-e", "--edges", help = "Input edges file")
parser.add_argument("-g", "--output_graph", help = "Output graph file")
parser.add_argument("-r", "--num_reads", help = "Number of shared reads", default=45, type=int)
args = parser.parse_args()

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


def visualize_graph(vertices_path, edges_path, num_reads, color_dict, output_path):
    # Read vertices and edges files
    df_vertices = build_vertices(vertices_path)
    df_edges = build_edges(edges_path, num_reads)

    # Build Graph
    spark = SparkSession.builder.appName("build_graph").getOrCreate()
    vertices = spark.createDataFrame(df_vertices)

    edges = spark.createDataFrame(df_edges)
    g = GraphFrame(vertices, edges)

    # Create iGraph from graphframes edges
    ig = Graph.TupleList(g.edges.collect(), directed=False)

    # Get label from graphframe vertices
    label = g.vertices.select("label").rdd.flatMap(lambda x: x).collect()

    # Create "color" attribute for vertex label
    ig.vs["color"] = [color_dict[label[int(name)]] for name in ig.vs["name"]]

    # Save plot as .png
    plot(ig, target=output_path, vertex_size=10, bbox=(0, 0, 500, 500))


palette = sns.color_palette(None, MAXIMUM_SPECIES).as_hex()
color_dict = {str(i):palette[i] for i in range(MAXIMUM_SPECIES)}
visualize_graph(args.vertices, args.edges, args.num_reads, color_dict, args.output_graph)