import re
import argparse
import os 
import itertools as it

parser = argparse.ArgumentParser()
parser.add_argument("-v", "--vertices", help = "Input vertices file")
parser.add_argument("-e", "--edges", help = "Input edges file")
parser.add_argument("-o", "--output", help = "Output file")
parser.add_argument("-g", "--output_graph", help = "Output graph file")
args = parser.parse_args()