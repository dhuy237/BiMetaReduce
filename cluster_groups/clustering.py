import pyspark
from pyspark.sql import SparkSession
from graphframes import GraphFrame

from sklearn.cluster import KMeans, SpectralClustering
from sklearn import preprocessing

import gensim
from gensim import corpora
import sys
import re
import json
import numpy as np

sys.path.append("../")  # Add "../" to utils folder path
from utils import globals

FILENAME_GL = globals.DATA_PATH + 'output_2_2.txt'
FILENAME_CORPUS = globals.DATA_PATH + 'output_1_3.txt'

DICTIONARY_PATH = globals.DATA_PATH + "dictionary.pkl"

NUM_OF_SPECIES = 2

def read_group(filename_gl):
    GL = []

    with open(filename_gl) as f:
        content_vertices = f.readlines()

    for line in content_vertices:
        clean_line = re.sub('[null\t\n\[\]\""]', '', line).replace(' ', '').split(',')
        GL.append(list(map(int, clean_line))) # Convert all strings in a list to int

    return GL


def load_dictionary(dictionary_path):
    dictionary = corpora.Dictionary.load(dictionary_path)
    return dictionary


def read_corpus(filename_corpus):
    corpus = []

    with open(filename_corpus) as f:
        content_corpus = f.readlines()

    for line in content_corpus:
        clean_line = json.loads(line.replace('null\t', '{"a":').replace("\n", "}"))["a"][1]
        corpus.append(clean_line)
    
    return corpus


def compute_dist(dist, groups, seeds, only_seed=True):
    res = []
    if only_seed:
        for seednodes in seeds:
            tmp = dist[seednodes, :]
            if globals.GROUP_AGGREGATION == "MEAN":
                res += [np.mean(tmp, axis=0)]
            elif globals.GROUP_AGGREGATION == "MEDIAN":
                res += [np.median(tmp, axis=0)]
    else:
        for groupnodes in groups:
            tmp = dist[groupnodes, :]
            if globals.GROUP_AGGREGATION == "MEAN":
                res += [np.mean(tmp, axis=0)]

            elif globals.GROUP_AGGREGATION == "MEDIAN":
                res += [np.median(tmp, axis=0)]
                
    return np.array(res)


def cluster_groups(group_dist):
    if globals.SCALING:
        scaler = preprocessing.StandardScaler()
        X_scaled = scaler.fit_transform(group_dist)
    else:
        X_scaled = group_dist

    if globals.CLUSTERING_METHOD == 'KMEANS':
        # clustering by k-means
        kmeans = KMeans(n_clusters=NUM_OF_SPECIES, init='k-means++')
        y_grp_cl = kmeans.fit_predict(X_scaled)
        
    elif globals.CLUSTERING_METHOD == 'SPECTRAL':
        spectral = SpectralClustering(n_clusters=NUM_OF_SPECIES, eigen_solver='arpack',
                                      affinity="nearest_neighbors")
        #spectral = SpectralClustering(n_clusters=NUM_OF_SPECIES, eigen_solver='arpack',
        #                              affinity="rbf")
        y_grp_cl = spectral.fit_predict(X_scaled)

    return y_grp_cl

def clustering(dictionary_path, filename_corpus, filename_gl):
    # Load dictionary, corpus and group list
    dictionary = load_dictionary(dictionary_path)
    corpus = read_corpus(filename_corpus)
    GL = read_group(filename_gl)

    corpus_m = gensim.matutils.corpus2dense(corpus, len(dictionary.keys())).T

    # Not implemented to get seed list yet
    SL = []

    # Training the clustering model
    kmer_group_dist = compute_dist(corpus_m, GL, SL, only_seed=False)

    return kmer_group_dist

kmer_clustering = clustering(DICTIONARY_PATH, FILENAME_CORPUS, FILENAME_GL)