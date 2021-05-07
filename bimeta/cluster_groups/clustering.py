from sklearn.cluster import KMeans, SpectralClustering
from sklearn import preprocessing

import gensim
from gensim import corpora
import sys
import re
import json
import numpy as np
from sklearn.metrics import confusion_matrix
import argparse

# sys.path.append("../")  # Add "../" to utils folder path
from bimeta.utils import globals

# FILENAME_GL = globals.DATA_PATH + 'output_2_2.txt'
# FILENAME_CORPUS = globals.DATA_PATH + 'output_1_3.txt'

# DICTIONARY_PATH = globals.DATA_PATH + "dictionary.pkl"

# NUM_OF_SPECIES = 2

NUM_OF_SPECIES = 2
GROUP_AGGREGATION = "MEAN"  # MEAN or MEDIAN
SCALING = True
CLUSTERING_METHOD = "KMEANS"

parser = argparse.ArgumentParser()
parser.add_argument("-g", "--group", help = "Input GL file")
parser.add_argument("-c", "--corpus", help = "Input corpus file")
parser.add_argument("-d", "--dictionary", help = "Input dictionary file")
parser.add_argument("-s", "--species", help = "Input number of species")
parser.add_argument("-l", "--labels", help = "Input labels file")
args = parser.parse_args()

def read_group(filename_gl):
    GL = []

    with open(filename_gl) as f:
        content_vertices = f.readlines()

    for line in content_vertices:
        clean_line = re.sub("[\t\n\[\]\'']", '', line).replace(' ', '').split(',')
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


def read_labels(filename_labels):
    labels = []
    
    with open(filename_labels) as f:
        content_labels = f.readlines()
    
    for line in content_labels:
        clean_line = int(re.sub('[null\t\n\[\]\""]', '', line).replace(' ', '').split(',')[2])
        labels.append(clean_line)
    
    return labels
    

def compute_dist(dist, groups, seeds, only_seed=True):
    res = []
    if only_seed:
        for seednodes in seeds:
            tmp = dist[seednodes, :]
            if GROUP_AGGREGATION == "MEAN":
                res += [np.mean(tmp, axis=0)]
            elif GROUP_AGGREGATION == "MEDIAN":
                res += [np.median(tmp, axis=0)]
    else:
        for groupnodes in groups:
            tmp = dist[groupnodes, :]
            if GROUP_AGGREGATION == "MEAN":
                res += [np.mean(tmp, axis=0)]

            elif GROUP_AGGREGATION == "MEDIAN":
                res += [np.median(tmp, axis=0)]
                
    return np.array(res)


def cluster_groups(group_dist, num_of_species):
    if SCALING:
        scaler = preprocessing.StandardScaler()
        X_scaled = scaler.fit_transform(group_dist)
    else:
        X_scaled = group_dist

    if CLUSTERING_METHOD == 'KMEANS':
        # clustering by k-means
        kmeans = KMeans(n_clusters=num_of_species, init='k-means++')
        y_grp_cl = kmeans.fit_predict(X_scaled)
        
    elif CLUSTERING_METHOD == 'SPECTRAL':
        spectral = SpectralClustering(n_clusters=num_of_species, eigen_solver='arpack',
                                      affinity="nearest_neighbors")
        #spectral = SpectralClustering(n_clusters=NUM_OF_SPECIES, eigen_solver='arpack',
        #                              affinity="rbf")
        y_grp_cl = spectral.fit_predict(X_scaled)

    return y_grp_cl


def clustering(dictionary_path, filename_corpus, filename_gl, num_of_species):
    # Load dictionary, corpus and group list
    dictionary = load_dictionary(dictionary_path)
    corpus = read_corpus(filename_corpus)
    GL = read_group(filename_gl)

    corpus_m = gensim.matutils.corpus2dense(corpus, len(dictionary.keys())).T

    # Not implemented to get seed list yet
    SL = []

    # Training the clustering model
    kmer_group_dist = compute_dist(corpus_m, GL, SL, only_seed=False)

    model = cluster_groups(kmer_group_dist, num_of_species)
    
    y_kmer_grp_cl = assign_cluster_2_reads(GL, model)

    return y_kmer_grp_cl


def assign_cluster_2_reads( groups, y_grp_cl ):
    label_cl_dict=dict()

    for idx, g in enumerate(groups):
        for r in g:
            label_cl_dict[r]=y_grp_cl[idx]
    
    y_cl=[]
    for i in sorted (label_cl_dict):
        y_cl.append(label_cl_dict[i])

    return y_cl


def evalQuality(y_true, y_pred, n_clusters=args.species):
    A = confusion_matrix(y_pred, y_true)
    if len(A) == 1:
      return 1, 1
    prec = sum([max(A[:,j]) for j in range(0,n_clusters)])/sum([sum(A[i,:]) for i in range(0,n_clusters)])
    rcal = sum([max(A[i,:]) for i in range(0,n_clusters)])/sum([sum(A[i,:]) for i in range(0,n_clusters)])

    return prec, rcal


kmer_clustering = clustering(args.dictionary, args.corpus, args.group, int(args.species))
labels = read_labels(args.labels)
prec, rcal = evalQuality(labels, kmer_clustering, n_clusters = int(args.species))
print( 'K-mer (group): Prec = %.4f, Recall = %.4f, F1 = %.4f' % (prec, rcal, 2.0/(1.0/prec+1.0/rcal)) )