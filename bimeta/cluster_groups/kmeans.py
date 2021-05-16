from gensim.models.tfidfmodel import precompute_idfs
from pyspark.ml.feature import VectorAssembler
from pyspark.ml.clustering import KMeans
from pyspark.sql import SparkSession

from sklearn import preprocessing

import gensim
from gensim import corpora
import sys
import re
import json
import numpy as np
from sklearn.metrics import confusion_matrix
import pandas as pd
import argparse
from datetime import datetime
import json

# FILENAME_GL = "/home/dhuy237/thesis/code/bimetaReduce/bimeta/data/test/output_2_2/part-00000"
# FILENAME_CORPUS = "/home/dhuy237/thesis/code/bimetaReduce/bimeta/data/test/output_1_3.txt"
# FILENAME_LABELS = "/home/dhuy237/thesis/code/bimetaReduce/bimeta/data/test/output_1_1/part-00000"
# DICTIONARY_PATH = "/home/dhuy237/thesis/code/bimetaReduce/bimeta/data/test/dictionary.pkl"

# NUM_OF_SPECIES = 2
GROUP_AGGREGATION = "MEAN"  # MEAN or MEDIAN

parser = argparse.ArgumentParser()
parser.add_argument("-g", "--group", help = "Input GL file")
parser.add_argument("-c", "--corpus", help = "Input corpus file")
parser.add_argument("-d", "--dictionary", help = "Input dictionary file")
parser.add_argument("-s", "--species", help = "Input number of species")
parser.add_argument("-l", "--labels", help = "Input labels file")
parser.add_argument("-t", "--time", help = "Output overview file")
args, unknown = parser.parse_known_args()

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

def assign_cluster_2_reads( groups, y_grp_cl ):
    label_cl_dict=dict()

    for idx, g in enumerate(groups):
        for r in g:
            label_cl_dict[r]=y_grp_cl[idx]
    
    y_cl=[]
    for i in sorted (label_cl_dict):
        y_cl.append(label_cl_dict[i])
    return y_cl

def evalQuality(y_true, y_pred, n_clusters):
    A = confusion_matrix(y_pred, y_true)
    if len(A) == 1:
      return 1, 1
    prec = sum([max(A[:,j]) for j in range(0,n_clusters)])/sum([sum(A[i,:]) for i in range(0,n_clusters)])
    rcal = sum([max(A[i,:]) for i in range(0,n_clusters)])/sum([sum(A[i,:]) for i in range(0,n_clusters)])

    return prec, rcal

def kmeans(dictionary_path, filename_corpus, filename_gl, filename_label, num_of_species):
    dictionary = load_dictionary(dictionary_path)
    corpus = read_corpus(filename_corpus)
    GL = read_group(filename_gl)

    corpus_m = gensim.matutils.corpus2dense(corpus, len(dictionary.keys())).T

    SL = []
    kmer_group_dist = compute_dist(corpus_m, GL, SL, only_seed=False)

    df = pd.DataFrame(kmer_group_dist)

    spark = SparkSession.builder.appName("kmeans").getOrCreate()
    group_dist_df = spark.createDataFrame(df)

    df_columns = group_dist_df.schema.names

    vecAssembler = VectorAssembler(inputCols=df_columns, outputCol="features")
    new_df = vecAssembler.transform(group_dist_df)

    kmeans = KMeans(k=2, seed=1)  # 2 clusters here
    model = kmeans.fit(new_df.select('features'))

    transformed = model.transform(new_df)
    transformed.select(["features", "prediction"]).show()  

    y_pred = transformed.select("prediction").rdd.flatMap(lambda x: x).collect()

    y_kmer_grp_cl = assign_cluster_2_reads(GL, y_pred)

    labels = read_labels(filename_label)

    prec, rcal = evalQuality(labels, y_kmer_grp_cl, n_clusters=num_of_species)

    return prec, rcal



start_time = datetime.now()
prec, rcal = kmeans(args.dictionary, args.corpus, args.group, args.labels, int(args.species))
execute_time = (datetime.now() - start_time).total_seconds()
print("Step 3:", execute_time)
print('K-mer (group): Prec = %.4f, Recall = %.4f, F1 = %.4f' % (prec, rcal, 2.0/(1.0/prec+1.0/rcal)))

F1 = 2 * (prec * rcal) / (prec + rcal)

data = {}
data["3"] = execute_time
data["Precision"] = prec
data["Recall"] = rcal
data["Fmeasure"] = F1

with open(args.time+'/overview.json', 'r+') as outfile:
    file = json.load(outfile)
    data["Execution"] = file["1.1"] + file["1.2"] + file["1.3"] + file["2.1"] + file["2.2"] + data["3"]
    file.update(data)
    outfile.seek(0)
    json.dump(file, outfile)