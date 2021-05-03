"""
Test countword program with cluster mode in PySpark

Run command:
spark-submit --master yarn --deploy-mode cluster \
  --conf "spark.yarn.access.hadoopFileSystems=hdfs://<folder_name>" \
  hello_spark_cluster.py hdfs://<folder_name>/<input_filename> 2
"""

import sys

from pyspark import SparkContext, SparkConf

if __name__ == "__main__":
    # create Spark context with Spark configuration
    conf = SparkConf().setAppName("SparkWordCount")
    sc = SparkContext(conf=conf)

    # get threshold
    threshold = int(sys.argv[2])

    # read in text file and split each document into words
    tokenized = sc.textFile(sys.argv[1]).flatMap(lambda line: line.split(" "))

    # count the occurrence of each word
    wordCounts = tokenized.map(lambda word: (word, 1)).reduceByKey(lambda v1,v2:v1 +v2)

    # filter out words with fewer than threshold occurrences
    filtered = wordCounts.filter(lambda pair:pair[1] >= threshold)

    # count characters
    charCounts = filtered.flatMap(lambda pair:pair[0]).map(lambda c: c).map(lambda c: (c, 1)).reduceByKey(lambda v1,v2:v1 +v2)

    result = charCounts.collect()

    print(repr(result)[1:-1])