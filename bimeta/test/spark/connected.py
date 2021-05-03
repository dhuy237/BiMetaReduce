import pyspark
from pyspark.sql import SparkSession

sc = pyspark.SparkContext.getOrCreate()

sc.addPyFile(
    "/home/dhuy237/.ivy2/jars/graphframes_graphframes-0.8.1-spark3.0-s_2.12.jar"
)

from graphframes.examples import Graphs
from graphframes import GraphFrame

sc.setCheckpointDir("/tmp/graphframes_cps")

g = Graphs(sqlContext).friends()  # Get example graph

result = g.connectedComponents()
result.select("id", "component").orderBy("component").show()
