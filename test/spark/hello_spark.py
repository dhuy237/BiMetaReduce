import pyspark

conf = pyspark.SparkConf().set("spark.driver.host", "127.0.0.1")
sc = pyspark.SparkContext(master="local", appName="myAppName", conf=conf)

txt = sc.textFile("file:////usr/share/doc/python3/copyright")
print(txt.count())

python_lines = txt.filter(lambda line: "python" in line.lower())
print(python_lines.count())