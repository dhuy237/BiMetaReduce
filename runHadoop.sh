#!/bin/bash

source $HOME/thesisEnv/bin/activate # Activate Virtual Environment

IN_FILE_NAME=$1
USR_SESSION=$2

CONF_FILE=$HOME/mrjob.conf
PERSONAL_LIB=$HOME/bimeta.zip
FOL_LCL_PATH=$HOME/Documents # Folder Local Path is used to save output of pure python program
OUT_FLR_WEB=$HOME/ServerWeb/BiMeta/userFolder/$USR_SESSION/output

RND_NO=$(echo $RANDOM % 100 + 1 | bc)
hdfs dfs -mkdir /user/session_$RND_NO

USR_HDFS=hdfs:///user
INP_HDFS=hdfs:///user/$IN_FILE_NAME
OUT_HDFS=hdfs:///user/session_$RND_NO

hdfs dfs -put $HOME/ServerWeb/BiMeta/userFolder/$USR_SESSION/input/$IN_FILE_NAME /user

# Step 1.1
python $HOME/ServerWeb/BiMeta/BimetaCode/bimeta/load_meta_reads/load_read.py $INP_HDFS \
--output $OUT_HDFS/output_1_1 \
-r hadoop \
--conf-path $CONF_FILE

hdfs dfs -get $OUT_HDFS/output_1_1/part-00000 $OUT_FLR_WEB
mv $OUT_FLR_WEB/part-00000 $OUT_FLR_WEB/OutStep_1_1

#  Step 1.2
python -m bimeta.parallel_create_document.create_dictionary --dictionary_path $FOL_LCL_PATH

python -m bimeta.parallel_create_document.load_create_document $OUT_HDFS/output_1_1/part-00000 \
--output $OUT_HDFS/output_1_2 \
-r hadoop \
--conf-path $CONF_FILE \
--py-files $PERSONAL_LIB

hdfs dfs -get $OUT_HDFS/output_1_2/part-00000 $OUT_FLR_WEB
mv $OUT_FLR_WEB/part-00000 $OUT_FLR_WEB/OutStep_1_2

#  Step 1.3 (pure python)
python -m bimeta.create_corpus.create_corpus \
--input $OUT_FLR_WEB/OutStep_1_2 \
--output $OUT_FLR_WEB/OutStep_1_3 \
--dictionary $FOL_LCL_PATH/dictionary.pkl

# Step 2.1
python -m bimeta.build_overlap_graph.build_overlap_graph $OUT_HDFS/output_1_1/part-00000 \
--output $OUT_HDFS/output_2_1 \
-r hadoop \
--conf-path $CONF_FILE \
--py-files $PERSONAL_LIB

hdfs dfs -get $OUT_HDFS/output_2_1/part-00000 $OUT_FLR_WEB
mv $OUT_FLR_WEB/part-00000 $OUT_FLR_WEB/OutStep_2_1

# Step 2.2
# spark-submit --packages graphframes:graphframes:0.8.1-spark3.0-s_2.12 \
# --py-files utils.zip \
# bimeta/build_overlap_graph/connected.py \
# --vertices $DATA_PATH_2/output_1_1_2.txt \
# --edges $DATA_PATH_2/output_2_1_2.txt \
# --checkpoint "/home/dhuy237/graphframes_cps/3" \
# --output "/home/dhuy237/graphframes_cps/3/4"

# Clean HDFS
hdfs dfs -rm /user/$IN_FILE_NAME
hdfs dfs -rm -r /user/session_$RND_NO

for i in {1..10..1}
  do 
    cp $OUT_FLR_WEB/OutStep_2_1 $OUT_FLR_WEB/OutStep_2_1_time_$i
 done

deactivate # Deactivate Virtual Environment