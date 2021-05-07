#!/bin/bash

# Activate Virtual Environment
source $HOME/thesisEnv/bin/activate

#export DATA_PATH_2=bimeta/data/R4_medium
INPUT_FILE=$1
DIC_PATH=$HOME/Documents
CONF_FILE=$HOME/mrjob.conf
UTILS_ZIP=$HOME/bimeta.zip
ID_OUT=$(echo $RANDOM % 100 + 1 | bc)
hdfs dfs -mkdir /user/session_$ID_OUT
DATA_PATH=hdfs:///user
DATA_PATH_2=hdfs:///user/session_$ID_OUT

# Step 1.1
python bimeta/load_meta_reads/load_read.py $DATA_PATH/$INPUT_FILE --output $DATA_PATH_2/output_1_1 -r hadoop --conf-path $CONF_FILE

#  Step 1.2
python -m bimeta.parallel_create_document.create_dictionary --dictionary_path $DIC_PATH 
python -m bimeta.parallel_create_document.load_create_document $DATA_PATH_2/output_1_1/part-00000 \
--output $DATA_PATH_2/output_1_2 -r hadoop --conf-path $CONF_FILE --py-files $UTILS_ZIP

hdfs dfs -get $DATA_PATH_2/output_1_2/part-00000 $DIC_PATH

#  Step 1.3
python -m bimeta.create_corpus.create_corpus --input $DIC_PATH/part-00000 --output $DIC_PATH/output_1_3.txt \
--dictionary $DIC_PATH/dictionary.pkl

# Step 2.1
python -m bimeta.build_overlap_graph.build_overlap_graph $DATA_PATH_2/output_1_1/part-00000 \
--output $DATA_PATH_2/output_2_1 -r hadoop \
--conf-path $CONF_FILE \
--py-files $UTILS_ZIP

#Step 2.2
# spark-submit --packages graphframes:graphframes:0.8.1-spark3.0-s_2.12 \
# --py-files utils.zip \
# bimeta/build_overlap_graph/connected.py \
# --vertices $DATA_PATH_2/output_1_1_2.txt \
# --edges $DATA_PATH_2/output_2_1_2.txt \
# --checkpoint "/home/dhuy237/graphframes_cps/3" \
# --output "/home/dhuy237/graphframes_cps/3/4"

# Deactivate Virtual Environment
deactivate