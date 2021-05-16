#!/bin/bash

source $HOME/thesisEnv/bin/activate # Activate Virtual Environment

IN_FILE_NAME=$1
USR_SESSION=$2

CONF_FILE=$HOME/mrjob.conf
PERSONAL_LIB=$HOME/bimeta.zip
FOL_LCL_PATH=$HOME/Documents # Folder Local Path is used to save output of pure python program
OUT_FLR_WEB=$HOME/ServerWeb/BiMeta/userFolder/$USR_SESSION/output
GRAPH_FLR_WEB=$HOME/ServerWeb/BiMeta/userFolder/$USR_SESSION/graph
JSON_FLR_WEB=$HOME/ServerWeb/BiMeta/jsonData
OVERVIEW=overview.json

RND_NO=$(echo $RANDOM % 100 + 1 | bc)
hdfs dfs -mkdir /user/session_$RND_NO

USR_HDFS=hdfs:///user
INP_HDFS=hdfs:///user/$IN_FILE_NAME
OUT_HDFS=hdfs:///user/session_$RND_NO

hdfs dfs -put $HOME/ServerWeb/BiMeta/userFolder/$USR_SESSION/input/$IN_FILE_NAME /user

# Python Code Variables START
LENGTHS_OF_K_MERS=4
LENGTH_OF_Q_MERS=30
NUM_SHARED_READS=45
NUM_OF_SPECIES=2
# Python Code Variables END

# Step 1.1
# Start 1.1----------------------------------------------------------------------
START_TIME=`date +%s%N`

python $HOME/ServerWeb/BiMeta/BimetaCode/bimeta/load_meta_reads/load_read.py $INP_HDFS \
--output $OUT_HDFS/output_1_1 \
-r hadoop \
--conf-path $CONF_FILE

END_TIME=`date +%s%N`

RUN_TIME=`expr $END_TIME - $START_TIME`
RUN_TIME_IN_S=$(echo "scale = 3; $RUN_TIME / 1000000000" | bc)
echo "{\"Step_1_1\":\"$RUN_TIME_IN_S\"," > $JSON_FLR_WEB/$OVERVIEW
# End----------------------------------------------------------------------

hdfs dfs -get $OUT_HDFS/output_1_1/part-00000 $OUT_FLR_WEB
mv $OUT_FLR_WEB/part-00000 $OUT_FLR_WEB/OutStep_1_1

#  Step 1.2 
python bimeta/parallel_create_document/create_dictionary.py \
--dictionary_path $FOL_LCL_PATH \
--k_mers $LENGTHS_OF_K_MERS

# Start 1.2----------------------------------------------------------------------
START_TIME=`date +%s%N`

python bimeta/parallel_create_document/load_create_document.py \
$OUT_HDFS/output_1_1/part-00000 \
--output $OUT_HDFS/output_1_2 \
-r hadoop \
--conf-path $CONF_FILE \
--k_mers $LENGTHS_OF_K_MERS

END_TIME=`date +%s%N`

RUN_TIME=`expr $END_TIME - $START_TIME`
RUN_TIME_IN_S=$(echo "scale = 3; $RUN_TIME / 1000000000" | bc)
echo "{\"Step_1_2\":\"$RUN_TIME_IN_S\"," > $JSON_FLR_WEB/$OVERVIEW
# End----------------------------------------------------------------------


hdfs dfs -get $OUT_HDFS/output_1_2/part-00000 $OUT_FLR_WEB
mv $OUT_FLR_WEB/part-00000 $OUT_FLR_WEB/OutStep_1_2

# Step 1.3 (pure python)
# Start 1.3----------------------------------------------------------------------
START_TIME=`date +%s%N`

python bimeta/create_corpus/create_corpus.py \
--input $OUT_FLR_WEB/OutStep_1_2 \
--output $OUT_FLR_WEB/OutStep_1_3 \
--dictionary $FOL_LCL_PATH/dictionary.pkl

END_TIME=`date +%s%N`

RUN_TIME=`expr $END_TIME - $START_TIME`
RUN_TIME_IN_S=$(echo "scale = 3; $RUN_TIME / 1000000000" | bc)
echo "{\"Step_1_3\":\"$RUN_TIME_IN_S\"," > $JSON_FLR_WEB/$OVERVIEW
# End----------------------------------------------------------------------


# Step 2.1
# Start 2.1----------------------------------------------------------------------
START_TIME=`date +%s%N`

python bimeta/build_overlap_graph/build_overlap_graph.py \
$OUT_HDFS/output_1_1/part-00000 \
--output $OUT_HDFS/output_2_1 \
-r hadoop \
--conf-path $CONF_FILE \
--q_mers $LENGTH_OF_Q_MERS

END_TIME=`date +%s%N`

RUN_TIME=`expr $END_TIME - $START_TIME`
RUN_TIME_IN_S=$(echo "scale = 3; $RUN_TIME / 1000000000" | bc)
echo "{\"Step_2_1\":\"$RUN_TIME_IN_S\"," > $JSON_FLR_WEB/$OVERVIEW
# End----------------------------------------------------------------------

hdfs dfs -get $OUT_HDFS/output_2_1/part-00000 $OUT_FLR_WEB
mv $OUT_FLR_WEB/part-00000 $OUT_FLR_WEB/OutStep_2_1

# Step 2.2
# spark-submit --packages graphframes:graphframes:0.8.1-spark3.0-s_2.12 \
# bimeta/build_overlap_graph/connected.py \
# --vertices $OUT_FLR_WEB/OutStep_1_1 \
# --edges $OUT_FLR_WEB/OutStep_2_1 \
# --checkpoint $OUT_HDFS/graphframes_cps \
# --output $OUT_HDFS/graphframes_cps/2 \
# --num_reads $NUM_SHARED_READS
# --conf "spark.yarn.access.hadoopFileSystems=${OUT_HDFS}" \
# --master yarn \
# --deploy-mode cluster \
# --files /home/tnhan/ServerWeb/BiMeta/userFolder/thesis1/output/OutStep_1_1,/home/tnhan/ServerWeb/BiMeta/userFolder/thesis1/output/OutStep_2_1 \

# Start 2.2----------------------------------------------------------------------
START_TIME=`date +%s%N`

spark-submit --packages graphframes:graphframes:0.8.1-spark3.0-s_2.12 \
bimeta/build_overlap_graph/connected.py \
--vertices $OUT_FLR_WEB/OutStep_1_1 \
--edges $OUT_FLR_WEB/OutStep_2_1 \
--checkpoint $OUT_HDFS/graphframes_cps \
--output $OUT_HDFS/graphframes_cps/2 \
--output_graph $GRAPH_FLR_WEB \
--num_reads $NUM_SHARED_READS

END_TIME=`date +%s%N`

RUN_TIME=`expr $END_TIME - $START_TIME`
RUN_TIME_IN_S=$(echo "scale = 3; $RUN_TIME / 1000000000" | bc)
echo "{\"Step_2_2\":\"$RUN_TIME_IN_S\"," > $JSON_FLR_WEB/$OVERVIEW
# End----------------------------------------------------------------------


hdfs dfs -get $OUT_HDFS/graphframes_cps/2/part-00000 $OUT_FLR_WEB
mv $OUT_FLR_WEB/part-00000 $OUT_FLR_WEB/OutStep_2_2

# Step 3
# Start 3----------------------------------------------------------------------
spark-submit bimeta/cluster_groups/kmeans.py \
--group $OUT_FLR_WEB/OutStep_2_2 \
--corpus $OUT_FLR_WEB/OutStep_1_3 \
--dictionary $FOL_LCL_PATH/dictionary.pkl \
--species $NUM_OF_SPECIES \
--labels $OUT_FLR_WEB/OutStep_1_1
--time $JSON_FLR_WEB/$OVERVIEW
# End----------------------------------------------------------------------

# # Clean HDFS
# # hdfs dfs -rm /user/$IN_FILE_NAME
# # hdfs dfs -rm -r /user/session_$RND_NO

# # for i in {1..10..1}
# # do 
# #   cp $OUT_FLR_WEB/OutStep_2_1 $OUT_FLR_WEB/OutStep_2_1_time_$i
# # done

deactivate # Deactivate Virtual Environment