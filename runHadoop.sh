#!/bin/bash

START_TIME=`date +%s`
echo "\n\n\n Start time is : \n\n"
echo $START_TIME
echo "\n\n\n"

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

# Python Code Variables START
LENGTHS_OF_K_MERS=4
LENGTH_OF_Q_MERS=30
NUM_SHARED_READS=45
NUM_OF_SPECIES=2
# Python Code Variables END

# Step 1.1
python $HOME/ServerWeb/BiMeta/BimetaCode/bimeta/load_meta_reads/load_read.py $INP_HDFS \
--output $OUT_HDFS/output_1_1 \
-r hadoop \
--conf-path $CONF_FILE

END_TIME=`date +%s`
echo "\n\n\n End time is : \n\n"
echo $END_TIME
echo "\n\n\n"

RUN_TIME=`expr $END_TIME - $START_TIME`
echo "\n\n\n Run time is : \n\n"
echo $RUN_TIME
echo "\n\n\n"

hdfs dfs -get $OUT_HDFS/output_1_1/part-00000 $OUT_FLR_WEB
mv $OUT_FLR_WEB/part-00000 $OUT_FLR_WEB/OutStep_1_1

#  Step 1.2
python bimeta/parallel_create_document/create_dictionary.py \
--dictionary_path $FOL_LCL_PATH \
--k_mers $LENGTHS_OF_K_MERS

python bimeta/parallel_create_document/load_create_document.py \
$OUT_HDFS/output_1_1/part-00000 \
--output $OUT_HDFS/output_1_2 \
-r hadoop \
--conf-path $CONF_FILE \
--k_mers $LENGTHS_OF_K_MERS

hdfs dfs -get $OUT_HDFS/output_1_2/part-00000 $OUT_FLR_WEB
mv $OUT_FLR_WEB/part-00000 $OUT_FLR_WEB/OutStep_1_2


#  Step 1.3 (pure python)
python -m bimeta.create_corpus.create_corpus \
--input $OUT_FLR_WEB/OutStep_1_2 \
--output $OUT_FLR_WEB/OutStep_1_3 \
--dictionary $FOL_LCL_PATH/dictionary.pkl


# Step 2.1
python bimeta/build_overlap_graph/build_overlap_graph.py \
$OUT_HDFS/output_1_1/part-00000 \
--output $OUT_HDFS/output_2_1 \
-r hadoop \
--conf-path $CONF_FILE \
--q_mers $LENGTH_OF_Q_MERS

hdfs dfs -get $OUT_HDFS/output_2_1/part-00000 $OUT_FLR_WEB
mv $OUT_FLR_WEB/part-00000 $OUT_FLR_WEB/OutStep_2_1


# Step 2.2
spark-submit --packages graphframes:graphframes:0.8.1-spark3.0-s_2.12 \
bimeta/build_overlap_graph/connected.py \
--vertices $OUT_FLR_WEB/OutStep_1_1 \
--edges $OUT_FLR_WEB/OutStep_2_1 \
--checkpoint $OUT_HDFS/graphframes_cps \
--output $OUT_HDFS/graphframes_cps/2 \
--num_reads $NUM_SHARED_READS
# --conf "spark.yarn.access.hadoopFileSystems=${OUT_HDFS}" \
# --master yarn \
# --deploy-mode cluster \
# --files /home/tnhan/ServerWeb/BiMeta/userFolder/thesis1/output/OutStep_1_1,/home/tnhan/ServerWeb/BiMeta/userFolder/thesis1/output/OutStep_2_1 \


hdfs dfs -get $OUT_HDFS/graphframes_cps/2/part-00000 $OUT_FLR_WEB
mv $OUT_FLR_WEB/part-00000 $OUT_FLR_WEB/OutStep_2_2


# # Step 3
python bimeta/cluster_groups/clustering.py \
--group $OUT_FLR_WEB/OutStep_2_2 \
--corpus $OUT_FLR_WEB/OutStep_1_3 \
--dictionary $FOL_LCL_PATH/dictionary.pkl \
--species $NUM_OF_SPECIES \
--labels $OUT_FLR_WEB/OutStep_1_1

# Clean HDFS
# hdfs dfs -rm /user/$IN_FILE_NAME
# hdfs dfs -rm -r /user/session_$RND_NO

# for i in {1..10..1}
# do 
#   cp $OUT_FLR_WEB/OutStep_2_1 $OUT_FLR_WEB/OutStep_2_1_time_$i
# done

deactivate # Deactivate Virtual Environment