source $HOME/thesisEnv/bin/activate # Activate Virtual Environment

PARAM_JSON=$1 # JSON file parameters
USR_SESSION=$2
IN_FILE_NAME=($(jq -r '.file' $PARAM_JSON))

CONF_FILE=$HOME/mrjob.conf
PERSONAL_LIB=$HOME/bimeta.zip

FOL_LCL_PATH=$HOME/Documents # Folder Local Path is used to save output of pure python program
OUT_FLR_WEB=$HOME/ServerWeb/BiMeta/userFolder/$USR_SESSION/output
GRAPH_FLR_WEB=$HOME/ServerWeb/BiMeta/userFolder/$USR_SESSION/graph
GRAPH_NAME=($(jq -r '.nodeGraph' $PARAM_JSON))
JSON_FLR_WEB=$HOME/ServerWeb/BiMeta/jsonData
OVERVIEW=overview.json

RND_NO=$(echo $RANDOM % 100 + 1 | bc)
hdfs dfs -mkdir /user/session_$RND_NO

USR_HDFS=hdfs:///user
INP_HDFS=hdfs:///user/$IN_FILE_NAME
OUT_HDFS=hdfs:///user/session_$RND_NO

hdfs dfs -put $HOME/ServerWeb/BiMeta/userFolder/$USR_SESSION/input/$IN_FILE_NAME /user

# Python Code Variables START
LENGTHS_OF_K_MERS=($(jq -r '.params.kmer' $PARAM_JSON))
LENGTH_OF_Q_MERS=($(jq -r '.params.lofqmer' $PARAM_JSON))
NUM_SHARED_READS=($(jq -r '.params.sharereads' $PARAM_JSON))
FLAG_NUM_OF_SPECIES=($(jq -r '.params.kNumber' $PARAM_JSON))
if [ "$FLAG_NUM_OF_SPECIES" = "false" ]
then
    NUM_OF_SPECIES=$3
else
    NUM_OF_SPECIES=$FLAG_NUM_OF_SPECIES
fi
# Python Code Variables END

STEP_1_1=($(jq -r '.steps.Step1' $PARAM_JSON))
STEP_1_2=($(jq -r '.steps.Step2' $PARAM_JSON))
STEP_1_3=($(jq -r '.steps.Step3' $PARAM_JSON))
STEP_2_1=($(jq -r '.steps.Step4' $PARAM_JSON))
STEP_2_2=($(jq -r '.steps.Step5' $PARAM_JSON))
STEP_3=($(jq -r '.steps.Step6' $PARAM_JSON))

# Step 1.1
# Start 1.1----------------------------------------------------------------------
echo "" >&2
echo "Step 1.1: Load Meta Read" >&2
echo "" >&2

START_TIME=`date +%s%N`

if [ "$STEP_1_1" = "true" ]
then
    python $HOME/ServerWeb/BiMeta/BimetaCode/bimeta/load_meta_reads/load_read_mr.py $INP_HDFS \
    --output $OUT_HDFS/output_1_1 \
    -r hadoop \
    --conf-path $CONF_FILE

    hdfs dfs -get $OUT_HDFS/output_1_1/part-00000 $OUT_FLR_WEB
    mv $OUT_FLR_WEB/part-00000 $OUT_FLR_WEB/OutStep_1_1
else
    python $HOME/ServerWeb/BiMeta/BimetaCode/bimeta/load_meta_reads/load_read.py \
    --input $HOME/ServerWeb/BiMeta/userFolder/$USR_SESSION/input/$IN_FILE_NAME \
    --output $OUT_FLR_WEB/OutStep_1_1

    hdfs dfs -mkdir $OUT_HDFS/output_1_1
    hdfs dfs -put $OUT_FLR_WEB/OutStep_1_1 $OUT_HDFS/output_1_1    
    hdfs dfs -mv $OUT_HDFS/output_1_1/OutStep_1_1 $OUT_HDFS/output_1_1/part-00000
fi

END_TIME=`date +%s%N`

RUN_TIME=`expr $END_TIME - $START_TIME`
RUN_TIME_IN_S=$(echo "scale = 3; $RUN_TIME / 1000000000" | bc)
echo "{\"Step_1_1\":\"$RUN_TIME_IN_S\"," > $JSON_FLR_WEB/$OVERVIEW
# End----------------------------------------------------------------------


#  Step 1.2 
echo "" >&2
echo "Step 1.2: Create Dictionary and Create Document" >&2
echo "" >&2

python bimeta/parallel_create_document/create_dictionary.py \
--dictionary_path $FOL_LCL_PATH \
--k_mers $LENGTHS_OF_K_MERS

# Start 1.2----------------------------------------------------------------------
START_TIME=`date +%s%N`

if [ "$STEP_1_2" = "true" ]
then
    python bimeta/parallel_create_document/create_document_mr.py \
    $OUT_HDFS/output_1_1/part-00000 \
    --output $OUT_HDFS/output_1_2 \
    -r hadoop \
    --conf-path $CONF_FILE \
    --k_mers $LENGTHS_OF_K_MERS

    hdfs dfs -get $OUT_HDFS/output_1_2/part-00000 $OUT_FLR_WEB
    mv $OUT_FLR_WEB/part-00000 $OUT_FLR_WEB/OutStep_1_2
else
    python bimeta/parallel_create_document/create_document.py \
    --input $OUT_FLR_WEB/OutStep_1_1 \
    --output $OUT_FLR_WEB/OutStep_1_2 \
    --k_mers $LENGTHS_OF_K_MERS

    hdfs dfs -mkdir $OUT_HDFS/output_1_2
    hdfs dfs -put $OUT_FLR_WEB/OutStep_1_2 $OUT_HDFS/output_1_2    
    hdfs dfs -mv $OUT_HDFS/output_1_2/OutStep_1_2 $OUT_HDFS/output_1_2/part-00000
fi

END_TIME=`date +%s%N`

RUN_TIME=`expr $END_TIME - $START_TIME`
RUN_TIME_IN_S=$(echo "scale = 3; $RUN_TIME / 1000000000" | bc)
echo "\"Step_1_2\":\"$RUN_TIME_IN_S\"," >> $JSON_FLR_WEB/$OVERVIEW
# End----------------------------------------------------------------------


# Step 1.3 (pure python)
# Start 1.3----------------------------------------------------------------------
echo "" >&2
echo "Step 1.3: Create Corpus" >&2
echo "" >&2
START_TIME=`date +%s%N`

if [ "$STEP_1_3" = "true" ]
then
    python bimeta/create_corpus/create_corpus_mr.py \
    $OUT_HDFS/output_1_2/part-00000 \
    --output $OUT_HDFS/output_1_3 \
    -r hadoop \
    --conf-path $CONF_FILE \
    --dictionary $FOL_LCL_PATH/dictionary.pkl

    hdfs dfs -get $OUT_HDFS/output_1_3/part-00000 $OUT_FLR_WEB
    mv $OUT_FLR_WEB/part-00000 $OUT_FLR_WEB/OutStep_1_3
else
    python bimeta/create_corpus/create_corpus.py \
    --input $OUT_FLR_WEB/OutStep_1_2 \
    --output $OUT_FLR_WEB/OutStep_1_3 \
    --dictionary $FOL_LCL_PATH/dictionary.pkl

    hdfs dfs -mkdir $OUT_HDFS/output_1_3
    hdfs dfs -put $OUT_FLR_WEB/OutStep_1_3 $OUT_HDFS/output_1_3    
    hdfs dfs -mv $OUT_HDFS/output_1_3/OutStep_1_3 $OUT_HDFS/output_1_3/part-00000
fi    

END_TIME=`date +%s%N`

RUN_TIME=`expr $END_TIME - $START_TIME`
RUN_TIME_IN_S=$(echo "scale = 3; $RUN_TIME / 1000000000" | bc)
echo "\"Step_1_3\":\"$RUN_TIME_IN_S\"," >> $JSON_FLR_WEB/$OVERVIEW
# End----------------------------------------------------------------------


# Step 2.1
# Start 2.1----------------------------------------------------------------------
echo "" >&2
echo "Step 2.1: Build Overlap Graph" >&2
echo "" >&2

START_TIME=`date +%s%N`

if [ "$STEP_2_1" = "true" ]
then
    python bimeta/build_overlap_graph/build_overlap_graph_mr.py \
    $OUT_HDFS/output_1_1/part-00000 \
    --output $OUT_HDFS/output_2_1 \
    -r hadoop \
    --conf-path $CONF_FILE \
    --q_mers $LENGTH_OF_Q_MERS

    hdfs dfs -get $OUT_HDFS/output_2_1/part-00000 $OUT_FLR_WEB
    mv $OUT_FLR_WEB/part-00000 $OUT_FLR_WEB/OutStep_2_1
else
    python bimeta/build_overlap_graph/build_overlap_graph.py \
    --input $OUT_FLR_WEB/OutStep_1_1 \
    --output $OUT_FLR_WEB/OutStep_2_1 \
    --q_mers $LENGTH_OF_Q_MERS \
    --num_reads $NUM_SHARED_READS

    hdfs dfs -mkdir $OUT_HDFS/output_2_1
    hdfs dfs -put $OUT_FLR_WEB/OutStep_2_1 $OUT_HDFS/output_2_1    
    hdfs dfs -mv $OUT_HDFS/output_2_1/OutStep_2_1 $OUT_HDFS/output_2_1/part-00000
fi

END_TIME=`date +%s%N`

RUN_TIME=`expr $END_TIME - $START_TIME`
RUN_TIME_IN_S=$(echo "scale = 3; $RUN_TIME / 1000000000" | bc)
echo "\"Step_2_1\":\"$RUN_TIME_IN_S\"," >> $JSON_FLR_WEB/$OVERVIEW
# End----------------------------------------------------------------------


# Start 2.2----------------------------------------------------------------------
echo "" >&2
echo "Step 2.2: Visualize Graph and Connect Nodes" >&2
echo "" >&2

START_TIME=`date +%s%N`

spark-submit --packages graphframes:graphframes:0.8.1-spark3.0-s_2.12 \
bimeta/build_overlap_graph/visualize_graph.py \
--vertices $OUT_FLR_WEB/OutStep_1_1 \
--edges $OUT_FLR_WEB/OutStep_2_1 \
--output_graph $GRAPH_FLR_WEB/$GRAPH_NAME \
--num_reads $NUM_SHARED_READS

if [ "$STEP_2_2" = "true" ]
then
    spark-submit --packages graphframes:graphframes:0.8.1-spark3.0-s_2.12 \
    bimeta/build_overlap_graph/connected.py \
    --vertices $OUT_FLR_WEB/OutStep_1_1 \
    --edges $OUT_FLR_WEB/OutStep_2_1 \
    --checkpoint $OUT_HDFS/graphframes_cps \
    --output $OUT_HDFS/graphframes_cps/2 \
    --num_reads $NUM_SHARED_READS

    hdfs dfs -get $OUT_HDFS/graphframes_cps/2/part-00000 $OUT_FLR_WEB
    mv $OUT_FLR_WEB/part-00000 $OUT_FLR_WEB/OutStep_2_2
else
    python bimeta/build_overlap_graph/build_connected.py \
    --vertices $OUT_FLR_WEB/OutStep_1_1 \
    --edges $OUT_FLR_WEB/OutStep_2_1 \
    --output $OUT_FLR_WEB/OutStep_2_2 \
    --output_graph $GRAPH_FLR_WEB/$GRAPH_NAME

    hdfs dfs -mkdir $OUT_HDFS/output_2_2
    hdfs dfs -put $OUT_FLR_WEB/OutStep_2_2 $OUT_HDFS/output_2_2    
    hdfs dfs -mv $OUT_HDFS/output_2_2/OutStep_2_2 $OUT_HDFS/output_2_2/part-00000
fi

END_TIME=`date +%s%N`

RUN_TIME=`expr $END_TIME - $START_TIME`
RUN_TIME_IN_S=$(echo "scale = 3; $RUN_TIME / 1000000000" | bc)
echo "\"Step_2_2\":\"$RUN_TIME_IN_S\"}" >> $JSON_FLR_WEB/$OVERVIEW
# End----------------------------------------------------------------------


# Step 3
# Start 3----------------------------------------------------------------------
echo "" >&2
echo "Step 3: Clustering" >&2
echo "" >&2

if [ "$STEP_3" = "true" ]
then
    spark-submit bimeta/cluster_groups/kmeans.py \
    --group $OUT_FLR_WEB/OutStep_2_2 \
    --corpus $OUT_FLR_WEB/OutStep_1_3 \
    --dictionary $FOL_LCL_PATH/dictionary.pkl \
    --species $NUM_OF_SPECIES \
    --labels $OUT_FLR_WEB/OutStep_1_1 \
    --time $JSON_FLR_WEB/$OVERVIEW
else
    python bimeta/cluster_groups/clustering.py \
    --group $OUT_FLR_WEB/OutStep_2_2 \
    --corpus $OUT_FLR_WEB/OutStep_1_3 \
    --dictionary $FOL_LCL_PATH/dictionary.pkl \
    --species $NUM_OF_SPECIES \
    --labels $OUT_FLR_WEB/OutStep_1_1 \
    --time $JSON_FLR_WEB/$OVERVIEW
fi
# End----------------------------------------------------------------------

# Clean HDFS
hdfs dfs -rm /user/$IN_FILE_NAME
# hdfs dfs -rm -r /user/session_$RND_NO

deactivate # Deactivate Virtual Environment