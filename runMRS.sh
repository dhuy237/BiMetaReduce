PARAM_JSON=/home/dhuy237/thesis/code/bimetaReduce/bimeta/data/testMRS/21_05_19_09_42_33.json


DATA_PATH=bimeta/data/testMRS
INPUT_FILE=($(jq -r '.file' $PARAM_JSON))

LENGTHS_OF_K_MERS=($(jq -r '.params.kmer' $PARAM_JSON))
LENGTH_OF_Q_MERS=($(jq -r '.params.lofqmer' $PARAM_JSON))
NUM_SHARED_READS=($(jq -r '.params.sharereads' $PARAM_JSON))

FLAG_NUM_OF_SPECIES=($(jq -r '.params.kNumber' $PARAM_JSON))

if [ "$FLAG_NUM_OF_SPECIES" = "false" ]
then
    NUM_OF_SPECIES=2
else
    NUM_OF_SPECIES=$FLAG_NUM_OF_SPECIES
fi

USR_HDFS=hdfs:///user/graphframes_cps/1
OVERVIEW=overview.json
OUTPUT_GRAPH=($(jq -r '.nodeGraph' $PARAM_JSON))

# Parameters for select Sequential mode or MR mode
# true: MR mode
# false: Sequential mode
STEP_1_1=($(jq -r '.steps.Step1' $PARAM_JSON))
STEP_1_2=($(jq -r '.steps.Step2' $PARAM_JSON))
STEP_1_3=($(jq -r '.steps.Step3' $PARAM_JSON))
STEP_2_1=($(jq -r '.steps.Step4' $PARAM_JSON))
STEP_2_2=($(jq -r '.steps.Step5' $PARAM_JSON))
STEP_3=($(jq -r '.steps.Step6' $PARAM_JSON))


# Step 1.1
# Start 1.1----------------------------------------------------------------------
START_TIME=`date +%s%N`

if [ "$STEP_1_1" = "true" ]
then
    python bimeta/load_meta_reads/load_read_mr.py \
    $DATA_PATH/$INPUT_FILE \
    --output $DATA_PATH/output_1_1
else
    python bimeta/load_meta_reads/load_read.py \
    --input $DATA_PATH/$INPUT_FILE \
    --output $DATA_PATH/output_1_1/part-00000
fi

END_TIME=`date +%s%N`

RUN_TIME=`expr $END_TIME - $START_TIME`
RUN_TIME_IN_S=$(echo "scale = 3; $RUN_TIME / 1000000000" | bc)
echo "{\"Step_1_1\":\"$RUN_TIME_IN_S\"," > $DATA_PATH/$OVERVIEW
# End----------------------------------------------------------------------


# Step 1.2
# Start 1.2----------------------------------------------------------------------
START_TIME=`date +%s%N`

python bimeta/parallel_create_document/create_dictionary.py \
--dictionary_path $DATA_PATH \
--k_mers $LENGTHS_OF_K_MERS


if [ "$STEP_1_2" = "true" ]
then
    python bimeta/parallel_create_document/create_document_mr.py \
    $DATA_PATH/output_1_1/part-00000 \
    --output $DATA_PATH/output_1_2 \
    --k_mers $LENGTHS_OF_K_MERS
else
    python bimeta/parallel_create_document/create_document.py \
    --input $DATA_PATH/output_1_1/part-00000 \
    --output $DATA_PATH/output_1_2/part-00000 \
    --k_mers $LENGTHS_OF_K_MERS
fi

END_TIME=`date +%s%N`

RUN_TIME=`expr $END_TIME - $START_TIME`
RUN_TIME_IN_S=$(echo "scale = 3; $RUN_TIME / 1000000000" | bc)
echo "\"Step_1_2\":\"$RUN_TIME_IN_S\"," >> $DATA_PATH/$OVERVIEW
# End----------------------------------------------------------------------


# Step 1.3
# Start 1.3----------------------------------------------------------------------
START_TIME=`date +%s%N`

if [ "$STEP_1_3" = "true" ]
then
    python bimeta/create_corpus/create_corpus_mr.py \
    $DATA_PATH/output_1_2/part-00000 \
    --output $DATA_PATH/output_1_3 \
    --dictionary $DATA_PATH/dictionary.pkl
else
    python bimeta/create_corpus/create_corpus.py \
    --input $DATA_PATH/output_1_2/part-00000 \
    --output $DATA_PATH/output_1_3/part-00000 \
    --dictionary $DATA_PATH/dictionary.pkl
fi

END_TIME=`date +%s%N`

RUN_TIME=`expr $END_TIME - $START_TIME`
RUN_TIME_IN_S=$(echo "scale = 3; $RUN_TIME / 1000000000" | bc)
echo "\"Step_1_3\":\"$RUN_TIME_IN_S\"," >> $DATA_PATH/$OVERVIEW
# End----------------------------------------------------------------------


# Step 2.1
# Start 2.1----------------------------------------------------------------------
START_TIME=`date +%s%N`

if [ "$STEP_2_1" = "true" ]
then
    python bimeta/build_overlap_graph/build_overlap_graph_mr.py \
    $DATA_PATH/output_1_1/part-00000 \
    --output $DATA_PATH/output_2_1 \
    --q_mers $LENGTH_OF_Q_MERS
else
    python bimeta/build_overlap_graph/build_overlap_graph.py \
    --input $DATA_PATH/output_1_1/part-00000 \
    --output $DATA_PATH/output_2_1/part-00000 \
    --q_mers $LENGTH_OF_Q_MERS \
    --num_reads $NUM_SHARED_READS
fi

END_TIME=`date +%s%N`

RUN_TIME=`expr $END_TIME - $START_TIME`
RUN_TIME_IN_S=$(echo "scale = 3; $RUN_TIME / 1000000000" | bc)
echo "\"Step_2_1\":\"$RUN_TIME_IN_S\"," >> $DATA_PATH/$OVERVIEW
# End----------------------------------------------------------------------


# Step 2.2
# Start 2.2----------------------------------------------------------------------
spark-submit --packages graphframes:graphframes:0.8.1-spark3.0-s_2.12 \
bimeta/build_overlap_graph/visualize_graph.py \
--vertices $DATA_PATH/output_1_1/part-00000 \
--edges $DATA_PATH/output_2_1/part-00000 \
--output_graph $DATA_PATH/$OUTPUT_GRAPH \
--num_reads $NUM_SHARED_READS

START_TIME=`date +%s%N`

if [ "$STEP_2_2" = "true" ]
then
    cat $DATA_PATH/output_2_1/* > $DATA_PATH/output_2_1.txt

    spark-submit --packages graphframes:graphframes:0.8.1-spark3.0-s_2.12 \
    bimeta/build_overlap_graph/connected.py \
    --vertices $DATA_PATH/output_1_1/part-00000 \
    --edges $DATA_PATH/output_2_1.txt \
    --checkpoint $USR_HDFS \
    --output $USR_HDFS/output \
    --num_reads $NUM_SHARED_READS

    mkdir $DATA_PATH/output_2_2/
    hdfs dfs -get $USR_HDFS/output/part-00000 $DATA_PATH/output_2_2/
else
    python bimeta/build_overlap_graph/build_connected.py \
    --vertices $DATA_PATH/output_1_1/part-00000 \
    --edges $DATA_PATH/output_2_1/part-00000 \
    --output $DATA_PATH/output_2_2/part-00000 \
    --output_graph $DATA_PATH
fi

END_TIME=`date +%s%N`

RUN_TIME=`expr $END_TIME - $START_TIME`
RUN_TIME_IN_S=$(echo "scale = 3; $RUN_TIME / 1000000000" | bc)
echo "\"Step_2_2\":\"$RUN_TIME_IN_S\"}" >> $DATA_PATH/$OVERVIEW
# End----------------------------------------------------------------------


# Step 3
# Start 3----------------------------------------------------------------------
if [ "$STEP_3" = "true" ]
then
    spark-submit bimeta/cluster_groups/kmeans.py \
    --group $DATA_PATH/output_2_2/part-00000 \
    --corpus $DATA_PATH/output_1_3/part-00000 \
    --dictionary $DATA_PATH/dictionary.pkl \
    --species $NUM_OF_SPECIES \
    --labels $DATA_PATH/output_1_1/part-00000 \
    --time $DATA_PATH/$OVERVIEW \
    --output $DATA_PATH/outputGroup
else
    python bimeta/cluster_groups/clustering.py \
    --group $DATA_PATH/output_2_2/part-00000 \
    --corpus $DATA_PATH/output_1_3/part-00000 \
    --dictionary $DATA_PATH/dictionary.pkl \
    --species $NUM_OF_SPECIES \
    --labels $DATA_PATH/output_1_1/part-00000 \
    --time $DATA_PATH/$OVERVIEW \
    --output $DATA_PATH/outputGroup
fi
# End----------------------------------------------------------------------
