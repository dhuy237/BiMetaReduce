conda activate pyenv

DATA_PATH=bimeta/data/testS
INPUT_FILE=R4_medium.fna
LENGTHS_OF_K_MERS=4
LENGTH_OF_Q_MERS=30
NUM_SHARED_READS=45
NUM_OF_SPECIES=2
USR_HDFS=hdfs:///user/graphframes_cps/4
OVERVIEW=overview.json

START_TIME=`date +%s%N`

# Step 1.1
python bimeta/load_meta_reads/load_read.py \
--input $DATA_PATH/$INPUT_FILE \
--output $DATA_PATH/output_1_1/part-00000

END_TIME=`date +%s%N`

RUN_TIME=`expr $END_TIME - $START_TIME`
RUN_TIME_IN_S=$(echo "scale = 3; $RUN_TIME / 1000000000" | bc)
echo "{\"Step_1_1\":\"$RUN_TIME_IN_S\"," > $DATA_PATH/$OVERVIEW


# #  Step 1.2
python bimeta/parallel_create_document/create_dictionary.py \
--dictionary_path $DATA_PATH \
--k_mers $LENGTHS_OF_K_MERS

START_TIME=`date +%s%N`

python bimeta/parallel_create_document/create_document.py \
--input $DATA_PATH/output_1_1/part-00000 \
--output $DATA_PATH/output_1_2/part-00000 \
--k_mers $LENGTHS_OF_K_MERS

END_TIME=`date +%s%N`

RUN_TIME=`expr $END_TIME - $START_TIME`
RUN_TIME_IN_S=$(echo "scale = 3; $RUN_TIME / 1000000000" | bc)
echo "\"Step_1_2\":\"$RUN_TIME_IN_S\"," >> $DATA_PATH/$OVERVIEW


#  Step 1.3
START_TIME=`date +%s%N`

python bimeta/create_corpus/create_corpus.py \
--input $DATA_PATH/output_1_2/part-00000 \
--output $DATA_PATH/output_1_3/part-00000 \
--dictionary $DATA_PATH/dictionary.pkl

END_TIME=`date +%s%N`

RUN_TIME=`expr $END_TIME - $START_TIME`
RUN_TIME_IN_S=$(echo "scale = 3; $RUN_TIME / 1000000000" | bc)
echo "\"Step_1_3\":\"$RUN_TIME_IN_S\"," >> $DATA_PATH/$OVERVIEW


# # Step 2.1
START_TIME=`date +%s%N`

python bimeta/build_overlap_graph/build_overlap_graph.py \
--input $DATA_PATH/output_1_1/part-00000 \
--output $DATA_PATH/output_2_1/part-00000 \
--q_mers $LENGTH_OF_Q_MERS \
--num_reads $NUM_SHARED_READS

END_TIME=`date +%s%N`

RUN_TIME=`expr $END_TIME - $START_TIME`
RUN_TIME_IN_S=$(echo "scale = 3; $RUN_TIME / 1000000000" | bc)
echo "\"Step_2_1\":\"$RUN_TIME_IN_S\"," >> $DATA_PATH/$OVERVIEW


# Step 2.2
START_TIME=`date +%s%N`

python bimeta/build_overlap_graph/build_connected.py \
--vertices $DATA_PATH/output_1_1/part-00000 \
--edges $DATA_PATH/output_2_1/part-00000 \
--output $DATA_PATH/output_2_2/part-00000 \
--output_graph $DATA_PATH

END_TIME=`date +%s%N`

RUN_TIME=`expr $END_TIME - $START_TIME`
RUN_TIME_IN_S=$(echo "scale = 3; $RUN_TIME / 1000000000" | bc)
echo "\"Step_2_2\":\"$RUN_TIME_IN_S\"}" >> $DATA_PATH/$OVERVIEW


# Step 3
python bimeta/cluster_groups/clustering.py \
--group $DATA_PATH/output_2_2/part-00000 \
--corpus $DATA_PATH/output_1_3/part-00000 \
--dictionary $DATA_PATH/dictionary.pkl \
--species $NUM_OF_SPECIES \
--labels $DATA_PATH/output_1_1/part-00000 \
--time $DATA_PATH/$OVERVIEW