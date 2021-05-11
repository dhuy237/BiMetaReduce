conda activate pyenv

export DATA_PATH=bimeta/data/test
export INPUT_FILE=R4_medium.fna
export LENGTHS_OF_K_MERS=4
export LENGTH_OF_Q_MERS=30
export NUM_SHARED_READS=45
export NUM_OF_SPECIES=2


# Step 1.1
python bimeta/load_meta_reads/load_read.py $DATA_PATH/$INPUT_FILE --output $DATA_PATH/output_1_1


#  Step 1.2
python bimeta/parallel_create_document/create_dictionary.py \
--dictionary_path $DATA_PATH \
--k_mers $LENGTHS_OF_K_MERS

python bimeta/parallel_create_document/load_create_document.py \
$DATA_PATH/output_1_1/part-00000 \
--output $DATA_PATH/output_1_2 \
--k_mers $LENGTHS_OF_K_MERS


#  Step 1.3
python bimeta/create_corpus/create_corpus.py \
--input $DATA_PATH/output_1_2/part-00000 \
--output $DATA_PATH/output_1_3.txt \
--dictionary $DATA_PATH/dictionary.pkl


# Step 2.1
python bimeta/build_overlap_graph/build_overlap_graph.py \
$DATA_PATH/output_1_1/part-00000 \
--output $DATA_PATH/output_2_1 \
--q_mers $LENGTH_OF_Q_MERS

# Somehow the output of MR job return many files, usually just part-00000
cat $DATA_PATH/output_2_1/* > $DATA_PATH/output_2_1.txt

# Step 2.2
spark-submit --packages graphframes:graphframes:0.8.1-spark3.0-s_2.12 \
bimeta/build_overlap_graph/connected.py \
--vertices $DATA_PATH/output_1_1/part-00000 \
--edges $DATA_PATH/output_2_1.txt \
--checkpoint "/home/dhuy237/graphframes_cps/4" \
--output "/home/dhuy237/graphframes_cps/4/5" \
--num_reads $NUM_SHARED_READS

mkdir $DATA_PATH/output_2_2/
hdfs dfs -get /home/dhuy237/graphframes_cps/4/5/part-00000 $DATA_PATH/output_2_2/


# Step 3
python bimeta/cluster_groups/clustering.py \
--group $DATA_PATH/output_2_2/part-00000 \
--corpus $DATA_PATH/output_1_3.txt \
--dictionary $DATA_PATH/dictionary.pkl \
--species $NUM_OF_SPECIES \
--labels $DATA_PATH/output_1_1/part-00000