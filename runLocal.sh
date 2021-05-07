conda activate pyenv

export DATA_PATH=bimeta/data/test
export INPUT_FILE=R4_medium.fna

# Step 1.1
python bimeta/load_meta_reads/load_read.py $DATA_PATH/$INPUT_FILE --output $DATA_PATH/output_1_1

#  Step 1.2
python -m bimeta.parallel_create_document.create_dictionary \
--dictionary_path $DATA_PATH

python -m bimeta.parallel_create_document.load_create_document $DATA_PATH/output_1_1/part-00000 \
--output $DATA_PATH/output_1_2

#  Step 1.3
python -m bimeta.create_corpus.create_corpus \
--input $DATA_PATH/output_1_2/part-00000 \
--output $DATA_PATH/output_1_3.txt \
--dictionary $DATA_PATH/dictionary.pkl

# Step 2.1
python -m bimeta.build_overlap_graph.build_overlap_graph $DATA_PATH/output_1_1/part-00000 \
--output $DATA_PATH/output_2_1

# Step 2.2
spark-submit --packages graphframes:graphframes:0.8.1-spark3.0-s_2.12 \
--py-files utils.zip \
bimeta/build_overlap_graph/connected.py \
--vertices $DATA_PATH/output_1_1/part-00000 \
--edges $DATA_PATH/output_2_1/part-00000 \
--checkpoint "/home/dhuy237/graphframes_cps/2" \
--output "/home/dhuy237/graphframes_cps/2/5"

hdfs dfs -get /home/dhuy237/graphframes_cps/2/5/part-00000 $DATA_PATH/output_2_2/

# Step 3
python -m bimeta.cluster_groups.clustering \
--group "/home/dhuy237/thesis/code/bimetaReduce/bimeta/data/test/output_2_2/part-00000" \
--corpus "/home/dhuy237/thesis/code/bimetaReduce/bimeta/data/test/output_1_3.txt" \
--dictionary "/home/dhuy237/thesis/code/bimetaReduce/bimeta/data/test/dictionary.pkl" \
--species 2 \
--labels "/home/dhuy237/thesis/code/bimetaReduce/bimeta/data/test/output_1_1/part-00000"