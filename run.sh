conda activate pyenv

export DATA_PATH=bimeta/data/test
export INPUT_FILE=R4_medium.fna

# Step 1.1
python bimeta/load_meta_reads/load_read.py $DATA_PATH/$INPUT_FILE > $DATA_PATH/output_1_1.txt

#  Step 1.2
python -m bimeta.parallel_create_document.create_dictionary --dictionary_path $DATA_PATH
python -m bimeta.parallel_create_document.load_create_document $DATA_PATH/output_1_1.txt > $DATA_PATH/output_1_2.txt

#  Step 1.3
python -m bimeta.create_corpus.create_corpus --input $DATA_PATH/output_1_2.txt --output $DATA_PATH/output_1_3.txt --dictionary $DATA_PATH/dictionary.pkl