conda activate pyenv

# Step 1.1
python load_meta_reads/load_read.py data/test/R4_medium.fna >> data/test/output_1_1.txt

#  Step 1.2
python parallel_create_document/load_create_document.py data/test/output_1_1.txt >> data/test/output_1_2.txt
