conda activate pyenv

# Step 1.1
python bimeta/load_meta_reads/load_read.py bimeta/data/test/R4_medium.fna >> bimeta/data/test/output_1_1.txt

#  Step 1.2
python -m bimeta.parallel_create_document.load_create_document bimeta/data/test/output_1_1.txt >> bimeta/data/test/output_1_2.txt