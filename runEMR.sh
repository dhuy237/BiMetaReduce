python bimeta/load_meta_reads/load_read_mr.py -r emr s3://mrjob-bucket-bimeta/input/R4.fna --conf-path mrjob.conf --output-dir=s3://mrjob-bucket-bimeta/output_1_1/

python bimeta/parallel_create_document/create_document_mr.py -r emr s3://mrjob-bucket-bimeta/output_1_1/part-00000 --conf-path mrjob.conf --output-dir=s3://mrjob-bucket-bimeta/output_1_2/