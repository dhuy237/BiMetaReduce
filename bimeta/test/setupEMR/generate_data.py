filenames = "/home/dhuy237/thesis/code/bimetaReduce/bimeta/data/test/R4.fna"

# Generate n copies from data file
# with open('/home/dhuy237/thesis/data/Simulated/generate/data-10gb.txt', 'w') as outfile:
    
#     for i in range(3):
#         with open(filenames) as infile:
#             for line in infile:
#                 outfile.write(line)

# Read first n line in text file
with open("/home/dhuy237/thesis/data/HMP/hmp-metadata.fasta", 'w') as output:

    with open("/home/dhuy237/thesis/data/HMP/hmp.fasta") as myfile:
        head = [next(myfile) for x in range(11)]
        for line in head:
            output.write(line)

