#!/usr/bin/env python

import sys
from Bio import SeqIO
from Bio.Seq import Seq
import re

def format_read(read):
    z = re.split('[|={,]+', read.description)
    return read.seq, z[3]
    
def load_meta_reads(filename, type='fasta'):
    try:
        seqs = list(SeqIO.parse(filename, type))
        #reads = []
        #labels = []
        is_paired_end = False
        if len(seqs) > 2 and seqs[0].id[-1:] != seqs[1].id[-1:]:
            is_paired_end = True

        label_list = dict()
        label_index = 0
        
        for i in range(0, len(seqs), 2 if is_paired_end else 1):
            read, label = format_read(seqs[i])
            if is_paired_end:
                read2, label2 = format_read(seqs[i + 1])
                read += read2
            print(str(i) + " " + str(read))
              
            if label not in label_list:
                label_list[label] = label_index
                label_index += 1
            print(str(i) + " " + str(label_list[label]))
            
        del seqs
        
    except:
        print('Error when loading file {} '.format(filename))


load_meta_reads(sys.stdin, type='fasta')