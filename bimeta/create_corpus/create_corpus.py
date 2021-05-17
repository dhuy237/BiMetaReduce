from gensim import corpora
from gensim.models.tfidfmodel import TfidfModel
from gensim.models import LogEntropyModel

import json
import re
import argparse
import json
import os 

# import sys
# sys.path.append("../")  # Add "../" to utils folder path
# from bimeta.utils import globals

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help = "Input file")
parser.add_argument("-o", "--output", help = "Output file")
parser.add_argument("-d", "--dictionary", help = "Dictionary file")
args = parser.parse_args()

# Not implemented yet in the Web UI
IS_TFIDF = False
SMARTIRS = None

def create_corpus(dictionary, documents, 
                  is_tfidf=False, 
                  smartirs=None, 
                  is_log_entropy=False, 
                  is_normalize=True):
    
    corpus = [dictionary.doc2bow(d, allow_update=False) for d in documents]
    if is_tfidf:
        tfidf = TfidfModel(corpus=corpus, smartirs=smartirs)
        corpus = tfidf[corpus]
    elif is_log_entropy:
        log_entropy_model = LogEntropyModel(corpus, normalize=is_normalize)
        corpus = log_entropy_model[corpus]

    # Will overwritten the existed file
    # Save new file because the dictionary allow to be updated
    # dictionary.save(args.dictionary)
    
    return corpus


def read_file(filename):
    """
    For reading output_1_2 file

    """
    documents = []

    with open(filename) as f:
        content = f.readlines()

    for line in content:
        clean_line = re.sub('[null\t\n\[\]\"]', '', line).replace(' ', '').split(',')[3:]
        documents.append(clean_line)
        
    return documents


def convert2json(corpus):
    """
    For saving to output_1_3 file

    """
    result = []
    for i, item in enumerate(corpus):
        item = [list(elem) for elem in item]
        result.append([i, item])
    return result


def save_file(result, output_path):
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, 'w+') as f:
        for item in result:
            f.write("null\t%s\n" % json.dumps(item))


documents = read_file(args.input)
dictionary = corpora.Dictionary.load(args.dictionary)
corpus = create_corpus(
            dictionary=dictionary,
            documents=documents,
            is_tfidf=IS_TFIDF,
            smartirs=SMARTIRS,
        )
result = convert2json(corpus)

save_file(result, args.output)
