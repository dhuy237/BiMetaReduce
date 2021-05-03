from gensim import corpora
from gensim.models.tfidfmodel import TfidfModel
from gensim.models import LogEntropyModel

import json
import re
import sys

sys.path.append("../")  # Add "../" to utils folder path
from utils import globals

# Don't know why cannot use this:
# DICTIONARY_PATH = globals.DATA_PATH + "dictionary.pkl"
# This path is used to save the updated dictionary.pkl file
DICTIONARY_PATH = "/home/dhuy237/thesis/code/bimetaReduce/data/R4_medium/dictionary.pkl"
FILENAME = globals.DATA_PATH + 'output_1_2.txt'

def create_corpus(dictionary, documents, 
                  is_tfidf=False, 
                  smartirs=None, 
                  is_log_entropy=False, 
                  is_normalize=True):
    
    corpus = [dictionary.doc2bow(d, allow_update=True) for d in documents]
    if is_tfidf:
        tfidf = TfidfModel(corpus=corpus, smartirs=smartirs)
        corpus = tfidf[corpus]
    elif is_log_entropy:
        log_entropy_model = LogEntropyModel(corpus, normalize=is_normalize)
        corpus = log_entropy_model[corpus]
    return corpus


def read_file(filename):
    documents = []

    with open(filename) as f:
        content = f.readlines()

    for line in content:
        clean_line = re.sub('[null\t\n\[\]\"]', '', line).replace(' ', '').split(',')[3:]
        documents.append(clean_line)
        
    return documents


def convert2json(corpus):
    result = []
    for i, item in enumerate(corpus):
        item = [list(elem) for elem in item]
        result.append([i, item])
    return result


def save_file(result):
    with open(globals.DATA_PATH+'output_1_3.txt', 'w+') as f:
        for item in result:
            f.write("null\t%s\n" % json.dumps(item))


documents = read_file(FILENAME)
dictionary = corpora.Dictionary.load(DICTIONARY_PATH)

corpus = create_corpus(
            dictionary=dictionary,
            documents=documents,
            is_tfidf=globals.IS_TFIDF,
            smartirs=globals.SMARTIRS,
        )

result = convert2json(corpus)

save_file(result)