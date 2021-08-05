from mrjob.job import MRJob
from mrjob.protocol import JSONProtocol

from gensim import corpora
from gensim.models.tfidfmodel import TfidfModel
from gensim.models import LogEntropyModel
import itertools as it
from Bio.Seq import Seq

def gen_kmers(klist):
    bases = ["A", "C", "G", "T"]
    kmers_list = []

    for k in klist:
        kmers_list += ["".join(p) for p in it.product(bases, repeat=k)]

    # reduce a half of k-mers due to symmetry
    kmers_dict = dict()
    for myk in kmers_list:
        k_reverse_complement = Seq(myk).reverse_complement()
        if not myk in kmers_dict and not str(k_reverse_complement) in kmers_dict:
            kmers_dict[myk] = 0

    return list(kmers_dict.keys())


def create_corpus(
    dictionary,
    documents,
    is_tfidf=False,
    smartirs=None,
    is_log_entropy=False,
    is_normalize=True
):
    corpus = dictionary.doc2bow(documents, allow_update=False)
    if is_tfidf:
        tfidf = TfidfModel(corpus=corpus, smartirs=smartirs)
        corpus = tfidf[corpus]
    elif is_log_entropy:
        log_entropy_model = LogEntropyModel(corpus, normalize=is_normalize)
        corpus = log_entropy_model[corpus]
    return corpus


class CreateCorpus(MRJob):

    INPUT_PROTOCOL = JSONProtocol
    
    def mapper(self, _, line):

        k_mers_set = [gen_kmers([4])]
        dictionary = corpora.Dictionary(k_mers_set)
        
        corpus = create_corpus(
            dictionary=dictionary,
            documents=line[3]
        )
        yield None, (line[0], corpus)

    def reducer(self, key, values):
        for value in values:
            yield key, value


CreateCorpus.run()
