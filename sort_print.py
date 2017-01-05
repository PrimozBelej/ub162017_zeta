import numpy as np
from collections import Counter


def sort_print(dict, reverse):
    for gene in (sorted(dict, key=dict.get, reverse=reverse)[:5]):
        print(gene, dict[gene])


def calc_similar_percentage(gene_score_dicts, all_gene_names):
    maxes = []
    for gene in all_gene_names:
        dists = [i[gene] for i in gene_score_dicts]
        maxes.append(np.argmax(dists))
    c = Counter(maxes)
    return [(i, c[i] / float(len(maxes)) * 100.0) for i in c]
