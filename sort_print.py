def sort_print(dict, reverse):
    for gene in (sorted(dict, key=dict.get, reverse=reverse)[:5]):
        print(gene, dict[gene])