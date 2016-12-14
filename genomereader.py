from Bio import SeqIO


INT_FIELDS = {0, 4, 5, 6, 7, 8, 16}
INT_LIST_FIELDS = {9, 10, 15}


class GenomeReader:

    def __init__(self, genome_path, genes_path):
        self.genome = None
        self.genes = None
        self.load_genome(genome_path)
        self.load_genes(genes_path)

    def load_genome(self, path):
        self.genome = SeqIO.read(path, 'fasta').seq

    def load_genes(self, path):
        self.genes = []
        with open(path, 'r') as in_file:
            labels = in_file.readline().strip().split()
            for row in in_file:
                values = row.strip().split('\t')
                for i, v in enumerate(values):
                    if i in INT_FIELDS:
                        values[i] = int(v)
                    elif i in INT_LIST_FIELDS:
                        values[i] = [int(n) for n in v.strip().split(',')[:-1]]
                self.genes.append(dict(zip(labels, values)))
