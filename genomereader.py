"""Parsing and storing of genome data along with gene annotations.

Attributes:
    INT_FIELDS ({int}): Set of indices of fields in gene annotation
        data that should be parsed as integers.
    INT_LIST_FIELDS ({int}): Set of indices of fields in gene
        annotation data that should be parsed as lists of integers. By
        convention, these fields are comma separated and terminated by
        an extra comma.
"""

from Bio import SeqIO


INT_FIELDS = {0, 4, 5, 6, 7, 8, 16}
INT_LIST_FIELDS = {9, 10, 15}


class GenomeReader:
    """Provides parsing and storing of genomes and gene annotations.

    Args:
        genome_path (str): Path to a file with genome data.
        genes_path (str): Path to a file with gene annotation data.

    Attributes:
        genome (str): Genome sequence.
        genes ([dict]): List of dictionaries containing gene
            annotations.
    """
    def __init__(self, genome_path, genes_path):
        self.genome = None
        self.genes = None
        self.load_genome(genome_path)
        self.load_genes(genes_path)

    def load_genome(self, path):
        """Loads genome data from uncompressed fasta file.

        Args:
            path (str): Path to an uncompressed fasta file.
        """
        self.genome = SeqIO.read(path, 'fasta').seq

    def load_genes(self, path):
        """Loads gene annotation data from a csv file.
        
        Args:
            path (str): Path to an uncompressed fasta file.
        """
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
