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
from Bio.Seq import translate, reverse_complement
import gzip

INT_FIELDS = {0, 4, 5, 6, 7, 8, 16}
INT_LIST_FIELDS = {9, 10, 15}

def load_compress_genome (path):
    """Loads genome data from compressed fasta file."""
    with gzip.open(path, 'rb') as f:
        f.readline()
        data = f.read().decode("utf-8").replace("\n", "")
    return data

class GenomeReader:
    """Provides parsing and storing of genomes and gene annotations.

    Args:
        genome_path (str): Path to a file with genome data.
        genes_path (str): Path to a file with gene annotation data.

    Attributes:
        genome (Bio.seq): Genome sequence.
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
        if path[-2:] == 'gz':
            self.genome = load_compress_genome(path)
        else:
            self.genome = SeqIO.read(path, 'fasta').seq

    def load_genes(self, path):
        """Loads gene annotation data from a csv file.

        Args:
            path (str): Path to an uncompressed fasta file.
        """
        self.genes = {}
        with open(path, 'r') as in_file:
            labels = in_file.readline().strip().split()
            for row in in_file:
                values = row.strip().split('\t')
                for i, v in enumerate(values):
                    if i in INT_FIELDS:
                        values[i] = int(v)
                    elif i in INT_LIST_FIELDS:
                        values[i] = [int(n) for n in v.strip().split(',')[:-1]]

                self.genes[values[1]] = dict(zip(labels, values))

    def gene_by_name(self, name):
        """Finds and returns gene anotation with a given name.

        Args:
            name (str): Name of a gene.

        Returns:
            str: Dictionary with gene annotations. If gene with a
                given name is not found the returned value is None.
        """
        return self.genes[name]

    def gene_seq(self, name):
        """Returns sequence of a specified gene.

        Returned sequence doesn not take into account the strand direction.

        Args:
            name (str): Name of a gene.

        Returns:
            str: Sequence of a specified gene. Empty string in case of
                non-existent gene.
        """
        gene = self.gene_by_name(name)
        return str(self.genome)[gene['txStart']:gene['txEnd']]

    def exons_seq(self, name, exon=-1):
        """Returns list of exon sequences.

        Returned sequences do not take into account the strand direction.

        Args:
            name (str)
            exon (int): Number of exon. Default value of -1 means all
                exons are required.

        Returns:
            [str] or str: List of exon sequences in order of
                occurrence.
                In case of specified exon number, single sequence is
                returned.
        """
        gene = self.gene_by_name(name)
        exons_pos = list(zip(gene['exonStarts'], gene['exonEnds']))
        if exon == -1:
            return [str(self.genome)[p[0]:p[1]] for p in exons_pos]

        if exon >= 0 and exon < gene['exonCount']:
            return str(self.genome)[exons_pos[exon][0]:exons_pos[exon][1]]

    def join_exons(self, name):
        """ Joins all exons of a given gene"""
        list_of_exons = self.exons_seq(name)
        
        r = "".join(list_of_exons)
            
        if self.genes[name]["strand"] == "+":
            #lower - case letters are mapped into upper case
            # (Source: wikipedia - FASTA_format)
            return r.upper()
        else:
            return reverse_complement(r).upper()


    def get_amino_acid_sequence(self, name):
        """return aminoacid sequence of a given gene"""
        
        if self.genes[name]["strand"] == "+":
            return translate(self.join_exons(name), stop_symbol="")
        else: #negative strand
            # for negative strand we take reverse complement
            return translate(reverse_complement(self.join_exons(name)), stop_symbol="")


