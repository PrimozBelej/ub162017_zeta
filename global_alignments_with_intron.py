import genomereader
import align
import dendrogram
from matplotlib import pyplot as plt

genome_pb2 = genomereader.GenomeReader(
    'data/genome/caePb2/chrUn.fa.gz',
    'data/genes/caePb2.augustusGene.txt')

genome_jap1 = genomereader.GenomeReader(
    'data/genome/caeJap1/chrUn.fa.gz',
    'data/genes/caeJap1.augustusGene.txt')


genome_rem3 = genomereader.GenomeReader(
    'data/genome/caeRem3/chrUn.fa.gz',
    'data/genes/caeRem3.augustusGene.txt')

# all genomes have the same genes
#all_gene_names = [g["name"] for g in genome_pb2.genes]
all_gene_names = genome_pb2.genes.keys()

a = align.Align("pam250", 5)  # eg. blosum62, pam250, pam30, rao ... or PAM250.txt (directly from file)

#sum of scores of global alignments of nucleotids over all genes with
# linear gap penalty 5

rem_jap_nucleotid_seq_pam = 0
pb_jap_nucleotid_seq_pam = 0
rem_pb_nucleotid_seq_pam = 0

#sum of scores of alignments of aminoacid sequnces over all genes with
#  pam250 scoring matrix
#  and linear gap penalty 5

rem_jap_aminoacid_seq_pam = 0
pb_jap_aminoacid_seq_pam = 0
rem_pb_aminoacid_seq_pam = 0

for gene in all_gene_names:
    #
    s1 = genome_rem3.gene_with_itron(gene)
    t1 = genome_pb2.gene_with_itron(gene)
    u1 = genome_jap1.gene_with_itron(gene)

    rem_pb_nucleotid_seq_pam += a.global_alignment(s1, t1)[0]
    rem_jap_nucleotid_seq_pam += a.global_alignment(s1, u1)[0]
    pb_jap_nucleotid_seq_pam += a.global_alignment(u1, t1)[0]

    #################################################################
    s2 = genome_rem3.get_amino_acid_sequence_with_itron(gene)
    t2 = genome_pb2.get_amino_acid_sequence_with_itron(gene)
    u2 = genome_jap1.get_amino_acid_sequence_with_itron(gene)

    rem_pb_aminoacid_seq_pam += a.global_alignment(s2, t2)[0]
    rem_jap_aminoacid_seq_pam += a.global_alignment(s2, u2)[0]
    pb_jap_aminoacid_seq_pam += a.global_alignment(u2, t2)[0]

print("Results when comparing nucleotid sequence (including introns)")
print(rem_jap_nucleotid_seq_pam)
print(pb_jap_nucleotid_seq_pam)
print(rem_pb_nucleotid_seq_pam)
print("Results when comparing aminoacid sequence (including introns)")
print(rem_jap_aminoacid_seq_pam)
print(pb_jap_aminoacid_seq_pam)
print(rem_pb_aminoacid_seq_pam)


dend = dendrogram.Dendrograms(1, 2)
# 1row, 2 columns for 2 dendrograms

d1 = [rem_jap_nucleotid_seq_pam, pb_jap_nucleotid_seq_pam, rem_pb_nucleotid_seq_pam]
dend.draw_dendrogram(d1, "nucleotids")
d2 = [rem_jap_aminoacid_seq_pam, pb_jap_aminoacid_seq_pam, rem_pb_aminoacid_seq_pam]
dend.draw_dendrogram(d2, "amino-acid")


plt.show()