import genomereader
import align
import dendrogram
from matplotlib import pyplot as plt
from sort_print import sort_print, calc_similar_percentage

# read genomes
genome_pb2 = genomereader.GenomeReader(
    '../../chC_alignment/chC_alignment_nematode/data/genome/caePb2/chrUn.fa.gz',
    '../../chC_alignment/chC_alignment_nematode/data/genes/caePb2.augustusGene.txt')

genome_jap1 = genomereader.GenomeReader(
    '../../chC_alignment/chC_alignment_nematode/data/genome/caeJap1/chrUn.fa.gz',
    '../../chC_alignment/chC_alignment_nematode/data/genes/caeJap1.augustusGene.txt')

genome_rem3 = genomereader.GenomeReader(
    '../../chC_alignment/chC_alignment_nematode/data/genome/caeRem3/chrUn.fa.gz',
    '../../chC_alignment/chC_alignment_nematode/data/genes/caeRem3.augustusGene.txt')

# all genomes have the same genes
all_gene_names = genome_pb2.genes.keys()

a = align.Align("pam250", 5)  # eg. blosum62, pam250, pam30, rao ... or PAM250.txt (directly from file)

#sum of scores of global alignments of nucleotids over all genes with
# linear gap penalty 5

rem_jap_nucleotid_seq_pam = 0
pb_jap_nucleotid_seq_pam = 0
rem_pb_nucleotid_seq_pam = 0

rem_jap_nucleotid_gene_distances = {}
pb_jap_nucleotid_gene_distances = {}
rem_pb_nucleotid_gene_distances = {}

#sum of scores of alignments of aminoacid sequnces over all genes with
#  pam250 scoring matrix
#  and linear gap penalty 5

rem_jap_aminoacid_seq_pam = 0
pb_jap_aminoacid_seq_pam = 0
rem_pb_aminoacid_seq_pam = 0

for gene in all_gene_names:

    s1 = genome_rem3.join_exons(gene)
    t1 = genome_pb2.join_exons(gene)
    u1 = genome_jap1.join_exons(gene)

    rem_pb_nucleotid_seq_pam += a.global_alignment(s1, t1)[0]
    rem_pb_nucleotid_gene_distances[gene] = a.global_alignment(s1, t1)[0]
    rem_jap_nucleotid_seq_pam += a.global_alignment(s1, u1)[0]
    rem_jap_nucleotid_gene_distances[gene] = a.global_alignment(s1, u1)[0]
    pb_jap_nucleotid_seq_pam += a.global_alignment(u1, t1)[0]
    pb_jap_nucleotid_gene_distances[gene] = a.global_alignment(u1, t1)[0]

    #################################################################
    s2 = genome_rem3.get_amino_acid_sequence(gene)
    t2 = genome_pb2.get_amino_acid_sequence(gene)
    u2 = genome_jap1.get_amino_acid_sequence(gene)

    rem_pb_aminoacid_seq_pam += a.global_alignment(s2, t2)[0]
    rem_jap_aminoacid_seq_pam += a.global_alignment(s2, u2)[0]
    pb_jap_aminoacid_seq_pam += a.global_alignment(u2, t2)[0]

print("Results when comparing nucleotid sequence")
print(rem_jap_nucleotid_seq_pam)
print(pb_jap_nucleotid_seq_pam)
print(rem_pb_nucleotid_seq_pam)
print("Genes with lowest alignment scores")
print("rem_jap:")
sort_print(rem_jap_nucleotid_gene_distances, False)
print("pb_jap:")
sort_print(pb_jap_nucleotid_gene_distances, False)
print("rem_pb:")
sort_print(rem_pb_nucleotid_gene_distances, False)
print("Genes with highest alingment scores")
print("rem_jap:")
sort_print(rem_jap_nucleotid_gene_distances, True)
print("pb_jap:")
sort_print(pb_jap_nucleotid_gene_distances, True)
print("rem_pb:")
sort_print(rem_pb_nucleotid_gene_distances, True)

print("Results when comparing aminoacid sequence")
print(rem_jap_aminoacid_seq_pam)
print(pb_jap_aminoacid_seq_pam)
print(rem_pb_aminoacid_seq_pam)

print("Results of percentages of similarities according to each gene in genome")
gene_percentages = calc_similar_percentage([rem_pb_nucleotid_gene_distances, rem_jap_nucleotid_gene_distances, pb_jap_nucleotid_gene_distances], all_gene_names)
print("rem_pb: " + str(gene_percentages[0][1]))
print("rem_jap: " + str(gene_percentages[1][1]))
print("pb_jap: " + str(gene_percentages[2][1]))

dend = dendrogram.Dendrograms(1, 2)
# 1row, 2 columns for 2 dendrograms

d1 = [rem_jap_nucleotid_seq_pam, pb_jap_nucleotid_seq_pam, rem_pb_nucleotid_seq_pam]
dend.draw_dendrogram(d1, "nucleotids")
d2 = [rem_jap_aminoacid_seq_pam, pb_jap_aminoacid_seq_pam, rem_pb_aminoacid_seq_pam]
dend.draw_dendrogram(d2, "amino-acid")

plt.show()
