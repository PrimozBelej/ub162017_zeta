import genomereader
import align
import matplotlib.pyplot as plt
import numpy as np


GENOMES_PATH = '../../chC_alignment/chC_alignment_nematode/data/genome/'
GENES_PATH = '../../chC_alignment/chC_alignment_nematode/data/genes/'

genome_pb2 = genomereader.GenomeReader(GENOMES_PATH+'caePb2/chrUn.fa',
                                       GENES_PATH+'caePb2.augustusGene.txt')
genome_jap1 = genomereader.GenomeReader(GENOMES_PATH+'caeJap1/chrUn.fa',
                                        GENES_PATH+'caeJap1.augustusGene.txt')
genome_rem3 = genomereader.GenomeReader(GENOMES_PATH+'caeRem3/chrUn.fa',
                                        GENES_PATH+'caeRem3.augustusGene.txt')

gene_names = [g['name'] for g in genome_pb2.genes]

alignment_scores_pb_jap = {}
alignment_scores_pb_rem = {}
alignment_scores_jap_rem = {}

amino_alignment_scores_pb_jap = {}
amino_alignment_scores_pb_rem = {}
amino_alignment_scores_jap_rem = {}

amino_a = align.Align('blosum62', 1)
dna_a = align.Align('dnafull.txt', 1)

for name in gene_names:
    n1 = genome_pb2.join_exons(name)
    n2 = genome_jap1.join_exons(name)
    n3 = genome_rem3.join_exons(name)

    alignment_scores_pb_jap[name] = dna_a.global_alignment(n1, n2)[0]
    alignment_scores_pb_rem[name] = dna_a.global_alignment(n1, n3)[0]
    alignment_scores_jap_rem[name] = dna_a.global_alignment(n2, n3)[0]

    a1 = genome_pb2.get_amino_acid_sequence(name)
    a2 = genome_jap1.get_amino_acid_sequence(name)
    a3 = genome_rem3.get_amino_acid_sequence(name)

    amino_alignment_scores_pb_jap[name] = amino_a.global_alignment(a1, a2)[0]
    amino_alignment_scores_pb_rem[name] = amino_a.global_alignment(a1, a3)[0]
    amino_alignment_scores_jap_rem[name] = amino_a.global_alignment(a2, a3)[0]


z_labels = ['Pb-Jap', 'Pb-Rem', 'Jap-Rem']

n = len(gene_names)
ind = np.arange(n)
width = 0.2

# Plot of nucleotide alignment scores
fig, ax = plt.subplots()
bars1 = ax.bar(ind-width, [alignment_scores_pb_jap[g] for g in gene_names],
               width, color='r')
bars2 = ax.bar(ind, [alignment_scores_pb_rem[g] for g in gene_names],
               width, color='g')
bars3 = ax.bar(ind + width, [alignment_scores_jap_rem[g] for g in
                               gene_names], width, color='b')
ax.set_ylabel('Score')
ax.set_title('Nucleotide alignment scores')
ax.set_xticks(ind+width)
ax.set_xticklabels(gene_names, rotation=45)
ax.legend((bars1[0], bars2[0], bars3[0]), ('Pb-Jap', 'Pb-Rem', 'Jap-Rem'))

# Plot of amino acid alignment scores
fig_amino, ax_amino = plt.subplots()
bars_a_1 = ax_amino.bar(ind-width, [amino_alignment_scores_pb_jap[g] for g in
                                gene_names], width, color='r')
bars_a_2 = ax_amino.bar(ind, [amino_alignment_scores_pb_rem[g] for g in
                                gene_names], width, color='g')
bars_a_3 = ax_amino.bar(ind+width, [amino_alignment_scores_jap_rem[g] for g in
                                gene_names], width, color='b')
ax_amino.set_ylabel('Score')
ax_amino.set_title('Amino acid alignment scores')
ax_amino.set_xticks(ind+width)
ax_amino.set_xticklabels(gene_names, rotation=45)
ax_amino.legend((bars_a_1[0], bars_a_2[0], bars_a_3[0]), ('Pb-Jap', 'Pb-Rem',
                                                          'Jap-Rem'))

plt.show()

