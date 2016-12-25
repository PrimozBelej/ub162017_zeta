import genomereader
import globalalignment

# read genomes
genome_pb2 = genomereader.GenomeReader(
    '../../chC_alignment/chC_alignment_nematode/data/genome/caePb2/chrUn.fa',
    '../../chC_alignment/chC_alignment_nematode/data/genes/caePb2.augustusGene.txt')

genome_jap1 = genomereader.GenomeReader(
    '../../chC_alignment/chC_alignment_nematode/data/genome/caeJap1/chrUn.fa',
    '../../chC_alignment/chC_alignment_nematode/data/genes/caeJap1.augustusGene.txt')


genome_rem3 = genomereader.GenomeReader(
    '../../chC_alignment/chC_alignment_nematode/data/genome/caeRem3/chrUn.fa',
    '../../chC_alignment/chC_alignment_nematode/data/genes/caeRem3.augustusGene.txt')

# all genomes have the same genes
all_gene_names = [g["name"] for g in genome_pb2.genes]

global_align = globalalignment.GlobalAlignment("blosum62.txt", 5)
# sum of scores over all genes
rem_jap = 0
pb_jap = 0
rem_pb = 0

for gene in all_gene_names:

    # ALERT!
    # small letters are treated the same as the one in caps
    # small letters denote repetitions, but this should be checked!
    s = genome_rem3.join_exons(gene).upper()
    t = genome_pb2.join_exons(gene).upper()
    u = genome_jap1.join_exons(gene).upper()

    m = len(t)
    n = len(s)
    o = len(u)

    M = global_align.dynamic_table(s, t)
    rem_pb += M[(m - 1, n - 1)]

    M = global_align.dynamic_table(u, s)
    rem_jap += M[n - 1, o - 1]

    M = global_align.dynamic_table(u, t)
    pb_jap += M[m - 1, o - 1]

print(rem_jap)
print(pb_jap)
print(rem_pb)



