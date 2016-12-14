import genomereader

genome_pb2 = genomereader.GenomeReader(
    '../../chC_alignment/chC_alignment_nematode/data/genome/caePb2/chrUn.fa',
    '../../chC_alignment/chC_alignment_nematode/data/genes/caePb2.augustusGene.txt')

genome_jap1=genomereader.GenomeReader(
    '../../chC_alignment/chC_alignment_nematode/data/genome/caeJap1/chrUn.fa',
    '../../chC_alignment/chC_alignment_nematode/data/genes/caeJap1.augustusGene.txt')


genome_rem3=genomereader.GenomeReader(
    '../../chC_alignment/chC_alignment_nematode/data/genome/caeRem3/chrUn.fa',
    '../../chC_alignment/chC_alignment_nematode/data/genes/caeRem3.augustusGene.txt')


# all genomes have the same genes
all_gene_names=[g["name"] for g in genome_pb2.genes]


def read_scoring_matrix(ime):
    scoring_matrix = {}
    with open(ime) as f:
        stolpci = f.readline().strip().split()
        for line in f:
            prvi, *vrstica = line.strip().split()
            for i, vrednost in enumerate(vrstica):
                scoring_matrix[(prvi, stolpci[i])] = int(vrednost)
    return scoring_matrix


def dinamicna_tabela(s, t, mat, lgp):
    M = {}
    M[(-1, -1)] = 0  # zgornji levi kot tabele
    n = len(s)
    m = len(t)
    # inicializacija
    for j in range(n):
        M[(-1, j)] = M[(-1, j - 1)] - lgp
    for i in range(m):
        M[(i, -1)] = M[(i - 1, -1)] - lgp
    # izracun vrednosti v tabeli
    for i in range(m):
        for j in range(n):
            M[(i, j)] = max(
                M[(i - 1, j)] - lgp,
                M[(i, j - 1)] - lgp,
                M[(i - 1, j - 1)] + mat[(t[i], s[j])]
            )

    return M

mat = read_scoring_matrix("blosum62.txt")
lgp = 5  # linear gap penalty

#sum of scores over all genes
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

    M = dinamicna_tabela(s, t, mat, lgp)
    rem_pb += M[(m - 1, n - 1)]

    M= dinamicna_tabela(u,s,mat,lgp)
    rem_jap += M[n - 1, o - 1]

    M =dinamicna_tabela(u,t, mat, lgp)
    pb_jap += M[m - 1, o - 1]

print(rem_jap)
print(pb_jap)
print(rem_pb)






