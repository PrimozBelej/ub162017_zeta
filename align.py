from Bio.SubsMat import MatrixInfo
from itertools import product
from itertools import combinations
from operator import sub
import numpy as np

def read_scoring_matrix(ime):
    """Loads scoring matrix from given name from Biopython MatrixInfo or from given path to file.
    Args:
        ime (str): Name of BioPythons MatrixInfo variable or name to a file with scoring matrix (eg. *.txt).
    Returns:
        dict: Dictionary scoring matrix.
    """
    scoring_matrix = {}

    try:
        scoring_matrix_temp = getattr(MatrixInfo, ime)  # use built in matrix if exists
        for key1, key2 in scoring_matrix_temp:  # biopython built in matrices are not symmetric -.-
            scoring_matrix[(key1, key2)] = scoring_matrix_temp[(key1, key2)]
            scoring_matrix[(key2, key1)] = scoring_matrix_temp[(key1, key2)]
    except AttributeError:  # else load from file
        with open(ime) as f:
            stolpci = f.readline().strip().split()
            for line in f:
                prvi, *vrstica = line.strip().split()
                for i, vrednost in enumerate(vrstica):
                    scoring_matrix[(prvi, stolpci[i])] = int(vrednost)
    return scoring_matrix


class Align:
    def __init__(self, scoring_matrix, linear_gap_penalty):

        self.sm = read_scoring_matrix(scoring_matrix)
        self.lgp = linear_gap_penalty

    def dinamicna_tabela_global(self, s, t):

        m = len(t)
        n = len(s)
        M = {}
        M[(-1, -1)] = 0  # zgornji levi kot tabele
        P = {}
        # inicializacija
        for j in range(n):
            M[(-1, j)], P[(-1, j)] = M[(-1, j - 1)] - self.lgp, (-1, j - 1)
        for i in range(m):
            M[(i, -1)], P[(i, -1)] = M[(i - 1, -1)] - self.lgp, (i - 1, -1)
        # izracun vrednosti v tabeli
        for i in range(m):
            for j in range(n):
                M[(i, j)], P[(i, j)] = max(
                    (M[(i - 1, j)] - self.lgp, (i - 1, j)),
                    (M[(i, j - 1)] - self.lgp, (i, j - 1)),
                    (M[(i - 1, j - 1)] + self.sm[(t[i], s[j])], (i - 1, j - 1))
                )
        return M, P

    def dinamicna_tabela_local(self, s, t):

        m = len(t)
        n = len(s)
        M = {}
        M[(-1, -1)] = 0  # zgornji levi kot tabele
        P = {}
        # inicializacija
        for j in range(n):
            M[(-1, j)], P[(-1, j)] = 0, (-2, -2)
        for i in range(m):
            M[(i, -1)], P[(i, -1)] = 0, (-2, -2)
        # izracun vrednosti v tabeli
        for i in range(m):
            for j in range(n):
                M[(i, j)], P[(i, j)] = max(
                    (M[(i - 1, j)] - self.lgp, (i - 1, j)),
                    (M[(i, j - 1)] - self.lgp, (i, j - 1)),
                    (M[(i - 1, j - 1)] + self.sm[(t[i], s[j])], (i - 1, j - 1)),
                    (0, (-2, -2))
                )
        return M, P

    def trace_back_local(self, s, t):
        M, P = self.dinamicna_tabela_local(s, t)
        razlika, (i2, j2) = max((vrednost, indeksa) for indeksa, vrednost in M.items())
        i1, j1 = i2, j2
        s1 = ""
        t1 = ""
        while M[(i1, j1)] != 0:
            x,y = P[(i1, j1)]
            if y == j1:
                s1 = "-" + s1
                t1 = t[i1] + t1
            elif x == i1:
                t1 = "-" + t1
                s1 = s[j1] + s1
            else:
                s1 = s[j1] + s1
                t1 = t[i1] + t1
            i1, j1 = P[(i1, j1)]
        return razlika, s1, t1

    def trace_back_global(self, s, t):
        m = len(t)
        n = len(s)
        M, P = self.dinamicna_tabela_global(s,t)
        s = iter(s)
        t = iter(t)
        el = (m - 1, n - 1)
        pot = []
        while True:
            pot.append(el)
            if el == (-1, -1):
                break
            el = P[el]
        pot = list(reversed(pot))
        s1 = []
        t1 = []
        for i in range(len(pot) - 1):
            prvi = pot[i]
            drugi = pot[i + 1]
            if prvi[0] == drugi[0]:
                t1.append("-")
            else:
                t1.append(next(t))
            if prvi[1] == drugi[1]:
                s1.append("-")
            else:
                s1.append(next(s))
        return M[(m - 1, n - 1)], "".join(s1), "".join(t1)

    def kazen(self, sez, indeksi, komb):
        """Punishment - auxiliary function for multiple alignment"""
        return sum(1 if not ((i * j == 1 and niz1[ind1] == niz2[ind2]) or (i == 0 and j == 0)) else 0 for
                   (niz1, ind1, i), (niz2, ind2, j) in combinations(zip(sez, indeksi, komb), 2))

    def dinamicna_tabela_multiple(self, sez):
        M = {}
        M[tuple(np.repeat(-1, len(sez)))] = 0  # zgornji levi kot tabele
        P = {}

        indexes = [tuple(np.subtract(dxs, np.ones(len(sez), dtype="int"))) for dxs in
                   list(np.ndindex(tuple(len(i) + 1 for i in sez)))]

        for index in indexes:
            if index != tuple(np.repeat(-1, len(sez))):
                kombinacije = product([0, 1], repeat=len(sez))
                next(kombinacije)  # same nicle izpustimo
                # ce kljuca ni, vzamemo zelo majhno vrednost (-100000), ki nima vpliva, saj iscemo max
                M[index], P[index] = max((M.get(tuple(map(sub, index, komb)), -100000) - self.kazen(sez, index, komb),
                                          tuple(map(sub, index, komb))) for komb in kombinacije)

        return M, P

    def trace_back_multiple(self, sez):
        M, P = self.dinamicna_tabela_multiple(sez)

        el = tuple(len(s) - 1 for s in sez)
        iters = [iter(s) for s in sez]
        pot = []
        while True:
            pot.append(el)
            if el == tuple(np.repeat(-1, len(sez))):
                break
            el = P[el]
        pot = list(reversed(pot))
        words = [[] for i in range(len(sez))]
        for i in range(len(pot) - 1):
            prvi = pot[i]
            drugi = pot[i + 1]
            for i in range(len(sez)):
                if prvi[i] == drugi[i]:
                    words[i].append("-")
                else:
                    words[i].append(next(iters[i]))
        return M[tuple(len(s)-1 for s in sez)], ["".join(s) for s in words]

    def global_alignment(self, s, t):
        """global alignment of strings s and t,
        returns score and aligned strings s1 and t1"""
        d, s1, t1 = self.trace_back_global(s, t)
        return d, s1, t1

    def local_alignment(self, s, t):
        """local alignment of strings s and t,
        returns score and aligned substrings s1 and t1"""
        d, s1, t1 = self.trace_back_local(s, t)
        return d, s1, t1

    def multiple_alignment(self, sez):
        """multiple alignment of all strings in array sez
        returns score and array of aligned substrings

        NOTE: IT IS NOT EFFICIENT and can take a lot of time
        """
        d, words = self.trace_back_multiple(sez)
        return d, words
