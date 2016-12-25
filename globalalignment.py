class GlobalAlignment:
    """Provides global alignment for specific scoring matrix and linear gap penalty for two genes

    Args:
        scoring_mat (str): Path to a file with scoring matrix.
        lgp (int): Linear gap penalty used in calculating dynamic table.
    """

    def __init__(self, scoring_mat, lgp):
        self.mat = self.read_scoring_matrix(scoring_mat)
        self.linear_gap_penalty = lgp

    def read_scoring_matrix(self, ime):
        """Loads scoring matrix from given path/name of the file.
        Args:
            ime (str): Path to a file with scoring matrix.
        Returns:
            dict: Dictionary scoring matrix.
        """
        scoring_matrix = {}
        with open(ime) as f:
            stolpci = f.readline().strip().split()
            for line in f:
                prvi, *vrstica = line.strip().split()
                for i, vrednost in enumerate(vrstica):
                    scoring_matrix[(prvi, stolpci[i])] = int(vrednost)
        return scoring_matrix

    def dynamic_table(self, s, t):
        """Calculates and return dynamic table for given genes

        Args:
            s (str): First gene.
            t (str): Second gene.

        Returns:
            dict: Dictionary dynamic table.
        """
        M = {}
        M[(-1, -1)] = 0  # zgornji levi kot tabele
        n = len(s)
        m = len(t)
        # inicializacija
        for j in range(n):
            M[(-1, j)] = M[(-1, j - 1)] - self.linear_gap_penalty
        for i in range(m):
            M[(i, -1)] = M[(i - 1, -1)] - self.linear_gap_penalty
        # izracun vrednosti v tabeli
        for i in range(m):
            for j in range(n):
                M[(i, j)] = max(
                    M[(i - 1, j)] - self.linear_gap_penalty,
                    M[(i, j - 1)] - self.linear_gap_penalty,
                    M[(i - 1, j - 1)] + self.mat[(t[i], s[j])]
                )

        return M
