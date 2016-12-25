from scipy.cluster.hierarchy import linkage, dendrogram
from matplotlib import pyplot as plt
import numpy

class Dendrograms:

    def __init__(self, rows, columns):
        """
        Provides drawing multiple dendrograms in table
        :param rows: number of rows in table
        :param columns: number of columns in table
        """
        self.rows = rows
        self.columns = columns
        self.labels = ["japonica", "remanei", "brenneri"]
        self.k = 1

    def draw_dendrogram(self, values, title):
        """
        draws one dendrogram in the table

        :param values: list of scores (three integers), either positive
        or negative in following order:
            1. score between remanei and japonica (rem_jap)
            2. score between brenneri and japonica (pb_jap)
            3. score between brenneri and remanei (rem_pb)
        !!! THIS METHOD WORKS ONLY WHEN HIGHER SCORE MEANS LOWER DISTANCE
        WHEN TURNING SCORE INTO DISTANCE, RESULT GETS TWISTED AND DOES NOT
        REPRESENT ACTUAL DISTANCE !!!

        TODO: CAN THIS BE DONE BETTER?

        :param title: title of the dendrogram
        :return:
        """
        values = numpy.array(values)
        if numpy.all(values > 0):
            # all positive
            # normalise (almost, to avoid 1 and then distance 0)
            #  and subtract from one
            y = 1 - (1 / (1.1*max(values))) * values
        elif numpy.all(values < 0):
            # all negative
            # normalise
            values = abs(values)
            y = (1 / (1.1*max(values))) * values
        else:
            # mixed, add offset to make all positive
            values = values - 2 * min(values)
            # minimum is negative
            y = 1 - (1 / (1.1 * max(values))) * values

        plt.subplot(self.rows, self.columns, self.k)
        self.k += 1
        plt.title(title)
        plt.xlabel('species')
        plt.ylabel('distance')
        Z = linkage(y)
        dendrogram(Z, labels=self.labels)