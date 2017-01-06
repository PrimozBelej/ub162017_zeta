# Evolution of genes - project C1
### Motivation

### Goals

The goal of the project was to implement global alignment algorithm. We scored each nucleotide in a given sequence by how many times it alings within the dataset. We experiment with different scoring functions and alingments. We also compared difference in alingment between aminoacid and nuqleotid sequence.

### Dataset

The dataset consist of genome of 3 different nematodes (round worms) - japanica, brenneri and remanei. It is common for this group of nematodes that their genome contains a vast ammount of redundancy. This means that not whole genes are translated into proteins. Parts that end up in protein are called exons nad parts that do not are called introns. Part of genes given in dataset contains at least two exons. 

### Methods

#### genomereader.py

Provides reading and parsing of the dataset. In addition to loading a whole genome we can access individual gene. Furthermore we can access one or all exons belonging to gene. For instructions and usage examples see [wiki](https://github.com/PrimozBelej/ub162017_zeta/wiki/Uporaba-modula-genomereader.py)

#### alignment.py

Accepts two strings and scores their alignment. Return score and alignment of strings. Local, global and multiple alignment are supported. Mutliple alignment is deprecated - it works only on short strings, in our case we have genes with more than 700 characters, which is computationaly unacceptable, so we did not use it in our project.

#### dendrogram.py

Visualization of differences between alignment of all genes between species. 

### Experiments and results

In the folowing section we present the experiments and their results. We test our implementation of gene alignment with different scoring matrices, with junk DNA and without etc. 

#### Global alignment of exons 

In the first experiment we perform global alignment on species pair by pair. For each two genes we extract exons and perform alignment on them. The sum of score is our metric to compare alignments of species pair by pair. After alignment using nucleotides we also experiment with aminoacid sequence. 

The resulting dendrogram is shown below:

![alt tag](http://shrani.si/f/2G/RH/454B7FSx/1/dendrogramglobal.png)

As we can see better results were obtained using nucleotides (note that scores were inverted into distances - higher score translate into lower distance). From dendrogram we can also conclude that remanei and japonica have higher similarity compared to brenner, no matter which sequence type was used. 

| Species pair | Nucleotide | Amino-acid |
|--------------|------------|------------|
| rem_jap      | 13848      | -2736      |
| pb_jap       | 9112       | -4108      |
| rem_pb       | 12144      | -3562      |

As we see from table above, the scores obtained using amino acid sequences are significantly lower. 

| Species pair | Nucleotide gene similarities (%)|
|--------------|------------|
| rem_jap      | 44.83      |
| pb_jap       | 31.03      |
| rem_pb       | 24.14      |

From table above, we see which species pair is the most similar to each other, according to scores of each gene seperately. We see that rem_jap is the most similar pair.

We also took a look at which genes have the highest and the lowest alignment score. The results are shown below, where 5 genes with lowest (first table) and highest score(second table) for each pair are presented with their name and alignment score.

| rem_jap      | pb_jap       | rem_pb      |
|--------------|--------------|-------------|
| g13065 -1241 | g8079 -663   | g8079 -1620 |
| g3039 -866   | g5838 -653   | g3149 -1020 |
| g1158 -524   | g12334 -592  | g5838 -463  |
| g14459 -423  | g13065 -492  | g5059 -319  |
| g10223 -319  | g10223 -344  | g992 -289   |

As we can see, some genes are appearing more than once (g13065, g5838, g10223). We can conclude that similar genes are causing lower alignment scores in different species pairs.

| rem_jap      | pb_jap       | rem_pb      |
|--------------|--------------|-------------|
| g5889 1214   | g5889 1675   | g13065 1378 |
| g832 1165    | g13779 1323  | g10223 1356 |
| g8079 1103   | g5621 1076   | g5889 1094  |
| g14404 1022  | g6414 1014   | g13779 867  |
| g13198 986   | g14404 938   | g12334 864  |

Similar conclusion can be drawn here. Genes that has the best aligment are common in more than one pair. 

#### Global alignment of exons with different scoring matrices

The goal here was to discover if different scoring matrices produce same results ie. if alignment between remanei and japonica will still be the highest scored alignment. We ran the alignment with the folowing matrices: blosum62, pam250, pam30, and rao. The results for nucleotide sequences are avaliable in the table below.

| Species pair | blosum62 | pam250 | pam30 | rao    |
|--------------|----------|--------|-------|--------|
| rem_jap      | 14655    | 13848  | 16267 | 114918 |
| pb_jap       | 10082    | 9112   | 12533 | 104259 |
| rem_pb       | 13630    | 12144  | 15967 | 117715 |

From the obtained results we see that remanei and japonica still have the best alignment score, except when using rao scoring matrix, where remanei and brenneri are scored better. In all cases the lowest scored alignment is between brenneri and japonica. We can conclude that even though choosing appropriate scoring matrix is important it rarely affects the results to such an extent that the order of species pairs is different. 

#### Global alignment if introns are included.

In this experiment we perform global alignment on species pair by pair. For each two genes we take a whole sequence and perform alignment on them. The sum of score is our metric to compare alignments of species pair by pair. After alignment using nucleotides we also experiment with aminoacid sequence. 

The resulting dendrogram is shown below:

![alt tag](http://shrani.si/f/3s/d8/1mWNyeFt/with-itron.png)

We can see results are different as in the first experiment. Now remanei and brenneri have higher similarity compared to japonica, no matter which sequence type was used. And also distances between them are slightly larger. 

| Species pair | Nucleotide | Amino-acid |
|--------------|------------|------------|
| rem_jap      | 14249      | -3584      |
| pb_jap       | 13593      | -3812      |
| rem_pb       | 20000      | -2548      |

As in the first experiment, the scores obtained using amino acid sequences are significantly lower.
