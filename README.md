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

Accepts two strings and scores their alignment. Return score and alignment of strings. Both local and global alignment are supported.  

#### dendrogram.py

Visualization of differences between alignment of all genes between species. 

### Experiments and results
