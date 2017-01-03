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

Accepts two strings and scores their alignment. Return score and alignment of strings. Local, global and multiple alignment are supported.  

#### dendrogram.py

Visualization of differences between alignment of all genes between species. 

### Experiments and results

In the folowing section we present the experiments and their results. We test our implementation of gene alignment with different scoring matrices, with junk DNA and without etc. 

#### Global alignment of exons with different scoring matrices

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

#### Global alignment if introns are included.

In this experiment we perform global alignment on species pair by pair. For each two genes we take a whole sequance and perform alignment on them. The sum of score is our metric to compare alignments of species pair by pair. After alignment using nucleotides we also experiment with aminoacid sequence. 

The resulting dendrogram is shown below:

![alt tag]( http://shrani.si/f/2H/bp/2sudstxP/with-itron.png)

We can see results are different as in the first experiment. Now remanei and brenneri have higher similarity compared to japonica, no matter which sequence type was used. And also distances between them are slightly larger. 

| Species pair | Nucleotide | Amino-acid |
|--------------|------------|------------|
| rem_jap      | 14249     | -3584  |
| pb_jap       | 13593      | -3812     |
| rem_pb       | 20000      | -2548  |

As in the first experiment, the scores obtained using amino acid sequences are significantly lower.
