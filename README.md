gene2csm.py
==================

## Description

`gene2csm.py` is a tool for efficient design of crRNA guide sequences for RNA mediated RNA processing enzymes. It is designed to satisfy the following features of the oligonucleotide:
1. Specificity (within target organism)

   Every sequence is the subject of BLAST[1] database search of coding and non-coding RNA sequences within the target organism to exclude off-traget effects of highly homologous sequences.
1. Isoform prevalence

   Considered sequence should span the region covered by maximal number of annotated isoforms for a complete target knock-down.
1. Absence of alternative, annotated start codons downstream of the cut site

   Considered sequence should not have any annotated start codons downstream of the target site on the protein coding isoforms.
1. Avoid exon-exon boundaries

   Sequence should not contain any annotated splice sites.
1. Not spanning know SNP sites

   Sequence should not contain any variable nucleotides contained in the dbSNP database and Ensembl resources.


Other features taken into consideration:
1. No low complexity regions

   Sequence should not contain any interspersed repeats and low complexity sequences, masked by RepeatMasker
1. Balanced GC content

   Sequences should have GC content within the given limits.
1. No target mRNA secondary structures within the binding region

   Target sequence should avoid highly structured regions of the transcript to assure the highest accessibility to the RNA strand. The RNA secondary structure modeling is performed with the ViennaRNA package [2].
1. No self-complementarity

   Considered sequence should not form stable homodimeric duplexes.
1. No crRNA secondary structures

   The sequence should not contain any local secondary structures.


## Output description:
The output table is sorted by the `score` column and contains 50 best scoring cRNA sequences characterized as follows:
* 1st column contains a unique index number; the numbering follows the genomic start position of the crRNA in descending order, starting from 0;
* `seq` column contains a sequence of the putative crRNA; reverse-complementary to the chosen target transcript;
* `GC` column contains the %GC content of the crRNA
* `ent` column contains the mean positional entropy value of the target mRNA sequence in the binding position of crRNA; it describes the structural well-definednes of the region; the higher the better;
* `dG_AA` column contains the change in Gibbs free energy of the homodimer duplex created by two crRNA oligos; the higher (less negative) the better- the less stable the homodimer complex is
* `G_A` column stores the Gibbs free energy of the monomeric crRNA, the higher (less negative) the better- the less structured the monomer is;
* `bitscore` column contains the bitscore value of the best alignment of the sequence to the sequences from the database as defined by the BLAST algorithm; the lower the better
* `nident` column contains the number of identical matching nucleotides in the best scoring blast alignment of the crRNA to the sequences from the database; the lower the better;
* `chr` column contains the name of the chromosome the target loci is on;
* `start` column contains the start position of the crRNA on the chromosome; 0-based;
* `end` column contains the stop position of the crRNA on the chromosome;
non-inclusive
* `score` column contains the cumulative rank score calculated from the entropy value and the bitscore value of the sequence only; no other characteristic is taken into consideration; the lower the better 

References:
1. Stephen F. Altschul, Thomas L. Madden, Alejandro A. Schäffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs", Nucleic Acids Res. 25:3389-3402.
1.  Lorenz, Ronny and Bernhart, Stephan H. and Höner zu Siederdissen, Christian and Tafer, Hakim and Flamm, Christoph and Stadler, Peter F. and Hofacker, Ivo L.; ViennaRNA Package 2.0; Algorithms for Molecular Biology, 6:1 26, 2011, doi:10.1186/1748-7188-6-26
