gene2csm.py
==================

## Example usage

See the included Jupyter [notebook](./gene2csm.ipynb).

## Dependencies

* `python3`
* `jupyter` (for notebook integration)
* `numpy`
* `matplotlib`
* `pandas`
* `biopython`
* `gffutils` [:link:](http://daler.github.io/gffutils/ "gffutils")
* `pybedtools` [:link:](https://daler.github.io/pybedtools/ "pybedtools")
* `pysam`
* `openpyxl` (for XLSX file generation)
* bedtools v2.26.0 [:link:](http://bedtools.readthedocs.io/en/latest/ "bedtools")
* NCBI BLAST 2.7.1 [:link:](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download "BLAST")
* ViennaRNA Package 2.4.3 [:link:](https://www.tbi.univie.ac.at/RNA/
  "ViennaRNA")

## Installation

#### Quick instalation with `pipenv` [:link:](http://pipenv.readthedocs.io/en/latest/ "pipenv")
To satisfy `python` dependencies install `pipenv` using `pip` or your prefered
package manager. Then run
```
pipenv install
```
from the program directory. This will install all the required modules with
their dependencies, as specified in the provided `Pipfile` 
and `Pipfile.lock` into the new virtual environment.

#### Quick installation with `conda` [:link:](https://conda.io/docs/ "conda")
For non-pythonic dependencies run
```
conda install -c bedtools blast viennarna
```

## Requirements

Ensembl GTF annotaion file, soft masked genomic sequence FASTA file and GVF
variation file are required, unless working only on user input (without any
database sequence lookup). 


Ensembl cDNA FASTA file and ncRNA FASTA file are required for BLAST database
searches.

Download the required files from [Ensembl FTP site](https://www.ensembl.org/info/data/ftp/index.html).

## Description

`gene2csm.py` is a tool for efficient design of crRNA guide sequences for RNA mediated RNA processing enzymes. It is designed to satisfy the following features of the oligonucleotide:

1. Specificity (within target organism)

   Every sequence is the subject of BLAST \[[1](#r1)\] database search of coding and non-coding RNA sequences within the target organism to exclude off-traget effects of highly homologous sequences.

2. Isoform prevalence

   Considered sequence should span the region covered by maximal number of annotated isoforms for a complete target knock-down.

3. Absence of alternative, annotated start codons downstream of the cut site

   Considered sequence should not have any annotated start codons downstream of the target site on the protein coding isoforms.

4. Avoid exon-exon boundaries

   Sequence should not contain any annotated splice sites.

5. Not spanning know SNP sites

   Sequence should not contain any variable nucleotides contained in the dbSNP database and Ensembl resources.


Other features taken into consideration:

1. No low complexity regions

   Sequence should not contain any interspersed repeats and low complexity sequences, masked by RepeatMasker

2. Balanced GC content

   Sequences should have GC content within the given limits.

3. No target mRNA secondary structures within the binding region

   Target sequence should avoid highly structured regions of the transcript to assure the highest accessibility to the RNA strand. The RNA secondary structure modeling is performed with the ViennaRNA package \[[2](#r2)\].

4. No self-complementarity

   Considered sequence should not form stable homodimeric duplexes.

5. No crRNA secondary structures

   The sequence should not contain any local secondary structures.


## Output description:

The output table is sorted by the `score` column and contains 50 best scoring cRNA sequences characterized as follows:

* 1st column contains a unique _index_ number; the numbering follows the genomic start position of the crRNA in descending order, starting from 0;
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

---

References:
1. <a name="r1"></a>Stephen F. Altschul, Thomas L. Madden, Alejandro A. Schäffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), _"Gapped BLAST and PSI-BLAST: a new generation of protein database search programs"_, Nucleic Acids Res. 25:3389-3402.
2. <a name="r2"></a>Lorenz, Ronny and Bernhart, Stephan H. and Höner zu Siederdissen, Christian and Tafer, Hakim and Flamm, Christoph and Stadler, Peter F. and Hofacker, Ivo L.; _ViennaRNA Package 2.0_; Algorithms for Molecular Biology, 6:1 26, 2011, doi:10.1186/1748-7188-6-26
