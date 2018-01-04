## `blastdb` directory
BLAST index directory. Put BLAST fasta index databases here.

`danRer_e91_allrna.nal` index alias file provided as example for search against combined cDNA and ncRNA fasta sequences of _Danio rerio_ GRCz10 genome assembly Ensembl release 91. 

Use the folowing example commands on Ensembl FASTA files:
```sh
makeblastdb -in '<Species_genome_assembly>.cdna.all.fa'\
            -parse_seqids\
            -dbtype nucl\
            -title '<Name>_cdna';
makeblastdb -in '<Species_genome_assemly>.GRCz10.ncrna.fa'\
            -parse_seqids\
            -dbtype nucl\
            -title '<Name>_ncrna';
blastdb_aliastool -dblist "./<Species_genome_assembly>.cdna.all.fa ./<Species_genome_assembly>.ncrna.fa"\
                  -dbtype nucl\
                  -out '<File_name_prefix>'\
                  -title "<Name>_allrna";
```
