## `blastdb` directory
BLAST index directory. Put BLAST fasta index databases here.

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
