Extract genomic coordinates of CpG probes on Illumina 450K and Illumina 27K arrays.
===

Prerequisites: Use "Download full table..." button to get full annotation table in text format for

- Illumina 450K: [http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL13534](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL13534)

- Illumina 27K: [http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL8490](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL8490)

Gzip the files. Run `make` to create BED files. Remember to [liftOver](http://genome.ucsc.edu/cgi-bin/hgLiftOver) Illumina 27K genomic coordinates from hg18 to hg19 human genome assembly coordinate system.