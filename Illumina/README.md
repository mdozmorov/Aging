Extract genomic coordinates of CpG probes on Illumina 850K, 450K and Illumina 27K arrays.
===

## 850K

- `Illumina850K.R` generates annotation file and a hg19 BED file

## 450K, 27K

Use "Download full table..." button to get full annotation table in text format for

- Illumina 450K: [http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL13534](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL13534)

- Illumina 27K: [http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL8490](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL8490)

Gzip the files. Run `make` to create BED files. Remember to [liftOver](http://genome.ucsc.edu/cgi-bin/hgLiftOver) Illumina 27K genomic coordinates from hg18 to hg19 human genome assembly coordinate system.

- `snpsites.txt` - file of SNP containing probes, from `wget --no-check-certificate https://www.rforge.net/IMA/snpsites.txt`

- `48639-non-specific-probes-Illumina450k.xlsx` and `48640-polymorphic-CpGs-Illumina450k.xlsx` - Data from Chen YA, et.al. "[Discovery of cross-reactive probes and polymorphic CpGs in the Illumina Infinium HumanMethylation450 microarray](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3592906/)." Epigenetics. [List of non-specific probes](http://www.sickkids.ca/MS-Office-Files/Research/Weksberg%20Lab/48639-non-specific-probes-Illumina450k.xlsx) - 29,233 non-specific 'cg' probes, 1,736 non-specific 'ch' probes; [List of polymorphic CpGs](http://www.sickkids.ca/MS-Office-Files/Research/Weksberg%20Lab/48640-polymorphic-CpGs-Illumina450k.xlsx) - 70,899 records (66,877 unique probes) about CpGs containing SNPs at or near single base extension (SBE) position, 316,034 records (220,582 unique probes) having SNPs in probe sequences.

- `NonspecificAndPolymorphic.xlsx` - a list of potential nonspecific probes and polymorphic probes of Illumina Human 27k Methylation Array, [BrainCloud](http://braincloud.jhmi.edu/downloads.htm), [source](http://braincloud.jhmi.edu/NonspecificAndPolymorphic.zip)

- "[raharris.Illumina_Infinium_450K_Array.pdf](http://genboree.org/theCommons/attachments/2296/raharris.Illumina_Infinium_450K_Array.pdf)", "[dataprocessing.pdf](www.bristol.ac.uk/caite/geocode/newcastleshortcourse/dataprocessing.pdf)" - data preprocessing recommendations