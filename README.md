# Analysis of age-associated genomic regions

Scripts and data files for the epigenomic enrichment analysis of age-associated differentially methylated regions, and genes changing their expression with age. The results are published in Dozmorov, Mikhail G. “Polycomb Repressive Complex 2 Epigenomic Signature Defines Age-Associated Hypermethylation and Gene Expression Changes.” Epigenetics 10, no. 6 (2015): 484–95. https://doi.org/10.1080/15592294.2015.1040619.

This README also collects papers and notes about DNA methylation and aging, aka epigenetic clock.

# Table of content

* [Files and folders description](#Files-and-folders-description)
* [Papers](#Papers)

## Files and folders description

- `data` - scripts to extract genomic coordinates of the age-associated genomic regions. The description of this folder contains the list of studies and the data available from them.

- `Illumina` - scripts to extract genomic coordinates of all CpG probes from Illumina Infinium 27K and Illumina Infinium 450K

- `GenomeRunner` - the results of the epigenomic enrichment analysis, and the scripts to process them

- `R.Aging` - R scripts for visualization of the results.

	- `genes_promoters_extract.Rmd` - extract promoters (+2000..-500bp arount TSS) of age-associated _genes_ (`data` folder).

## Papers

- Horvath, Steve, and Kenneth Raj. “DNA Methylation-Based Biomarkers and the Epigenetic Clock Theory of Ageing.” Nature Reviews Genetics, April 11, 2018. https://doi.org/10.1038/s41576-018-0004-3. - Epigenetic clock review. Horvath clock, Hannum, PhenoAge, other developments. Age-related conditions linked to epigenetic age. Epigenetic age is 40% heritable.

