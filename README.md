# Analysis of age-associated genomic regions

Scripts and data files for the epigenomic enrichment analysis of age-associated differentially methylated regions, and genes changing their expression with age. The results are published in Dozmorov, Mikhail G. “Polycomb Repressive Complex 2 Epigenomic Signature Defines Age-Associated Hypermethylation and Gene Expression Changes.” Epigenetics 10, no. 6 (2015): 484–95. https://doi.org/10.1080/15592294.2015.1040619.

This README also collects papers and notes about DNA methylation and aging, aka epigenetic clock.

# Table of content

* [Papers](#Papers)
  * [Human epigenetic clock](#Human-epigenetic-clock)
  * [Mouse epigenetic clock](#Mouse-epigenetic-clock)
  * [Reviews](#Reviews)
* [Files and folders description](#Files-and-folders-description)
  * [data](#data)
    * [Mouse](#mouse)

## Papers

### Human epigenetic clock

- Weidner, Carola Ingrid, Qiong Lin, Carmen Maike Koch, Lewin Eisele, Fabian Beier, Patrick Ziegler, Dirk Olaf Bauerschlag, et al. “Aging of Blood Can Be Tracked by DNA Methylation Changes at Just Three CpG Sites.” Genome Biology 15, no. 2 (2014): R24. https://doi.org/10.1186/gb-2014-15-2-r24. - 3-CpG-age predictive model, selected from 102 CpGs. Pearson correlation, then LM, formula. Hypermethylated are enriched in CGIs, bivalent histone modifications. Non-tissue specific (refs to tissue specific). Sex influences age. Epigenomic enrichment using Fisher's exact. Intro about CpG-age epigenomic enrichments. Prediction studies.

- Bocklandt, Sven, Wen Lin, Mary E. Sehl, Francisco J. Sánchez, Janet S. Sinsheimer, Steve Horvath, and Eric Vilain. “Epigenetic Predictor of Age.” PloS One 6, no. 6 (2011): e14821. https://doi.org/10.1371/journal.pone.0014821. - Methylation age predictor. Illumina 27K,  88 sites in or near 80 genes. 34 male twins, GSE28746. Three sites, EDARADD, TOM1L1, NPTX2 explain >73% of age. ComBat-corrected data for prediction. [Table S1.](https://journals.plos.org/plosone/article/file?type=supplementary&id=info:doi/10.1371/journal.pone.0014821.s002) - 88 loci significantly correlated with age TargetID represents the exact Illumina probe on the array, Chr: chromosome number, Gene_ID: NCBI Gene database locator, Symbol: gene name, r: correlation coefficient, p-value: significance of correlation, q-value: significance corrected for multiple comparisons.

### Mouse epigenetic clock

- BI Ageing Clock Team, Thomas M. Stubbs, Marc Jan Bonder, Anne-Katrien Stark, Felix Krueger, Ferdinand von Meyenn, Oliver Stegle, and Wolf Reik. “Multi-Tissue DNA Methylation Age Predictor in Mouse.” Genome Biology 18, no. 1 (December 2017). https://doi.org/10.1186/s13059-017-1203-5. - Mouse epigenetic age clock, 329 CpGs. Elastic net regression. https://github.com/EpigenomeClock/MouseEpigeneticClock


### Reviews

- Horvath, Steve, and Kenneth Raj. “DNA Methylation-Based Biomarkers and the Epigenetic Clock Theory of Ageing.” Nature Reviews Genetics, April 11, 2018. https://doi.org/10.1038/s41576-018-0004-3. - Epigenetic clock review. Horvath clock, Hannum, PhenoAge, other developments. Age-related conditions linked to epigenetic age. Epigenetic age is 40% heritable.

- Field, Adam E., Neil A. Robertson, Tina Wang, Aaron Havas, Trey Ideker, and Peter D. Adams. “DNA Methylation Clocks in Aging: Categories, Causes, and Consequences.” Molecular Cell 71, no. 6 (September 2018): 882–95. https://doi.org/10.1016/j.molcel.2018.08.008. - Review of DNA methylation age clocks. The dynamic nature of DNA methylation. References to methylation clock studies in humans, mice, other organisms. Clocks derived from multiple tissue data. Poor overlap among studies. Latest specialized clocks like EpiTOC, PhenoAge,  Box 1 - penalized regression framework for epigenetic clock.


## Files and folders description

### data

Scripts to extract genomic coordinates of the age-associated genomic regions. The description of this folder contains the list of studies and the data available from them.

#### Mouse

- `Wang_2017.xlsx` - Mouse epigenetic clock, CpG sites and genes. From Wang, Tina, Brian Tsui, Jason F. Kreisberg, Neil A. Robertson, Andrew M. Gross, Michael Ku Yu, Hannah Carter, Holly M. Brown-Borg, Peter D. Adams, and Trey Ideker. “Epigenetic Aging Signatures in Mice Livers Are Slowed by Dwarfism, Calorie Restriction and Rapamycin Treatment.” Genome Biology 18, no. 1 (28 2017): 57. https://doi.org/10.1186/s13059-017-1186-2. [Source: Additional file 3](https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-017-1186-2/MediaObjects/13059_2017_1186_MOESM3_ESM.xlsx)



- `Illumina` - scripts to extract genomic coordinates of all CpG probes from Illumina Infinium 27K and Illumina Infinium 450K

- `GenomeRunner` - the results of the epigenomic enrichment analysis, and the scripts to process them

- `R.Aging` - R scripts for visualization of the results.

	- `genes_promoters_extract.Rmd` - extract promoters (+2000..-500bp arount TSS) of age-associated _genes_ (`data` folder).

