Makefiles for creating BED files of genomic coordinates for individual studies.
===

### `Alisch`

Alisch R. et. al. ["Age-associated DNA methylation in pediatric populations"](http://genome.cshlp.org/content/22/4/623.full)

Illumina 450K. BED coordinates: 460 up, 1,505 dn, 1,965 total.

- [Alisch_et_al_Sup_Table2.*](http://genome.cshlp.org/content/suppl/2012/02/01/gr.125187.111.DC1/Alisch_et_al_Sup_Table2.xls) - 2,078 regions correlated with age by methylation in peripheral blood; 479 of them showed positive correlation with age by methylation, 1,599 were found to be negatively correlated with age ("SSC.age.cor", 13th column).


### `Cory.age.summary`

- `AGE Atlas gender effects, part 2.csv` - Data from [AGE Atlas: gender effects, part 2](http://blog.corygil.es/age-atlas/gender-effects-part-2.html)


### `deMagalhaes`

De Magalhães, J. P., et.al. ["Meta-analysis of age-related gene expression profiles identifies common signatures of aging"](http://bioinformatics.oxfordjournals.org/content/25/7/875.short)

- [signatures_supplement.zip](http://genomics.senescence.info/gene_expression/signatures_supplement.zip) A Conserved Gene Expression Signature of Mammalian Aging on [Human Ageing Genomic Resources](http://genomics.senescence.info/gene_expression/signatures.html).

	- `supplementary_material.pdf` - description of the supplementary tables
	- `supplementary_tables*` - tables S3 and S4 contain 232 overespressed/146 underexpressed with age genes, consensus across all tissues. EntrezIDs and gene symbols, mixture of human and mouse IDs.
	- `genes_deMagalhaes_over/under.txt/bed` - gene symbols of 232 overespressed/146 underexpressed with age genes. BED files contain promoters (+2000..-500bp) of these genes, extracted with `R.Aging/genes_promoters_extract.Rmd`

### `Fernandez`

Illumina 450K. BED coordinates: 18,735 up, 45,407 dn, 64,142 total.

Supplementary data from Fernández AF, et.al. ["H3K4me1 marks DNA regions hypomethylated during aging in human stem and differentiated cells"](http://genome.cshlp.org/content/25/1/27/suppl/DC1) Genome Res 2015

- `make.sh` - download and pre-process all supplementary data.

- `Supplemental_TableLegends.txt` - headers of the supplemental tables

The most useful at the moment are Supplemental Tables 2-3 (MSCs CpGs) and 4-5 (blood CpGs), processed into BED files with `make`. But there's more, about twin age-methylation differences.

### `Florath`

Illumina 450K. BED coordinates: 119 up, 43 dn, 162 total.

Table 2 from Florath et al., [“Cross-Sectional and Longitudinal Changes in DNA Methylation with Age.”](http://hmg.oxfordjournals.org/content/23/5/1186.long), 162 differentially methylated CpG sites for age with Bonferroni-corrected statistical significance (P < 2.5 × 10-4)


### `genes/GenAge`

Data from [Human Ageing Genomic Resources](http://genomics.senescence.info/download.html)

- `[human_genes.zip](http://genomics.senescence.info/genes/human_genes.zip)` and `genage_human.csv` - 305 genes associated with age. Non-directional
- `genes_GenAge.bed` - promoters (+2000..-500bp) of these genes, extracted with `R.Aging/genes_promoters_extract.Rmd`


### `genes/JenAge

Data from [JenAge/AgeFactDB](http://agefactdb.jenage.de/cgi-bin/jaDB.cgi?VIEW=download). Use 'Ageing Factor Type': 'gene', 'Ageing Relevance Evidence': 'Any', 'File Type': 'TSV'


### `genes/LongevityMap`

Data from [Human Ageing Genomic Resources](http://genomics.senescence.info/download.html)

- `[longevity_genes.zip](http://genomics.senescence.info/longevity/longevity_genes.zip)` and `longevity.csv` - population-specific genes and rsids associated with longevity. 

### `genes/Peters`

Peters MJ, et. al. ["The transcriptional landscape of age in human peripheral blood"](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4639797). Nat Commun 2015

Supplementary Data 1 provides gene lists, statistics and, importantly, directionality of correlations. The directionality is clearly defined for the whole blood only. For other tissues, there are association statistics (p-values), with Z-scores not present for brain tissues.

Supplementary Data 2A lists tissue-specific genes not/associated with age as 0/1. Use directionality!

- `ncomms9570.pdf` - manuscript
- `ncomms9570-s1.pdf` - Supplementary Information: Supplementary Figures 1-10, Supplementary Tables 1-21, Supplementary Notes 1-2, Supplementary Acknowledgements and Supplementary References
- `ncomms9570-s2.xlsx` - Supplementary Data 1: Discovery, replication and generalization results. Y = yes, N = no; - is negative direction, + = positive direction; weight = sample size, P = p-value; WB = whole blood.
- `ncomms9570-s3.xlsx` - Supplementary Data 2: Number of genes that associate with age across all tissues. 2A: A gene-based overview for age-association across all tissues. 2B: A gene-based overview for significant expression across all tissues.
- `ncomms9570-s4.xlsx` - Supplementary Data 3: Up/down regulated gene set analyses in WEBGESTALT and GENENETWORK
- `ncomms9570-s5.xlsx` - Supplementary Data 4: Methylation results: Sobel test for methylation-based mediation of age-expression association
- `ncomms9570-s6.xlsx` - Supplementary Data 5: Transcriptomic age prediction formulas. 5A: Prediction formula weights for all participating cohorts (based on Leave-One-Out Meta-Analyses). 5B: The average "transcriptomic predictor formula" - the GENERAL predictor - for external cohorts

### `genes/Sood`

Sood et al., [“A Novel Multi-Tissue RNA Diagnostic of Healthy Ageing Relates to Cognitive Health Status.”](http://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0750-x)

- `13059_2015_750_MOESM1_ESM.xlsx` - 150 multi-tissue biological age classifier. [Source](https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-015-0750-x/MediaObjects/13059_2015_750_MOESM1_ESM.xlsx). An Excel spreadsheet containing data related to our study with six tabs: 1) analysis of the top 150 genes from age prototype classifier in PUBMED; 2) the top 670 genes from the first stage of the project; 3) phenodata for the training data set and validation data sets in our study; 4) a list of prior markers in the literature for Alzheimer’s disease; 5) positional gene enrichment analysis; 6) sample information for the BrainEac study


### `Hannum`

Hannum G. et. al. ["Genome-wide methylation profiles reveal quantitative views of human aging rates"](http://www.sciencedirect.com/science/article/pii/S1097276512008933).

Illumina 450K. BED coordinates: 89 all, 71 primary.

[Table S3. Aging Model Markers, Related to Figure 2.](http://www.sciencedirect.com/science/MiamiMultiMediaURL/1-s2.0-S1097276512008933/1-s2.0-S1097276512008933-mmc2.xlsx/272198/FULL/S1097276512008933/91cff6863693f2d294890d6fd28662e0/mmc2.xlsx) - A table of the methylation markers included in the primary aging model. The coefficient listed for each marker is its regression coefficient within the model. A second table is provided for the model based on all samples (primary and validation).

- `mmc2_Model_Primary_data.*` - The 71 age-associated methylation regions included in the primary aging model.
- `mmc2_Model_All_data.*` - The 89 age-associated methylation regions identified using all data. 60 of them are the same as those identified in primary data (see mmc2.png).

[Table S5. Aging Model Markers for TCGA Data, Related to Figure 4.](http://www.sciencedirect.com/science/MiamiMultiMediaURL/1-s2.0-S1097276512008933/1-s2.0-S1097276512008933-mmc3.xlsx/272198/FULL/S1097276512008933/c280262b4e861fc11051c187f3388ef1/mmc3.xlsx) - To investigate the similarities and differences between the tissues, we built an age model for breast, kidney, and lung tissues. The skin cohort did not have enough samples to build a model. The markers and coefficients of each model are listed here.

- `mmc3_Model_Breast/Kidney/Lung.bed` - tissue-specific age-associated markers. See mmc3 info about their overlap.

[Table S6. Genes Associated with Aging in Both the Methylome and the Transcriptome, Related to Figure 7.](http://www.sciencedirect.com/science/MiamiMultiMediaURL/1-s2.0-S1097276512008933/1-s2.0-S1097276512008933-mmc4.xlsx/272198/FULL/S1097276512008933/60c29e9569bef2044ecead91da601f47/mmc4.xlsx) - A list of genes that mapped to age-associated methylation markers and showed age-associated changes the transcriptome.

- `mmc4_Genes_S6.*` - 326 genes correlated with age by expression, whole blood, data not related to methylation data. From [Table S6. Genes Associated with Aging in Both the Methylome and the Transcriptome, Related to Figure 7.](http://www.sciencedirect.com/science/MiamiMultiMediaURL/1-s2.0-S1097276512008933/1-s2.0-S1097276512008933-mmc4.xlsx/272198/html/S1097276512008933/a609fb76dfb564b4e3f57887c2df17dc/mmc4.xlsx)

[Table S7. Transcriptome Aging Model, Related to Figure 7.](http://www.sciencedirect.com/science/MiamiMultiMediaURL/1-s2.0-S1097276512008933/1-s2.0-S1097276512008933-mmc5.xlsx/272198/FULL/S1097276512008933/bdeaadf32df47e709bf1e98e9f7405e2/mmc5.xlsx)

 - `mmc5_Genes_all/neg/pos_S7.*` - 55 genes predictive of age based by their expression, whole blood, data not related to methylation data. 24 of these genes showed positive correlation with age by expression, 31 were found to be negatively correlated with age. From [Table S7. Transcriptome Aging Model, Related to Figure 7.](http://www.sciencedirect.com/science/MiamiMultiMediaURL/1-s2.0-S1097276512008933/1-s2.0-S1097276512008933-mmc5.xlsx/272198/html/S1097276512008933/a0cf53be370e97f7d3fabcb80fdb74ec/mmc5.xlsx)


### `Heyn`

Illumina 450K. BED coordinates: 1,219 up, 1,986 dn, 3,205 total.

Heyn, H., et.al. ["Distinct DNA methylomes of newborns and centenarians"](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3387108/)

- [Dataset_S02 (XLSX)](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3387108/bin/1120658109_sd02.xlsx) - 3,205 age-associated methylation markers differentially methylated in cord blood of newborns and CD4+ T cells of centenarians; 1,219 of them showed hypermethylation with age, 1,986 were found to be hypomethylated with age


### `Horvath`

Illumina 450K. BED coordinates: 156 up, 196 dn, 352 total.

Classical epigenetic age clock predictor from Horvath S: "[DNA methylation age of human tissues and cell types](http://genomebiology.biomedcentral.com/articles/10.1186/gb-2013-14-10-r115#MOESM1)". Genome Biol 2013

- "[13059_2013_3156_MOESM3_ESM.csv](https://static-content.springer.com/esm/art%3A10.1186%2Fgb-2013-14-10-r115/MediaObjects/13059_2013_3156_MOESM3_ESM.csv)" - Additional file 3: Coefficient values for the DNAm age predictor. This Excel file provides detailed information on the multi-tissue age predictor defined using the training set data. The multi-tissue age predictor uses 353 CpGs, of which 193 and 160 have positive and negative correlations with age, respectively. The table also represents the coefficient values for the shrunken age predictor that is based on a subset of 110 CpGs (a subset of the 353 CpGs).

- `make` - extracts genomic coordinates of CpG cites positively/negatively correlating with age as BED files.

Software at [DNA methylation age and the epigenetic clock](https://labs.genetics.ucla.edu/horvath/dnamage/)


### `Marttila`

Illumina 450K. BED coordinates: 3,925 up, 4,615 dn, 8,540 total.

Data from Marttila S, et.al. "[Ageing-associated changes in the human DNA methylome: genomic locations and effects on gene expression](http://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-015-1381-z)". BMC Genomics 2015

- `12864_2015_1381_MOESM1_ESM.xlsx` - data from [Additional file 1:](https://static-content.springer.com/esm/art%3A10.1186%2Fs12864-015-1381-z/MediaObjects/12864_2015_1381_MOESM1_ESM.xlsx). Leukocyte-specific age-associated CpGs, 3,925 CpG probes positively, and 4,615 probes negatively correlated with age (8,540 total).


### `Rakyan`

Illumina 27K. BED coordinates: 123 total

Rakyan V. et. al. ["Human aging-associated DNA hypermethylation occurs preferentially at bivalent chromatin domains"](http://genome.cshlp.org/content/early/2010/03/09/gr.103101.109)

- [Supp_Table_3.*](http://genome.cshlp.org/content/suppl/2010/03/11/gr.103101.109.DC1/Supp_Table_3.xls) - 131 regions correlated with age by methylation in whole blood, and show the same directional age-associated DNA methylation change in CD4+ T-cells and CD14+ monocytes.


### `Reynolds`

Illumina 450K. BED coordinates: Monocytes - 355 up, 1,439 dn, 1,794 total; T-cells - 2,049 up, 546 dn, 2,595 total.

Reynolds, L. M. et. al. ["Age-related variations in the methylome associated with gene expression in human monocytes and T cells."](http://www.nature.com/ncomms/2014/141118/ncomms6366/full/ncomms6366.html#supplementary-information)

- [Supplementary Data 1](http://www.nature.com/ncomms/2014/141118/ncomms6366/extref/ncomms6366-s2.xlsx) - 2,595 regions correlated with age by methylation in CD4+ T-cells. 2,049 of them showed positive correlation with age by methylation, 546 were found to be negatively correlated with age.

- [Supplementary Data 2](http://www.nature.com/ncomms/2014/141118/ncomms6366/extref/ncomms6366-s3.xlsx) - 2,259 regions correlated with age by methylation on CD14+ monocytes. 468 of them showed positive correlation with age by methylation, 1,791 were found to be negatively correlated with age.


### `SNPs`

Age-associated SNPs

Pilling et al., [“Human Longevity Is Influenced by Many Genetic Variants.”](http://www.impactaging.com/papers/v8/n3/full/100930.html)

- `SupTable1.*` - Lists the genetic variants included for each genetic risk score. SNPs associated with different diseases that have been tested for association with longevity
- `SupTable2.*` - Includes the 1,000 most-associated variants from each of the four GWAS performed

### `Steegenga`

Illumina 450K. BED coordinates: Age-CpG-gene expression - 726; Age-CpG-not gene expression - 4,554; CpG meta analysis - 7,401. 

Steegenga, W. T., et. al. ["Genome-wide age-related changes in DNA methylation and gene expression in human PBMCs."](http://link.springer.com/content/pdf/10.1007%2Fs11357-014-9648-x.pdf)

- [Supplementary Data](http://static-content.springer.com/esm/art%3A10.1007%2Fs11357-014-9648-x/MediaObjects/11357_2014_9648_MOESM2_ESM.xls) contains several supplementary tables:

Supplemental table 2S: 726 methylation markers correlated with age- and gene expression changes in peripheral blood (cpg_and_gene)

Supplemental table 3S: 4,552 methylation markers correlated with age- but not with gene expression changes in peripheral blood (cpg_not_gene)

Supplemental table 4S: 7,477 age-associated methylation markers identified in multiple studies (meta)

### `Weidner`

Illumina 450K. BED coordinates: 42 up, 57 dn, 99 total.

Additional file 2 from Weidner et al., [“Aging of Blood Can Be Tracked by DNA Methylation Changes at Just Three CpG Sites.”](http://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-2-r24) - Beta values for 102 AR-GpGs from 575 samples
