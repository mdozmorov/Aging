Makefiles for creating BED files of genomic coordinates for individual studies.
===

### `Alisch`

Alisch R. et. al. ["Age-associated DNA methylation in pediatric populations"](http://genome.cshlp.org/content/22/4/623.full)

Illumina 450K

- [Alisch_et_al_Sup_Table2.*](http://genome.cshlp.org/content/suppl/2012/02/01/gr.125187.111.DC1/Alisch_et_al_Sup_Table2.xls) - 2,078 regions correlated with age by methylation in peripheral blood; 479 of them showed positive correlation with age by methylation, 1,599 were found to be negatively correlated with age ("SSC.age.cor", 13th column). 


### `Cory.age.summary`

- `AGE Atlas gender effects, part 2.csv` - Data from [AGE Atlas: gender effects, part 2](http://blog.corygil.es/age-atlas/gender-effects-part-2.html)


### `deMagalhaes`

De Magalhães, J. P., et.al. ["Meta-analysis of age-related gene expression profiles identifies common signatures of aging"](http://bioinformatics.oxfordjournals.org/content/25/7/875.short)

- [signatures_supplement.zip](http://genomics.senescence.info/gene_expression/signatures_supplement.zip) A Conserved Gene Expression Signature of Mammalian Aging on [Human Ageing Genomic Resources](http://genomics.senescence.info/gene_expression/signatures.html).

	- `supplementary_material.pdf` - description of the supplementary tables
	- `supplementary_tables*` - tables S3 and S4 contain 232 overespressed/146 underexpressed with age genes, consensus across all tissues. EntrezIDs and gene symbols, mixture of human and mouse IDs.
	- `genes_deMagalhaes_over/under.txt/bed` - gene symbols of 232 overespressed/146 underexpressed with age genes. BED files contain promoters (+2000..-500bp) of these genes, extracted with `R.Aging/genes_promoters_extract.Rmd`


### `genes/GenAge`

Data from [Human Ageing Genomic Resources](http://genomics.senescence.info/download.html)

- `[human_genes.zip](http://genomics.senescence.info/genes/human_genes.zip)` and `genage_human.csv` - 305 genes associated with age. Non-directional
- `genes_GenAge.bed` - promoters (+2000..-500bp) of these genes, extracted with `R.Aging/genes_promoters_extract.Rmd`


### `genes/JenAge

Data from [JenAge/AgeFactDB](http://agefactdb.jenage.de/cgi-bin/jaDB.cgi?VIEW=download). Use 'Ageing Factor Type': 'gene', 'Ageing Relevance Evidence': 'Any', 'File Type': 'TSV'


### `genes/LongevityMap`

Data from [Human Ageing Genomic Resources](http://genomics.senescence.info/download.html)

- `[longevity_genes.zip](http://genomics.senescence.info/longevity/longevity_genes.zip)` and `longevity.csv` - population-specific genes and rsids associated with longevity. 


### `Hannum`

Hannum G. et. al. ["Genome-wide methylation profiles reveal quantitative views of human aging rates"](http://www.sciencedirect.com/science/article/pii/S1097276512008933).

Illumina 450K

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

Heyn, H., et.al. ["Distinct DNA methylomes of newborns and centenarians"](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3387108/)

Illumina 450K

- [Dataset_S02 (XLSX)](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3387108/bin/1120658109_sd02.xlsx) - 3,205 age-associated methylation markers differentially methylated in cord blood of newborns and CD4+ T cells of centenarians; 1,219 of them showed hypermethylation with age, 1,986 were found to be hypomethylated with age


### `Peters`

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

### `Rakyan`

Rakyan V. et. al. ["Human aging-associated DNA hypermethylation occurs preferentially at bivalent chromatin domains"](http://genome.cshlp.org/content/early/2010/03/09/gr.103101.109)

Illumina 27K

- [Supp_Table_3.*](http://genome.cshlp.org/content/suppl/2010/03/11/gr.103101.109.DC1/Supp_Table_3.xls) - 131 regions correlated with age by methylation in whole blood, and show the same directional age-associated DNA methylation change in CD4+ T-cells and CD14+ monocytes.


### `Reynolds`

Reynolds, L. M. et. al. ["Age-related variations in the methylome associated with gene expression in human monocytes and T cells."](http://www.nature.com/ncomms/2014/141118/ncomms6366/full/ncomms6366.html#supplementary-information)

Illumina 450K

- [Supplementary Data 1](http://www.nature.com/ncomms/2014/141118/ncomms6366/extref/ncomms6366-s2.xlsx) - 2,595 regions correlated with age by methylation in CD4+ T-cells. 2,049 of them showed positive correlation with age by methylation, 546 were found to be negatively correlated with age.

- [Supplementary Data 2](http://www.nature.com/ncomms/2014/141118/ncomms6366/extref/ncomms6366-s3.xlsx) - 2,259 regions correlated with age by methylation on CD14+ monocytes. 468 of them showed positive correlation with age by methylation, 1,791 were found to be negatively correlated with age.


### `Steegenga`

Steegenga, W. T., et. al. ["Genome-wide age-related changes in DNA methylation and gene expression in human PBMCs."](http://link.springer.com/content/pdf/10.1007%2Fs11357-014-9648-x.pdf)

Illumina 450K

- [Supplementary Data](http://static-content.springer.com/esm/art%3A10.1007%2Fs11357-014-9648-x/MediaObjects/11357_2014_9648_MOESM2_ESM.xls) contains several supplementary tables:

Supplemental table 2S: 726 methylation markers correlated with age- and gene expression changes in peripheral blood (cpg_and_gene)

Supplemental table 3S: 4,552 methylation markers correlated with age- but not with gene expression changes in peripheral blood (cpg_not_gene)

Supplemental table 4S: 7,477 age-associated methylation markers identified in multiple studies (meta)

