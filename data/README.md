Makefiles for creating BED files of genomic coordinates for individual studies.
===

Hannum G. et. al. ["Genome-wide methylation profiles reveal quantitative views of human aging rates"](http://www.sciencedirect.com/science/article/pii/S1097276512008933).
---

Illumina 450K

[Table S3. Aging Model Markers, Related to Figure 2.](http://www.sciencedirect.com/science/MiamiMultiMediaURL/1-s2.0-S1097276512008933/1-s2.0-S1097276512008933-mmc2.xlsx/272198/FULL/S1097276512008933/91cff6863693f2d294890d6fd28662e0/mmc2.xlsx)

A table of the methylation markers included in the primary aging model. The coefficient listed for each marker is its regression coefficient within the model. A second table is provided for the model based on all samples (primary and validation).

- mmc2_Model_Primary_data.* - The 71 age-associated methylation regions included in the primary aging model.
- mmc2_Model_All_data.* - The 89 age-associated methylation regions identified using all data. 60 of them are the same as those identified in primary data (see mmc2.png).

[Table S5. Aging Model Markers for TCGA Data, Related to Figure 4.](http://www.sciencedirect.com/science/MiamiMultiMediaURL/1-s2.0-S1097276512008933/1-s2.0-S1097276512008933-mmc3.xlsx/272198/FULL/S1097276512008933/c280262b4e861fc11051c187f3388ef1/mmc3.xlsx)

To investigate the similarities and differences between the tissues, we built an age model for breast, kidney, and lung tissues. The skin cohort did not have enough samples to build a model. The markers and coefficients of each model are listed here.

- mmc3_Model_Breast/Kidney/Lung.bed - tissue-specific age-associated markers. See mmc3 info about their overlap.

[Table S6. Genes Associated with Aging in Both the Methylome and the Transcriptome, Related to Figure 7.](http://www.sciencedirect.com/science/MiamiMultiMediaURL/1-s2.0-S1097276512008933/1-s2.0-S1097276512008933-mmc4.xlsx/272198/FULL/S1097276512008933/60c29e9569bef2044ecead91da601f47/mmc4.xlsx)

A list of genes that mapped to age-associated methylation markers and showed age-associated changes the transcriptome.

- mmc4_Genes_S6.* - 326 genes correlated with age by expression, whole blood, data not related to methylation data 

[Table S7. Transcriptome Aging Model, Related to Figure 7.](http://www.sciencedirect.com/science/MiamiMultiMediaURL/1-s2.0-S1097276512008933/1-s2.0-S1097276512008933-mmc5.xlsx/272198/FULL/S1097276512008933/bdeaadf32df47e709bf1e98e9f7405e2/mmc5.xlsx)

 - mmc5_Genes_all/neg/pos_S7.* - 55 genes predictive of age based by their expression, whole blood, data not related to methylation data. 24 of these genes showed positive correlation with age by expression, 31 were found to be negatively correlated with age


Alisch R. et. al. ["Age-associated DNA methylation in pediatric populations"](http://genome.cshlp.org/content/22/4/623.full)
---

Illumina 450K

- [Alisch_et_al_Sup_Table2.*](http://genome.cshlp.org/content/suppl/2012/02/01/gr.125187.111.DC1/Alisch_et_al_Sup_Table2.xls) - 2,078 regions correlated with age by methylation in peripheral blood; 479 of them showed positive correlation with age by methylation, 1,599 were found to be negatively correlated with age ("SSC.age.cor", 13th column). 


Rakyan V. et. al. ["Human aging-associated DNA hypermethylation occurs preferentially at bivalent chromatin domains"](http://genome.cshlp.org/content/early/2010/03/09/gr.103101.109)
---

Illumina 27K

- [Supp_Table_3.*](http://genome.cshlp.org/content/suppl/2010/03/11/gr.103101.109.DC1/Supp_Table_3.xls) - 131 regions correlated with age by methylation in whole blood, and show the same directional age-associated DNA methylation change in CD4+ T-cells and CD14+ monocytes.


Reynolds, L. M. et. al. ["Age-related variations in the methylome associated with gene expression in human monocytes and T cells."](http://www.nature.com/ncomms/2014/141118/ncomms6366/full/ncomms6366.html#supplementary-information)
---

Illumina 450K

- [Supplementary Data 1](http://www.nature.com/ncomms/2014/141118/ncomms6366/extref/ncomms6366-s2.xlsx) - 2,595 regions correlated with age by methylation in CD4+ T-cells. 2,049 of them showed positive correlation with age by methylation, 546 were found to be negatively correlated with age.

- [Supplementary Data 2](http://www.nature.com/ncomms/2014/141118/ncomms6366/extref/ncomms6366-s3.xlsx) - 2,259 regions correlated with age by methylation on CD14+ monocytes. 468 of them showed positive correlation with age by methylation, 1,791 were found to be negatively correlated with age.


Steegenga, W. T., et. al. ["Genome-wide age-related changes in DNA methylation and gene expression in human PBMCs."](http://link.springer.com/content/pdf/10.1007%2Fs11357-014-9648-x.pdf)
---

Illumina 450K

- [Supplementary Data](http://static-content.springer.com/esm/art%3A10.1007%2Fs11357-014-9648-x/MediaObjects/11357_2014_9648_MOESM2_ESM.xls) contains several supplementary tables:

1) Supplemental table 2S: 726 methylation markers correlated with age- and gene expression changes in peripheral blood

2) Supplemental table 3S: 4,552 methylation markers correlated with age- but not with gene expression changes in peripheral blood

3) Supplemental table 4S: 7,477 age-associated methylation markers identified in multiple studies


Heyn, H., et.al. ["Distinct DNA methylomes of newborns and centenarians"](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3387108/)
---

Illumina 450K

- [Dataset_S02 (XLSX)](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3387108/bin/1120658109_sd02.xlsx) - 3,205 age-associated methylation markers differentially methylated in cord blood of newborns and CD4+ T cells of centenarians; 1,219 of them showed hypermethylation with age, 1,986 were found to be hypomethylated with age

De Magalhães, J. P., et.al. ["Meta-analysis of age-related gene expression profiles identifies common signatures of aging"](http://bioinformatics.oxfordjournals.org/content/25/7/875.short)

- [signatures_supplement.zip](http://genomics.senescence.info/gene_expression/signatures_supplement.zip) A Conserved Gene Expression Signature of Mammalian Aging on [Human Ageing Genomic Resources](http://genomics.senescence.info/gene_expression/signatures.html) - 73 age-associated genes. 56 genes overexpressed with age; 17 genes underexpressed with age.