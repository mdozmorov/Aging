Analysis of age-related genomic regions
----------------------------------------

Hannum G. et. al. paper ["Genome-wide methylation profiles reveal quantitative views of human aging rates"](http://www.sciencedirect.com/science/article/pii/S1097276512008933).
===

[Table S3. Aging Model Markers, Related to Figure 2.](http://www.sciencedirect.com/science/MiamiMultiMediaURL/1-s2.0-S1097276512008933/1-s2.0-S1097276512008933-mmc2.xlsx/272198/FULL/S1097276512008933/91cff6863693f2d294890d6fd28662e0/mmc2.xlsx)

A table of the methylation markers included in the primary aging model. The coefficient listed for each marker is its regression coefficient within the model. A second table is provided for the model based on all samples (primary and validation).

    - `mmc2_Model_Primary_data.*` - The 71 age-associated methylation regions included in the primary aging model.
    - `mmc2_Model_All_data.*` - The 89 age-associated methylation regions identified using all data. 60 of them are the same as those identified in primary data (see mmc2.png).

[Table S5. Aging Model Markers for TCGA Data, Related to Figure 4.](http://www.sciencedirect.com/science/MiamiMultiMediaURL/1-s2.0-S1097276512008933/1-s2.0-S1097276512008933-mmc3.xlsx/272198/FULL/S1097276512008933/c280262b4e861fc11051c187f3388ef1/mmc3.xlsx)

To investigate the similarities and differences between the tissues, we built an age model for breast, kidney, and lung tissues. The skin cohort did not have enough samples to build a model. The markers and coefficients of each model are listed here.

    - `mmc3_Model_Breadt/Kidney/Lung.bed` - tissue-specific age-associated markers. See mmc3 info about their overlap.

[Table S6. Genes Associated with Aging in Both the Methylome and the Transcriptome, Related to Figure 7.](http://www.sciencedirect.com/science/MiamiMultiMediaURL/1-s2.0-S1097276512008933/1-s2.0-S1097276512008933-mmc4.xlsx/272198/FULL/S1097276512008933/60c29e9569bef2044ecead91da601f47/mmc4.xlsx)

A list of genes that mapped to age-associated methylation markers and showed age-associated changes the transcriptome.

    - `mmc4_Genes_S6.*` - 326 genes correlated with age by expression 

[Table S7. Transcriptome Aging Model, Related to Figure 7.](http://www.sciencedirect.com/science/MiamiMultiMediaURL/1-s2.0-S1097276512008933/1-s2.0-S1097276512008933-mmc5.xlsx/272198/FULL/S1097276512008933/bdeaadf32df47e709bf1e98e9f7405e2/mmc5.xlsx)

The list of genes and coefficients used for predicting age based on transcriptome data.

    - `mmc5_Genes_all/neg/pos_S7.*` 54 genes correlated with age by expression and by the presence of age-related DRMs. Split in negatively/positively correlated

Use 450K for background

Alisch R. et. al. paper ["Age-associated DNA methylation in pediatric populations"](http://genome.cshlp.org/content/22/4/623.full)
===
- [Alisch_et_al_Sup_Table2.*](http://genome.cshlp.org/content/suppl/2012/02/01/gr.125187.111.DC1/Alisch_et_al_Sup_Table2.xls) - a total of 2078 age-associated DMRs, 1601 demethylated, 477 methylated ("SSC.age.cor", 13th column). 

Use 450K for background

Rakyan V. et. al. paper ["Human aging-associated DNA hypermethylation occurs preferentially at bivalent chromatin domains"](http://genome.cshlp.org/content/early/2010/03/09/gr.103101.109)
===
- [Supp_Table_3.*](http://genome.cshlp.org/content/suppl/2010/03/11/gr.103101.109.DC1/Supp_Table_3.xls) - List of statistically significant whole blood hyper-aDMRs that show the same directional age-associated DNA methylation change in CD4+ and CD14+ cells.

130 age-related hypermethylated regions in CD4+ and CD14+ cells. Use 27K for background

Reynolds, L. M. et. al. paper ["Age-related variations in the methylome associated with gene expression in human monocytes and T cells."](http://www.nature.com/ncomms/2014/141118/ncomms6366/full/ncomms6366.html#supplementary-information)
===
- [Supplementary Data 1](http://www.nature.com/ncomms/2014/141118/ncomms6366/extref/ncomms6366-s2.xlsx) - CpG sites with age-associated methylation (age-dMS) detected in 227 T-cell samples (FDR<0.001)

- [Supplementary Data 2](http://www.nature.com/ncomms/2014/141118/ncomms6366/extref/ncomms6366-s3.xlsx) - Age and cis-gene expression associated methylation sites (age-eMS) detected in 1,264 monocyte samples (FDR<0.001)

Coordinates are in hg19. Manually extract them from .xls files converting chr17:27369892 into chr17   27369892   27369893.

Miscellaneous
--------------
`R.Aging` - Scripts to create .BED files for Illumina Infinium 27K and 450K arrays. Also, extract genomic coordinates of subsets of probes. Note that point coordinates are outputted. Ideally, one need to convert them to 50bp probes depending on strand. Currently, simply add +/-50bp flanking regions around each site.

`GenomeRunner` - Results of GenomeRunner analysis

`make` - Misc. helpful commands