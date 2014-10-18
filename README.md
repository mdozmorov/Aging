Analysis of age-related genomic regions
----------------------------------------

##Hannum G. et. al. paper ["Genome-wide methylation profiles reveal quantitative views of human aging rates"](http://www.sciencedirect.com/science/article/pii/S1097276512008933).**

- `mmc2_Model_Primary_data.*` - The 71 age-associated methylation regions included in the primary aging model.
- `mmc2_Model_All_data.*` - The 89 age-associated methylation regions identified using all data. 60 of them are the same as those identified in primary data (see mmc2.png).
- `mmc3_Model_Breadt/Kidney/Lung.bed` - tissue-specific age-associated markers. See mmc3 info about their overlap.
- `mmc4_Genes_S6.*` - 334 genes correlated with age by expression 
- `mmc5_Genes_all/neg/pos_S7.*` 54 genes correlated with age by expression and by the presence of age-related DRMs. Split in negatively/positively correlated

Use 450K for background

##Alisch R. et. al. paper ["Age-associated DNA methylation in pediatric populations"](http://genome.cshlp.org/content/22/4/623.full)**

- `Alisch_et_al_Sup_Table2.*` - a total of 2078 age-associated DMRs, 1601 demethylated, 477 methylated ("SSC.age.cor" column). Use 450K for background

##Rakyan V. et. al. paper ["Human aging-associated DNA hypermethylation occurs preferentially at bivalent chromatin domains"](http://genome.cshlp.org/content/early/2010/03/09/gr.103101.109)

- `Supp_Table_3.*` - 130 age-related hypermethylated regions in CD4+ and CD14+ cells. Use 27K for background

Miscellaneous
--------------
`R.Aging` - Scripts to create .BED files for Illumina Infinium 27K and 450K arrays. Also, extract genomic coordinates of subsets of probes. Note that point coordinates are outputted. Ideally, one need to convert them to 50bp probes depending on strand. Currently, simply add +/-50bp flanking regions around each site.

`GenomeRunner` - Results of GenomeRunner analysis

`make` - Misc. helpful commands