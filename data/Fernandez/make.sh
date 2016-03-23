# Download all 27 XLSX supplemental tables
for file in `cat Supplementary_URLs.txt`; do wget $file; done
# Convert XLSX to txt, and remove XLSX
for file in *.xlsx; do echo $file; Rscript make_xlsx2txt.R $file; done
# Extract table headers
for file in `ls Supplemental_Table_*.txt | sort`; do echo $file; head -n 1 $file | sed 's/NA//g' >> Supplemental_TableLegends.txt; done




