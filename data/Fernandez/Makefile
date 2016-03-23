all:	Fernandez_MSCs_cpg_pos.bed Fernandez_MSCs_cpg_neg.bed Fernandez_blood_cpg_pos.bed Fernandez_blood_cpg_neg.bed


Fernandez_MSCs_cpg_pos.bed:	Supplemental_Table_2.txt
	sed '1,2d' $< | cut -f1 | grep "^cg" | sort | uniq > tmp.txt
	join -t "	" -1 1 -2 4 tmp.txt ../../Illumina/Illumina450K.bed | awk 'BEGIN {OFS="\t"} {print $$2,$$3,$$4,$$1,$$5,$$6}' > $@ && \
	echo "chr2	1226245	1226368	ch.3.1226245F" >> $@ ; \
	rm tmp.txt

Fernandez_MSCs_cpg_neg.bed:	Supplemental_Table_3.txt
	sed '1,2d' $< | cut -f1 | grep "^cg" | sort | uniq > tmp.txt
	join -t "	" -1 1 -2 4 tmp.txt ../../Illumina/Illumina450K.bed | awk 'BEGIN {OFS="\t"} {print $$2,$$3,$$4,$$1,$$5,$$6}' > $@ && \
	rm tmp.txt

Fernandez_blood_cpg_pos.bed:	Supplemental_Table_4.txt
	sed '1,2d' $< | cut -f1 | grep "^cg" | sort | uniq > tmp.txt
	join -t "	" -1 1 -2 4 tmp.txt ../../Illumina/Illumina450K.bed | awk 'BEGIN {OFS="\t"} {print $$2,$$3,$$4,$$1,$$5,$$6}' > $@ && \
	sed '1,2d' $< | cut -f1 | grep -v "^cg" | tr "." "\t" | tr -d "F|R" | awk '{OFS="\t"} {print $$1"r"$$2, $$3, $$3+123}' >> $@ ; \
	rm tmp.txt

Fernandez_blood_cpg_neg.bed:	Supplemental_Table_5.txt
	sed '1,2d' $< | cut -f1 | grep "^cg" | sort | uniq > tmp.txt
	join -t "	" -1 1 -2 4 tmp.txt ../../Illumina/Illumina450K.bed | awk 'BEGIN {OFS="\t"} {print $$2,$$3,$$4,$$1,$$5,$$6}' > $@ && \
	sed '1,2d' $< | cut -f1 | grep -v "^cg" | tr "." "\t" | tr -d "F|R" | awk '{OFS="\t"} {print $$1"r"$$2, $$3, $$3+123}' >> $@ ; \
	rm tmp.txt


