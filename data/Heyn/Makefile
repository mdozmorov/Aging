all:	Heyn_markers_pos.bed Heyn_markers_neg.bed

Heyn_markers_pos.bed:	sd02.txt
	sed '1,4d' $< | awk '$$4 == "Hypermethylated" {print $$5}' | sort | uniq > tmp.txt && \
	join -t "	" -1 1 -2 4 tmp.txt ../../Illumina/Illumina450K.bed | awk 'BEGIN {OFS="\t"} {print $$2,$$3,$$4,$$1,$$5,$$6}' > $@ && \
	rm tmp.txt

Heyn_markers_neg.bed:	sd02.txt
	sed '1,4d' $< | awk '$$4 == "Hypomethylated" {print $$5}' | sort | uniq > tmp.txt && \
	join -t "	" -1 1 -2 4 tmp.txt ../../Illumina/Illumina450K.bed | awk 'BEGIN {OFS="\t"} {print $$2,$$3,$$4,$$1,$$5,$$6}' > $@ && \
	rm tmp.txt