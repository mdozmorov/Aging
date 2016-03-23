all:	Reynolds_cpg_tcells_all.bed Reynolds_cpg_monocytes_all.bed Reynolds_cpg_tcells_pos.bed Reynolds_cpg_monocytes_pos.bed Reynolds_cpg_tcells_neg.bed Reynolds_cpg_monocytes_neg.bed

SHELL=/bin/bash

Reynolds_cpg_tcells_all.bed:	ncomms6366-s2.txt
	cat $< | sed '1,4d' | cut -f1-2 | tr ":" "\t" | awk 'BEGIN {OFS="\t"} {print $$2,$$3,$$3+123,$$1}' > $@

Reynolds_cpg_monocytes_all.bed:	ncomms6366-s3.txt
	cat $< | sed '1,4d' | cut -f1-2 | tr ":" "\t" | awk 'BEGIN {OFS="\t"} {print $$2,$$3,$$3+123,$$1}' > $@

Reynolds_cpg_tcells_pos.bed:	ncomms6366-s2.txt
	cat $< | sed '1,4d' | cut -f1-2,6 | awk '{if ( $$3 > 0 ) { print $0 } }' | tr ":" "\t" | awk 'BEGIN {OFS="\t"} {print $$2,$$3,$$3+123,$$1}' > $@

Reynolds_cpg_monocytes_pos.bed:	ncomms6366-s3.txt
	cat $< | sed '1,4d' | cut -f1-2,9 | awk '{if ( $$3 > 0 ) { print $0 } }' | tr ":" "\t" | awk 'BEGIN {OFS="\t"} {print $$2,$$3,$$3+123,$$1}' > $@

Reynolds_cpg_tcells_neg.bed:	ncomms6366-s2.txt
	cat $< | sed '1,4d' | cut -f1-2,6 | awk '{if ( $$3 < 0 ) { print $0 } }' | tr ":" "\t" | awk 'BEGIN {OFS="\t"} {print $$2,$$3,$$3+123,$$1}' > $@

Reynolds_cpg_monocytes_neg.bed:	ncomms6366-s3.txt
	cat $< | sed '1,4d' | cut -f1-2,9 | awk '{if ( $$3 < 0 ) { print $0 } }' | tr ":" "\t" | awk 'BEGIN {OFS="\t"} {print $$2,$$3,$$3+123,$$1}' > $@
