#!/usr/bin/env python
#
# Given a list of gene names in a file, get genomic coordinates of the promoters
# A promoter is defined as a region 2000bp upstream of transcription start site
# Strand specificity matters. HG19 is currently hard-coded
#
# Example: pyton refgene.py txt/all.refGene.txt | bedtools sort | mergeBed -s -nms -i > all.refGene.bed
# Another output is "notfound" file, containing genes that weren't found in the database
#
import csv, sys, mysql.connector, argparse, os

db = mysql.connector.connect(user="genome", 
    password="", database="hg19", host="genome-mysql.cse.ucsc.edu")
c = db.cursor()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("snp_file", nargs=1)
    args = parser.parse_args()

    genelist = [row.strip() for row in open(args.snp_file[0])]
    
    goodchrom = "'chr1','chr2','chr3','chr4','chr5','chr6','chr7','chrX','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr20','chrY','chr19','chr22','chr21','chrM'";
    q = """SELECT chrom, txStart, txEnd, name2, strand from refGene WHERE name2 IN (%s) AND chrom IN (%s);""" % ("'"+"','".join(genelist)+"'", goodchrom)
    c.execute(q)
    d=c.fetchall()
    genesout = []
    for row in d:
    	genesout.append(row[3]) # Collect genes that were obtained from the database
    	if row[4] == "+":
    		print "\t".join([row[0], str(row[1]-2000), str(row[1]), row[3], "0", row[4]])
    	else:
    		print "\t".join([row[0], str(row[2]-2000), str(row[2]), row[3], "0", row[4]])
    genesdiff = list(set(genelist) - set(genesout)) # Get the difference between input and output gene lists
    h2 = open(args.snp_file[0][:-4]+".notfound.txt", "w") # Save them
    h2.write("\n".join(genesdiff))
    h2.close()
db.close()
