#!/bin/bash


#create the project on the desktop
python geneysis.py create --wd ~/Desktop/geneysis

#load all the phages
for GENBANK in $( ls genbank/*.gb ); do
    python geneysis.py import --wd ~/Desktop/geneysis --file $GENBANK
done

#create fasta
python geneysis.py create_fasta --wd ~/Desktop/geneysis

#create blast database and blast all proteins
python geneysis.py blast --wd ~/Desktop/geneysis

#run clustalo
python geneysis.py clustalo --wd ~/Desktop/geneysis

#cluster genes
python geneysis.py cluster --wd ~/Desktop/geneysis