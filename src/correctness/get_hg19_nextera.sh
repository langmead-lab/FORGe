#!/bin/sh

if [ ! -f "nexterarapidcapture_exome_targetedregions_v1.2.bed" ] ; then
    wget http://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/samplepreps_nextera/nexterarapidcapture/nexterarapidcapture_exome_targetedregions_v1.2.bed
fi

if [ ! -f "exome_chr9.bed" ] ; then
    cat nexterarapidcapture_exome_targetedregions_v1.2.bed | awk '$1 == "chr9"' >  exome_chr9.bed
fi
