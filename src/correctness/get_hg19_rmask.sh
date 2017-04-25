#!/bin/sh

if [ ! -f "hg19.fa.out.gz" ] ; then
    wget http://repeatmasker.org/genomes/hg19/RepeatMasker-rm330-db20120124/hg19.fa.out.gz
fi

if [ ! -f "hg19.chr1.fa.out" ] ; then
    gzip -dc hg19.fa.out.gz | awk '$5 == "chr1"' | sed 's/ chr1 / 1 /' > hg19.chr1.fa.out
fi
