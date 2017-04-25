Scripts for:

* Annotating SAM files according to the genomic features the reads were simulated from
    * Stored in the ZZ:Z:... extra field
* Determining whether each alignment is correct or incorrect
    * Stored in the ZC:i:... extra field.  -1 means unaligned, 0 means aligned incorrectly, 1 means aligned correctly

For example, if `hapA_fnr_40.sam` is in the current directory:

    sh get_hg19_rmask.sh  # if you haven't already; gets repeatmasker file for GRCh37
    python rep.py --rep hg19.chr1.fa.out \
                  --sam-in hapA_fnr30_first100k.sam \
                  --exome-bed exome_chr1.bed \
                  --mason-source | python correctness.py

Because this uses `bx-tools` (for interval tree) and `pysam` (for SAM/BAM parsing and writing), you can't use `pypy`.

One fairly simply strategy for analyzing the results in R broken out by repeat family/subfamily is to do:

    echo -e "rep\texome\tcorrect" > hapA_fnr_40.rep.tsv
    python rep.py --rep hg19.chr1.fa.out \
                  --sam-in hapA_fnr30_first100k.sam \
                  --exome-bed exome_chr1.bed \
                  --mason-source | \
                  python correctness.py | \
                  tee hapA_fnr30_first100k.extra.sam | \
                  sed 's/.*ZZ:Z://' | \
                  sed 's/ZX:Z://' | \
                  sed 's/ZC:i://' >> hapA_fnr_40.rep.tsv

Then, in R do:

    tab <- read.table('hapA_fnr_40.rep.tsv', sep='\t', header=T)
    table(tab$V1, tab$V2)
