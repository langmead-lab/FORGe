# FORGe tool for ranking variants and building an optimal graph genome

FORGe consists of 2 primary scripts -- 'rank.py' and 'build.py' -- as well as a helper script 'vcf_to_1ksnp.py' for generating input files in the proper format.

## vcf_to_1ksnp.py ##

Rather than reading variant from a VCF file, FORGe takes a input a 1ksnp file containing variant information and an optional phasing file. The 1ksnp file format simplifies a complex VCF file. Each alternate allele is stored as a separate line with the following columns:
Chromosome	Position	Reference Allele	Alternate Allele	Population Frequency	??	# Alternates	Variant Name

The phasing file should store a comma-separated table with a column for each individual and row for each variant, in the same order they appear in the 1ksnp file. Each cell contains the allele that the given individual contains for the given variant, with 0 indicating the reference allele and 'k' indicating the kth alternate allele, in order, appearing in the 1ksnp file.

FORGe includes a helper script 'vcf_to_1ksnp.py' to facilitate convertsion from a VCF file and set of ingroup individuals to FORGe variant and phasing files. Ingroup individuals can be specified as a list of names to either include or exclude (the first should use the --ingroup argument, the second the --outgroup argument). Sample usage of this script would look like this:

## rank.py ##

rank.py takes as input a linear reference genome fasta, a variant file in 1ksnp format, and an optional file containing variant phasing information. The user should also specify a ranking method (popcov, popcov-blowup, or hybrid) and a window size.

Optionally the user can provide a specific chromosome to process, otherwise the full genome is processed. When running the hybrid ranking method for window sizes over ~35, we recommend adding the --prune argument limit the amount of blowup in dense regions of variants. A limit of 15 variants per window seems to perform well in practice.

## build.py ##

build.py takes a input a set of ranked variants (produced by rank.py) and a percentage of variants to include in the graph. It produces the necessary input files to build an index with HISAT2 or with Bowtie (ERG).

## Full pipeline

From beginning to end, running the FORGe pipeline might look like this:

> ./vcf_to_1ksnp.py --reference ref.fa --vcf variants.vcf --ingroup names.txt --out variants.1ksnp --individuals phasing.txt
> ./rank.py --method popcov --reference ref.fa --vars variants.1ksnp --window-size 100 --phasing phasing.txt --output ordered.txt
> ./build.py --reference ref.fa --vars variants.1ksnp --window-size 100 --hisat variants.snp --sorted ordered.txt --pct 10
> $HISAT_HOME/hisat2-build --snp variants.snp ref.fa index


