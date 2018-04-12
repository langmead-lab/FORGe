# FORGe tool for ranking variants and building an optimal graph genome

FORGe consists of 2 primary scripts -- `rank.py` and `build.py` -- as well as a helper script `vcf_to_1ksnp.py` for generating input files in the proper format.

## vcf_to_1ksnp.py ##

FORGe takes information about genetic variants in a `1ksnp` and an optional phasing file. The `1ksnp` format is less complex than VCF, but you can convert VCF to using the included `vcf_to_1ksnp.py`.

In a `1ksnp` file, each alternate allele is stored as a separate line with the following columns:

```
Chromosome	Position	Reference Allele	Alternate Allele	Population Frequency	??	# Alternates	Variant Name
```

The phasing file is a comma-separated table with a column for each individual and row for each variant, in the same order they appear in the 1ksnp file. An element contains the allele for the corresponding individual (column) and variant (row), with `0` indicating the reference allele and `k` indicating the kth alternate allele, in order, appearing in the 1ksnp file.

FORGe includes a helper script `vcf_to_1ksnp.py` to facilitate convertsion from a VCF file and set of ingroup individuals to FORGe variant and phasing files.  Ingroup individuals can be specified as a list of names to either include (`--ingroup`) or exclude (`--outgroup`). Sample usage of this script would look like this:

```
TODO
```

## rank.py ##

`rank.py` takes a linear reference genome fasta, a `1ksnp` variant file and an optional file containing phasing information. The user also specifies a model (`popcov`, `popcov-blowup`, or `hybrid`) and window size.  Note that the `hybrid` model produces scores both with and without "blowup avoidance," whereas `popcov` and `popcov-blowup` do this separately.

The user can indicate a specific chromosome to process using `--chrom`.  Otherwise the full genome is used. When running the hybrid ranking method for window sizes over ~35, we recommend adding the `--prune` which limits blowup in regions dense with genetic variants. A limit of 15 variants per window performed well in practice.

## build.py ##

`build.py` takes as input a set of ranked variants (output by `rank.py`) and a percentage of variants to include in the graph. It produces the necessary input files to build an index with HISAT2 or with Bowtie (ERG).

## Full pipeline

From beginning to end, running the FORGe pipeline with HISAT2 might look like this:

```
./vcf_to_1ksnp.py --reference ref.fa --vcf variants.vcf --ingroup names.txt --out variants.1ksnp --individuals phasing.txt
./rank.py --method popcov --reference ref.fa --vars variants.1ksnp --window-size 100 --phasing phasing.txt --output ordered.txt
./build.py --reference ref.fa --vars variants.1ksnp --window-size 100 --hisat variants.snp --sorted ordered.txt --pct 10
$HISAT_HOME/hisat2-build --snp variants.snp ref.fa index
```



