#!gmake

# Start with a relatively dense collection of variants, then try progressively more inclusive SNP files

HSA=$(HOME)/software/hisat2-2.0.5/hisat2
HSB=$(HOME)/software/hisat2-2.0.5/hisat2-build
MASON=$(HOME)/software/mason/mason
FA=$(HOME)/fasta/hg38.fa
GENOME_LEN=1000000
READLEN=75
NREADS=100000
NTHREADS=8
CORRECT_PY=../../src/correctness/correctness.py

.SECONDARY:

.PHONY: cors reads

overall.csv: cors
	pypy tabulate_cor.py > $@

define type_stride_pct
$(1)_base_stride_$(2)_pct_$(3).cor: $(1)_base_stride_$(2)_pct_$(3).sam
	pypy $$(CORRECT_PY) < $$< | awk -v FS='\t' '{print $$$$2,$$$$NF}' | sed 's/ZC:i://' > $$@

$(1)_base_stride_$(2)_pct_$(3)_alts.cor: $(1)_base_stride_$(2)_pct_$(3)_alts.sam
	pypy $$(CORRECT_PY) < $$< | awk -v FS='\t' '{print $$$$2,$$$$NF}' | sed 's/ZC:i://' > $$@

$(1)_base_stride_$(2)_pct_$(3)_alts.sam: $(1)_base_stride_$(2).fq $(1)_base_stride_$(2)_pct_$(3)_alts.1.ht2
	$$(HSA) --no-spliced-alignment --sam-no-qname-trunc -p $$(NTHREADS) -x $(1)_base_stride_$(2)_pct_$(3)_alts -U $$< -S $$@

$(1)_base_stride_$(2)_pct_$(3)_alts.1.ht2: $(1)_base_stride_$(2).fa
	python subsample_snps.py $(3) < $$(<:%.fa=%.snp) > $$(<:%.fa=%_pct_$(3).snp)
	$$(HSB) $$< $$(<:%.fa=%_pct_$(3)_alts) --snp $$(<:%.fa=%_pct_$(3).snp)

cors: $(1)_base_stride_$(2)_pct_$(3)_alts.cor
endef

define type_stride
$(1)_base_stride_$(2).sam: $(1)_base_stride_$(2).fq $(1)_base_stride_$(2).1.ht2
	$$(HSA) --no-spliced-alignment --sam-no-qname-trunc -p $$(NTHREADS) -x $(1)_base_stride_$(2) -U $$< -S $$@

$(1)_base_stride_$(2).fq: $(1)_alt_stride_$(2).fa
	$$(MASON) illumina -hn 2 -i -s 3466 -sq -n $(READLEN) -N $(NREADS) -o $$@ $$<

$(1)_base_stride_$(2).1.ht2: $(1)_base_stride_$(2).fa
	$$(HSB) $$< $$(<:%.fa=%)
endef

TARGETS=snp indel
PCTS=10 20 30 40 60 80 100
STRIDES=900 450

$(foreach target,$(TARGETS),$(foreach pct,$(PCTS),$(foreach stride,$(STRIDES),$(eval $(call type_stride_pct,$(target),$(stride),$(pct))))))
$(foreach target,$(TARGETS),$(foreach stride,$(STRIDES),$(eval $(call type_stride,$(target),$(stride)))))

define stride_group
indel_base_stride_$(1).fa indel_alt_stride_$(1).fa:
	pypy ../util/make_ref.py indel $$(FA) $(1) $$(GENOME_LEN)

snp_base_stride_$(1).fa snp_alt_stride_$(1).fa:
	pypy ../util/make_ref.py snp $$(FA) $(1) $$(GENOME_LEN)

cors: snp_base_stride_$(1).cor indel_base_stride_$(1).cor

reads: snp_base_stride_$(1).fq indel_base_stride_$(1).fq
endef

$(foreach stride,$(STRIDES),$(eval $(call stride_group,$(stride))))

.PHONY: clean
clean:
	rm -f *.fa *.ht2 *.fq *.snp *.sam *.cor
