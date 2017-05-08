#!/bin/bash

# Simulate a random genome
./generate_genome.py ref.fa

# Add simulated indels to the genome
./simulate_indels.py ref.fa indels.fa

# Simulate reads from genome with indels
./mason.sh indels.fa reads.fastq

# Build HISAT index for ref genome
/scratch0/langmead-fs1/user/jacob/enhanced_reference_genome/software/hisat2/hisat2-build ref.fa index

# Align reads to HISAT index
/scratch0/langmead-fs1/user/jacob/enhanced_reference_genome/software/hisat2/hisat2 --sam-no-qname-trunc -t -k 1 -x index -U reads.fastq -S aligned.sam

# map reads from ref genome coordinates to modified chromosome coordinates
/scratch0/langmead-fs1/user/jacob/vis/src/remap_reads.py aligned.sam indels.txt aligned_remapped.sam

# compute accuracy
echo ''
echo 'Accuracy against reference coordinate system'
cat aligned.sam | /usr/local/bin/python2.7 /scratch0/langmead-fs1/user/jacob/vis/src/correctness/correctness.py

echo ''
echo 'Accuracy against new coordinate system'
cat aligned_remapped.sam | /usr/local/bin/python2.7 /scratch0/langmead-fs1/user/jacob/vis/src/correctness/correctness.py

