#!/bin/bash

export NUM_READS=10000
export LEN=35

/scratch0/langmead-fs1/shared/mason/bin/mason illumina -N ${NUM_READS} -n $LEN -hn 1 -i -sq -o $2 $1
