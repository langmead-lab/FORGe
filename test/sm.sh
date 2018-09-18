#!/bin/sh

INTERP=$1
[ -z "${INTERP}" ] && echo "Must specify interpreter as arg" && exit 1
shift

set -ex

echo "*** POP COV ***"
${INTERP} src/rank.py \
    --method popcov \
    --reference test/sm.fa \
    --vars test/sm.1ksnp \
    --window-size 40 \
    --prune 5 \
    --output popcov.txt \
    $*

echo "*** HYBRID w/ phasing ***"
${INTERP} src/rank.py \
    --method hybrid \
    --reference test/sm.fa \
    --vars test/sm.1ksnp \
    --window-size 40 \
    --prune 5 \
    --output hybrid-haplo.txt \
    --phasing test/sm.phasing \
    $*

echo "*** HYBRID w/o phasing ***"
${INTERP} src/rank.py \
    --method hybrid \
    --reference test/sm.fa \
    --vars test/sm.1ksnp \
    --window-size 40 \
    --prune 5 \
    --output hybrid.txt \
    $*
