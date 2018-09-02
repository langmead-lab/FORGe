#!/bin/sh

cat >ingroup.txt <<EOF
HG00133
HG00136
HG00137
HG00138
HG00139
HG00140
HG00141
HG00142
HG00143
HG00145
HG00146
HG00148
HG00149
HG00150
HG00151
HG00154
HG00155
HG00157
EOF

python ../src/vcf_to_1ksnp.py \
    --reference sm.fa \
    --vcf sm.vcf \
    --out small.1ksnp \
    --ingroup ingroup.txt \
    --individuals sm.phasing
