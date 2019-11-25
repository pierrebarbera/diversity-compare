#!/bin/bash

BASE=$(cd `dirname "${BASH_SOURCE[0]}"` && pwd)
DATA=${BASE}/data

SEQ=${DATA}/chunks/chunk_0.fasta
TREE=${BASE}/ref/Bact_Constr_Backbone.newick
REF=${BASE}/ref/tax_cons_border.phylip
OUT=${BASE}/query

cd $OUT
papara -t $TREE -s $REF -q $SEQ -r -n align -j 40 "$@"

epa-ng --split ../ref/tax_cons_border.phylip papara_alignment.align

gzip query.fasta
