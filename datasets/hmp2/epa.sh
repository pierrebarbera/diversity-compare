#!/bin/bash

BASE=$(cd `dirname "${BASH_SOURCE[0]}"` && pwd)

SEQ=${BASE}/query/query.fasta
TREE=${BASE}/ref/Bact_Constr_Backbone.newick
REF=${BASE}/query/reference.fasta
MOD=${BASE}/ref/RAxML_info.infofile
OUT=${BASE}/place

cd $OUT
epa-ng --tree $TREE --msa $REF --query $SEQ --threads 40 --model $MOD "$@"
cd -