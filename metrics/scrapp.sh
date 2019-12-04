#!/bin/bash

BASE=$(cd `dirname "${BASH_SOURCE[0]}"`/.. && pwd)

## BV data
# DATA=${BASE}/datasets/bv
# WORKDIR=${BASE}/workdir/scrapp

# JPDIR=${DATA}/03_place/rare
# SEQDIR=${DATA}/03_place/queries

## HMP2 data
DATA=${BASE}/datasets/hmp2
WORKDIR=${BASE}/workdir/hmp2/scrapp_boot_1000

mkdir -p ${WORKDIR}

JPDIR=${DATA}/place/samples
SEQDIR=${JPDIR}

# set -e
shopt -s nullglob

N=10
(
for jplace in ${JPDIR}/*.jplace; do

  name=${jplace##*/}
  name=${name%.*}

  OUT=${WORKDIR}/${name}
  QRY=${SEQDIR}/${name}.fasta

  ((i=i%N)); ((i++==0)) && wait
  ~/scrapp/scrapp.py --cluster-above 1000 --bootstrap --jplace ${jplace} --alignment ${QRY} --work-dir ${OUT} --parallel threads --num-threads 4 --min-weight 0.5 &

done
)
