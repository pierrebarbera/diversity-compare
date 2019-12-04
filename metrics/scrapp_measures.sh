#!/bin/bash

BASE=$(cd `dirname "${BASH_SOURCE[0]}"`/.. && pwd)

DATA=${BASE}/datasets/bv
WORKDIR=${BASE}/workdir/hmp2/scrapp
SF=${WORKDIR}/*/summary.newick

#JPDIR=${DATA}/03_place/queries

set -e
shopt -s nullglob

mkdir -p ${WORKDIR}
cd ${WORKDIR}

# export OMP_NUM_THREADS=10

scrapp-diversity ${SF} > result.csv
