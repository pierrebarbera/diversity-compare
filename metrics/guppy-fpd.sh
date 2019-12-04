#!/bin/bash

BASE=$(cd `dirname "${BASH_SOURCE[0]}"`/.. && pwd)

DS=hmp2
#DS=bv

DATA=${BASE}/datasets/${DS}
WORKDIR=${BASE}/workdir/${DS}/guppy-fpd

#JPDIR=${DATA}/03_place/rare
JPDIR=${DATA}/place/samples

set -e
shopt -s nullglob

mkdir -p ${WORKDIR}
cd ${WORKDIR}

guppy fpd --theta 0.0,0.25,0.5,0.75,1.0 --out-dir ${WORKDIR} --csv -o result.csv ${JPDIR}/*.jplace
