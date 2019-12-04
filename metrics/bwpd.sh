#!/bin/bash

BASE=$(cd `dirname "${BASH_SOURCE[0]}"`/.. && pwd)

#DS=bv
DS=hmp2

DATA=${BASE}/datasets/${DS}
WORKDIR=${BASE}/workdir/${DS}/bwpd

#JPDIR=${DATA}/03_place/rare
JPDIR=${DATA}/place/samples

set -e
shopt -s nullglob

mkdir -p ${WORKDIR}
cd ${WORKDIR}

jplace-diversity ${JPDIR}/*.jplace > result.csv
