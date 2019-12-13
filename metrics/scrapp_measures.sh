#!/bin/bash
BASE=$(cd `dirname "${BASH_SOURCE[0]}"`/.. && pwd)

die () {
    echo >&2 "$@"
    exit 1
}

[ "$#" -eq 1 ] || die "specify dataset!"

DS=$1
DATA=${BASE}/datasets/${DS}

[ ! -d "${DATA}" ] && die "No such dataset."

WORKDIR=${BASE}/workdir/${DS}/scrapp_boot_1000
SF=${WORKDIR}/*/summary.newick

#JPDIR=${DATA}/03_place/queries

set -e
shopt -s nullglob

mkdir -p ${WORKDIR}
cd ${WORKDIR}

# export OMP_NUM_THREADS=10

scrapp-diversity ${SF} > result.csv
