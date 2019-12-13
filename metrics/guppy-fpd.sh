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


WORKDIR=${BASE}/workdir/${DS}/guppy-fpd

JPDIR=${DATA}/place/samples

set -e
shopt -s nullglob

mkdir -p ${WORKDIR}
cd ${WORKDIR}

guppy fpd --theta 0.0,0.25,0.5,0.75,1.0 --out-dir ${WORKDIR} --csv -o result.csv ${JPDIR}/*.jplace
