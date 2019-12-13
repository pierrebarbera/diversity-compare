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

WORKDIR=${BASE}/workdir/${DS}/bwpd

JPDIR=${DATA}/place/samples

set -e
shopt -s nullglob

mkdir -p ${WORKDIR}
cd ${WORKDIR}

jplace-diversity ${JPDIR}/*.jplace > result.csv
